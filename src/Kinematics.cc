/*
 * RBDL - Rigid Body Dynamics Library
 * Copyright (c) 2011-2012 Martin Felis <martin.felis@iwr.uni-heidelberg.de>
 *
 * Licensed under the zlib license. See LICENSE for more details.
 */

#include <iostream>
#include <limits>
#include <cstring>
#include <assert.h>

#include "rbdl/rbdl_mathutils.h"
#include "rbdl/Logging.h"

#include "rbdl/Model.h"
#include "rbdl/Kinematics.h"

namespace RigidBodyDynamics {

  using namespace Math;

  void UpdateKinematics (Model &model,
                         const VectorNd &Q,
                         const VectorNd &QDot,
                         const VectorNd &QDDot
                         ) {
    LOG << "-------- " << __func__ << " --------" << std::endl;

    unsigned int i;

    SpatialVector spatial_gravity (0., 0., 0., model.gravity[0], model.gravity[1], model.gravity[2]);

    model.a[0].setZero();
    //model.a[0] = spatial_gravity;

    for (i = 1; i < model.mBodies.size(); i++) {
      SpatialTransform X_J;
      SpatialVector v_J;
      SpatialVector c_J;
      Joint joint = model.mJoints[i];
      unsigned int lambda = model.lambda[i];

      jcalc (model, i, X_J, model.S[i], v_J, c_J, Q[i - 1], QDot[i - 1]);

      model.X_lambda[i] = X_J * model.X_T[i];

      if (lambda != 0) {
        model.X_base[i] = model.X_lambda[i] * model.X_base.at(lambda);
        model.v[i] = model.X_lambda[i].apply(model.v[lambda]) + v_J;
        model.c[i] = c_J + crossm(model.v[i],v_J);
        model.a_bias[i] = model.X_lambda[i].apply(model.a_bias[lambda]) + model.c[i];
      }	else {
        model.X_base[i] = model.X_lambda[i];
        model.v[i] = v_J;
        model.c[i].setZero();
        model.a_bias[i].setZero();
      }
      /// @todo optimize this, there are repeated computations
      model.a[i] = model.X_lambda[i].apply(model.a[lambda]) + model.c[i];
      model.a[i] = model.a[i] + model.S[i] * QDDot[i - 1];
    }

    for (i = 1; i < model.mBodies.size(); i++) {
      LOG << "a[" << i << "] = " << model.a[i].transpose() << std::endl;
    }

  }

  void UpdateKinematicsCustom (Model &model,
                               const VectorNd *Q,
                               const VectorNd *QDot,
                               const VectorNd *QDDot
                               ) {
    LOG << "-------- " << __func__ << " --------" << std::endl;

    unsigned int i;

    if (Q) {
      for (i = 1; i < model.mBodies.size(); ++i) {
        SpatialVector v_J;
        SpatialVector c_J;
        SpatialTransform X_J;
        unsigned int lambda = model.lambda[i];

        jcalc (model, i, X_J, model.S[i], v_J, c_J, (*Q)[i - 1], 0.);

        model.X_lambda[i] = X_J * model.X_T[i];

        if (lambda != 0) {
          // using at you get a safer data acces the vector
          model.X_base[i] = model.X_lambda[i] * model.X_base.at(lambda);
        }	else {
          model.X_base[i] = model.X_lambda[i];
        }
      }
    }

    if (QDot) {
      for (i = 1; i < model.mBodies.size(); i++) {
        SpatialVector v_J;
        SpatialVector c_J;
        SpatialTransform X_J;
        unsigned int lambda = model.lambda[i];

        jcalc (model, i, X_J, model.S[i], v_J, c_J, (*Q)[i - 1], (*QDot)[i - 1]);

        if (lambda != 0) {
          model.v[i] = model.X_lambda[i].apply(model.v[lambda]) + v_J;
          model.c[i] = c_J + crossm(model.v[i],v_J);
          model.a_bias[i] = model.X_lambda[i].apply(model.a_bias[lambda]) + model.c[i];
        }	else {
          model.v[i] = v_J;
          model.c[i].setZero();
          model.a_bias[i].setZero();
        }
      }
    }

    if (QDDot) {
      for (i = 1; i < model.mBodies.size(); i++) {
        unsigned int lambda = model.lambda[i];

        if (lambda != 0) {
          model.a[i] = model.X_lambda[i].apply(model.a[lambda]) + model.c[i];
        }	else {
          model.a[i].setZero();
        }

        model.a[i] = model.a[i] + model.S[i] * (*QDDot)[i - 1];
      }
    }
  }

  Vector3d CalcBodyToBaseCoordinates (
      Model &model,
      const VectorNd &Q,
      unsigned int body_id,
      const Vector3d &point_body_coordinates,
      bool update_kinematics) {
    // update the Kinematics if necessary
    if (update_kinematics) {
      UpdateKinematicsCustom (model, &Q, NULL, NULL);
    }

    if (body_id >= model.fixed_body_discriminator) {
      unsigned int fbody_id = body_id - model.fixed_body_discriminator;
      unsigned int parent_id = model.mFixedBodies[fbody_id].mMovableParent;

      Matrix3d fixed_rotation = model.mFixedBodies[fbody_id].mParentTransform.E.transpose();
      Vector3d fixed_position = model.mFixedBodies[fbody_id].mParentTransform.r;

      Matrix3d parent_body_rotation = model.X_base[parent_id].E.transpose();
      Vector3d parent_body_position = model.X_base[parent_id].r;
      return parent_body_position + parent_body_rotation * (fixed_position + fixed_rotation * (point_body_coordinates));
    }

    Matrix3d body_rotation = model.X_base[body_id].E.transpose();
    Vector3d body_position = model.X_base[body_id].r;

    return body_position + body_rotation * point_body_coordinates;
  }

  Vector3d CalcBaseToBodyCoordinates (
      Model &model,
      const VectorNd &Q,
      unsigned int body_id,
      const Vector3d &point_base_coordinates,
      bool update_kinematics) {
    if (update_kinematics) {
      UpdateKinematicsCustom (model, &Q, NULL, NULL);
    }

    if (body_id >= model.fixed_body_discriminator) {
      unsigned int fbody_id = body_id - model.fixed_body_discriminator;
      unsigned int parent_id = model.mFixedBodies[fbody_id].mMovableParent;

      Matrix3d fixed_rotation = model.mFixedBodies[fbody_id].mParentTransform.E;
      Vector3d fixed_position = model.mFixedBodies[fbody_id].mParentTransform.r;

      Matrix3d parent_body_rotation = model.X_base[parent_id].E;
      Vector3d parent_body_position = model.X_base[parent_id].r;

      return fixed_rotation * ( - fixed_position - parent_body_rotation * (parent_body_position - point_base_coordinates));
    }

    Matrix3d body_rotation = model.X_base[body_id].E;
    Vector3d body_position = model.X_base[body_id].r;

    return body_rotation * (point_base_coordinates - body_position);
  }

  Matrix3d CalcBodyWorldOrientation (
      Model &model,
      const VectorNd &Q,
      const unsigned int body_id,
      bool update_kinematics)
  {
    // update the Kinematics if necessary
    if (update_kinematics) {
      UpdateKinematicsCustom (model, &Q, NULL, NULL);
    }

    if (body_id >= model.fixed_body_discriminator) {
      unsigned int fbody_id = body_id - model.fixed_body_discriminator;
      model.mFixedBodies[fbody_id].mBaseTransform = model.mFixedBodies[fbody_id].mParentTransform*model.X_base[model.mFixedBodies[fbody_id].mMovableParent];

      return model.mFixedBodies[fbody_id].mBaseTransform.E;
    }

    return model.X_base[body_id].E;
  }

  void CalcPointJacobian (
      Model &model,
      const VectorNd &Q,
      unsigned int body_id,
      const Vector3d &point_position,
      MatrixNd &G,
      bool update_kinematics
      ) {
    LOG << "-------- " << __func__ << " --------" << std::endl;

    // update the Kinematics if necessary
    if (update_kinematics) {
      UpdateKinematicsCustom (model, &Q, NULL, NULL);
    }

    Vector3d point_base_pos = CalcBodyToBaseCoordinates (model, Q, body_id, point_position, false);
    SpatialMatrix point_trans = Xtrans_mat (point_base_pos);

    assert (G.rows() == 3 && G.cols() == model.dof_count );

    G.setZero();

    // we have to make sure that only the joints that contribute to the
    // bodies motion also get non-zero columns in the jacobian.
    // VectorNd e = VectorNd::Zero(Q.size() + 1);
    char *e = new char[Q.size() + 1];
    if (e == NULL) {
      std::cerr << "Error: allocating memory." << std::endl;
      abort();
    }
    memset (&e[0], 0, Q.size() + 1);

    unsigned int reference_body_id = body_id;

    if (model.IsFixedBodyId(body_id)) {
      unsigned int fbody_id = body_id - model.fixed_body_discriminator;
      reference_body_id = model.mFixedBodies[fbody_id].mMovableParent;
    }

    unsigned int j = reference_body_id;

    // e[j] is set to 1 if joint j contributes to the jacobian that we are
    // computing. For all other joints the column will be zero.
    while (j != 0) {
      e[j] = 1;
      j = model.lambda[j];
    }

    for (j = 1; j < model.mBodies.size(); j++) {
      if (e[j] == 1) {
        SpatialVector S_base;
        S_base = point_trans * spatial_inverse(model.X_base[j].toMatrix()) * model.S[j];

        G(0, j - 1) = S_base[3];
        G(1, j - 1) = S_base[4];
        G(2, j - 1) = S_base[5];
      }
    }

    delete[] e;
  }

  void CalcOrientationJacobian (
      Model &model,
      const VectorNd &Q,
      unsigned int body_id,
      const Vector3d &point_position,
      MatrixNd &G,
      bool update_kinematics
      ) {
    LOG << "-------- " << __func__ << " --------" << std::endl;

    // update the Kinematics if necessary
    if (update_kinematics) {
      UpdateKinematicsCustom (model, &Q, NULL, NULL);
    }

    Vector3d point_base_pos = CalcBodyToBaseCoordinates (model, Q, body_id, point_position, false);
    SpatialMatrix point_trans = Xtrans_mat (point_base_pos);

    assert (G.rows() == 3 && G.cols() == model.dof_count );

    G.setZero();

    // we have to make sure that only the joints that contribute to the
    // bodies motion also get non-zero columns in the jacobian.
    // VectorNd e = VectorNd::Zero(Q.size() + 1);
    char *e = new char[Q.size() + 1];
    if (e == NULL) {
      std::cerr << "Error: allocating memory." << std::endl;
      abort();
    }
    memset (&e[0], 0, Q.size() + 1);

    unsigned int reference_body_id = body_id;

    if (model.IsFixedBodyId(body_id)) {
      unsigned int fbody_id = body_id - model.fixed_body_discriminator;
      reference_body_id = model.mFixedBodies[fbody_id].mMovableParent;
    }

    unsigned int j = reference_body_id;

    // e[j] is set to 1 if joint j contributes to the jacobian that we are
    // computing. For all other joints the column will be zero.
    while (j != 0) {
      e[j] = 1;
      j = model.lambda[j];
    }

    for (j = 1; j < model.mBodies.size(); j++) {
      if (e[j] == 1) {
        SpatialVector S_base;
        S_base = point_trans * spatial_inverse(model.X_base[j].toMatrix()) * model.S[j];

        G(0, j - 1) = S_base[0];
        G(1, j - 1) = S_base[1];
        G(2, j - 1) = S_base[2];
      }
    }
    delete[] e;
  }

  void CalcPointJacobianGlobalFrame (
      Model &model,
      const VectorNd &Q,
      unsigned int body_id,
      const Vector3d &point_position,
      MatrixNd &G,
      bool update_kinematics
      ) {
    LOG << "-------- " << __func__ << " --------" << std::endl;

    // update the Kinematics if necessary
    if (update_kinematics) {
      UpdateKinematicsCustom (model, &Q, NULL, NULL);
    }

    assert (G.rows() == 3 && G.cols() == model.dof_count );

    G.setZero();

    // we have to make sure that only the joints that contribute to the
    // bodies motion also get non-zero columns in the jacobian.
    // VectorNd e = VectorNd::Zero(Q.size() + 1);
    char *e = new char[Q.size() + 1];
    if (e == NULL) {
      std::cerr << "Error: allocating memory." << std::endl;
      abort();
    }
    memset (&e[0], 0, Q.size() + 1);

    unsigned int reference_body_id = body_id;

    if (model.IsFixedBodyId(body_id)) {
      unsigned int fbody_id = body_id - model.fixed_body_discriminator;
      reference_body_id = model.mFixedBodies[fbody_id].mMovableParent;
    }

    unsigned int j = reference_body_id;

    // e[j] is set to 1 if joint j contributes to the jacobian that we are
    // computing. For all other joints the column will be zero.
    while (j != 0) {
      e[j] = 1;
      j = model.lambda[j];
    }

    for (j = 1; j < model.mBodies.size(); j++) {
      if (e[j] == 1) {
        SpatialVector S_base;
        S_base = spatial_inverse(model.X_base[j].toMatrix()) * model.S[j];

        G(0, j - 1) = S_base[3];
        G(1, j - 1) = S_base[4];
        G(2, j - 1) = S_base[5];
      }
    }

    delete[] e;
  }


  void CalcPoseJacobian (
      Model &model,
      const VectorNd &Q,
      unsigned int body_id,
      const Vector3d &point_position,
      MatrixNd &G,
      bool update_kinematics
      ) {
    LOG << "-------- " << __func__ << " --------" << std::endl;

    // update the Kinematics if necessary
    if (update_kinematics) {
      UpdateKinematicsCustom (model, &Q, NULL, NULL);
    }

    Vector3d point_base_pos = CalcBodyToBaseCoordinates (model, Q, body_id, point_position, false);
    SpatialMatrix point_trans = Xtrans_mat (point_base_pos);

    assert (G.rows() == 6 && G.cols() == model.dof_count );

    G.setZero();

    // we have to make sure that only the joints that contribute to the
    // bodies motion also get non-zero columns in the jacobian.
    // VectorNd e = VectorNd::Zero(Q.size() + 1);
    char *e = new char[Q.size() + 1];
    if (e == NULL) {
      std::cerr << "Error: allocating memory." << std::endl;
      abort();
    }
    memset (&e[0], 0, Q.size() + 1);

    unsigned int reference_body_id = body_id;

    if (model.IsFixedBodyId(body_id)) {
      unsigned int fbody_id = body_id - model.fixed_body_discriminator;
      reference_body_id = model.mFixedBodies[fbody_id].mMovableParent;
    }

    unsigned int j = reference_body_id;

    // e[j] is set to 1 if joint j contributes to the jacobian that we are
    // computing. For all other joints the column will be zero.
    while (j != 0) {
      e[j] = 1;
      j = model.lambda[j];
    }

    for (j = 1; j < model.mBodies.size(); j++) {
      if (e[j] == 1) {
        SpatialVector S_base;
        S_base =  spatial_inverse(model.X_base[j].toMatrix()) * model.S[j];

        G(0, j - 1) = S_base[0];
        G(1, j - 1) = S_base[1];
        G(2, j - 1) = S_base[2];
        G(3, j - 1) = S_base[3];
        G(4, j - 1) = S_base[4];
        G(5, j - 1) = S_base[5];
      }
    }

    //Optimized transformation of the spatial jacobian to be expressed in the coordinate frame of the tip the jacobian (Same as tranlation from the global coordinate
    //frame to the tip
    G = point_trans *G;
    delete[] e;
  }

  void CalcPoseJacobianGlobalFrame (
      Model &model,
      const VectorNd &Q,
      unsigned int body_id,
      MatrixNd &G,
      bool update_kinematics
      ) {
    LOG << "-------- " << __func__ << " --------" << std::endl;

    // update the Kinematics if necessary
    if (update_kinematics) {
      UpdateKinematicsCustom (model, &Q, NULL, NULL);
    }

    assert (G.rows() == 6 && G.cols() == model.dof_count );

    G.setZero();

    // we have to make sure that only the joints that contribute to the
    // bodies motion also get non-zero columns in the jacobian.
    // VectorNd e = VectorNd::Zero(Q.size() + 1);
    char *e = new char[Q.size() + 1];
    if (e == NULL) {
      std::cerr << "Error: allocating memory." << std::endl;
      abort();
    }
    memset (&e[0], 0, Q.size() + 1);

    unsigned int reference_body_id = body_id;

    if (model.IsFixedBodyId(body_id)) {
      unsigned int fbody_id = body_id - model.fixed_body_discriminator;
      reference_body_id = model.mFixedBodies[fbody_id].mMovableParent;
    }

    unsigned int j = reference_body_id;

    // e[j] is set to 1 if joint j contributes to the jacobian that we are
    // computing. For all other joints the column will be zero.
    while (j != 0) {
      e[j] = 1;
      j = model.lambda[j];
    }

    for (j = 1; j < model.mBodies.size(); j++) {
      if (e[j] == 1) {
        SpatialVector S_base;
        S_base =  spatial_inverse(model.X_base[j].toMatrix()) * model.S[j];

        G(0, j - 1) = S_base[0];
        G(1, j - 1) = S_base[1];
        G(2, j - 1) = S_base[2];
        G(3, j - 1) = S_base[3];
        G(4, j - 1) = S_base[4];
        G(5, j - 1) = S_base[5];
      }
    }

    delete[] e;
  }

  void CalcPoseJacobianBodyFrame (
      Model &model,
      const VectorNd &Q,
      unsigned int body_id,
      const Vector3d &point_position,
      MatrixNd &G,
      bool update_kinematics
      ) {
    LOG << "-------- " << __func__ << " --------" << std::endl;

    // update the Kinematics if necessary
    if (update_kinematics) {
      UpdateKinematicsCustom (model, &Q, NULL, NULL);
    }

    /// @todo This has to be finished if its desired to have a body jacobian offsetted by some amount
    //Vector3d point_base_pos = CalcBaseToBodyCoordinates(model, Q, body_id, point_position, false);
    Vector3d point_base_pos; point_base_pos.setZero();
    SpatialMatrix point_trans = Xtrans_mat (-point_base_pos);

    Eigen::Vector3d zero; zero.setZero();
    Vector3d body_r = CalcBodyToBaseCoordinates(model, Q, body_id, zero, false);
    Matrix3d body_E = CalcBodyWorldOrientation(model, Q, body_id, false);
    SpatialTransform body_trans(body_E, body_r);

    assert (G.rows() == 6 && G.cols() == model.dof_count );

    G.setZero();

    // we have to make sure that only the joints that contribute to the
    // bodies motion also get non-zero columns in the jacobian.
    // VectorNd e = VectorNd::Zero(Q.size() + 1);
    char *e = new char[Q.size() + 1];
    if (e == NULL) {
      std::cerr << "Error: allocating memory." << std::endl;
      abort();
    }
    memset (&e[0], 0, Q.size() + 1);

    unsigned int reference_body_id = body_id;

    if (model.IsFixedBodyId(body_id)) {
      unsigned int fbody_id = body_id - model.fixed_body_discriminator;
      reference_body_id = model.mFixedBodies[fbody_id].mMovableParent;
    }

    unsigned int j = reference_body_id;

    // e[j] is set to 1 if joint j contributes to the jacobian that we are
    // computing. For all other joints the column will be zero.
    while (j != 0) {
      e[j] = 1;
      j = model.lambda[j];
    }

    for (j = 1; j < model.mBodies.size(); j++) {
      if (e[j] == 1) {
        SpatialVector S_base;
        S_base =  spatial_inverse(model.X_base[j].toMatrix()) * model.S[j];

        G(0, j - 1) = S_base[0];
        G(1, j - 1) = S_base[1];
        G(2, j - 1) = S_base[2];
        G(3, j - 1) = S_base[3];
        G(4, j - 1) = S_base[4];
        G(5, j - 1) = S_base[5];
      }
    }

    //Optimized transformation of the spatial jacobian to be expressed in the coordinate frame of the tip the jacobian (Same as tranlation from the global coordinate
    //frame to the tip
    G = point_trans*(body_trans.toMatrix())*G;
    //G = point_trans*(model.X_base[body_id].toMatrix())*G;

    delete[] e;
  }

  /* THIS DOES NOT WORK DO TO THAT SPATIAL ACCELERATIONS DO NOT MAP TO 3D ACCELERATIONS AS EASY AS SPATIAL VELOCITIES DO
  void CalcPoseJacobianDot (
      Model &model,
      const VectorNd &Q,
      const VectorNd &QDot,
      unsigned int body_id,
      const Vector3d &point_position,
      MatrixNd &G,
      bool update_kinematics
      ) {

    // update the Kinematics if necessary
    if (update_kinematics) {
      UpdateKinematicsCustom (model, &Q, &QDot, NULL);
    }

    Vector3d point_base_pos = CalcBodyToBaseCoordinates (model, Q, body_id, point_position, false);
    SpatialMatrix point_trans = Xtrans_mat (point_base_pos);

    assert (G.rows() == 6 && G.cols() == model.dof_count );

    G.setZero();

    // we have to make sure that only the joints that contribute to the
    // bodies motion also get non-zero columns in the jacobian.
    // VectorNd e = VectorNd::Zero(Q.size() + 1);
    char *e = new char[Q.size() + 1];
    if (e == NULL) {
      std::cerr << "Error: allocating memory." << std::endl;
      abort();
    }
    memset (&e[0], 0, Q.size() + 1);

    unsigned int reference_body_id = body_id;

    if (model.IsFixedBodyId(body_id)) {
      unsigned int fbody_id = body_id - model.fixed_body_discriminator;
      reference_body_id = model.mFixedBodies[fbody_id].mMovableParent;
    }

    unsigned int j = reference_body_id;

    // e[j] is set to 1 if joint j contributes to the jacobian that we are
    // computing. For all other joints the column will be zero.
    while (j != 0) {
      e[j] = 1;
      j = model.lambda[j];
    }

    for (j = 1; j < model.mBodies.size(); j++) {
      if (e[j] == 1) {

        /// @todo optimize this
        SpatialVector body_global_velocity (spatial_inverse(model.X_base[j].toMatrix()) * model.v[j]);
        SpatialVector S_dot_base;
        S_dot_base =   crossm(body_global_velocity,  spatial_inverse(model.X_base[j].toMatrix()) * model.S[j] );
        //S_dot_base =   crossm(body_global_velocity,  spatial_inverse(model.X_base[j].toMatrix()) * model.S[j] );

        G(0, j - 1) = S_dot_base[0];
        G(1, j - 1) = S_dot_base[1];
        G(2, j - 1) = S_dot_base[2];
        G(3, j - 1) = S_dot_base[3];
        G(4, j - 1) = S_dot_base[4];
        G(5, j - 1) = S_dot_base[5];
      }
    }

    //Optimized transformation of the spatial jacobian to be expressed in the coordinate frame of the tip the jacobian (Same as tranlation from the global coordinate
    //frame to the tip
    G = point_trans *G;
    delete[] e;
  }
*/

  void CalcCOMJacobian_ineficient (
      Model &model,
      const VectorNd &Q,
      MatrixNd &COMJ,
      bool update_kinematics
      ) {
    
    assert (COMJ.rows() == 3 && COMJ.cols() == model.dof_count );

    // update the Kinematics if necessary
    if (update_kinematics) {
      UpdateKinematicsCustom (model, &Q, NULL, NULL);
    }

    double total_mass = 0;
    for (unsigned int j = 1; j < model.mBodies.size(); j++) {
      double mass_link = model.mBodies[j].mMass;
      total_mass += mass_link;
    }

    COMJ.setZero();
    Vector3d zero;
    zero.setZero();
    Eigen::MatrixXd jacobian_link(3, model.dof_count);
    for (unsigned int j = 1; j < model.mBodies.size(); ++j) {

      double mass_link = model.mBodies[j].mMass;
      /// @todo: make optinional to not especify extra displacement of the tip
      CalcPointJacobian(model, Q, j,
                        model.mBodies[j].mCenterOfMass, jacobian_link, false);
      COMJ += mass_link*jacobian_link;
    }
    COMJ = COMJ/total_mass;
  }

  // We have to acumulate the spatial transforms, it seems that the cross product of acumulationg the
  // com displacements is a wrong assumption
  void CalcAcumulatedMass(Model &model, const VectorNd &Q){
    assert(model.acumulated_mass.size() == model.mBodies.size());

    for(unsigned int i = 1; i<model.mBodies.size(); ++i){
      Vector3d comi = CalcBodyToBaseCoordinates(model, Q, i, model.mBodies[i].mCenterOfMass, false);
      model.acumulated_mass[i] = model.mBodies[i].mMass*Xtrans_mat(comi);
    }

    for (unsigned int i = model.mBodies.size() - 1; i > 0; i--) {
      unsigned int lambda = model.lambda[i];
      model.acumulated_mass[lambda] += model.acumulated_mass[i];
    }
  }


  void CalcCOMJacobian (
      Model &model,
      const VectorNd &Q,
      MatrixNd &COMJ,
      bool update_kinematics
      ){

    assert (COMJ.rows() == 3 && COMJ.cols() == model.dof_count );

    COMJ.setZero();
    Vector3d zero;
    zero.setZero();

    // update the Kinematics if necessary
    if (update_kinematics) {
      UpdateKinematicsCustom (model, &Q, NULL, NULL);
    }

    double total_mass = 0;
    for (unsigned int j = 1; j < model.mBodies.size(); j++) {
      double mass_link = model.mBodies[j].mMass;
      total_mass += mass_link;
    }

    CalcAcumulatedMass(model, Q);

   
    //Compute the total mass in each link of the acumulated subtree links
    
    
    //Compute all the screws of the system
    
    
    //Combine the screws 
   
    for (unsigned int j = 1; j < model.mBodies.size(); j++) {
      SpatialVector S_base;

      S_base =  model.acumulated_mass[j]*spatial_inverse(model.X_base[j].toMatrix()) * model.S[j];

      //In the spactial vector angular comes first then linear part
      COMJ(0, j - 1) = S_base[3];
      COMJ(1, j - 1) = S_base[4];
      COMJ(2, j - 1) = S_base[5];
    }

    COMJ = COMJ/total_mass;

  }

  Vector3d CalCOM(Model &model,
                  const VectorNd &Q,
                  bool update_kinematics){

    // update the Kinematics if necessary
    if (update_kinematics) {
      UpdateKinematicsCustom (model, &Q, NULL, NULL);
    }
    Eigen::Vector3d com;
    com.setZero();
    double total_mass = 0;

    for(unsigned int i=1; i<model.mBodies.size(); ++i){
      int body_id = i;//rbdl_model_.GetBodyId(link_names_[i].c_str());
      Eigen::Vector3d link_com;
      link_com = CalcBodyToBaseCoordinates(model, Q,
                                           body_id, model.mBodies[body_id].mCenterOfMass, false);

      com += model.mBodies[body_id].mMass*link_com;
      total_mass +=  model.mBodies[body_id].mMass;
    }
    return  com/total_mass;
  }

  Vector3d CalCOMVelocity(Model &model,
                          const VectorNd &Q,
                          const VectorNd &QDot,
                          bool update_kinematics){

    // update the Kinematics if necessary
    if (update_kinematics) {
      UpdateKinematicsCustom (model, &Q, &QDot, NULL);
    }
    Eigen::Vector3d com_vel;
    com_vel.setZero();
    double total_mass = 0;

    for(unsigned int i=1; i<model.mBodies.size(); ++i){
      int body_id = i;//rbdl_model_.GetBodyId(link_names_[i].c_str());
      Eigen::Vector3d link_com_vel;
      link_com_vel = CalcPointVelocity(model, Q, QDot,
                                       body_id, model.mBodies[body_id].mCenterOfMass, false);

      com_vel += model.mBodies[body_id].mMass*link_com_vel;
      total_mass +=  model.mBodies[body_id].mMass;
    }
    return  com_vel/total_mass;
  }

  Math::Vector3d CalCOMAcceleartion(Model &model,
                                    const Math::VectorNd &Q,
                                    const Math::VectorNd &QDot,
                                    const Math::VectorNd &QDDot,
                                    bool update_kinematics){

    // update the Kinematics if necessary
    if (update_kinematics) {
      UpdateKinematics(model, Q, QDot, QDDot);
    }
    Eigen::Vector3d com_acc;
    com_acc.setZero();
    double total_mass = 0;

    for(unsigned int i=1; i<model.mBodies.size(); ++i){
      int body_id = i;
      Eigen::Vector3d link_com_acc;
      link_com_acc = CalcPointAcceleration(model, Q, QDot, QDDot,
                                           body_id, model.mBodies[body_id].mCenterOfMass, false);

      com_acc += model.mBodies[body_id].mMass*link_com_acc;
      total_mass +=  model.mBodies[body_id].mMass;
    }
    return  com_acc/total_mass;

  }

  Math::Vector3d CalCOMAccelerationBias(Model &model,
                                                const Math::VectorNd &Q,
                                                const Math::VectorNd &QDot,
                                                bool update_kinematics){

    // update the Kinematics if necessary
    if (update_kinematics) {
      UpdateKinematicsCustom(model, &Q, &QDot, NULL);
    }
    Eigen::Vector3d com_acc_bias;
    com_acc_bias.setZero();
    double total_mass = 0;

    for(unsigned int i=1; i<model.mBodies.size(); ++i){
      int body_id = i;
      Eigen::Vector3d link_com_acc_bias;
      link_com_acc_bias = CalcPointAccelerationBias(model, Q, QDot,
                                                    body_id, model.mBodies[body_id].mCenterOfMass, false);

      com_acc_bias += model.mBodies[body_id].mMass*link_com_acc_bias;
      total_mass +=  model.mBodies[body_id].mMass;
    }
    return  com_acc_bias/total_mass;
  }

  Vector3d CalcPointVelocity (
      Model &model,
      const VectorNd &Q,
      const VectorNd &QDot,
      unsigned int body_id,
      const Vector3d &point_position,
      bool update_kinematics
      ) {

    // Reset the velocity of the root body
    model.v[0].setZero();

    // update the Kinematics with zero acceleration
    if (update_kinematics) {
      UpdateKinematicsCustom (model, &Q, &QDot, NULL);
    }

    Vector3d point_abs_pos = CalcBodyToBaseCoordinates (model, Q, body_id, point_position, false);

    unsigned int reference_body_id = body_id;

    if (model.IsFixedBodyId(body_id)) {
      unsigned int fbody_id = body_id - model.fixed_body_discriminator;
      reference_body_id = model.mFixedBodies[fbody_id].mMovableParent;
    }

    // Now we can compute the spatial velocity at the given point
    //	SpatialVector body_global_velocity (global_velocities.at(body_id));
    SpatialVector point_spatial_velocity = Xtrans_mat (point_abs_pos) * spatial_inverse(model.X_base[reference_body_id].toMatrix()) * model.v[reference_body_id];

    return Vector3d (
          point_spatial_velocity[3],
          point_spatial_velocity[4],
          point_spatial_velocity[5]
          );
  }

  Vector3d CalcPointAngularVelocity (
      Model &model,
      const VectorNd &Q,
      const VectorNd &QDot,
      unsigned int body_id,
      const Vector3d &point_position,
      bool update_kinematics
      ) {
    LOG << "-------- " << __func__ << " --------" << std::endl;
    assert (model.IsBodyId(body_id));
    assert (model.mBodies.size() == Q.size() + 1);
    assert (model.mBodies.size() == QDot.size() + 1);

    // Reset the velocity of the root body
    model.v[0].setZero();

    // update the Kinematics with zero acceleration
    if (update_kinematics) {
      UpdateKinematicsCustom (model, &Q, &QDot, NULL);
    }

    Vector3d point_abs_pos = CalcBodyToBaseCoordinates (model, Q, body_id, point_position, false);

    unsigned int reference_body_id = body_id;

    if (model.IsFixedBodyId(body_id)) {
      unsigned int fbody_id = body_id - model.fixed_body_discriminator;
      reference_body_id = model.mFixedBodies[fbody_id].mMovableParent;
    }

    LOG << "body_index     = " << body_id << std::endl;
    LOG << "point_pos      = " << point_position.transpose() << std::endl;
    //	LOG << "global_velo    = " << global_velocities.at(body_id) << std::endl;
    LOG << "body_transf    = " << std::endl << model.X_base[reference_body_id].toMatrix() << std::endl;
    LOG << "point_abs_ps   = " << point_abs_pos.transpose() << std::endl;
    LOG << "X   = " << std::endl << Xtrans_mat (point_abs_pos) * spatial_inverse(model.X_base[reference_body_id].toMatrix()) << std::endl;
    LOG << "v   = " << model.v[reference_body_id].transpose() << std::endl;

    // Now we can compute the spatial velocity at the given point
    //	SpatialVector body_global_velocity (global_velocities.at(body_id));
    SpatialVector point_spatial_velocity = Xtrans_mat (point_abs_pos) * spatial_inverse(model.X_base[reference_body_id].toMatrix()) * model.v[reference_body_id];

    LOG << "point_angular_velocity = " <<	Vector3d (
             point_spatial_velocity[0],
             point_spatial_velocity[1],
             point_spatial_velocity[2]
             ) << std::endl;

    return Vector3d (
          point_spatial_velocity[0],
          point_spatial_velocity[1],
          point_spatial_velocity[2]
          );
  }


  Vector3d CalcPointAcceleration (
      Model &model,
      const VectorNd &Q,
      const VectorNd &QDot,
      const VectorNd &QDDot,
      unsigned int body_id,
      const Vector3d &point_position,
      bool update_kinematics
      )
  {
    LOG << "-------- " << __func__ << " --------" << std::endl;

    // Reset the velocity of the root body
    model.v[0].setZero();
    model.a[0].setZero();

    if (update_kinematics)
      UpdateKinematics (model, Q, QDot, QDDot);

    LOG << std::endl;

    unsigned int reference_body_id = body_id;
    Vector3d reference_point = point_position;

    if (model.IsFixedBodyId(body_id)) {
      unsigned int fbody_id = body_id - model.fixed_body_discriminator;
      reference_body_id = model.mFixedBodies[fbody_id].mMovableParent;
      Vector3d base_coords = CalcBodyToBaseCoordinates (model, Q, body_id, point_position, false);
      reference_point = CalcBaseToBodyCoordinates (model, Q, reference_body_id, base_coords, false);
    }

    SpatialVector body_global_velocity (spatial_inverse(model.X_base[reference_body_id].toMatrix()) * model.v[reference_body_id]);

    LOG << " orientation " << std::endl << CalcBodyWorldOrientation (model, Q, reference_body_id, false) << std::endl;
    LOG << " orientationT " << std::endl <<  CalcBodyWorldOrientation (model, Q, reference_body_id, false).transpose() << std::endl;

    Matrix3d global_body_orientation_inv = CalcBodyWorldOrientation (model, Q, reference_body_id, false).transpose();

    SpatialTransform p_X_i (global_body_orientation_inv, reference_point);

    LOG << "p_X_i = " << std::endl << p_X_i << std::endl;

    SpatialVector p_v_i = p_X_i.apply(model.v[reference_body_id]);
    SpatialVector p_a_i = p_X_i.apply(model.a[reference_body_id]);

    SpatialVector frame_acceleration =
        crossm( SpatialVector(0., 0., 0., p_v_i[3], p_v_i[4], p_v_i[5]), (body_global_velocity));

    /// @todo all loggin statements should be removed, they are inneficient like this
    LOG << "v_i                = " << model.v[reference_body_id].transpose() << std::endl;
    LOG << "a_i                = " << model.a[reference_body_id].transpose() << std::endl;
    LOG << "p_X_i              = " << std::endl << p_X_i << std::endl;
    LOG << "p_v_i              = " << p_v_i.transpose() << std::endl;
    LOG << "p_a_i              = " << p_a_i.transpose() << std::endl;
    LOG << "body_global_vel    = " << body_global_velocity.transpose() << std::endl;
    LOG << "frame_acceleration = " << frame_acceleration.transpose() << std::endl;

    SpatialVector p_a_i_dash = p_a_i - frame_acceleration;

    LOG << "point_acceleration = " <<	Vector3d (
             p_a_i_dash[3],
             p_a_i_dash[4],
             p_a_i_dash[5]
             ).transpose() << std::endl;

    return Vector3d (
          p_a_i_dash[3],
          p_a_i_dash[4],
          p_a_i_dash[5]
          );
  }

  Vector3d CalcPointAccelerationBias (
      Model &model,
      const VectorNd &Q,
      const VectorNd &QDot,
      unsigned int body_id,
      const Vector3d &point_position,
      bool update_kinematics
      )
  {
    LOG << "-------- " << __func__ << " --------" << std::endl;

    // Reset the velocity of the root body
    /// @todo why are we reseting this?
    model.v[0].setZero();
    model.a[0].setZero();
    model.a_bias[0].setZero();

    if (update_kinematics)
      UpdateKinematicsCustom(model, &Q, &QDot, NULL);

    LOG << std::endl;

    unsigned int reference_body_id = body_id;
    Vector3d reference_point = point_position;

    if (model.IsFixedBodyId(body_id)) {
      unsigned int fbody_id = body_id - model.fixed_body_discriminator;
      reference_body_id = model.mFixedBodies[fbody_id].mMovableParent;
      Vector3d base_coords = CalcBodyToBaseCoordinates (model, Q, body_id, point_position, false);
      reference_point = CalcBaseToBodyCoordinates (model, Q, reference_body_id, base_coords, false);
    }

    SpatialVector body_global_velocity (spatial_inverse(model.X_base[reference_body_id].toMatrix()) * model.v[reference_body_id]);

    LOG << " orientation " << std::endl << CalcBodyWorldOrientation (model, Q, reference_body_id, false) << std::endl;
    LOG << " orientationT " << std::endl <<  CalcBodyWorldOrientation (model, Q, reference_body_id, false).transpose() << std::endl;

    Matrix3d global_body_orientation_inv = CalcBodyWorldOrientation (model, Q, reference_body_id, false).transpose();

    SpatialTransform p_X_i (global_body_orientation_inv, reference_point);

    LOG << "p_X_i = " << std::endl << p_X_i << std::endl;

    SpatialVector p_v_i = p_X_i.apply(model.v[reference_body_id]);
    SpatialVector p_a_i = p_X_i.apply(model.a_bias[reference_body_id]);

    SpatialVector frame_acceleration =
        crossm( SpatialVector(0., 0., 0., p_v_i[3], p_v_i[4], p_v_i[5]), (body_global_velocity));

    SpatialVector p_a_i_dash = p_a_i - frame_acceleration;

    LOG << "point_acceleration = " <<	Vector3d (
             p_a_i_dash[3],
             p_a_i_dash[4],
             p_a_i_dash[5]
             ).transpose() << std::endl;

    return Vector3d (
          p_a_i_dash[3],
          p_a_i_dash[4],
          p_a_i_dash[5]
          );
  }

  SpatialVector CalcPoseAccelerationBias (
      Model &model,
      const VectorNd &Q,
      const VectorNd &QDot,
      unsigned int body_id,
      const Vector3d &point_position,
      bool update_kinematics
      )
  {
    LOG << "-------- " << __func__ << " --------" << std::endl;

    // Reset the velocity of the root body
    /// @todo why are we reseting this?
    model.v[0].setZero();
    model.a[0].setZero();
    model.a_bias[0].setZero();

    if (update_kinematics)
      UpdateKinematicsCustom(model, &Q, &QDot, NULL);

    LOG << std::endl;

    unsigned int reference_body_id = body_id;
    Vector3d reference_point = point_position;

    if (model.IsFixedBodyId(body_id)) {
      unsigned int fbody_id = body_id - model.fixed_body_discriminator;
      reference_body_id = model.mFixedBodies[fbody_id].mMovableParent;
      Vector3d base_coords = CalcBodyToBaseCoordinates (model, Q, body_id, point_position, false);
      reference_point = CalcBaseToBodyCoordinates (model, Q, reference_body_id, base_coords, false);
    }

    SpatialVector body_global_velocity (spatial_inverse(model.X_base[reference_body_id].toMatrix()) * model.v[reference_body_id]);

    Matrix3d global_body_orientation_inv = CalcBodyWorldOrientation (model, Q, reference_body_id, false).transpose();

    SpatialTransform p_X_i (global_body_orientation_inv, reference_point);

    SpatialVector p_v_i = p_X_i.apply(model.v[reference_body_id]);
    SpatialVector p_a_i = p_X_i.apply(model.a_bias[reference_body_id]);

    SpatialVector frame_acceleration =
        crossm( SpatialVector(0., 0., 0., p_v_i[3], p_v_i[4], p_v_i[5]), (body_global_velocity));

    SpatialVector p_a_i_dash = p_a_i - frame_acceleration;

    Eigen::Matrix<double, 6, 1> res = p_a_i_dash;
    return res;
  }

  SpatialVector CalcPoseAccelerationBiasGlobalFrame (
      Model &model,
      const VectorNd &Q,
      const VectorNd &QDot,
      unsigned int body_id,
      const Vector3d &point_position,
      bool update_kinematics
      )
  {
    LOG << "-------- " << __func__ << " --------" << std::endl;

    // Reset the velocity of the root body
    /// @todo why are we reseting this?
    model.v[0].setZero();
    model.a[0].setZero();
    model.a_bias[0].setZero();

    if (update_kinematics)
      UpdateKinematicsCustom(model, &Q, &QDot, NULL);

    LOG << std::endl;

    unsigned int reference_body_id = body_id;
    Vector3d reference_point = point_position;

    if (model.IsFixedBodyId(body_id)) {
      unsigned int fbody_id = body_id - model.fixed_body_discriminator;
      reference_body_id = model.mFixedBodies[fbody_id].mMovableParent;
      Vector3d base_coords = CalcBodyToBaseCoordinates (model, Q, body_id, point_position, false);
      reference_point = CalcBaseToBodyCoordinates (model, Q, reference_body_id, base_coords, false);
    }

    Matrix3d global_body_orientation_inv = CalcBodyWorldOrientation (model, Q, reference_body_id, false).transpose();

    SpatialTransform p_X_i (global_body_orientation_inv, reference_point);

    SpatialVector p_a_i = p_X_i.apply(model.a_bias[reference_body_id]);

    return p_a_i;
  }

  bool InverseKinematics (
      Model &model,
      const VectorNd &Qinit,
      const std::vector<unsigned int>& body_id,
      const std::vector<Vector3d>& body_point,
      const std::vector<Vector3d>& target_pos,
      VectorNd &Qres,
      double step_tol,
      double lambda,
      unsigned int max_iter
      ) {

    assert (Qinit.size() == model.dof_count);
    assert (body_id.size() == body_point.size());
    assert (body_id.size() == target_pos.size());

    MatrixNd J = MatrixNd::Zero(3 * body_id.size(), model.dof_count);
    VectorNd e = VectorNd::Zero(3 * body_id.size());

    Qres = Qinit;

    for (unsigned int ik_iter = 0; ik_iter < max_iter; ik_iter++) {
      UpdateKinematicsCustom (model, &Qres, NULL, NULL);

      for (unsigned int k = 0; k < body_id.size(); k++) {
        MatrixNd G (3, model.dof_count);
        CalcPointJacobian (model, Qres, body_id[k], body_point[k], G, false);
        Vector3d point_base = CalcBodyToBaseCoordinates (model, Qres, body_id[k], body_point[k], false);
        LOG << "current_pos = " << point_base.transpose() << std::endl;

        for (unsigned int i = 0; i < 3; i++) {
          for (unsigned int j = 0; j < model.dof_count; j++) {
            unsigned int row = k * 3 + i;
            LOG << "i = " << i << " j = " << j << " k = " << k << " row = " << row << " col = " << j << std::endl;
            J(row, j) = G (i,j);
          }

          e[k * 3 + i] = target_pos[k][i] - point_base[i];
        }

        LOG << J << std::endl;

        // abort if we are getting "close"
        if (e.norm() < step_tol) {
          LOG << "Reached target close enough after " << ik_iter << " steps" << std::endl;
          return true;
        }
      }

      LOG << "J = " << J << std::endl;
      LOG << "e = " << e.transpose() << std::endl;

      MatrixNd JJTe_lambda2_I = J * J.transpose() + lambda*lambda * MatrixNd::Identity(e.size(), e.size());

      VectorNd z (body_id.size() * 3);
#ifndef RBDL_USE_SIMPLE_MATH
      z = JJTe_lambda2_I.colPivHouseholderQr().solve (e);
#else
      bool solve_successful = LinSolveGaussElimPivot (JJTe_lambda2_I, e, z);
      assert (solve_successful);
#endif

      LOG << "z = " << z << std::endl;

      VectorNd delta_theta = J.transpose() * z;
      LOG << "change = " << delta_theta << std::endl;

      Qres = Qres + delta_theta;
      LOG << "Qres = " << Qres.transpose() << std::endl;

      if (delta_theta.norm() < step_tol) {
        LOG << "reached convergence after " << ik_iter << " steps" << std::endl;
        return true;
      }

      VectorNd test_1 (z.size());
      VectorNd test_res (z.size());

      test_1.setZero();

      for (unsigned int i = 0; i < z.size(); i++) {
        test_1[i] = 1.;

        VectorNd test_delta = J.transpose() * test_1;

        test_res[i] = test_delta.squaredNorm();

        test_1[i] = 0.;
      }

      LOG << "test_res = " << test_res.transpose() << std::endl;
    }

    return false;
  }

}
