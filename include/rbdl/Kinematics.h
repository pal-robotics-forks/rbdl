/*
 * RBDL - Rigid Body Dynamics Library
 * Copyright (c) 2011-2016 Martin Felis <martin.felis@iwr.uni-heidelberg.de>
 *
 * Licensed under the zlib license. See LICENSE for more details.
 */

#ifndef RBDL_KINEMATICS_H
#define RBDL_KINEMATICS_H

#include "rbdl/Logging.h"
#include "rbdl/Model.h"
#include "rbdl/rbdl_math.h"
#include <assert.h>
#include <iostream>
#include <rbdl/Joint.h>

namespace RigidBodyDynamics
{
using namespace Math;
/** \page kinematics_page Kinematics
 * All functions related to kinematics are specified in the \ref
 * kinematics_group "Kinematics Module".
 *
 * \note Please note that in the Rigid %Body Dynamics Library all angles
 * are specified in radians.
 *
 * \defgroup kinematics_group Kinematics
 * @{
 *
 * \note Please note that in the Rigid %Body Dynamics Library all angles
 * are specified in radians.
 */

/** \brief Updates and computes velocities and accelerations of the bodies
 *
 * This function updates the kinematic variables such as body velocities
 * and accelerations in the model to reflect the variables passed to this
 * function.
 *
 * \param model the model
 * \param Q     the positional variables of the model
 * \param QDot  the generalized velocities of the joints
 * \param QDDot the generalized accelerations of the joints
 */
template <typename T>
RBDL_DLLAPI void UpdateKinematics(const Model &model, ModelData<T> &model_data,
                                  const VectorN<T> &Q, const VectorN<T> &QDot,
                                  const VectorN<T> &QDDot)
{
  //  LOG << "-------- " << __func__ << " --------" << std::endl;

  assert(Q.rows() == model.q_size);
  assert(QDot.rows() == model.qdot_size);
  assert(QDDot.rows() == model.qdot_size);

  unsigned int i;

  SpatialVectord spatial_gravity(0., 0., 0., model.gravity[0], model.gravity[1],
                                 model.gravity[2]);

  model_data.a[0].setZero();
  // model_data.a[0] = spatial_gravity;

  for (i = 1; i < model.mBodies.size(); i++)
  {
    unsigned int q_index = model.mJoints[i].q_index;

    Joint joint = model.mJoints[i];
    unsigned int lambda = model.lambda[i];

    jcalc(model, model_data, i, Q, QDot);

    model_data.X_lambda[i] = model_data.X_J[i] * model.X_T[i];

    if (lambda != 0)
    {
      model_data.X_base[i] = model_data.X_lambda[i] * model_data.X_base[lambda];
      model_data.v[i] = model_data.X_lambda[i].apply(model_data.v[lambda]) + model_data.v_J[i];
      model_data.a_bias[i] =
          model_data.X_lambda[i].apply(model_data.a_bias[lambda]) + model_data.c[i];
    }
    else
    {
      model_data.X_base[i] = model_data.X_lambda[i];
      model_data.v[i] = model_data.v_J[i];
      model_data.a_bias[i].setZero();
    }

    model_data.c[i] = model_data.c_J[i] + crossm(model_data.v[i], model_data.v_J[i]);
    model_data.a[i] = model_data.X_lambda[i].apply(model_data.a[lambda]) + model_data.c[i];

    // if(model.mJoints[i].mJointType != JointTypeCustom){
    if (model.mJoints[i].mDoFCount == 1)
    {
      model_data.a[i] = model_data.a[i] + model_data.S[i] * QDDot[q_index];
    }
    else if (model.mJoints[i].mDoFCount == 3)
    {
      Vector3d omegadot_temp(QDDot[q_index], QDDot[q_index + 1], QDDot[q_index + 2]);
      model_data.a[i] = model_data.a[i] + model_data.multdof3_S[i] * omegadot_temp;
    }
    //    } else {
    //      unsigned int custom_index = model.mJoints[i].custom_joint_index;
    //      const CustomJoint* custom_joint = model.mCustomJoints[custom_index];
    //      unsigned int joint_dof_count = custom_joint->mDoFCount;

    //      model_data.a[i] = model_data.a[i]
    //        + ( model.mCustomJoints[custom_index]->S
    //            * QDDot.block(q_index, 0, joint_dof_count, 1));
    //    }
  }

  //  for (i = 1; i < model.mBodies.size(); i++) {
  //    //LOG << "a[" << i << "] = " << model_data.a[i].transpose() <<
  //    std::endl;
  //  }
}

void UpdateKinematics(Model &model, const VectorNd &Q, const VectorNd &QDot,
                      const VectorNd &QDDot);

/** \brief Selectively updates model internal states of body positions,
 velocities and/or accelerations.
 *
 * This function updates the kinematic variables such as body velocities and
 * accelerations in the model to reflect the variables passed to this function.
 *
 * In contrast to UpdateKinematics() this function allows to update the model
 * state with values one is interested and thus reduce computations (e.g. only
 * positions, only positions + accelerations, only velocities, etc.).

 * \param model the model
 * \param Q     the positional variables of the model
 * \param QDot  the generalized velocities of the joints
 * \param QDDot the generalized accelerations of the joints
 */
template <typename T>
void UpdateKinematicsCustom(const Model &model, ModelData<T> &model_data,
                            const Math::VectorN<T> *Q, const Math::VectorN<T> *QDot,
                            const Math::VectorN<T> *QDDot)
{
  unsigned int i;

  if (Q)
  {
    assert(Q->rows() == model.q_size);
    for (i = 1; i < model.mBodies.size(); i++)
    {
      unsigned int lambda = model.lambda[i];

      VectorN<T> QDot_zero(VectorN<T>::Zero(model.q_size));

      jcalc<T>(model, model_data, i, (*Q), QDot_zero);

      model_data.X_lambda[i] = model_data.X_J[i] * model.X_T[i];

      if (lambda != 0)
      {
        model_data.X_base[i] = model_data.X_lambda[i] * model_data.X_base[lambda];
      }
      else
      {
        model_data.X_base[i] = model_data.X_lambda[i];
      }
    }
  }

  if (QDot)
  {
    assert(Q->rows() == model.q_size);
    assert(QDot->rows() == model.qdot_size);

    for (i = 1; i < model.mBodies.size(); i++)
    {
      unsigned int lambda = model.lambda[i];

      jcalc(model, model_data, i, *Q, *QDot);

      if (lambda != 0)
      {
        model_data.v[i] =
            model_data.X_lambda[i].apply(model_data.v[lambda]) + model_data.v_J[i];
        model_data.c[i] = model_data.c_J[i] + crossm(model_data.v[i], model_data.v_J[i]);
        model_data.a_bias[i] =
            model_data.X_lambda[i].apply(model_data.a_bias[lambda]) + model_data.c[i];
      }
      else
      {
        model_data.v[i] = model_data.v_J[i];
        model_data.c[i] = model_data.c_J[i] + crossm(model_data.v[i], model_data.v_J[i]);
        model_data.a_bias[i] =
            model_data.X_lambda[i].apply(model_data.a_bias[lambda]) + model_data.c[i];
      }
      // LOG << "v[" << i << "] = " << model_data.v[i].transpose() << std::endl;
    }
  }

  if (QDDot)
  {
    assert(Q->rows() == model.q_size);
    assert(QDot->rows() == model.qdot_size);
    assert(QDDot->rows() == model.qdot_size);

    for (i = 1; i < model.mBodies.size(); i++)
    {
      unsigned int q_index = model.mJoints[i].q_index;

      unsigned int lambda = model.lambda[i];

      if (lambda != 0)
      {
        model_data.a[i] = model_data.X_lambda[i].apply(model_data.a[lambda]) + model_data.c[i];
      }
      else
      {
        model_data.a[i] = model_data.c[i];
      }

      //      if( model.mJoints[i].mJointType != JointTypeCustom){
      if (model.mJoints[i].mDoFCount == 1)
      {
        model_data.a[i] = model_data.a[i] + model_data.S[i] * (*QDDot)[q_index];
      }
      else if (model.mJoints[i].mDoFCount == 3)
      {
        Vector3<T> omegadot_temp((*QDDot)[q_index], (*QDDot)[q_index + 1],
                                 (*QDDot)[q_index + 2]);
        model_data.a[i] = model_data.a[i] + model_data.multdof3_S[i] * omegadot_temp;
      }
      //      } else {
      //        unsigned int k = model.mJoints[i].custom_joint_index;

      //        const CustomJoint* custom_joint = model.mCustomJoints[k];
      //        unsigned int joint_dof_count = custom_joint->mDoFCount;

      //        model_data.a[i] = model_data.a[i]
      //                          + (  (model.mCustomJoints[k]->S). template
      //                          cast<T>()
      //                               *(QDDot->block(q_index, 0,
      //                               joint_dof_count, 1)));
      //      }
    }
  }
}

void UpdateKinematicsCustom(Model &model, const Math::VectorNd *Q,
                            const Math::VectorNd *QDot, const Math::VectorNd *QDDot);

/** \brief Returns the orientation of a given body as 3x3 matrix
 *
 * \param model the rigid body model
 * \param Q the curent genereralized positions
 * \param body_id id of the body for which the point coordinates are expressed
 * \param update_kinematics whether UpdateKinematics() should be called or not
 * (default: true).
 *
 * \returns An orthonormal 3x3 matrix that rotates vectors from base coordinates
 * to body coordinates.
 */

template <typename T>
RBDL_DLLAPI Math::Matrix3<T> CalcBodyWorldOrientation(const Model &model, ModelData<T> &model_data,
                                                      const Math::VectorN<T> &Q,
                                                      const unsigned int body_id,
                                                      bool update_kinematics = true)
{
  return CalcBodyWorldOrientation<T>(model, model_data, Q, body_id,
                                     Eigen::Matrix<T, 3, 3>::Identity(), update_kinematics);
}

template <typename T>
RBDL_DLLAPI Math::Matrix3<T> CalcBodyWorldOrientation(const Model &model, ModelData<T> &model_data,
                                                      const Math::VectorN<T> &Q,
                                                      const unsigned int body_id,
                                                      const Eigen::Matrix<T, 3, 3> &rot,
                                                      bool update_kinematics = true)
{
  // update the Kinematics if necessary
  if (update_kinematics)
  {
    UpdateKinematicsCustom<T>(model, model_data, &Q, NULL, NULL);
  }

  if (body_id >= model.fixed_body_discriminator)
  {
    unsigned int fbody_id = body_id - model.fixed_body_discriminator;
    return rot.transpose() *
           (model.mFixedBodies[fbody_id].mParentTransform.cast<T>() *
            model_data.X_base[model.mFixedBodies[fbody_id].mMovableParent])
               .E;
  }
  return rot.transpose() * model_data.X_base[body_id].E;
}

Math::Matrix3d CalcBodyWorldOrientation(Model &model, const Math::VectorNd &Q,
                                        const unsigned int body_id,
                                        bool update_kinematics = true);

Math::Matrix3d CalcBodyWorldOrientation(Model &model, const Math::VectorNd &Q,
                                        const unsigned int body_id, const Math::Matrix3d &rot,
                                        bool update_kinematics = true);


/** \brief Returns the base coordinates of a point given in body coordinates.
 *
 * \param model the rigid body model
 * \param Q the curent genereralized positions
 * \param body_id id of the body for which the point coordinates are expressed
 * \param body_point_position coordinates of the point in body coordinates
 * \param update_kinematics whether UpdateKinematics() should be called
 * or not (default: true)
 *
 * \returns a 3-D vector with coordinates of the point in base coordinates
 */
Vector3d CalcBodyToBaseCoordinates(Model &model, const VectorNd &Q, unsigned int body_id,
                                   const Vector3d &point_body_coordinates,
                                   bool update_kinematics = true);

template <typename T>
Vector3<T> CalcBodyToBaseCoordinates(const Model &model, ModelData<T> &model_data,
                                     const VectorN<T> &Q, unsigned int body_id,
                                     const Vector3_t<T> &point_body_coordinates,
                                     bool update_kinematics = true)
{
  // update the Kinematics if necessary
  if (update_kinematics)
  {
    UpdateKinematicsCustom<T>(model, model_data, &Q, NULL, NULL);
  }

  if (body_id >= model.fixed_body_discriminator)
  {
    unsigned int fbody_id = body_id - model.fixed_body_discriminator;
    unsigned int parent_id = model.mFixedBodies[fbody_id].mMovableParent;

    Matrix3<T> fixed_rotation =
        model.mFixedBodies[fbody_id].mParentTransform.E.transpose().cast<T>();
    Vector3<T> fixed_position = model.mFixedBodies[fbody_id].mParentTransform.r.cast<T>();

    Matrix3<T> parent_body_rotation = model_data.X_base[parent_id].E.transpose();
    Vector3<T> parent_body_position = model_data.X_base[parent_id].r;

    return (parent_body_position +
            (parent_body_rotation * (fixed_position + fixed_rotation * point_body_coordinates)));
  }

  Matrix3<T> body_rotation = model_data.X_base[body_id].E.transpose();
  Vector3<T> body_position = model_data.X_base[body_id].r;

  return body_position + body_rotation * point_body_coordinates;
}

/** \brief Returns the body coordinates of a point given in base coordinates.
 *
 * \param model the rigid body model
 * \param Q the curent genereralized positions
 * \param body_id id of the body for which the point coordinates are expressed
 * \param base_point_position coordinates of the point in base coordinates
 * \param update_kinematics whether UpdateKinematics() should be called or not
 * (default: true).
 *
 * \returns a 3-D vector with coordinates of the point in body coordinates
 */
template <typename T>
RBDL_DLLAPI Math::Vector3<T> CalcBaseToBodyCoordinates(
    const Model &model, ModelData<T> &model_data, const Math::VectorN<T> &Q, unsigned int body_id,
    const Math::Vector3<T> &base_point_position, bool update_kinematics = true)
{
  if (update_kinematics)
  {
    UpdateKinematicsCustom<T>(model, model_data, &Q, NULL, NULL);
  }

  if (body_id >= model.fixed_body_discriminator)
  {
    unsigned int fbody_id = body_id - model.fixed_body_discriminator;
    unsigned int parent_id = model.mFixedBodies[fbody_id].mMovableParent;

    Matrix3<T> fixed_rotation = model.mFixedBodies[fbody_id].mParentTransform.E;
    Vector3<T> fixed_position = model.mFixedBodies[fbody_id].mParentTransform.r;

    Matrix3<T> parent_body_rotation = model_data.X_base[parent_id].E;
    Vector3<T> parent_body_position = model_data.X_base[parent_id].r;

    return (fixed_rotation *
            (-fixed_position -
             parent_body_rotation * (parent_body_position - base_point_position)));
  }

  Matrix3<T> body_rotation = model_data.X_base[body_id].E;
  Vector3<T> body_position = model_data.X_base[body_id].r;

  Vector3<T> position = body_rotation * (base_point_position - body_position);
  return position;
}

Math::Vector3d CalcBaseToBodyCoordinates(Model &model, const Math::VectorNd &Q,
                                         unsigned int body_id,
                                         const Math::Vector3d &base_point_position,
                                         bool update_kinematics = true);


/** \brief Computes the point jacobian for a point on a body
 *
 * If a position of a point is computed by a function \f$g(q(t))\f$ for which
 * its
 * time derivative is \f$\frac{d}{dt} g(q(t)) = G(q)\dot{q}\f$ then this
 * function computes the jacobian matrix \f$G(q)\f$.
 *
 * \param model   rigid body model
 * \param Q       state vector of the internal joints
 * \param body_id the id of the body
 * \param point_position the position of the point in body-local data
 * \param G       a matrix of dimensions 3 x \#qdot_size where the result will
 * be stored in
 * \param update_kinematics whether UpdateKinematics() should be called or not
 * (default: true)
 *
 * The result will be returned via the G argument.
 *
 * \note This function only evaluates the entries of G that are non-zero. One
 * Before calling this function one has to ensure that all other values
 * have been set to zero, e.g. by calling G.setZero().
 *
 */
RBDL_DLLAPI void CalcPointJacobian(const Model &model, ModelDatad &model_data,
                                   const Math::VectorNd &Q, unsigned int body_id,
                                   const Math::Vector3d &point_position,
                                   Math::MatrixNd &G, bool update_kinematics = true);

void CalcPointJacobian(Model &model, const Math::VectorNd &Q, unsigned int body_id,
                       const Math::Vector3d &point_position, Math::MatrixNd &G,
                       bool update_kinematics = true);

RBDL_DLLAPI void CalcOrientationJacobian(const Model &model, ModelDatad &model_data,
                                         const Math::VectorNd &Q, unsigned int body_id,
                                         Math::MatrixNd &G, bool update_kinematics = true);

void CalcOrientationJacobian(Model &model, const Math::VectorNd &Q, unsigned int body_id,
                             Math::MatrixNd &G, bool update_kinematics = true);

/** \brief Computes a 6-D Jacobian for a point on a body
 *
 * Computes the 6-D Jacobian \f$G(q)\f$ that when multiplied with
 * \f$\dot{q}\f$ gives a 6-D vector that has the angular velocity as the
 * first three entries and the linear velocity as the last three entries.
 *
 * \param model   rigid body model
 * \param Q       state vector of the internal joints
 * \param body_id the id of the body
 * \param point_position the position of the point in body-local data
 * \param G       a matrix of dimensions 6 x \#qdot_size where the result will
 * be stored in
 * \param update_kinematics whether UpdateKinematics() should be called or not
 * (default: true)
 *
 * The result will be returned via the G argument.
 *
 * \note This function only evaluates the entries of G that are non-zero. One
 * Before calling this function one has to ensure that all other values
 * have been set to zero, e.g. by calling G.setZero().
 *
 */

RBDL_DLLAPI void CalcPointJacobian6D(const Model &model, ModelDatad &model_data,
                                     const Math::VectorNd &Q, unsigned int body_id,
                                     const Math::Vector3d &point_position,
                                     Math::MatrixNd &G, bool update_kinematics = true);

RBDL_DLLAPI void CalcPointJacobian6D(Model &model, const Math::VectorNd &Q,
                                     unsigned int body_id, const Math::Vector3d &point_position,
                                     Math::MatrixNd &G, bool update_kinematics = true);

RBDL_DLLAPI void CalcPointJacobian6DRelative(const Model &model, ModelDatad &model_data,
                                             const Math::VectorNd &Q, unsigned int body_id,
                                             unsigned int respect_body_id,
                                             const Math::Vector3d &point_position,
                                             Math::MatrixNd &G, bool update_kinematics = true);

RBDL_DLLAPI void CalcPointJacobian6DRelative(Model &model, const Math::VectorNd &Q,
                                             unsigned int body_id, unsigned int respect_body_id,
                                             const Math::Vector3d &point_position,
                                             Math::MatrixNd &G, bool update_kinematics = true);


RBDL_DLLAPI void CalcPointJacobian6D(Model &model, const Math::VectorNd &Q,
                                     unsigned int body_id, const Math::Vector3d &point_position,
                                     unsigned int respect_body_id, Math::MatrixNd &G,
                                     bool update_kinematics = true);

RBDL_DLLAPI void CalcPointJacobian6D(const Model &model, ModelDatad &model_data,
                                     const Math::VectorNd &Q, unsigned int body_id,
                                     const Math::Isometry3d &pose, Math::MatrixNd &G,
                                     bool update_kinematics = true);

RBDL_DLLAPI void CalcPointJacobian6D(Model &model, const Math::VectorNd &Q,
                                     unsigned int body_id, const Math::Isometry3d &pose,
                                     Math::MatrixNd &G, bool update_kinematics = true);

RBDL_DLLAPI void CalcPointJacobian6DBodyFrame(const Model &model, ModelDatad &model_data,
                                              const Math::VectorNd &Q, unsigned int body_id,
                                              const Math::Vector3d &point_position,
                                              Math::MatrixNd &G, bool update_kinematics = true);

RBDL_DLLAPI void CalcPointJacobian6DBodyFrame(Model &model, const Math::VectorNd &Q,
                                              unsigned int body_id,
                                              const Math::Vector3d &point_position,
                                              Math::MatrixNd &G, bool update_kinematics = true);

/** \brief Computes the spatial jacobian for a body
 *
 * The spatial velocity of a body at the origin of the base coordinate
 * system can be expressed as \f${}^0 \hat{v}_i = G(q) * \dot{q}\f$. The
 * matrix \f$G(q)\f$ is called the spatial body jacobian of the body and
 * can be computed using this function.
 *
 * \param model   rigid body model
 * \param Q       state vector of the internal joints
 * \param body_id the id of the body
 * \param G       a matrix of size 6 x \#qdot_size where the result will be
 * stored in
 * \param update_kinematics whether UpdateKinematics() should be called or not
 * (default: true)
 *
 * The result will be returned via the G argument and represents the
 * body Jacobian expressed at the origin of the body.
 *
 * \note This function only evaluates the entries of G that are non-zero. One
 * Before calling this function one has to ensure that all other values
 * have been set to zero, e.g. by calling G.setZero().
 */
RBDL_DLLAPI void CalcBodySpatialJacobian(const Model &model, ModelDatad &model_data,
                                         const Math::VectorNd &Q, unsigned int body_id,
                                         Math::MatrixNd &G, bool update_kinematics = true);

/** \brief Computes the velocity of a point on a body
 *
 * \param model   rigid body model
 * \param Q       state vector of the internal joints
 * \param QDot    velocity vector of the internal joints
 * \param body_id the id of the body
 * \param point_position the position of the point in body-local data
 * \param update_kinematics whether UpdateKinematics() should be called or not
 * (default: true)
 *
 * \returns The cartesian velocity of the point in global frame (output)
 */

RBDL_DLLAPI Math::Vector3d CalcPointVelocity(const Model &model, ModelDatad &model_data,
                                             const Math::VectorNd &Q,
                                             const Math::VectorNd &QDot, unsigned int body_id,
                                             const Math::Vector3d &point_position,
                                             bool update_kinematics = true);

/// @todo point_position is not used, it has to be deleted from params
RBDL_DLLAPI Math::Vector3d CalcPointVelocity(Model &model, const Math::VectorNd &Q,
                                             const Math::VectorNd &QDot, unsigned int body_id,
                                             const Math::Vector3d &point_position,
                                             bool update_kinematics = true);


RBDL_DLLAPI Math::Vector3d CalcPointVelocity(const Model &model, ModelDatad &model_data,
                                             const Math::VectorNd &Q, const Math::VectorNd &QDot,
                                             unsigned int body_id, const Math::Isometry3d &pose,
                                             bool update_kinematics = true);

RBDL_DLLAPI Math::Vector3d CalcPointVelocity(Model &model, const Math::VectorNd &Q,
                                             const Math::VectorNd &QDot, unsigned int body_id,
                                             const Math::Isometry3d &pose,
                                             bool update_kinematics = true);

RBDL_DLLAPI Math::Vector3d CalcPointAngularVelocity(const Model &model, ModelDatad &model_data,
                                                    const Math::VectorNd &Q,
                                                    const Math::VectorNd &QDot,
                                                    unsigned int body_id,
                                                    const Math::Vector3d &point_position,
                                                    bool update_kinematics = true);

RBDL_DLLAPI Math::Vector3d CalcPointAngularVelocity(Model &model, const Math::VectorNd &Q,
                                                    const Math::VectorNd &QDot,
                                                    unsigned int body_id,
                                                    const Math::Vector3d &point_position,
                                                    bool update_kinematics = true);

RBDL_DLLAPI Math::Vector3d CalcPointAngularVelocityRelative(
    const Model &model, ModelDatad &model_data, const Math::VectorNd &Q,
    const Math::VectorNd &QDot, unsigned int body_id, unsigned int respect_body_id,
    const Math::Vector3d &point_position, bool update_kinematics = true);

RBDL_DLLAPI Math::Vector3d CalcPointAngularVelocityRelative(
    Model &model, const Math::VectorNd &Q, const Math::VectorNd &QDot,
    unsigned int body_id, unsigned int respect_body_id,
    const Math::Vector3d &point_position, bool update_kinematics = true);

RBDL_DLLAPI Math::Vector3d CalcPointAngularVelocity(const Model &model, ModelDatad &model_data,
                                                    const Math::VectorNd &Q,
                                                    const Math::VectorNd &QDot,
                                                    unsigned int body_id,
                                                    const Math::Isometry3d &pose,
                                                    bool update_kinematics = true);

RBDL_DLLAPI Math::Vector3d CalcPointAngularVelocity(Model &model, const Math::VectorNd &Q,
                                                    const Math::VectorNd &QDot,
                                                    unsigned int body_id,
                                                    const Math::Isometry3d &pose,
                                                    bool update_kinematics = true);

/** \brief Computes angular and linear velocity of a point that is fixed on a body
 *
 * \param model   rigid body model
 * \param Q       state vector of the internal joints
 * \param QDot    velocity vector of the internal joints
 * \param body_id the id of the body
 * \param point_position the position of the point in body-local data
 * \param update_kinematics whether UpdateKinematics() should be called or not
 * (default: true)
 *
 * \returns The a 6-D vector for which the first three elements are the
 * angular velocity and the last three elements the linear velocity in the
 * global reference system.
 */
RBDL_DLLAPI
Math::SpatialVectord CalcPointVelocity6D(const Model &model, ModelDatad &model_data,
                                         const Math::VectorNd &Q,
                                         const Math::VectorNd &QDot, unsigned int body_id,
                                         const Math::Vector3d &point_position,
                                         bool update_kinematics = true);

Math::SpatialVectord CalcPointVelocity6D(Model &model, const Math::VectorNd &Q,
                                         const Math::VectorNd &QDot, unsigned int body_id,
                                         const Math::Vector3d &point_position,
                                         bool update_kinematics = true);

/** \brief Computes the linear acceleration of a point on a body
 *
 * \param model   rigid body model
 * \param Q       state vector of the internal joints
 * \param QDot    velocity vector of the internal joints
 * \param QDDot    velocity vector of the internal joints
 * \param body_id the id of the body
 * \param point_position the position of the point in body-local data
 * \param update_kinematics whether UpdateKinematics() should be called or not
 * (default: true)
 *
 * \returns The cartesian acceleration of the point in global frame (output)
 *
 * The kinematic state of the model has to be updated before valid
 * values can be obtained. This can either be done by calling
 * UpdateKinematics() or setting the last parameter update_kinematics to
 * true (default).
 *
 * \warning  If this function is called after ForwardDynamics() without
 * an update of the kinematic state one has to add the gravity
 * acceleration has to be added to the result.
 */
RBDL_DLLAPI
Math::Vector3d CalcPointAcceleration(const Model &model, ModelDatad &model_data,
                                     const Math::VectorNd &Q, const Math::VectorNd &QDot,
                                     const Math::VectorNd &QDDot, unsigned int body_id,
                                     const Math::Vector3d &point_position,
                                     bool update_kinematics = true);

Math::Vector3d CalcPointAcceleration(Model &model, const Math::VectorNd &Q,
                                     const Math::VectorNd &QDot, const Math::VectorNd &QDDot,
                                     unsigned int body_id, const Math::Vector3d &point_position,
                                     bool update_kinematics = true);

Vector3d CalcPointAngularAcceleration(const Model &model, ModelDatad &model_data,
                                      const VectorNd &Q, const VectorNd &QDot,
                                      const VectorNd &QDDot, unsigned int body_id,
                                      const Vector3d &point_position,
                                      bool update_kinematics = true);

Vector3d CalcPointAngularAcceleration(Model &model, const VectorNd &Q,
                                      const VectorNd &QDot, const VectorNd &QDDot,
                                      unsigned int body_id, const Vector3d &point_position,
                                      bool update_kinematics = true);

RBDL_DLLAPI
Math::Vector3d CalcPointAccelerationBias(const Model &model, ModelDatad &model_data,
                                         const Math::VectorNd &Q,
                                         const Math::VectorNd &QDot, unsigned int body_id,
                                         const Math::Vector3d &point_position,
                                         bool update_kinematics = true);

Math::Vector3d CalcPointAccelerationBias(Model &model, const Math::VectorNd &Q,
                                         const Math::VectorNd &QDot, unsigned int body_id,
                                         const Math::Vector3d &point_position,
                                         bool update_kinematics = true);

/** \brief Computes linear and angular acceleration of a point on a body
 *
 * \param model   rigid body model
 * \param Q       state vector of the internal joints
 * \param QDot    velocity vector of the internal joints
 * \param QDDot    velocity vector of the internal joints
 * \param body_id the id of the body
 * \param point_position the position of the point in body-local data
 * \param update_kinematics whether UpdateKinematics() should be called or not
 * (default: true)
 *
 * \returns A 6-D vector where the first three elements are the angular
 * acceleration and the last three elements the linear accelerations of
 * the point.
 *
 * The kinematic state of the model has to be updated before valid
 * values can be obtained. This can either be done by calling
 * UpdateKinematics() or setting the last parameter update_kinematics to
 * true (default).
 *
 * \warning  If this function is called after ForwardDynamics() without
 * an update of the kinematic state one has to add the gravity
 * acceleration has to be added to the result.
 */
RBDL_DLLAPI
Math::SpatialVectord CalcPointAcceleration6D(const Model &model, ModelDatad &model_data,
                                             const Math::VectorNd &Q, const Math::VectorNd &QDot,
                                             const Math::VectorNd &QDDot, unsigned int body_id,
                                             const Math::Vector3d &point_position,
                                             bool update_kinematics = true);

Math::SpatialVectord CalcPointAcceleration6D(Model &model, const Math::VectorNd &Q,
                                             const Math::VectorNd &QDot,
                                             const Math::VectorNd &QDDot, unsigned int body_id,
                                             const Math::Vector3d &point_position,
                                             bool update_kinematics = true);

RBDL_DLLAPI
Math::SpatialVectord CalcPointAcceleration6DBias(const Model &model, ModelDatad &model_data,
                                                 const Math::VectorNd &Q,
                                                 const Math::VectorNd &QDot, unsigned int body_id,
                                                 const Math::Vector3d &point_position,
                                                 bool update_kinematics = true);

Math::SpatialVectord CalcPointAcceleration6DBias(Model &model, const Math::VectorNd &Q,
                                                 const Math::VectorNd &QDot, unsigned int body_id,
                                                 const Math::Vector3d &point_position,
                                                 bool update_kinematics = true);
RBDL_DLLAPI
Math::SpatialVectord CalcPointAcceleration6DBias(const Model &model, ModelDatad &model_data,
                                                 const Math::VectorNd &Q,
                                                 const Math::VectorNd &QDot, unsigned int body_id,
                                                 const Math::Isometry3d &pose,
                                                 bool update_kinematics = true);

Math::SpatialVectord CalcPointAcceleration6DBias(Model &model, const Math::VectorNd &Q,
                                                 const Math::VectorNd &QDot, unsigned int body_id,
                                                 const Math::Isometry3d &pose,
                                                 bool update_kinematics = true);

/** \brief Computes the inverse kinematics iteratively using a damped Levenberg-Marquardt
 * method (also known as Damped Least Squares method)
 *
 * \param model rigid body model
 * \param Qinit initial guess for the state
 * \param body_id a vector of all bodies for which we we have kinematic target positions
 * \param body_point a vector of points in body local coordinates that are
 * to be matched to target positions
 * \param target_pos a vector of target positions
 * \param Qres output of the computed inverse kinematics
 * \param step_tol tolerance used for convergence detection
 * \param lambda damping factor for the least squares function
 * \param max_iter maximum number of steps that should be performed
 * \returns true on success, false otherwise
 *
 * This function repeatedly computes
 *   \f[ Qres = Qres + \Delta \theta\f]
 *   \f[ \Delta \theta = G^T (G^T G + \lambda^2 I)^{-1} e \f]
 * where \f$G = G(q) = \frac{d}{dt} g(q(t))\f$ and \f$e\f$ is the
 * correction of the body points so that they coincide with the target
 * positions. The function returns true when \f$||\Delta \theta||_2 \le\f$
 * step_tol or if the error between body points and target gets smaller
 * than step_tol. Otherwise it returns false.
 *
 * The parameter \f$\lambda\f$ is the damping factor that has to
 * be chosen carefully. In case of unreachable positions higher values (e.g
 * 0.9) can be helpful. Otherwise values of 0.0001, 0.001, 0.01, 0.1 might
 * yield good results. See the literature for best practices.
 *
 * \warning The actual accuracy might be rather low (~1.0e-2)! Use this function with a
 * grain of suspicion.
 */
RBDL_DLLAPI
bool InverseKinematics(Model &model, ModelDatad &model_data, const Math::VectorNd &Qinit,
                       const std::vector<unsigned int> &body_id,
                       const std::vector<Math::Vector3d> &body_point,
                       const std::vector<Math::Vector3d> &target_pos,
                       Math::VectorNd &Qres, double step_tol = 1.0e-12,
                       double lambda = 0.01, unsigned int max_iter = 50);

RBDL_DLLAPI Math::Vector3d CalcAngularVelocityfromMatrix(const Math::Matrix3d &RotMat);

struct RBDL_DLLAPI InverseKinematicsConstraintSet
{
  enum ConstraintType
  {
    ConstraintTypePosition = 0,
    ConstraintTypeOrientation,
    ConstraintTypeFull
  };

  InverseKinematicsConstraintSet();

  Math::MatrixNd J;  /// the Jacobian of all constraints
  Math::MatrixNd G;  /// temporary storage of a single body Jacobian
  Math::VectorNd e;  /// Vector with all the constraint residuals.

  unsigned int num_constraints;  // size of all constraints
  double lambda;  /// Damping factor, the default value of 1.0e-6 is reasonable for most
                  /// problems
  unsigned int num_steps;  // The number of iterations performed
  unsigned int max_steps;  // Maximum number of steps (default 300), abort if more steps
                           // are performed.
  double step_tol;  // Step tolerance (default = 1.0e-12). If the computed step length is
                    // smaller than this value the algorithm terminates successfully (i.e.
                    // returns true). If error_norm is still larger than constraint_tol
                    // then this usually means that the target is unreachable.
  double constraint_tol;  // Constraint tolerance (default = 1.0e-12). If error_norm is
                          // smaller than this value the algorithm terminates
                          // successfully, i.e. all constraints are satisfied.
  double error_norm;      // Norm of the constraint residual vector.

  // everything to define a IKin constraint
  std::vector<ConstraintType> constraint_type;
  std::vector<unsigned int> body_ids;
  std::vector<Math::Vector3d> body_points;
  std::vector<Math::Vector3d> target_positions;
  std::vector<Math::Matrix3d> target_orientations;
  std::vector<unsigned int> constraint_row_index;

  // Adds a point constraint that tries to get a body point close to a
  // point described in base coordinates.
  unsigned int AddPointConstraint(unsigned int body_id, const Math::Vector3d &body_point,
                                  const Math::Vector3d &target_pos);
  // Adds an orientation constraint that tries to align a body to the
  // orientation specified as a rotation matrix expressed in base
  // coordinates.
  unsigned int AddOrientationConstraint(unsigned int body_id,
                                        const Math::Matrix3d &target_orientation);
  // Adds a constraint on both location and orientation of a body.
  unsigned int AddFullConstraint(unsigned int body_id, const Math::Vector3d &body_point,
                                 const Math::Vector3d &target_pos,
                                 const Math::Matrix3d &target_orientation);
  // Clears all entries of the constraint setting
  unsigned int ClearConstraints();
};

RBDL_DLLAPI bool InverseKinematics(Model &model, const Math::VectorNd &Qinit,
                                   InverseKinematicsConstraintSet &CS, Math::VectorNd &Qres);

/** @} */
}

/* RBDL_KINEMATICS_H */
#endif
