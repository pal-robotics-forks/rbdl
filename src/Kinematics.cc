/*
 * RBDL - Rigid Body Dynamics Library
 * Copyright (c) 2011-2016 Martin Felis <martin.felis@iwr.uni-heidelberg.de>
 *
 * Licensed under the zlib license. See LICENSE for more details.
 */

#include <assert.h>
#include <cstring>
#include <iostream>
#include <limits>

#include "rbdl/Logging.h"
#include "rbdl/rbdl_mathutils.h"

#include "rbdl/Kinematics.h"
#include "rbdl/Model.h"

namespace RigidBodyDynamics
{
using namespace Math;

void UpdateKinematics(Model &model, const VectorNd &Q, const VectorNd &QDot, const VectorNd &QDDot)
{
  UpdateKinematics<double>(model, *model.getModelData(), Q, QDot, QDDot);
}

void UpdateKinematicsCustom(Model &model, const Math::VectorNd *Q,
                            const Math::VectorNd *QDot, const Math::VectorNd *QDDot)
{
  UpdateKinematicsCustom<double>(model, *model.getModelData(), Q, QDot, QDDot);
}

Vector3d CalcBodyToBaseCoordinates(Model &model, const VectorNd &Q, unsigned int body_id,
                                   const Vector3d &point_body_coordinates, bool update_kinematics)
{
  return CalcBodyToBaseCoordinates<double>(model, *model.getModelData(), Q, body_id,
                                           point_body_coordinates, update_kinematics);
}

Math::Matrix3d CalcBodyWorldOrientation(Model &model, const Math::VectorNd &Q,
                                        const unsigned int body_id, bool update_kinematics)
{
  return CalcBodyWorldOrientation<double>(model, *model.getModelData(), Q, body_id,
                                          Eigen::Matrix3d::Identity(), update_kinematics);
}

Math::Matrix3d CalcBodyWorldOrientation(Model &model, const Math::VectorNd &Q,
                                        const unsigned int body_id,
                                        const Math::Matrix3d &rot, bool update_kinematics)
{
  return CalcBodyWorldOrientation<double>(model, *model.getModelData(), Q, body_id, rot,
                                          update_kinematics);
}

Math::Vector3d CalcBaseToBodyCoordinates(Model &model, const VectorNd &Q, unsigned int body_id,
                                         const Math::Vector3d &base_point_position,
                                         bool update_kinematics)
{
  return CalcBaseToBodyCoordinates<double>(model, *model.getModelData(), Q, body_id,
                                           base_point_position, update_kinematics);
}

RBDL_DLLAPI void CalcPointJacobian(const Model &model, ModelDatad &model_data,
                                   const VectorNd &Q, unsigned int body_id,
                                   const Vector3d &point_position, MatrixNd &G,
                                   bool update_kinematics)
{
  LOG << "-------- " << __func__ << " --------" << std::endl;

  G.setZero();
  // update the Kinematics if necessary
  if (update_kinematics)
  {
    UpdateKinematicsCustom<double>(model, model_data, &Q, NULL, NULL);
  }

  SpatialTransformd point_trans = SpatialTransformd(
      Matrix3d::Identity(),
      CalcBodyToBaseCoordinates(model, model_data, Q, body_id, point_position, false));

  assert(G.rows() == 3 && G.cols() == model.qdot_size);

  unsigned int reference_body_id = body_id;

  if (model.IsFixedBodyId(body_id))
  {
    unsigned int fbody_id = body_id - model.fixed_body_discriminator;
    reference_body_id = model.mFixedBodies[fbody_id].mMovableParent;
  }

  unsigned int j = reference_body_id;

  // e[j] is set to 1 if joint j contributes to the jacobian that we are
  // computing. For all other joints the column will be zero.
  while (j != 0)
  {
    unsigned int q_index = model.mJoints[j].q_index;

    // if(model.mJoints[j].mJointType != JointTypeCustom){
    if (model.mJoints[j].mDoFCount == 1)
    {
      G.block(0, q_index, 3, 1) =
          point_trans.apply(model_data.X_base[j].inverse().apply(model_data.S[j])).block(3, 0, 3, 1);
    }
    else if (model.mJoints[j].mDoFCount == 3)
    {
      G.block(0, q_index, 3, 3) =
          ((point_trans * model_data.X_base[j].inverse()).toMatrix() * model_data.multdof3_S[j])
              .block(3, 0, 3, 3);
    }
    //    } else {
    //      unsigned int k = model.mJoints[j].custom_joint_index;

    //      G.block(0, q_index, 3, model.mCustomJoints[k]->mDoFCount) =
    //        ((point_trans
    //          * model_data.X_base[j].inverse()).toMatrix()
    //         * model.mCustomJoints[k]->S).block(
    //           3,0,3,model.mCustomJoints[k]->mDoFCount);
    //    }

    j = model.lambda[j];
  }
}

void CalcPointJacobian(Model &model, const VectorNd &Q, unsigned int body_id,
                       const Vector3d &point_position, MatrixNd &G, bool update_kinematics)
{
  CalcPointJacobian(model, *model.getModelData(), Q, body_id, point_position, G, update_kinematics);
}

RBDL_DLLAPI void CalcOrientationJacobian(const Model &model, ModelDatad &model_data,
                                         const VectorNd &Q, unsigned int body_id,
                                         MatrixNd &G, bool update_kinematics)
{
  LOG << "-------- " << __func__ << " --------" << std::endl;

  G.setZero();
  Math::Vector3d point_position = Math::Vector3d::Zero();

  // update the Kinematics if necessary
  if (update_kinematics)
  {
    UpdateKinematicsCustom<double>(model, model_data, &Q, NULL, NULL);
  }

  SpatialTransformd point_trans = SpatialTransformd(
      Matrix3d::Identity(),
      CalcBodyToBaseCoordinates(model, model_data, Q, body_id, point_position, false));

  assert(G.rows() == 3 && G.cols() == model.qdot_size);

  unsigned int reference_body_id = body_id;

  if (model.IsFixedBodyId(body_id))
  {
    unsigned int fbody_id = body_id - model.fixed_body_discriminator;
    reference_body_id = model.mFixedBodies[fbody_id].mMovableParent;
  }

  unsigned int j = reference_body_id;

  // e[j] is set to 1 if joint j contributes to the jacobian that we are
  // computing. For all other joints the column will be zero.
  while (j != 0)
  {
    unsigned int q_index = model.mJoints[j].q_index;

    //    if(model.mJoints[j].mJointType != JointTypeCustom){
    if (model.mJoints[j].mDoFCount == 1)
    {
      G.block(0, q_index, 3, 1) =
          point_trans.apply(model_data.X_base[j].inverse().apply(model_data.S[j])).block(0, 0, 3, 1);
    }
    else if (model.mJoints[j].mDoFCount == 3)
    {
      G.block(0, q_index, 3, 3) =
          ((point_trans * model_data.X_base[j].inverse()).toMatrix() * model_data.multdof3_S[j])
              .block(0, 0, 3, 3);
    }
    //    } else {
    //      unsigned int k = model.mJoints[j].custom_joint_index;

    //      G.block(0, q_index, 3, model.mCustomJoints[k]->mDoFCount) =
    //        ((point_trans
    //          * model_data.X_base[j].inverse()).toMatrix()
    //         * model.mCustomJoints[k]->S).block(
    //           0,0,3,model.mCustomJoints[k]->mDoFCount);
    //    }

    j = model.lambda[j];
  }
}

void CalcOrientationJacobian(Model &model, const Math::VectorNd &Q, unsigned int body_id,
                             Math::MatrixNd &G, bool update_kinematics)
{
  CalcOrientationJacobian(model, *model.getModelData(), Q, body_id, G, update_kinematics);
}

RBDL_DLLAPI void CalcPointJacobian6D(const Model &model, ModelDatad &model_data,
                                     const Math::VectorNd &Q, unsigned int body_id,
                                     const Math::Isometry3d &pose, Math::MatrixNd &G,
                                     bool update_kinematics)
{
  return CalcPointJacobian6D(model, model_data, Q, body_id, Vector3d(pose.translation()),
                             G, update_kinematics);
}

RBDL_DLLAPI void CalcPointJacobian6D(Model &model, const Math::VectorNd &Q,
                                     unsigned int body_id, const Math::Isometry3d &pose,
                                     Math::MatrixNd &G, bool update_kinematics)
{
  return CalcPointJacobian6D(model, Q, body_id, Vector3d(pose.translation()), G, update_kinematics);
}

RBDL_DLLAPI void CalcPointJacobian6D(const Model &model, ModelDatad &model_data,
                                     const VectorNd &Q, unsigned int body_id,
                                     const Vector3d &point_position, MatrixNd &G,
                                     bool update_kinematics)
{
  LOG << "-------- " << __func__ << " --------" << std::endl;

  G.setZero();

  // update the Kinematics if necessary
  if (update_kinematics)
  {
    UpdateKinematicsCustom<double>(model, model_data, &Q, NULL, NULL);
  }

  SpatialTransformd point_trans = SpatialTransformd(
      Matrix3d::Identity(),
      CalcBodyToBaseCoordinates(model, model_data, Q, body_id, point_position, false));

  assert(G.rows() == 6 && G.cols() == model.qdot_size);

  unsigned int reference_body_id = body_id;

  if (model.IsFixedBodyId(body_id))
  {
    unsigned int fbody_id = body_id - model.fixed_body_discriminator;
    reference_body_id = model.mFixedBodies[fbody_id].mMovableParent;
  }

  unsigned int j = reference_body_id;

  while (j != 0)
  {
    unsigned int q_index = model.mJoints[j].q_index;

    //    if(model.mJoints[j].mJointType != JointTypeCustom){
    if (model.mJoints[j].mDoFCount == 1)
    {
      G.block(0, q_index, 6, 1) =
          point_trans.apply(model_data.X_base[j].inverse().apply(model_data.S[j])).block(0, 0, 6, 1);
    }
    else if (model.mJoints[j].mDoFCount == 3)
    {
      G.block(0, q_index, 6, 3) =
          ((point_trans * model_data.X_base[j].inverse()).toMatrix() * model_data.multdof3_S[j])
              .block(0, 0, 6, 3);
    }

    j = model.lambda[j];
  }
}

void CalcPointJacobian6D(Model &model, const Math::VectorNd &Q, unsigned int body_id,
                         const Math::Vector3d &point_position, Math::MatrixNd &G,
                         bool update_kinematics)
{
  CalcPointJacobian6D(model, *model.getModelData(), Q, body_id, point_position, G,
                      update_kinematics);
}

void CalcSpatialJacobian(const Model &model, ModelDatad &model_data, const VectorNd &Q,
                         unsigned int body_id, MatrixNd &G, bool update_kinematics)
{
  LOG << "-------- " << __func__ << " --------" << std::endl;

  // update the Kinematics if necessary
  if (update_kinematics)
  {
    UpdateKinematicsCustom<double>(model, model_data, &Q, NULL, NULL);
  }

  assert(G.rows() == 6 && G.cols() == model.qdot_size);

  unsigned int reference_body_id = body_id;

  if (model.IsFixedBodyId(body_id))
  {
    unsigned int fbody_id = body_id - model.fixed_body_discriminator;

    reference_body_id = model.mFixedBodies[fbody_id].mMovableParent;
  }

  unsigned int j = reference_body_id;

  while (j != 0)
  {
    unsigned int q_index = model.mJoints[j].q_index;

    // if(model.mJoints[j].mJointType != JointTypeCustom){
    if (model.mJoints[j].mDoFCount == 1)
    {
      G.block(0, q_index, 6, 1) = (model_data.X_base[j].inverse().apply(model_data.S[j]));
    }
    else if (model.mJoints[j].mDoFCount == 3)
    {
      G.block(0, q_index, 6, 3) =
          (model_data.X_base[j].inverse()).toMatrix() * model_data.multdof3_S[j];
    }

    j = model.lambda[j];
  }
}

RBDL_DLLAPI void CalcRelativeBodySpatialJacobian(const Model &model, ModelDatad &model_data,
                                                 const VectorNd &Q, unsigned int body_id,
                                                 unsigned int respect_body_id,
                                                 MatrixNd &G, bool update_kinematics)
{
  LOG << "-------- " << __func__ << " --------" << std::endl;

  // update the Kinematics if necessary
  if (update_kinematics)
  {
    UpdateKinematicsCustom<double>(model, model_data, &Q, NULL, NULL);
  }

  assert(G.rows() == 6 && G.cols() == model.qdot_size);

  unsigned int reference_body_id = respect_body_id;

  SpatialTransformd base_to_body;

  if (model.IsFixedBodyId(body_id))
  {
    unsigned int fbody_id = body_id - model.fixed_body_discriminator;

    reference_body_id = model.mFixedBodies[fbody_id].mMovableParent;

    base_to_body =
        model.mFixedBodies[fbody_id].mParentTransform * model_data.X_base[respect_body_id];
  }
  else
  {
    base_to_body = model_data.X_base[respect_body_id];
  }

  unsigned int j = body_id;

  while (j != 0)
  {
    unsigned int q_index = model.mJoints[j].q_index;

    if (model.mJoints[j].mDoFCount == 1)
    {
      G.block(0, q_index, 6, 1) =
          base_to_body.apply(model_data.X_base[j].inverse().apply(model_data.S[j]));
    }
    else if (model.mJoints[j].mDoFCount == 3)
    {
      G.block(0, q_index, 6, 3) = (base_to_body * model_data.X_base[j].inverse()).toMatrix() *
                                  model_data.multdof3_S[j];
    }

    j = model.lambda[j];
  }
}

RBDL_DLLAPI void CalcPointJacobian6DRelative(const Model &model, ModelDatad &model_data,
                                             const Math::VectorNd &Q, unsigned int body_id,
                                             unsigned int respect_body_id,
                                             const Math::Vector3d &point_position,
                                             Math::MatrixNd &G, bool update_kinematics)
{
  if (update_kinematics)
  {
    RigidBodyDynamics::UpdateKinematicsCustom<double>(model, model_data, &Q, NULL, NULL);
  }

  Eigen::Vector3d p =
      CalcBodyToBaseCoordinates(model, model_data, Q, body_id, point_position, false);

  Eigen::Vector3d p_local =
      CalcBaseToBodyCoordinates<double>(model, model_data, Q, respect_body_id, p, false);

  SpatialTransformd point_trans = SpatialTransformd(Matrix3d::Identity(), p_local);


  Eigen::MatrixXd body_jac(6, model.dof_count);
  body_jac.setZero();
  Eigen::MatrixXd respect_body_jac(6, model.dof_count);
  respect_body_jac.setZero();

  CalcRelativeBodySpatialJacobian(model, model_data, Q, body_id, respect_body_id, body_jac, false);
  CalcRelativeBodySpatialJacobian(model, model_data, Q, respect_body_id, respect_body_id,
                                  respect_body_jac, false);

  G = point_trans.toMatrix() * (body_jac - respect_body_jac);
}

RBDL_DLLAPI void CalcPointJacobian6DRelative(Model &model, const Math::VectorNd &Q,
                                             unsigned int body_id, unsigned int respect_body_id,
                                             const Math::Vector3d &point_position,
                                             Math::MatrixNd &G, bool update_kinematics)
{
  CalcPointJacobian6DRelative(model, *model.getModelData(), Q, body_id, respect_body_id,
                              point_position, G, update_kinematics);
}


RBDL_DLLAPI void CalcPointJacobian6DBodyFrame(const Model &model, ModelDatad &model_data,
                                              const VectorNd &Q, unsigned int body_id,
                                              const Vector3d &point_position, MatrixNd &G,
                                              bool update_kinematics)
{
  LOG << "-------- " << __func__ << " --------" << std::endl;

  // update the Kinematics if necessary
  if (update_kinematics)
  {
    UpdateKinematicsCustom<double>(model, model_data, &Q, NULL, NULL);
  }

  SpatialTransformd point_trans = SpatialTransformd(
      Matrix3d::Identity(),
      CalcBodyToBaseCoordinates(model, model_data, Q, body_id, point_position, false));

  assert(G.rows() == 6 && G.cols() == model.qdot_size);

  unsigned int reference_body_id = body_id;

  if (model.IsFixedBodyId(body_id))
  {
    unsigned int fbody_id = body_id - model.fixed_body_discriminator;
    reference_body_id = model.mFixedBodies[fbody_id].mMovableParent;
  }

  unsigned int j = reference_body_id;

  while (j != 0)
  {
    unsigned int q_index = model.mJoints[j].q_index;

    // if(model.mJoints[j].mJointType != JointTypeCustom){
    if (model.mJoints[j].mDoFCount == 1)
    {
      G.block(0, q_index, 6, 1) =
          point_trans.apply(model_data.X_base[j].inverse().apply(model_data.S[j])).block(0, 0, 6, 1);
    }
    else if (model.mJoints[j].mDoFCount == 3)
    {
      G.block(0, q_index, 6, 3) =
          ((point_trans * model_data.X_base[j].inverse()).toMatrix() * model_data.multdof3_S[j])
              .block(0, 0, 6, 3);
    }

    j = model.lambda[j];
  }


  Vector3d body_r = CalcBodyToBaseCoordinates<double>(model, model_data, Q, body_id,
                                                      Eigen::Vector3d::Zero(), false);
  Matrix3d body_E = CalcBodyWorldOrientation<double>(model, model_data, Q, body_id, false);
  SpatialTransformd body_trans(body_E, body_r);

  G = point_trans.toMatrix() * (body_trans.toMatrix()) * G;
}

RBDL_DLLAPI void CalcPointJacobian6DBodyFrame(Model &model, const Math::VectorNd &Q,
                                              unsigned int body_id,
                                              const Math::Vector3d &point_position,
                                              Math::MatrixNd &G, bool update_kinematics)
{
  return CalcPointJacobian6DBodyFrame(model, *model.getModelData(), Q, body_id,
                                      point_position, G, update_kinematics);
}

RBDL_DLLAPI void CalcBodySpatialJacobian(const Model &model, ModelDatad &model_data,
                                         const VectorNd &Q, unsigned int body_id,
                                         MatrixNd &G, bool update_kinematics)
{
  LOG << "-------- " << __func__ << " --------" << std::endl;

  // update the Kinematics if necessary
  if (update_kinematics)
  {
    UpdateKinematicsCustom<double>(model, model_data, &Q, NULL, NULL);
  }

  assert(G.rows() == 6 && G.cols() == model.qdot_size);

  unsigned int reference_body_id = body_id;

  SpatialTransformd base_to_body;

  if (model.IsFixedBodyId(body_id))
  {
    unsigned int fbody_id = body_id - model.fixed_body_discriminator;

    reference_body_id = model.mFixedBodies[fbody_id].mMovableParent;

    base_to_body = model.mFixedBodies[fbody_id].mParentTransform *
                   model_data.X_base[reference_body_id];
  }
  else
  {
    base_to_body = model_data.X_base[reference_body_id];
  }

  unsigned int j = reference_body_id;

  while (j != 0)
  {
    unsigned int q_index = model.mJoints[j].q_index;

    // if(model.mJoints[j].mJointType != JointTypeCustom){
    if (model.mJoints[j].mDoFCount == 1)
    {
      G.block(0, q_index, 6, 1) =
          base_to_body.apply(model_data.X_base[j].inverse().apply(model_data.S[j]));
    }
    else if (model.mJoints[j].mDoFCount == 3)
    {
      G.block(0, q_index, 6, 3) = (base_to_body * model_data.X_base[j].inverse()).toMatrix() *
                                  model_data.multdof3_S[j];
    }
    //    }else if(model.mJoints[j].mJointType == JointTypeCustom) {
    //      unsigned int k = model.mJoints[j].custom_joint_index;

    //      G.block(0,q_index,6,model.mCustomJoints[k]->mDoFCount ) =
    //        (base_to_body * model_data.X_base[j].inverse()
    //        ).toMatrix() * model.mCustomJoints[k]->S;
    //    }

    j = model.lambda[j];
  }
}

RBDL_DLLAPI Vector3d CalcPointVelocity(const Model &model, ModelDatad &model_data,
                                       const VectorNd &Q, const VectorNd &QDot,
                                       unsigned int body_id, const Vector3d &point_position,
                                       bool update_kinematics)
{
  LOG << "-------- " << __func__ << " --------" << std::endl;
  assert(model.IsBodyId(body_id));
  assert(model.q_size == Q.size());
  assert(model.qdot_size == QDot.size());

  // Reset the velocity of the root body
  model_data.v[0].setZero();

  // update the Kinematics with zero acceleration
  if (update_kinematics)
  {
    UpdateKinematicsCustom<double>(model, model_data, &Q, &QDot, NULL);
  }

  unsigned int reference_body_id = body_id;
  Vector3d reference_point = point_position;

  if (model.IsFixedBodyId(body_id))
  {
    unsigned int fbody_id = body_id - model.fixed_body_discriminator;
    reference_body_id = model.mFixedBodies[fbody_id].mMovableParent;
    Vector3d base_coords =
        CalcBodyToBaseCoordinates(model, model_data, Q, body_id, point_position, false);
    reference_point = CalcBaseToBodyCoordinates(model, model_data, Q, reference_body_id,
                                                base_coords, false);
  }

  SpatialVectord point_spatial_velocity =
      SpatialTransformd(
          CalcBodyWorldOrientation(model, model_data, Q, reference_body_id, false).transpose(),
          reference_point)
          .apply(model_data.v[reference_body_id]);

  return Vector3d(point_spatial_velocity[3], point_spatial_velocity[4],
                  point_spatial_velocity[5]);
}

RBDL_DLLAPI Vector3d CalcPointVelocity(Model &model, const VectorNd &Q,
                                       const VectorNd &QDot, unsigned int body_id,
                                       const Vector3d &point_position, bool update_kinematics)
{
  return CalcPointVelocity(model, *model.getModelData(), Q, QDot, body_id, point_position,
                           update_kinematics);
}

RBDL_DLLAPI Vector3d CalcPointVelocityRelative(const Model &model, ModelDatad &model_data,
                                               const VectorNd &Q, const VectorNd &QDot,
                                               unsigned int body_id, unsigned int respect_body_id,
                                               const Vector3d &point_position,
                                               bool update_kinematics)
{
  LOG << "-------- " << __func__ << " --------" << std::endl;
  assert(model.IsBodyId(body_id));
  assert(model.q_size == Q.size());
  assert(model.qdot_size == QDot.size());

  // Reset the velocity of the root body
  model_data.v[0].setZero();

  // update the Kinematics with zero acceleration
  if (update_kinematics)
  {
    UpdateKinematicsCustom<double>(model, model_data, &Q, &QDot, NULL);
  }

  unsigned int reference_body_id = body_id;
  Vector3d reference_point = point_position;

  if (model.IsFixedBodyId(body_id))
  {
    unsigned int fbody_id = body_id - model.fixed_body_discriminator;
    reference_body_id = model.mFixedBodies[fbody_id].mMovableParent;
    Vector3d base_coords =
        CalcBodyToBaseCoordinates(model, model_data, Q, body_id, point_position, false);
    reference_point = CalcBaseToBodyCoordinates(model, model_data, Q, reference_body_id,
                                                base_coords, false);
  }

  SpatialTransformd base_to_body;
  if (model.IsFixedBodyId(body_id))
  {
    unsigned int fbody_id = body_id - model.fixed_body_discriminator;
    reference_body_id = model.mFixedBodies[fbody_id].mMovableParent;
    base_to_body =
        model.mFixedBodies[fbody_id].mParentTransform * model_data.X_base[respect_body_id];
  }
  else
  {
    base_to_body = model_data.X_base[respect_body_id];
  }

  SpatialVectord body_vel =
      base_to_body.apply(model_data.X_base[body_id].inverse().apply(model_data.v[body_id]));
  SpatialVectord respect_body_vel = base_to_body.apply(
      model_data.X_base[respect_body_id].inverse().apply(model_data.v[respect_body_id]));

  SpatialVectord relative_vel = body_vel - respect_body_vel;

  Eigen::Vector3d p =
      CalcBodyToBaseCoordinates(model, model_data, Q, body_id, point_position, false);
  Eigen::Vector3d p_local =
      CalcBaseToBodyCoordinates<double>(model, model_data, Q, respect_body_id, p, false);
  SpatialTransformd point_trans = SpatialTransformd(Matrix3d::Identity(), p_local);
  SpatialVectord point_spatial_velocity = point_trans.apply(relative_vel);


  return Vector3d(point_spatial_velocity[3], point_spatial_velocity[4],
                  point_spatial_velocity[5]);
}

RBDL_DLLAPI SpatialVectord CalcPointVelocityRelative6D(
    const Model &model, ModelDatad &model_data, const VectorNd &Q, const VectorNd &QDot,
    unsigned int body_id, unsigned int respect_body_id, const Vector3d &point_position,
    bool update_kinematics)
{
  LOG << "-------- " << __func__ << " --------" << std::endl;
  assert(model.IsBodyId(body_id));
  assert(model.q_size == Q.size());
  assert(model.qdot_size == QDot.size());

  // Reset the velocity of the root body
  model_data.v[0].setZero();

  // update the Kinematics with zero acceleration
  if (update_kinematics)
  {
    UpdateKinematicsCustom<double>(model, model_data, &Q, &QDot, NULL);
  }

  unsigned int reference_body_id = body_id;
  Vector3d reference_point = point_position;

  if (model.IsFixedBodyId(body_id))
  {
    unsigned int fbody_id = body_id - model.fixed_body_discriminator;
    reference_body_id = model.mFixedBodies[fbody_id].mMovableParent;
    Vector3d base_coords =
        CalcBodyToBaseCoordinates(model, model_data, Q, body_id, point_position, false);
    reference_point = CalcBaseToBodyCoordinates(model, model_data, Q, reference_body_id,
                                                base_coords, false);
  }

  SpatialTransformd base_to_body;
  if (model.IsFixedBodyId(body_id))
  {
    unsigned int fbody_id = body_id - model.fixed_body_discriminator;
    reference_body_id = model.mFixedBodies[fbody_id].mMovableParent;
    base_to_body =
        model.mFixedBodies[fbody_id].mParentTransform * model_data.X_base[respect_body_id];
  }
  else
  {
    base_to_body = model_data.X_base[respect_body_id];
  }

  SpatialVectord body_vel =
      base_to_body.apply(model_data.X_base[body_id].inverse().apply(model_data.v[body_id]));
  SpatialVectord respect_body_vel = base_to_body.apply(
      model_data.X_base[respect_body_id].inverse().apply(model_data.v[respect_body_id]));

  SpatialVectord relative_vel = body_vel - respect_body_vel;

  Eigen::Vector3d p =
      CalcBodyToBaseCoordinates(model, model_data, Q, body_id, point_position, false);
  Eigen::Vector3d p_local =
      CalcBaseToBodyCoordinates<double>(model, model_data, Q, respect_body_id, p, false);
  SpatialTransformd point_trans = SpatialTransformd(Matrix3d::Identity(), p_local);
  SpatialVectord point_spatial_velocity = point_trans.apply(relative_vel);


  return point_spatial_velocity;
}

RBDL_DLLAPI Vector3d CalcPointAngularVelocity(const Model &model, ModelDatad &model_data,
                                              const VectorNd &Q, const VectorNd &QDot,
                                              unsigned int body_id, const Vector3d &point_position,
                                              bool update_kinematics)
{
  LOG << "-------- " << __func__ << " --------" << std::endl;
  assert(model.IsBodyId(body_id));
  assert(model.q_size == Q.size());
  assert(model.qdot_size == QDot.size());

  // Reset the velocity of the root body
  model_data.v[0].setZero();

  // update the Kinematics with zero acceleration
  if (update_kinematics)
  {
    UpdateKinematicsCustom<double>(model, model_data, &Q, &QDot, NULL);
  }

  unsigned int reference_body_id = body_id;
  Vector3d reference_point = point_position;

  if (model.IsFixedBodyId(body_id))
  {
    unsigned int fbody_id = body_id - model.fixed_body_discriminator;
    reference_body_id = model.mFixedBodies[fbody_id].mMovableParent;
    Vector3d base_coords =
        CalcBodyToBaseCoordinates(model, model_data, Q, body_id, point_position, false);
    reference_point = CalcBaseToBodyCoordinates(model, model_data, Q, reference_body_id,
                                                base_coords, false);
  }

  SpatialVectord point_spatial_velocity =
      SpatialTransformd(
          CalcBodyWorldOrientation(model, model_data, Q, reference_body_id, false).transpose(),
          reference_point)
          .apply(model_data.v[reference_body_id]);

  return Vector3d(point_spatial_velocity[0], point_spatial_velocity[1],
                  point_spatial_velocity[2]);
}

RBDL_DLLAPI Vector3d CalcPointAngularVelocity(Model &model, const VectorNd &Q,
                                              const VectorNd &QDot, unsigned int body_id,
                                              const Vector3d &point_position,
                                              bool update_kinematics)
{
  return CalcPointAngularVelocity(model, *model.getModelData(), Q, QDot, body_id,
                                  point_position, update_kinematics);
}

RBDL_DLLAPI Math::Vector3d CalcPointVelocity(const Model &model, ModelDatad &model_data,
                                             const Math::VectorNd &Q,
                                             const Math::VectorNd &QDot, unsigned int body_id,
                                             const Math::Isometry3d &pose, bool update_kinematics)
{
  return CalcPointVelocity(model, model_data, Q, QDot, body_id,
                           Vector3d(pose.translation()), update_kinematics);
}

RBDL_DLLAPI Math::Vector3d CalcPointVelocity(Model &model, const Math::VectorNd &Q,
                                             const Math::VectorNd &QDot, unsigned int tip_id,
                                             const Math::Isometry3d &pose, bool update_kinematics)
{
  return CalcPointVelocity(model, Q, QDot, tip_id, Vector3d(pose.translation()), update_kinematics);
}

RBDL_DLLAPI Math::Vector3d CalcPointAngularVelocity(const Model &model, ModelDatad &model_data,
                                                    const Math::VectorNd &Q,
                                                    const Math::VectorNd &QDot,
                                                    unsigned int body_id,
                                                    const Math::Isometry3d &pose,
                                                    bool update_kinematics)
{
  return CalcPointAngularVelocity(model, model_data, Q, QDot, body_id,
                                  Vector3d(pose.translation()), update_kinematics);
}

RBDL_DLLAPI Math::Vector3d CalcPointAngularVelocity(Model &model, const Math::VectorNd &Q,
                                                    const Math::VectorNd &QDot,
                                                    unsigned int tip_id,
                                                    const Math::Isometry3d &pose,
                                                    bool update_kinematics)
{
  return CalcPointAngularVelocity(model, Q, QDot, tip_id, Vector3d(pose.translation()),
                                  update_kinematics);
}

RBDL_DLLAPI Math::SpatialVectord CalcPointVelocity6D(const Model &model, ModelDatad &model_data,
                                                     const Math::VectorNd &Q,
                                                     const Math::VectorNd &QDot,
                                                     unsigned int body_id,
                                                     const Math::Vector3d &point_position,
                                                     bool update_kinematics)
{
  LOG << "-------- " << __func__ << " --------" << std::endl;
  assert(model.IsBodyId(body_id));
  assert(model.q_size == Q.size());
  assert(model.qdot_size == QDot.size());

  // Reset the velocity of the root body
  model_data.v[0].setZero();

  // update the Kinematics with zero acceleration
  if (update_kinematics)
  {
    UpdateKinematicsCustom<double>(model, model_data, &Q, &QDot, NULL);
  }

  unsigned int reference_body_id = body_id;
  Vector3d reference_point = point_position;

  if (model.IsFixedBodyId(body_id))
  {
    unsigned int fbody_id = body_id - model.fixed_body_discriminator;
    reference_body_id = model.mFixedBodies[fbody_id].mMovableParent;
    Vector3d base_coords =
        CalcBodyToBaseCoordinates(model, model_data, Q, body_id, point_position, false);
    reference_point = CalcBaseToBodyCoordinates(model, model_data, Q, reference_body_id,
                                                base_coords, false);
  }

  return SpatialTransformd(
             CalcBodyWorldOrientation(model, model_data, Q, reference_body_id, false).transpose(),
             reference_point)
      .apply(model_data.v[reference_body_id]);
}

RBDL_DLLAPI Math::SpatialVectord CalcPointVelocity6D(Model &model, const Math::VectorNd &Q,
                                                     const Math::VectorNd &QDot,
                                                     unsigned int body_id,
                                                     const Math::Vector3d &point_position,
                                                     bool update_kinematics)
{
  return CalcPointVelocity6D(model, *model.getModelData(), Q, QDot, body_id,
                             point_position, update_kinematics);
}

RBDL_DLLAPI Vector3d CalcPointAcceleration(const Model &model, ModelDatad &model_data,
                                           const VectorNd &Q, const VectorNd &QDot,
                                           const VectorNd &QDDot, unsigned int body_id,
                                           const Vector3d &point_position, bool update_kinematics)
{
  LOG << "-------- " << __func__ << " --------" << std::endl;

  // Reset the velocity of the root body
  model_data.v[0].setZero();
  model_data.a[0].setZero();

  if (update_kinematics)
    UpdateKinematics(model, model_data, Q, QDot, QDDot);

  LOG << std::endl;

  unsigned int reference_body_id = body_id;
  Vector3d reference_point = point_position;

  if (model.IsFixedBodyId(body_id))
  {
    unsigned int fbody_id = body_id - model.fixed_body_discriminator;
    reference_body_id = model.mFixedBodies[fbody_id].mMovableParent;
    Vector3d base_coords =
        CalcBodyToBaseCoordinates(model, model_data, Q, body_id, point_position, false);
    reference_point = CalcBaseToBodyCoordinates(model, model_data, Q, reference_body_id,
                                                base_coords, false);
  }

  SpatialTransformd p_X_i(
      CalcBodyWorldOrientation(model, model_data, Q, reference_body_id, false).transpose(),
      reference_point);

  SpatialVectord p_v_i = p_X_i.apply(model_data.v[reference_body_id]);
  Vector3d a_dash =
      Vector3d(p_v_i[0], p_v_i[1], p_v_i[2]).cross(Vector3d(p_v_i[3], p_v_i[4], p_v_i[5]));
  SpatialVectord p_a_i = p_X_i.apply(model_data.a[reference_body_id]);

  return Vector3d(p_a_i[3] + a_dash[0], p_a_i[4] + a_dash[1], p_a_i[5] + a_dash[2]);
}

Vector3d CalcPointAcceleration(Model &model, const VectorNd &Q, const VectorNd &QDot,
                               const VectorNd &QDDot, unsigned int body_id,
                               const Vector3d &point_position, bool update_kinematics)
{
  return CalcPointAcceleration(model, *model.getModelData(), Q, QDot, QDDot, body_id,
                               point_position, update_kinematics);
}

Vector3d CalcPointAngularAcceleration(const Model &model, ModelDatad &model_data,
                                      const VectorNd &Q, const VectorNd &QDot,
                                      const VectorNd &QDDot, unsigned int body_id,
                                      const Vector3d &point_position, bool update_kinematics)
{
  LOG << "-------- " << __func__ << " --------" << std::endl;

  // Reset the velocity of the root body
  model_data.v[0].setZero();
  model_data.a[0].setZero();

  if (update_kinematics)
    UpdateKinematics(model, model_data, Q, QDot, QDDot);

  LOG << std::endl;

  unsigned int reference_body_id = body_id;
  Vector3d reference_point = point_position;

  if (model.IsFixedBodyId(body_id))
  {
    unsigned int fbody_id = body_id - model.fixed_body_discriminator;
    reference_body_id = model.mFixedBodies[fbody_id].mMovableParent;
    Vector3d base_coords =
        CalcBodyToBaseCoordinates(model, model_data, Q, body_id, point_position, false);
    reference_point = CalcBaseToBodyCoordinates(model, model_data, Q, reference_body_id,
                                                base_coords, false);
  }

  SpatialTransformd p_X_i(
      CalcBodyWorldOrientation(model, model_data, Q, reference_body_id, false).transpose(),
      reference_point);

  return p_X_i.apply(model_data.a[reference_body_id]).segment(0, 3);
}

Vector3d CalcPointAngularAcceleration(Model &model, const VectorNd &Q, const VectorNd &QDot,
                                      const VectorNd &QDDot, unsigned int body_id,
                                      const Vector3d &point_position, bool update_kinematics)
{
  return CalcPointAngularAcceleration(model, *model.getModelData(), Q, QDot, QDDot,
                                      body_id, point_position, update_kinematics);
}

RBDL_DLLAPI Vector3d CalcPointAccelerationBias(const Model &model, ModelDatad &model_data,
                                               const VectorNd &Q, const VectorNd &QDot,
                                               unsigned int body_id, const Vector3d &point_position,
                                               bool update_kinematics)
{
  LOG << "-------- " << __func__ << " --------" << std::endl;

  // Reset the velocity of the root body
  model_data.v[0].setZero();
  model_data.a[0].setZero();

  if (update_kinematics)
    UpdateKinematicsCustom<double>(model, model_data, &Q, &QDot, NULL);

  LOG << std::endl;

  unsigned int reference_body_id = body_id;
  Vector3d reference_point = point_position;

  if (model.IsFixedBodyId(body_id))
  {
    unsigned int fbody_id = body_id - model.fixed_body_discriminator;
    reference_body_id = model.mFixedBodies[fbody_id].mMovableParent;
    Vector3d base_coords =
        CalcBodyToBaseCoordinates(model, model_data, Q, body_id, point_position, false);
    reference_point = CalcBaseToBodyCoordinates(model, model_data, Q, reference_body_id,
                                                base_coords, false);
  }

  SpatialTransformd p_X_i(
      CalcBodyWorldOrientation(model, model_data, Q, reference_body_id, false).transpose(),
      reference_point);

  SpatialVectord p_v_i = p_X_i.apply(model_data.v[reference_body_id]);
  Vector3d a_dash =
      Vector3d(p_v_i[0], p_v_i[1], p_v_i[2]).cross(Vector3d(p_v_i[3], p_v_i[4], p_v_i[5]));
  SpatialVectord p_a_i = p_X_i.apply(model_data.a_bias[reference_body_id]);

  return Vector3d(p_a_i[3] + a_dash[0], p_a_i[4] + a_dash[1], p_a_i[5] + a_dash[2]);
}

Vector3d CalcPointAccelerationBias(Model &model, const VectorNd &Q, const VectorNd &QDot,
                                   unsigned int body_id, const Vector3d &point_position,
                                   bool update_kinematics)
{
  return CalcPointAccelerationBias(model, *model.getModelData(), Q, QDot, body_id,
                                   point_position, update_kinematics);
}

RBDL_DLLAPI SpatialVectord CalcPointAcceleration6D(const Model &model, ModelDatad &model_data,
                                                   const VectorNd &Q, const VectorNd &QDot,
                                                   const VectorNd &QDDot, unsigned int body_id,
                                                   const Vector3d &point_position,
                                                   bool update_kinematics)
{
  LOG << "-------- " << __func__ << " --------" << std::endl;

  // Reset the velocity of the root body
  model_data.v[0].setZero();
  model_data.a[0].setZero();

  if (update_kinematics)
    UpdateKinematics(model, model_data, Q, QDot, QDDot);

  LOG << std::endl;

  unsigned int reference_body_id = body_id;
  Vector3d reference_point = point_position;

  if (model.IsFixedBodyId(body_id))
  {
    unsigned int fbody_id = body_id - model.fixed_body_discriminator;
    reference_body_id = model.mFixedBodies[fbody_id].mMovableParent;
    Vector3d base_coords =
        CalcBodyToBaseCoordinates(model, model_data, Q, body_id, point_position, false);
    reference_point = CalcBaseToBodyCoordinates(model, model_data, Q, reference_body_id,
                                                base_coords, false);
  }

  SpatialTransformd p_X_i(
      CalcBodyWorldOrientation(model, model_data, Q, reference_body_id, false).transpose(),
      reference_point);

  SpatialVectord p_v_i = p_X_i.apply(model_data.v[reference_body_id]);
  Vector3d a_dash =
      Vector3d(p_v_i[0], p_v_i[1], p_v_i[2]).cross(Vector3d(p_v_i[3], p_v_i[4], p_v_i[5]));
  return (p_X_i.apply(model_data.a[reference_body_id]) +
          SpatialVectord(0, 0, 0, a_dash[0], a_dash[1], a_dash[2]));
}

SpatialVectord CalcPointAcceleration6D(Model &model, const VectorNd &Q, const VectorNd &QDot,
                                       const VectorNd &QDDot, unsigned int body_id,
                                       const Vector3d &point_position, bool update_kinematics)
{
  return CalcPointAcceleration6D(model, *model.getModelData(), Q, QDot, QDDot, body_id,
                                 point_position, update_kinematics);
}

RBDL_DLLAPI SpatialVectord CalcPointAcceleration6DBias(
    const Model &model, ModelDatad &model_data, const VectorNd &Q, const VectorNd &QDot,
    unsigned int body_id, const Vector3d &point_position, bool update_kinematics)
{
  LOG << "-------- " << __func__ << " --------" << std::endl;

  // Reset the velocity of the root body
  model_data.v[0].setZero();
  model_data.a[0].setZero();

  if (update_kinematics)
    UpdateKinematicsCustom<double>(model, model_data, &Q, &QDot, NULL);

  LOG << std::endl;

  unsigned int reference_body_id = body_id;
  Vector3d reference_point = point_position;

  if (model.IsFixedBodyId(body_id))
  {
    unsigned int fbody_id = body_id - model.fixed_body_discriminator;
    reference_body_id = model.mFixedBodies[fbody_id].mMovableParent;
    Vector3d base_coords =
        CalcBodyToBaseCoordinates(model, model_data, Q, body_id, point_position, false);
    reference_point = CalcBaseToBodyCoordinates(model, model_data, Q, reference_body_id,
                                                base_coords, false);
  }

  SpatialTransformd p_X_i(
      CalcBodyWorldOrientation(model, model_data, Q, reference_body_id, false).transpose(),
      reference_point);

  SpatialVectord p_v_i = p_X_i.apply(model_data.v[reference_body_id]);
  Vector3d a_dash =
      Vector3d(p_v_i[0], p_v_i[1], p_v_i[2]).cross(Vector3d(p_v_i[3], p_v_i[4], p_v_i[5]));
  return (p_X_i.apply(model_data.a_bias[reference_body_id]) +
          SpatialVectord(0, 0, 0, a_dash[0], a_dash[1], a_dash[2]));
}

SpatialVectord CalcPointAcceleration6DBias(Model &model, const VectorNd &Q,
                                           const VectorNd &QDot, unsigned int body_id,
                                           const Vector3d &point_position, bool update_kinematics)
{
  return CalcPointAcceleration6DBias(model, *model.getModelData(), Q, QDot, body_id,
                                     point_position, update_kinematics);
}

RBDL_DLLAPI
Math::SpatialVectord CalcPointAcceleration6DBias(const Model &model, ModelDatad &model_data,
                                                 const Math::VectorNd &Q,
                                                 const Math::VectorNd &QDot, unsigned int body_id,
                                                 const Math::Isometry3d &pose,
                                                 bool update_kinematics)
{
  return CalcPointAcceleration6DBias(model, model_data, Q, QDot, body_id,
                                     Vector3d(pose.translation()), update_kinematics);
}

RBDL_DLLAPI
Math::SpatialVectord CalcPointAcceleration6DBias(Model &model, const Math::VectorNd &Q,
                                                 const Math::VectorNd &QDot, unsigned int body_id,
                                                 const Math::Isometry3d &pose,
                                                 bool update_kinematics)
{
  return CalcPointAcceleration6DBias(model, Q, QDot, body_id,
                                     Vector3d(pose.translation()), update_kinematics);
}

RBDL_DLLAPI bool InverseKinematics(Model &model, ModelDatad &model_data, const VectorNd &Qinit,
                                   const std::vector<unsigned int> &body_id,
                                   const std::vector<Vector3d> &body_point,
                                   const std::vector<Vector3d> &target_pos, VectorNd &Qres,
                                   double step_tol, double lambda, unsigned int max_iter)
{
  assert(Qinit.size() == model.q_size);
  assert(body_id.size() == body_point.size());
  assert(body_id.size() == target_pos.size());

  MatrixNd J = MatrixNd::Zero(3 * body_id.size(), model.qdot_size);
  VectorNd e = VectorNd::Zero(3 * body_id.size());

  Qres = Qinit;

  for (unsigned int ik_iter = 0; ik_iter < max_iter; ik_iter++)
  {
    UpdateKinematicsCustom<double>(model, model_data, &Qres, NULL, NULL);

    for (unsigned int k = 0; k < body_id.size(); k++)
    {
      MatrixNd G(MatrixNd::Zero(3, model.qdot_size));
      CalcPointJacobian(model, model_data, Qres, body_id[k], body_point[k], G, false);
      Vector3d point_base = CalcBodyToBaseCoordinates(model, model_data, Qres, body_id[k],
                                                      body_point[k], false);
      LOG << "current_pos = " << point_base.transpose() << std::endl;

      for (unsigned int i = 0; i < 3; i++)
      {
        for (unsigned int j = 0; j < model.qdot_size; j++)
        {
          unsigned int row = k * 3 + i;
          LOG << "i = " << i << " j = " << j << " k = " << k << " row = " << row
              << " col = " << j << std::endl;
          J(row, j) = G(i, j);
        }

        e[k * 3 + i] = target_pos[k][i] - point_base[i];
      }
    }

    LOG << "J = " << J << std::endl;
    LOG << "e = " << e.transpose() << std::endl;

    // abort if we are getting "close"
    if (e.norm() < step_tol)
    {
      LOG << "Reached target close enough after " << ik_iter << " steps" << std::endl;
      return true;
    }

    MatrixNd JJTe_lambda2_I =
        J * J.transpose() + lambda * lambda * MatrixNd::Identity(e.size(), e.size());

    VectorNd z(body_id.size() * 3);
#ifndef RBDL_USE_SIMPLE_MATH
    z = JJTe_lambda2_I.colPivHouseholderQr().solve(e);
#else
    bool solve_successful = LinSolveGaussElimPivot(JJTe_lambda2_I, e, z);
    assert(solve_successful);
#endif

    LOG << "z = " << z << std::endl;

    VectorNd delta_theta = J.transpose() * z;
    LOG << "change = " << delta_theta << std::endl;

    Qres = Qres + delta_theta;
    LOG << "Qres = " << Qres.transpose() << std::endl;

    if (delta_theta.norm() < step_tol)
    {
      LOG << "reached convergence after " << ik_iter << " steps" << std::endl;
      return true;
    }

    VectorNd test_1(z.size());
    VectorNd test_res(z.size());

    test_1.setZero();

    for (unsigned int i = 0; i < z.size(); i++)
    {
      test_1[i] = 1.;

      VectorNd test_delta = J.transpose() * test_1;

      test_res[i] = test_delta.squaredNorm();

      test_1[i] = 0.;
    }

    LOG << "test_res = " << test_res.transpose() << std::endl;
  }

  return false;
}

RBDL_DLLAPI
Vector3d CalcAngularVelocityfromMatrix(const Matrix3d &RotMat)
{
  double tol = 1e-12;

  Vector3d l = Vector3d(RotMat(2, 1) - RotMat(1, 2), RotMat(0, 2) - RotMat(2, 0),
                        RotMat(1, 0) - RotMat(0, 1));
  if (l.norm() > tol)
  {
    double preFactor = atan2(l.norm(), (RotMat.trace() - 1.0)) / l.norm();
    return preFactor * l;
  }
  else if ((RotMat(0, 0) > 0 && RotMat(1, 1) > 0 && RotMat(2, 2) > 0) || l.norm() < tol)
  {
    return Vector3d(0.0, 0.0, 0.0);
  }
  else
  {
    double PI = atan(1) * 4.0;
    return Vector3d(PI / 2 * (RotMat(0, 0) + 1.0), PI / 2 * (RotMat(1, 1) + 1.0),
                    PI / 2 * (RotMat(2, 2) + 1.0));
  }
}

RBDL_DLLAPI
InverseKinematicsConstraintSet::InverseKinematicsConstraintSet()
{
  lambda = 1e-9;
  num_steps = 0;
  max_steps = 300;
  step_tol = 1e-12;
  constraint_tol = 1e-12;
  num_constraints = 0;
}

RBDL_DLLAPI
unsigned int InverseKinematicsConstraintSet::AddPointConstraint(unsigned int body_id,
                                                                const Vector3d &body_point,
                                                                const Vector3d &target_pos)
{
  constraint_type.push_back(ConstraintTypePosition);
  body_ids.push_back(body_id);
  body_points.push_back(body_point);
  target_positions.push_back(target_pos);
  target_orientations.push_back(Matrix3d::Zero(3, 3));
  constraint_row_index.push_back(num_constraints);
  num_constraints = num_constraints + 3;
  return constraint_type.size() - 1;
}

RBDL_DLLAPI
unsigned int InverseKinematicsConstraintSet::AddOrientationConstraint(unsigned int body_id,
                                                                      const Matrix3d &target_orientation)
{
  constraint_type.push_back(ConstraintTypeOrientation);
  body_ids.push_back(body_id);
  body_points.push_back(Vector3d::Zero());
  target_positions.push_back(Vector3d::Zero());
  target_orientations.push_back(target_orientation);
  constraint_row_index.push_back(num_constraints);
  num_constraints = num_constraints + 3;
  return constraint_type.size() - 1;
}

RBDL_DLLAPI
unsigned int InverseKinematicsConstraintSet::AddFullConstraint(unsigned int body_id,
                                                               const Vector3d &body_point,
                                                               const Vector3d &target_pos,
                                                               const Matrix3d &target_orientation)
{
  constraint_type.push_back(ConstraintTypeFull);
  body_ids.push_back(body_id);
  body_points.push_back(body_point);
  target_positions.push_back(target_pos);
  target_orientations.push_back(target_orientation);
  constraint_row_index.push_back(num_constraints);
  num_constraints = num_constraints + 6;
  return constraint_type.size() - 1;
}

RBDL_DLLAPI
unsigned int InverseKinematicsConstraintSet::ClearConstraints()
{
  for (unsigned int i = 0; i < constraint_type.size(); i++)
  {
    constraint_type.pop_back();
    body_ids.pop_back();
    body_points.pop_back();
    target_positions.pop_back();
    target_orientations.pop_back();
    num_constraints = 0;
  }
  return constraint_type.size();
}


RBDL_DLLAPI
bool InverseKinematics(Model &model, const Math::VectorNd &Qinit,
                       InverseKinematicsConstraintSet &CS, Math::VectorNd &Qres)
{
  assert(Qinit.size() == model.q_size);
  assert(Qres.size() == Qinit.size());

  CS.J = MatrixNd::Zero(CS.num_constraints, model.qdot_size);
  CS.e = VectorNd::Zero(CS.num_constraints);

  Qres = Qinit;

  for (CS.num_steps = 0; CS.num_steps < CS.max_steps; CS.num_steps++)
  {
    UpdateKinematicsCustom(model, &Qres, NULL, NULL);

    for (unsigned int k = 0; k < CS.body_ids.size(); k++)
    {
      CS.G = MatrixNd::Zero(6, model.qdot_size);
      CalcPointJacobian6D(model, Qres, CS.body_ids[k], CS.body_points[k], CS.G, false);
      Vector3d point_base =
          CalcBodyToBaseCoordinates(model, Qres, CS.body_ids[k], CS.body_points[k], false);
      Matrix3d R = CalcBodyWorldOrientation(model, Qres, CS.body_ids[k], false);
      Vector3d angular_velocity =
          R.transpose() *
          CalcAngularVelocityfromMatrix(R * CS.target_orientations[k].transpose());

      // assign offsets and Jacobians
      if (CS.constraint_type[k] == InverseKinematicsConstraintSet::ConstraintTypeFull)
      {
        for (unsigned int i = 0; i < 3; i++)
        {
          unsigned int row = CS.constraint_row_index[k] + i;
          CS.e[row + 3] = CS.target_positions[k][i] - point_base[i];
          CS.e[row] = angular_velocity[i];
          for (unsigned int j = 0; j < model.qdot_size; j++)
          {
            CS.J(row + 3, j) = CS.G(i + 3, j);
            CS.J(row, j) = CS.G(i, j);
          }
        }
      }
      else if (CS.constraint_type[k] == InverseKinematicsConstraintSet::ConstraintTypeOrientation)
      {
        for (unsigned int i = 0; i < 3; i++)
        {
          unsigned int row = CS.constraint_row_index[k] + i;
          CS.e[row] = angular_velocity[i];
          for (unsigned int j = 0; j < model.qdot_size; j++)
          {
            CS.J(row, j) = CS.G(i, j);
          }
        }
      }
      else if (CS.constraint_type[k] == InverseKinematicsConstraintSet::ConstraintTypePosition)
      {
        for (unsigned int i = 0; i < 3; i++)
        {
          unsigned int row = CS.constraint_row_index[k] + i;
          CS.e[row] = CS.target_positions[k][i] - point_base[i];
          for (unsigned int j = 0; j < model.qdot_size; j++)
          {
            CS.J(row, j) = CS.G(i + 3, j);
          }
        }
      }
      else
      {
        assert(false && !"Invalid inverse kinematics constraint");
      }
    }

    LOG << "J = " << CS.J << std::endl;
    LOG << "e = " << CS.e.transpose() << std::endl;
    CS.error_norm = CS.e.norm();

    // abort if we are getting "close"
    if (CS.error_norm < CS.step_tol)
    {
      LOG << "Reached target close enough after " << CS.num_steps << " steps" << std::endl;
      return true;
    }

    //     // "task space" from puppeteer
    //     MatrixNd Ek = MatrixNd::Zero (CS.e.size(), CS.e.size());
    //
    //     for (size_t ei = 0; ei < CS.e.size(); ei ++) {
    // //      Ek(ei,ei) = CS.error_norm * CS.error_norm * 0.5 + CS.lambda;
    //       Ek(ei,ei) = CS.e[ei]*CS.e[ei] * 0.5 + CS.lambda;
    //     }
    //
    //     MatrixNd JJT_Ek_wnI = CS.J * CS.J.transpose() + Ek;
    //
    //     VectorNd delta_theta = CS.J.transpose() *
    //     JJT_Ek_wnI.colPivHouseholderQr().solve (CS.e);
    //
    //     LOG << "change = " << delta_theta << std::endl;


    // "joint space" from puppeteer

    double Ek = 0.;

    for (unsigned int ei = 0; ei < CS.e.size(); ei++)
    {
      Ek += CS.e[ei] * CS.e[ei] * 0.5;
    }

    VectorNd ek = CS.J.transpose() * CS.e;
    MatrixNd Wn = MatrixNd::Zero(Qres.size(), Qres.size());

    assert(ek.size() == Qres.size());

    for (unsigned int wi = 0; wi < Qres.size(); wi++)
    {
      Wn(wi, wi) = ek[wi] * ek[wi] * 0.5 + CS.lambda;
      //      Wn(wi, wi) = Ek + 1.0e-3;
    }

    MatrixNd A = CS.J.transpose() * CS.J + Wn;
    VectorNd delta_theta = A.colPivHouseholderQr().solve(CS.J.transpose() * CS.e);

    Qres = Qres + delta_theta;
    if (delta_theta.norm() < CS.step_tol)
    {
      LOG << "reached convergence after " << CS.num_steps << " steps" << std::endl;
      return true;
    }
  }

  return false;
}
}
