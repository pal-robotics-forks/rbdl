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


Vector3d CalcBodyToBaseCoordinates(Model &model, ModelDatad &model_data, const VectorNd &Q,
                                   unsigned int body_id, const Vector3d &point_body_coordinates,
                                   bool update_kinematics)
{
  return CalcBodyToBaseCoordinates<double>(model, model_data, Q, body_id,
                                             point_body_coordinates, update_kinematics);
}

Math::Matrix3d CalcBodyWorldOrientation(Model &model, const Math::VectorNd &Q,
                                        const unsigned int body_id, bool update_kinematics)
{
  return CalcBodyWorldOrientation<double>(model, *model.getModelData(), Q, body_id,
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
                                   bool update_kinematics) {
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
                                     const Math::Isometry3d &pose,
                                     Math::MatrixNd &G, bool update_kinematics)
{
  return CalcPointJacobian6D(model, model_data, Q, body_id, Vector3d(pose.translation()), G, update_kinematics);
}

RBDL_DLLAPI void CalcPointJacobian6D (
    Model &model,
    const Math::VectorNd &Q,
    unsigned int body_id,
    const Math::Isometry3d &pose,
    Math::MatrixNd &G,
    bool update_kinematics)
{
  return CalcPointJacobian6D(model, Q, body_id, Vector3d(pose.translation()), G, update_kinematics);
}

RBDL_DLLAPI void CalcPointJacobian6D(const Model &model, ModelDatad &model_data,
                                     const VectorNd &Q, unsigned int body_id,
                                     const Vector3d &point_position, MatrixNd &G,
                                     bool update_kinematics) {
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
    //    } else {
    //      unsigned int k = model.mJoints[j].custom_joint_index;

    //      G.block(0, q_index, 6, model.mCustomJoints[k]->mDoFCount)
    //        = ((point_trans
    //              * model_data.X_base[j].inverse()).toMatrix()
    //            * model.mCustomJoints[k]->S).block(
    //              0,0,6,model.mCustomJoints[k]->mDoFCount);
    //    }

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
    //    } else {
    //      unsigned int k = model.mJoints[j].custom_joint_index;

    //      G.block(0, q_index, 6, model.mCustomJoints[k]->mDoFCount)
    //        = ((point_trans
    //              * model_data.X_base[j].inverse()).toMatrix()
    //            * model.mCustomJoints[k]->S).block(
    //              0,0,6,model.mCustomJoints[k]->mDoFCount);
    //    }

    j = model.lambda[j];
  }

  Eigen::Vector3d zero;
  zero.setZero();
  Vector3d body_r = CalcBodyToBaseCoordinates<double>(model, model_data, Q, body_id, zero, false);
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
                                             const Math::Isometry3d &pose,
                                             bool update_kinematics)
{
  return CalcPointVelocity(model, model_data, Q, QDot, body_id, Vector3d(pose.translation()), update_kinematics);
}

RBDL_DLLAPI Math::Vector3d CalcPointVelocity(
    Model &model,
    const Math::VectorNd &Q,
    const Math::VectorNd &QDot,
    unsigned int tip_id,
    const Math::Isometry3d &pose,
    bool update_kinematics)
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
  return CalcPointAngularVelocity(model, model_data, Q, QDot, body_id, Vector3d(pose.translation()), update_kinematics);
}

RBDL_DLLAPI Math::Vector3d CalcPointAngularVelocity(
    Model &model,
    const Math::VectorNd &Q,
    const Math::VectorNd &QDot,
    unsigned int tip_id,
    const Math::Isometry3d &pose,
    bool update_kinematics)
{
  return CalcPointAngularVelocity(model, Q, QDot, tip_id, Vector3d(pose.translation()), update_kinematics);
}

RBDL_DLLAPI Math::SpatialVectord CalcPointVelocity6D(const Model &model, ModelDatad &model_data,
                                                     const Math::VectorNd &Q,
                                                     const Math::VectorNd &QDot,
                                                     unsigned int body_id,
                                                     const Math::Vector3d &point_position,
                                                     bool update_kinematics) {
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
  return CalcPointAcceleration6DBias(model, model_data, Q, QDot, body_id, Vector3d(pose.translation()), update_kinematics);
}

RBDL_DLLAPI
  Math::SpatialVectord CalcPointAcceleration6DBias (
      Model &model,
      const Math::VectorNd &Q,
      const Math::VectorNd &QDot,
      unsigned int body_id,
      const Math::Isometry3d &pose,
      bool update_kinematics)
{
  return CalcPointAcceleration6DBias(model, Q, QDot, body_id, Vector3d(pose.translation()), update_kinematics);
}

RBDL_DLLAPI bool InverseKinematics(Model &model, ModelDatad &model_data, const VectorNd &Qinit,
                                   const std::vector<unsigned int> &body_id,
                                   const std::vector<Vector3d> &body_point,
                                   const std::vector<Vector3d> &target_pos, VectorNd &Qres,
                                   double step_tol, double lambda, unsigned int max_iter) {
  assert (Qinit.size() == model.q_size);
  assert (body_id.size() == body_point.size());
  assert (body_id.size() == target_pos.size());

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
}
