/*
 * RBDL - Rigid Body Dynamics Library
 * Copyright (c) 2011-2016 Martin Felis <martin.felis@iwr.uni-heidelberg.de>
 *
 * Licensed under the zlib license. See LICENSE for more details.
 */

#include <iostream>
#include <limits>
#include <assert.h>

#include "rbdl/Logging.h"

#include "rbdl/Model.h"
#include "rbdl/Joint.h"

namespace RigidBodyDynamics {

using namespace Math;

RBDL_DLLAPI void jcalc (
    const Model &model,
    ModelData &model_data,
    unsigned int joint_id,
    const VectorNd &q,
    const VectorNd &qdot
    ) {
  // exception if we calculate it for the root body
  assert (joint_id > 0);

  if (model.mJoints[joint_id].mJointType == JointTypeRevoluteX) {
    model_data.X_J[joint_id] = Xrotx (q[model.mJoints[joint_id].q_index]);
    model_data.v_J[joint_id][0] = qdot[model.mJoints[joint_id].q_index];
  } else if (model.mJoints[joint_id].mJointType == JointTypeRevoluteY) {
    model_data.X_J[joint_id] = Xroty (q[model.mJoints[joint_id].q_index]);
    model_data.v_J[joint_id][1] = qdot[model.mJoints[joint_id].q_index];
  } else if (model.mJoints[joint_id].mJointType == JointTypeRevoluteZ) {
    model_data.X_J[joint_id] = Xrotz (q[model.mJoints[joint_id].q_index]);
    model_data.v_J[joint_id][2] = qdot[model.mJoints[joint_id].q_index];
  } else if (model.mJoints[joint_id].mDoFCount == 1 &&
      model.mJoints[joint_id].mJointType != JointTypeCustom) {
    model_data.X_J[joint_id] = jcalc_XJ (model, model_data, joint_id, q);

    model_data.v_J[joint_id] =
      model_data.S[joint_id] * qdot[model.mJoints[joint_id].q_index];
  } else if (model.mJoints[joint_id].mJointType == JointTypeSpherical) {
    model_data.X_J[joint_id] =
      SpatialTransformd (model.GetQuaternion(joint_id, q).toMatrix(),
          Vector3d (0., 0., 0.));

    model_data.multdof3_S[joint_id](0,0) = 1.;
    model_data.multdof3_S[joint_id](1,1) = 1.;
    model_data.multdof3_S[joint_id](2,2) = 1.;

    Vector3d omega (qdot[model.mJoints[joint_id].q_index],
        qdot[model.mJoints[joint_id].q_index+1],
        qdot[model.mJoints[joint_id].q_index+2]);

    model_data.v_J[joint_id] = SpatialVectord(
        omega[0], omega[1], omega[2],
        0., 0., 0.);
  } else if (model.mJoints[joint_id].mJointType == JointTypeEulerZYX) {
    double q0 = q[model.mJoints[joint_id].q_index];
    double q1 = q[model.mJoints[joint_id].q_index + 1];
    double q2 = q[model.mJoints[joint_id].q_index + 2];

    double s0 = sin (q0);
    double c0 = cos (q0);
    double s1 = sin (q1);
    double c1 = cos (q1);
    double s2 = sin (q2);
    double c2 = cos (q2);

    model_data.X_J[joint_id].E = Matrix3d(
        c0 * c1, s0 * c1, -s1,
        c0 * s1 * s2 - s0 * c2, s0 * s1 * s2 + c0 * c2, c1 * s2,
        c0 * s1 * c2 + s0 * s2, s0 * s1 * c2 - c0 * s2, c1 * c2
        );

    model_data.multdof3_S[joint_id](0,0) = -s1;
    model_data.multdof3_S[joint_id](0,2) = 1.;

    model_data.multdof3_S[joint_id](1,0) = c1 * s2;
    model_data.multdof3_S[joint_id](1,1) = c2;

    model_data.multdof3_S[joint_id](2,0) = c1 * c2;
    model_data.multdof3_S[joint_id](2,1) = - s2;

    double qdot0 = qdot[model.mJoints[joint_id].q_index];
    double qdot1 = qdot[model.mJoints[joint_id].q_index + 1];
    double qdot2 = qdot[model.mJoints[joint_id].q_index + 2];

    model_data.v_J[joint_id] =
      model_data.multdof3_S[joint_id] * Vector3d (qdot0, qdot1, qdot2);

    model_data.c_J[joint_id].set(
        -c1*qdot0*qdot1,
        -s1*s2*qdot0*qdot1 + c1*c2*qdot0*qdot2 - s2*qdot1*qdot2,
        -s1*c2*qdot0*qdot1 - c1*s2*qdot0*qdot2 - c2*qdot1*qdot2,
        0.,0., 0.);
  } else if (model.mJoints[joint_id].mJointType == JointTypeEulerXYZ) {
    double q0 = q[model.mJoints[joint_id].q_index];
    double q1 = q[model.mJoints[joint_id].q_index + 1];
    double q2 = q[model.mJoints[joint_id].q_index + 2];

    double s0 = sin (q0);
    double c0 = cos (q0);
    double s1 = sin (q1);
    double c1 = cos (q1);
    double s2 = sin (q2);
    double c2 = cos (q2);

    model_data.X_J[joint_id].E = Matrix3d(
        c2 * c1, s2 * c0 + c2 * s1 * s0, s2 * s0 - c2 * s1 * c0,
        -s2 * c1, c2 * c0 - s2 * s1 * s0, c2 * s0 + s2 * s1 * c0,
        s1, -c1 * s0, c1 * c0
        );

    model_data.multdof3_S[joint_id](0,0) = c2 * c1;
    model_data.multdof3_S[joint_id](0,1) = s2;

    model_data.multdof3_S[joint_id](1,0) = -s2 * c1;
    model_data.multdof3_S[joint_id](1,1) = c2;

    model_data.multdof3_S[joint_id](2,0) = s1;
    model_data.multdof3_S[joint_id](2,2) = 1.;

    double qdot0 = qdot[model.mJoints[joint_id].q_index];
    double qdot1 = qdot[model.mJoints[joint_id].q_index + 1];
    double qdot2 = qdot[model.mJoints[joint_id].q_index + 2];

    model_data.v_J[joint_id] =
      model_data.multdof3_S[joint_id] * Vector3d (qdot0, qdot1, qdot2);

    model_data.c_J[joint_id].set(
        -s2*c1*qdot2*qdot0 - c2*s1*qdot1*qdot0 + c2*qdot2*qdot1,
        -c2*c1*qdot2*qdot0 + s2*s1*qdot1*qdot0 - s2*qdot2*qdot1,
        c1*qdot1*qdot0,
        0., 0., 0.
        );
  } else if (model.mJoints[joint_id].mJointType == JointTypeEulerYXZ) {
    double q0 = q[model.mJoints[joint_id].q_index];
    double q1 = q[model.mJoints[joint_id].q_index + 1];
    double q2 = q[model.mJoints[joint_id].q_index + 2];

    double s0 = sin (q0);
    double c0 = cos (q0);
    double s1 = sin (q1);
    double c1 = cos (q1);
    double s2 = sin (q2);
    double c2 = cos (q2);

    model_data.X_J[joint_id].E = Matrix3d(
        c2 * c0 + s2 * s1 * s0, s2 * c1, -c2 * s0 + s2 * s1 * c0,
        -s2 * c0 + c2 * s1 * s0, c2 * c1,  s2 * s0 + c2 * s1 * c0,
        c1 * s0,    - s1,                 c1 * c0);

    model_data.multdof3_S[joint_id](0,0) = s2 * c1;
    model_data.multdof3_S[joint_id](0,1) = c2;

    model_data.multdof3_S[joint_id](1,0) = c2 * c1;
    model_data.multdof3_S[joint_id](1,1) = -s2;

    model_data.multdof3_S[joint_id](2,0) = -s1;
    model_data.multdof3_S[joint_id](2,2) = 1.;

    double qdot0 = qdot[model.mJoints[joint_id].q_index];
    double qdot1 = qdot[model.mJoints[joint_id].q_index + 1];
    double qdot2 = qdot[model.mJoints[joint_id].q_index + 2];

    model_data.v_J[joint_id] =
      model_data.multdof3_S[joint_id] * Vector3d (qdot0, qdot1, qdot2);

    model_data.c_J[joint_id].set(
        c2*c1*qdot2*qdot0 - s2*s1*qdot1*qdot0 - s2*qdot2*qdot1,
        -s2*c1*qdot2*qdot0 - c2*s1*qdot1*qdot0 - c2*qdot2*qdot1,
        -c1*qdot1*qdot0,
        0., 0., 0.
        );
  } else if(model.mJoints[joint_id].mJointType == JointTypeTranslationXYZ){
    double q0 = q[model.mJoints[joint_id].q_index];
    double q1 = q[model.mJoints[joint_id].q_index + 1];
    double q2 = q[model.mJoints[joint_id].q_index + 2];

    model_data.X_J[joint_id].E = Matrix3d::Identity();
    model_data.X_J[joint_id].r = Vector3d (q0, q1, q2);

    model_data.multdof3_S[joint_id](3,0) = 1.;
    model_data.multdof3_S[joint_id](4,1) = 1.;
    model_data.multdof3_S[joint_id](5,2) = 1.;

    double qdot0 = qdot[model.mJoints[joint_id].q_index];
    double qdot1 = qdot[model.mJoints[joint_id].q_index + 1];
    double qdot2 = qdot[model.mJoints[joint_id].q_index + 2];

    model_data.v_J[joint_id] =
      model_data.multdof3_S[joint_id] * Vector3d (qdot0, qdot1, qdot2);

    model_data.c_J[joint_id].set(0., 0., 0., 0., 0., 0.);
  } else if (model.mJoints[joint_id].mJointType == JointTypeCustom) {
    const Joint &joint = model.mJoints[joint_id];
    CustomJoint *custom_joint = 
      model.mCustomJoints[joint.custom_joint_index];
    custom_joint->jcalc (model, model_data, joint_id, q, qdot);
  } else {
    std::cerr << "Error: invalid joint type " << model.mJoints[joint_id].mJointType << " at id " << joint_id << std::endl;
    abort();
  }

  model_data.X_lambda[joint_id] = model_data.X_J[joint_id] * model.X_T[joint_id];
}

RBDL_DLLAPI Math::SpatialTransformd jcalc_XJ (
    const Model &model,
    ModelData &model_data,
    unsigned int joint_id,
    const Math::VectorNd &q) {
  // exception if we calculate it for the root body
  assert (joint_id > 0);

  if (model.mJoints[joint_id].mDoFCount == 1
      && model.mJoints[joint_id].mJointType != JointTypeCustom) {
    if (model.mJoints[joint_id].mJointType == JointTypeRevolute) {
      return Xrot (q[model.mJoints[joint_id].q_index], Vector3d (
            model.mJoints[joint_id].mJointAxes[0][0],
            model.mJoints[joint_id].mJointAxes[0][1],
            model.mJoints[joint_id].mJointAxes[0][2]
            ));
    } else if (model.mJoints[joint_id].mJointType == JointTypePrismatic) {
      return Xtrans ( Vector3d (
            model.mJoints[joint_id].mJointAxes[0][3]
            * q[model.mJoints[joint_id].q_index],
            model.mJoints[joint_id].mJointAxes[0][4] 
            * q[model.mJoints[joint_id].q_index],
            model.mJoints[joint_id].mJointAxes[0][5] 
            * q[model.mJoints[joint_id].q_index]
            )
          );
    }
  }
  std::cerr << "Error: invalid joint type!" << std::endl;
  abort();
  return SpatialTransformd();
}

RBDL_DLLAPI void jcalc_X_lambda_S (
    const Model &model,
    ModelData &model_data,
    unsigned int joint_id,
    const VectorNd &q
    ) {
  // exception if we calculate it for the root body
  assert (joint_id > 0);

  if (model.mJoints[joint_id].mJointType == JointTypeRevoluteX) {
    model_data.X_lambda[joint_id] =
      Xrotx (q[model.mJoints[joint_id].q_index]) * model.X_T[joint_id];
    model_data.S[joint_id] = model.mJoints[joint_id].mJointAxes[0];
  } else if (model.mJoints[joint_id].mJointType == JointTypeRevoluteY) {
    model_data.X_lambda[joint_id] =
      Xroty (q[model.mJoints[joint_id].q_index]) * model.X_T[joint_id];
    model_data.S[joint_id] = model.mJoints[joint_id].mJointAxes[0];
  } else if (model.mJoints[joint_id].mJointType == JointTypeRevoluteZ) {
    model_data.X_lambda[joint_id] =
      Xrotz (q[model.mJoints[joint_id].q_index]) * model.X_T[joint_id];
    model_data.S[joint_id] = model.mJoints[joint_id].mJointAxes[0];
  } else if (model.mJoints[joint_id].mDoFCount == 1
      && model.mJoints[joint_id].mJointType != JointTypeCustom){
    model_data.X_lambda[joint_id] =
      jcalc_XJ (model, model_data, joint_id, q) * model.X_T[joint_id];
    // Set the joint axis
    model_data.S[joint_id] = model.mJoints[joint_id].mJointAxes[0];
  } else if (model.mJoints[joint_id].mJointType == JointTypeSpherical) {
    model_data.X_lambda[joint_id] = SpatialTransformd (
        model.GetQuaternion (joint_id, q).toMatrix(),
        Vector3d (0., 0., 0.))
      * model.X_T[joint_id];

    model_data.multdof3_S[joint_id].setZero();

    model_data.multdof3_S[joint_id](0,0) = 1.;
    model_data.multdof3_S[joint_id](1,1) = 1.;
    model_data.multdof3_S[joint_id](2,2) = 1.;
  } else if (model.mJoints[joint_id].mJointType == JointTypeEulerZYX) {
    double q0 = q[model.mJoints[joint_id].q_index];
    double q1 = q[model.mJoints[joint_id].q_index + 1];
    double q2 = q[model.mJoints[joint_id].q_index + 2];

    double s0 = sin (q0);
    double c0 = cos (q0);
    double s1 = sin (q1);
    double c1 = cos (q1);
    double s2 = sin (q2);
    double c2 = cos (q2);

    model_data.X_lambda[joint_id] = SpatialTransformd (
        Matrix3d(
          c0 * c1, s0 * c1, -s1,
          c0 * s1 * s2 - s0 * c2, s0 * s1 * s2 + c0 * c2, c1 * s2,
          c0 * s1 * c2 + s0 * s2, s0 * s1 * c2 - c0 * s2, c1 * c2
          ),
        Vector3d (0., 0., 0.))
      * model.X_T[joint_id];

    model_data.multdof3_S[joint_id].setZero();

    model_data.multdof3_S[joint_id](0,0) = -s1;
    model_data.multdof3_S[joint_id](0,2) = 1.;

    model_data.multdof3_S[joint_id](1,0) = c1 * s2;
    model_data.multdof3_S[joint_id](1,1) = c2;

    model_data.multdof3_S[joint_id](2,0) = c1 * c2;
    model_data.multdof3_S[joint_id](2,1) = - s2;
  } else if (model.mJoints[joint_id].mJointType == JointTypeEulerXYZ) {
    double q0 = q[model.mJoints[joint_id].q_index];
    double q1 = q[model.mJoints[joint_id].q_index + 1];
    double q2 = q[model.mJoints[joint_id].q_index + 2];

    double s0 = sin (q0);
    double c0 = cos (q0);
    double s1 = sin (q1);
    double c1 = cos (q1);
    double s2 = sin (q2);
    double c2 = cos (q2);

    model_data.X_lambda[joint_id] = SpatialTransformd (
        Matrix3d(
          c2 * c1, s2 * c0 + c2 * s1 * s0, s2 * s0 - c2 * s1 * c0,
          -s2 * c1, c2 * c0 - s2 * s1 * s0, c2 * s0 + s2 * s1 * c0,
          s1, -c1 * s0, c1 * c0
          ),
        Vector3d (0., 0., 0.))
      * model.X_T[joint_id];

    model_data.multdof3_S[joint_id].setZero();

    model_data.multdof3_S[joint_id](0,0) = c2 * c1;
    model_data.multdof3_S[joint_id](0,1) = s2;

    model_data.multdof3_S[joint_id](1,0) = -s2 * c1;
    model_data.multdof3_S[joint_id](1,1) = c2;

    model_data.multdof3_S[joint_id](2,0) = s1;
    model_data.multdof3_S[joint_id](2,2) = 1.;
  } else if (model.mJoints[joint_id].mJointType == JointTypeEulerYXZ ) {
    double q0 = q[model.mJoints[joint_id].q_index];
    double q1 = q[model.mJoints[joint_id].q_index + 1];
    double q2 = q[model.mJoints[joint_id].q_index + 2];

    double s0 = sin (q0);
    double c0 = cos (q0);
    double s1 = sin (q1);
    double c1 = cos (q1);
    double s2 = sin (q2);
    double c2 = cos (q2);

    model_data.X_lambda[joint_id] = SpatialTransformd (
        Matrix3d(
          c2 * c0 + s2 * s1 * s0, s2 * c1, -c2 * s0 + s2 * s1 * c0,
          -s2 * c0 + c2 * s1 * s0, c2 * c1, s2 * s0 + c2 * s1 * c0,
          c1 * s0, - s1, c1 * c0
          ),
        Vector3d (0., 0., 0.))
      * model.X_T[joint_id];

    model_data.multdof3_S[joint_id].setZero();

    model_data.multdof3_S[joint_id](0,0) = s2 * c1;
    model_data.multdof3_S[joint_id](0,1) = c2;

    model_data.multdof3_S[joint_id](1,0) = c2 * c1;
    model_data.multdof3_S[joint_id](1,1) = -s2;

    model_data.multdof3_S[joint_id](2,0) = -s1;
    model_data.multdof3_S[joint_id](2,2) = 1.;
  } else if (model.mJoints[joint_id].mJointType == JointTypeTranslationXYZ) {
    double q0 = q[model.mJoints[joint_id].q_index];
    double q1 = q[model.mJoints[joint_id].q_index + 1];
    double q2 = q[model.mJoints[joint_id].q_index + 2];

    model_data.X_lambda[joint_id] = SpatialTransformd (
        Matrix3d::Identity (3,3),
        Vector3d (q0, q1, q2))
      * model.X_T[joint_id];

    model_data.multdof3_S[joint_id].setZero();

    model_data.multdof3_S[joint_id](3,0) = 1.;
    model_data.multdof3_S[joint_id](4,1) = 1.;
    model_data.multdof3_S[joint_id](5,2) = 1.;
  } else if (model.mJoints[joint_id].mJointType == JointTypeCustom) {
    const Joint &joint = model.mJoints[joint_id];
    CustomJoint *custom_joint 
      = model.mCustomJoints[joint.custom_joint_index];

    custom_joint->jcalc_X_lambda_S (model, model_data, joint_id, q);
  } else {
    std::cerr << "Error: invalid joint type!" << std::endl;
    abort();
  }
}
}
