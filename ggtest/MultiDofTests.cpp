#include <eigen_checks/gtest.h>

#include <iostream>

#include "FixturesTests.h"
#include "Human36FixtureGTest.h"
#include "rbdl/rbdl_mathutils.h"
#include "rbdl/Logging.h"

#include "rbdl/Model.h"
#include "rbdl/Kinematics.h"
#include "rbdl/Dynamics.h"
#include "rbdl/Contacts.h"

using namespace std;
using namespace RigidBodyDynamics;
using namespace RigidBodyDynamics::Math;

const double TEST_PREC = 1.0e-12;

class SphericalJoint : public ::testing::Test {
protected:
    virtual void SetUp () {
    ClearLogOutput();

    emulated_model.gravity = Vector3d (0., 0., -9.81);
    multdof3_model.gravity = Vector3d (0., 0., -9.81);
    eulerzyx_model.gravity = Vector3d (0., 0., -9.81);

    body = Body (1., Vector3d (1., 0., 0.), Vector3d (1., 1., 1.));

    joint_rot_zyx = Joint (
        SpatialVectord (0., 0., 1., 0., 0., 0.),
        SpatialVectord (0., 1., 0., 0., 0., 0.),
        SpatialVectord (1., 0., 0., 0., 0., 0.)
        );
    joint_spherical = Joint (JointTypeSpherical);
    joint_eulerzyx = Joint (JointTypeEulerZYX);

    joint_rot_y = Joint (SpatialVectord (0., 1., 0., 0., 0., 0.));

    emulated_model.AppendBody (Xtrans(Vector3d (0., 0., 0.)), joint_rot_y, body);
    emu_body_id = emulated_model.AppendBody (Xtrans (Vector3d (1., 0., 0.)), joint_rot_zyx, body);
    emu_child_id = emulated_model.AppendBody (Xtrans (Vector3d (1., 0., 0.)), joint_rot_y, body);

    multdof3_model.AppendBody (Xtrans(Vector3d (0., 0., 0.)), joint_rot_y, body);
    sph_body_id = multdof3_model.AppendBody (Xtrans (Vector3d (1., 0., 0.)), joint_spherical, body);
    sph_child_id = multdof3_model.AppendBody (Xtrans (Vector3d (1., 0., 0.)), joint_rot_y, body);

    eulerzyx_model.AppendBody (Xtrans(Vector3d (0., 0., 0.)), joint_rot_y, body);
    eulerzyx_body_id = eulerzyx_model.AppendBody (Xtrans (Vector3d (1., 0., 0.)), joint_eulerzyx, body);
    eulerzyx_child_id = eulerzyx_model.AppendBody (Xtrans (Vector3d (1., 0., 0.)), joint_rot_y, body);

    emuQ = VectorNd::Zero ((size_t) emulated_model.q_size);
    emuQDot = VectorNd::Zero ((size_t) emulated_model.qdot_size);
    emuQDDot = VectorNd::Zero ((size_t) emulated_model.qdot_size);
    emuTau = VectorNd::Zero ((size_t) emulated_model.qdot_size);

    sphQ = VectorNd::Zero ((size_t) multdof3_model.q_size);
    sphQDot = VectorNd::Zero ((size_t) multdof3_model.qdot_size);
    sphQDDot = VectorNd::Zero ((size_t) multdof3_model.qdot_size);
    sphTau = VectorNd::Zero ((size_t) multdof3_model.qdot_size);

    eulerzyxQ = VectorNd::Zero ((size_t) eulerzyx_model.q_size);
    eulerzyxQDot = VectorNd::Zero ((size_t) eulerzyx_model.qdot_size);
    eulerzyxQDDot = VectorNd::Zero ((size_t) eulerzyx_model.qdot_size);
    eulerzyxTau = VectorNd::Zero ((size_t) eulerzyx_model.qdot_size);
  }

  virtual void TearDown() {}

  Joint joint_rot_zyx;
  Joint joint_spherical;
  Joint joint_eulerzyx;
  Joint joint_rot_y;
  Body body;

  unsigned int emu_body_id, emu_child_id,
               sph_body_id, sph_child_id,
               eulerzyx_body_id, eulerzyx_child_id;

  Model emulated_model;
  Model multdof3_model;
  Model eulerzyx_model;

  VectorNd emuQ;
  VectorNd emuQDot;
  VectorNd emuQDDot;
  VectorNd emuTau;

  VectorNd sphQ;
  VectorNd sphQDot;
  VectorNd sphQDDot;
  VectorNd sphTau;

  VectorNd eulerzyxQ;
  VectorNd eulerzyxQDot;
  VectorNd eulerzyxQDDot;
  VectorNd eulerzyxTau;
};

void ConvertQAndQDotFromEmulated (
    const Model &emulated_model, const VectorNd &q_emulated, const VectorNd &qdot_emulated,
    const Model &multdof3_model, VectorNd *q_spherical, VectorNd *qdot_spherical) {
  for (unsigned int i = 1; i < multdof3_model.mJoints.size(); i++) {
    unsigned int q_index = multdof3_model.mJoints[i].q_index;

    if (multdof3_model.mJoints[i].mJointType == JointTypeSpherical) {
      Quaterniond quat = Quaterniond::fromZYXAngles ( Vector3d (
            q_emulated[q_index + 0], q_emulated[q_index + 1], q_emulated[q_index + 2]));
      multdof3_model.SetQuaternion (i, quat, (*q_spherical));

      Vector3d omega = angular_velocity_from_angle_rates (
          Vector3d (q_emulated[q_index], q_emulated[q_index + 1], q_emulated[q_index + 2]),
          Vector3d (qdot_emulated[q_index], qdot_emulated[q_index + 1], qdot_emulated[q_index + 2])
          );

      (*qdot_spherical)[q_index] = omega[0];
      (*qdot_spherical)[q_index + 1] = omega[1];
      (*qdot_spherical)[q_index + 2] = omega[2];
    } else {
      (*q_spherical)[q_index] = q_emulated[q_index];
      (*qdot_spherical)[q_index] = qdot_emulated[q_index];
    }
  }
}

TEST(MultiDofTests, TestQuaternionIntegration ) {
  double timestep = 0.001;

  Vector3d zyx_angles_t0 (0.1, 0.2, 0.3);
  Vector3d zyx_rates (3., 5., 2.);
  Vector3d zyx_angles_t1 = zyx_angles_t0 + timestep * zyx_rates;
  Quaterniond q_zyx_t1 = Quaterniond::fromZYXAngles (zyx_angles_t1);

  Quaterniond q_t0 = Quaterniond::fromZYXAngles (zyx_angles_t0);
  Vector3d w_base = global_angular_velocity_from_rates (zyx_angles_t0, zyx_rates);
  Quaterniond q_t1 = q_t0.timeStep (w_base, timestep);

  // Note: we test with a rather crude precision. My guess for the error is
  // that we compare two different things:
  //   A) integration under the assumption that the euler rates are
  //   constant
  //   B) integration under the assumption that the angular velocity is
  //   constant
  // However I am not entirely sure about this...
  EXPECT_TRUE(EIGEN_MATRIX_NEAR (q_zyx_t1, q_t1, 1.0e-5));
}

TEST_F (SphericalJoint, TestQIndices) {
  EXPECT_EQ (0u, multdof3_model.mJoints[1].q_index);
  EXPECT_EQ (1u, multdof3_model.mJoints[2].q_index);
  EXPECT_EQ (4u, multdof3_model.mJoints[3].q_index);

  EXPECT_EQ (5u, emulated_model.q_size);
  EXPECT_EQ (5u, emulated_model.qdot_size);

  EXPECT_EQ (6u, multdof3_model.q_size);
  EXPECT_EQ (5u, multdof3_model.qdot_size);
  EXPECT_EQ (5u, multdof3_model.multdof3_w_index[2]);
}

TEST_F (SphericalJoint, TestGetQuaternion) {
  multdof3_model.AppendBody (Xtrans (Vector3d (1., 0., 0.)), joint_spherical, body);

  sphQ = VectorNd::Zero ((size_t) multdof3_model.q_size);
  sphQDot = VectorNd::Zero ((size_t) multdof3_model.qdot_size);
  sphQDDot = VectorNd::Zero ((size_t) multdof3_model.qdot_size);
  sphTau = VectorNd::Zero ((size_t) multdof3_model.qdot_size);

  EXPECT_EQ (10u, multdof3_model.q_size);
  EXPECT_EQ (8u, multdof3_model.qdot_size);

  EXPECT_EQ (0u, multdof3_model.mJoints[1].q_index);
  EXPECT_EQ (1u, multdof3_model.mJoints[2].q_index);
  EXPECT_EQ (4u, multdof3_model.mJoints[3].q_index);
  EXPECT_EQ (5u, multdof3_model.mJoints[4].q_index);

  EXPECT_EQ (8u, multdof3_model.multdof3_w_index[2]);
  EXPECT_EQ (9u, multdof3_model.multdof3_w_index[4]);

  sphQ[0] = 100.;
  sphQ[1] = 0.;
  sphQ[2] = 1.;
  sphQ[3] = 2.;
  sphQ[4] = 100.;
  sphQ[5] = -6.;
  sphQ[6] = -7.;
  sphQ[7] = -8;
  sphQ[8] = 4.;
  sphQ[9] = -9.;

  Quaterniond reference_1 (0., 1., 2., 4.);
  Quaterniond quat_1 = multdof3_model.GetQuaternion (2, sphQ);
  EXPECT_TRUE(EIGEN_MATRIX_EQUAL_DOUBLE (reference_1, quat_1));

  Quaterniond reference_3 (-6., -7., -8., -9.);
  Quaterniond quat_3 = multdof3_model.GetQuaternion (4, sphQ);
  EXPECT_TRUE(EIGEN_MATRIX_EQUAL_DOUBLE (reference_3, quat_3));
}

TEST_F (SphericalJoint, TestSetQuaternion) {
  multdof3_model.AppendBody (Xtrans (Vector3d (1., 0., 0.)), joint_spherical, body);

  sphQ = VectorNd::Zero ((size_t) multdof3_model.q_size);
  sphQDot = VectorNd::Zero ((size_t) multdof3_model.qdot_size);
  sphQDDot = VectorNd::Zero ((size_t) multdof3_model.qdot_size);
  sphTau = VectorNd::Zero ((size_t) multdof3_model.qdot_size);

  Quaterniond reference_1 (0., 1., 2., 3.);
  multdof3_model.SetQuaternion (2, reference_1, sphQ);
  Quaterniond test = multdof3_model.GetQuaternion (2, sphQ);
  EXPECT_TRUE(EIGEN_MATRIX_EQUAL_DOUBLE (reference_1, test));

  Quaterniond reference_2 (11., 22., 33., 44.);
  multdof3_model.SetQuaternion (4, reference_2, sphQ);
  test = multdof3_model.GetQuaternion (4, sphQ);
  EXPECT_TRUE(EIGEN_MATRIX_EQUAL_DOUBLE (reference_2, test));
}

TEST_F (SphericalJoint, TestOrientation) {
  emuQ[0] = 1.1;
  emuQ[1] = 1.1;
  emuQ[2] = 1.1;
  emuQ[3] = 1.1;

  for (unsigned int i = 0; i < emuQ.size(); i++) {
    sphQ[i] = emuQ[i];
  }

  Quaterniond quat =  Quaterniond::fromAxisAngle (Vector3d (1., 0., 0.), emuQ[2])
    * Quaterniond::fromAxisAngle (Vector3d (0., 1., 0.), emuQ[1])
    * Quaterniond::fromAxisAngle (Vector3d (0., 0., 1.), emuQ[0]);
  multdof3_model.SetQuaternion (2, quat, sphQ);

  Matrix3d emu_orientation = CalcBodyWorldOrientation (emulated_model, emuQ, emu_child_id);
  Matrix3d sph_orientation = CalcBodyWorldOrientation (multdof3_model, sphQ, sph_child_id);

  EXPECT_TRUE(EIGEN_MATRIX_NEAR (emu_orientation, sph_orientation, TEST_PREC));
}

TEST_F (SphericalJoint, TestUpdateKinematics) {
  emuQ[0] = 1.;
  emuQ[1] = 1.;
  emuQ[2] = 1.;
  emuQ[3] = 1.;
  emuQ[4] = 1.;

  emuQDot[0] = 1.;
  emuQDot[1] = 1.;
  emuQDot[2] = 1.;
  emuQDot[3] = 1.;
  emuQDot[4] = 1.;

  emuQDDot[0] = 1.;
  emuQDDot[1] = 1.;
  emuQDDot[2] = 1.;
  emuQDDot[3] = 1.;
  emuQDDot[4] = 1.;

  ConvertQAndQDotFromEmulated (emulated_model, emuQ, emuQDot, multdof3_model, &sphQ, &sphQDot);
  ConvertQAndQDotFromEmulated (emulated_model, emuQ, emuQDDot, multdof3_model, &sphQ, &sphQDDot);

  Vector3d a = angular_acceleration_from_angle_rates (
      Vector3d (emuQ[3], emuQ[2], emuQ[1]),
      Vector3d (emuQDot[3], emuQDot[2], emuQDot[1]),
      Vector3d (emuQDDot[3], emuQDDot[2], emuQDDot[1])
      );

  sphQDDot[0] = emuQDDot[0];
  sphQDDot[1] = a[0];
  sphQDDot[2] = a[1];
  sphQDDot[3] = a[2];
  sphQDDot[4] = emuQDDot[4];

  UpdateKinematicsCustom (emulated_model, &emuQ, &emuQDot, &emuQDDot);
  UpdateKinematicsCustom (multdof3_model, &sphQ, &sphQDot, &sphQDDot);

  EXPECT_TRUE(EIGEN_MATRIX_NEAR (emulated_model_data.v[emu_body_id], multdof3_model_data.v[sph_body_id], TEST_PREC));
  EXPECT_TRUE(EIGEN_MATRIX_NEAR (emulated_model_data.a[emu_body_id], multdof3_model_data.a[sph_body_id],  TEST_PREC));

  UpdateKinematics (multdof3_model, sphQ, sphQDot, sphQDDot);

  EXPECT_TRUE(EIGEN_MATRIX_NEAR (emulated_model_data.v[emu_child_id], multdof3_model_data.v[sph_child_id],  TEST_PREC));
  EXPECT_TRUE(EIGEN_MATRIX_NEAR (emulated_model_data.a[emu_child_id], multdof3_model_data.a[sph_child_id],  TEST_PREC));
}

TEST_F (SphericalJoint, TestSpatialVelocities) {
  emuQ[0] = 1.;
  emuQ[1] = 2.;
  emuQ[2] = 3.;
  emuQ[3] = 4.;

  emuQDot[0] = 4.;
  emuQDot[1] = 2.;
  emuQDot[2] = 3.;
  emuQDot[3] = 6.;

  ConvertQAndQDotFromEmulated (emulated_model, emuQ, emuQDot, multdof3_model, &sphQ, &sphQDot);

  UpdateKinematicsCustom (emulated_model, &emuQ, &emuQDot, NULL);
  UpdateKinematicsCustom (multdof3_model, &sphQ, &sphQDot, NULL);

  EXPECT_TRUE(EIGEN_MATRIX_NEAR (emulated_model_data.v[emu_child_id], multdof3_model_data.v[sph_child_id],  TEST_PREC));
}

TEST_F (SphericalJoint, TestForwardDynamicsQAndQDot) {
  emuQ[0] = 1.1;
  emuQ[1] = 1.2;
  emuQ[2] = 1.3;
  emuQ[3] = 1.4;

  emuQDot[0] = 2.2;
  emuQDot[1] = 2.3;
  emuQDot[2] = 2.4;
  emuQDot[3] = 2.5;

  ConvertQAndQDotFromEmulated (emulated_model, emuQ, emuQDot, multdof3_model, &sphQ, &sphQDot);

  ForwardDynamics (emulated_model, emuQ, emuQDot, emuTau, emuQDDot);
  ForwardDynamics (multdof3_model, sphQ, sphQDot, sphTau, sphQDDot);

  EXPECT_TRUE(EIGEN_MATRIX_NEAR (emulated_model_data.a[emu_child_id], multdof3_model_data.a[sph_child_id], TEST_PREC));
}

TEST_F (SphericalJoint, TestDynamicsConsistencyRNEA_ABA ) {
  emuQ[0] = 1.1;
  emuQ[1] = 1.2;
  emuQ[2] = 1.3;
  emuQ[3] = 1.4;
  emuQ[4] = 1.5;

  emuQDot[0] = 1.;
  emuQDot[1] = 2.;
  emuQDot[2] = 3.;
  emuQDot[3] = 4.;
  emuQDot[4] = 5.;

  sphTau[0] = 5.;
  sphTau[1] = 4.;
  sphTau[2] = 7.;
  sphTau[3] = 3.;
  sphTau[4] = 2.;

  ConvertQAndQDotFromEmulated (emulated_model, emuQ, emuQDot, multdof3_model, &sphQ, &sphQDot);

  ForwardDynamics (multdof3_model, sphQ, sphQDot, sphTau, sphQDDot);

  VectorNd tau_id (VectorNd::Zero (multdof3_model.qdot_size));
  InverseDynamics (multdof3_model, sphQ, sphQDot, sphQDDot, tau_id);

  EXPECT_TRUE(EIGEN_MATRIX_NEAR (sphTau, tau_id,  TEST_PREC));
}

TEST_F (SphericalJoint, TestCRBA ) {
  emuQ[0] = 1.1;
  emuQ[1] = 1.2;
  emuQ[2] = 1.3;
  emuQ[3] = 1.4;
  emuQ[4] = 1.5;

  emuQDot[0] = 1.;
  emuQDot[1] = 2.;
  emuQDot[2] = 3.;
  emuQDot[3] = 4.;
  emuQDot[4] = 5.;

  sphTau[0] = 5.;
  sphTau[1] = 4.;
  sphTau[2] = 7.;
  sphTau[3] = 3.;
  sphTau[4] = 2.;

  ConvertQAndQDotFromEmulated (emulated_model, emuQ, emuQDot, multdof3_model, &sphQ, &sphQDot);

  MatrixNd H_crba (MatrixNd::Zero (multdof3_model.qdot_size, multdof3_model.qdot_size));

  UpdateKinematicsCustom (multdof3_model, &sphQ, NULL, NULL);
  CompositeRigidBodyAlgorithm (multdof3_model, sphQ, H_crba, false);

  MatrixNd H_id (MatrixNd::Zero (multdof3_model.qdot_size, multdof3_model.qdot_size));
  VectorNd H_col = VectorNd::Zero (multdof3_model.qdot_size);
  VectorNd QDDot_zero = VectorNd::Zero (multdof3_model.qdot_size);

  for (unsigned int i = 0; i < multdof3_model.qdot_size; i++) {
    // compute each column
    VectorNd delta_a = VectorNd::Zero (multdof3_model.qdot_size);
    delta_a[i] = 1.;

    // compute ID (model, q, qdot, delta_a)
    VectorNd id_delta = VectorNd::Zero (multdof3_model.qdot_size);
    InverseDynamics (multdof3_model, sphQ, sphQDot, delta_a, id_delta);

    // compute ID (model, q, qdot, zero)
    VectorNd id_zero = VectorNd::Zero (multdof3_model.qdot_size);
    InverseDynamics (multdof3_model, sphQ, sphQDot, QDDot_zero, id_zero);

    H_col = id_delta - id_zero;
    H_id.block(0, i, multdof3_model.qdot_size, 1) = H_col;
  }

  EXPECT_TRUE(EIGEN_MATRIX_NEAR (H_id, H_crba, TEST_PREC));
}

TEST_F (SphericalJoint, TestForwardDynamicsLagrangianVsABA ) {
  emuQ[0] = 1.1;
  emuQ[1] = 1.2;
  emuQ[2] = 1.3;
  emuQ[3] = 1.4;
  emuQ[4] = 1.5;

  emuQDot[0] = 1.;
  emuQDot[1] = 2.;
  emuQDot[2] = 3.;
  emuQDot[3] = 4.;
  emuQDot[4] = 5.;

  sphTau[0] = 5.;
  sphTau[1] = 4.;
  sphTau[2] = 7.;
  sphTau[3] = 3.;
  sphTau[4] = 2.;

  ConvertQAndQDotFromEmulated (emulated_model, emuQ, emuQDot, multdof3_model, &sphQ, &sphQDot);

  VectorNd QDDot_aba = VectorNd::Zero (multdof3_model.qdot_size);
  VectorNd QDDot_lag = VectorNd::Zero (multdof3_model.qdot_size);

  ForwardDynamicsLagrangian (multdof3_model, sphQ, sphQDot, sphTau, QDDot_lag);
  ForwardDynamics (multdof3_model, sphQ, sphQDot, sphTau, QDDot_aba);

  EXPECT_TRUE(EIGEN_MATRIX_NEAR (QDDot_lag, QDDot_aba, TEST_PREC));
}

TEST_F (SphericalJoint, TestContactsLagrangian) {
  ConstraintSet constraint_set_emu;

  constraint_set_emu.AddConstraint (emu_child_id, Vector3d (0., 0., -1.), Vector3d (1., 0., 0.));
  constraint_set_emu.AddConstraint (emu_child_id, Vector3d (0., 0., -1.), Vector3d (0., 1., 0.));
  constraint_set_emu.AddConstraint (emu_child_id, Vector3d (0., 0., -1.), Vector3d (0., 0., 1.));

  constraint_set_emu.Bind(emulated_model);

  ConstraintSet constraint_set_sph;

  constraint_set_sph.AddConstraint (sph_child_id, Vector3d (0., 0., -1.), Vector3d (1., 0., 0.));
  constraint_set_sph.AddConstraint (sph_child_id, Vector3d (0., 0., -1.), Vector3d (0., 1., 0.));
  constraint_set_sph.AddConstraint (sph_child_id, Vector3d (0., 0., -1.), Vector3d (0., 0., 1.));

  constraint_set_sph.Bind(multdof3_model);

  ForwardDynamicsContactsDirect (emulated_model, emuQ, emuQDot, emuTau, constraint_set_emu, emuQDDot);
  VectorNd emu_force_lagrangian = constraint_set_emu.force;
  ForwardDynamicsContactsDirect (multdof3_model, sphQ, sphQDot, sphTau, constraint_set_sph, sphQDDot);
  VectorNd sph_force_lagrangian = constraint_set_sph.force;

  EXPECT_TRUE(EIGEN_MATRIX_NEAR (emu_force_lagrangian, sph_force_lagrangian, TEST_PREC));
}

TEST_F (SphericalJoint, TestContacts) {
  ConstraintSet constraint_set_emu;

  constraint_set_emu.AddConstraint (emu_child_id, Vector3d (0., 0., -1.), Vector3d (1., 0., 0.));
  constraint_set_emu.AddConstraint (emu_child_id, Vector3d (0., 0., -1.), Vector3d (0., 1., 0.));
  constraint_set_emu.AddConstraint (emu_child_id, Vector3d (0., 0., -1.), Vector3d (0., 0., 1.));

  constraint_set_emu.Bind(emulated_model);

  ConstraintSet constraint_set_sph;

  constraint_set_sph.AddConstraint (sph_child_id, Vector3d (0., 0., -1.), Vector3d (1., 0., 0.));
  constraint_set_sph.AddConstraint (sph_child_id, Vector3d (0., 0., -1.), Vector3d (0., 1., 0.));
  constraint_set_sph.AddConstraint (sph_child_id, Vector3d (0., 0., -1.), Vector3d (0., 0., 1.));

  constraint_set_sph.Bind(multdof3_model);

  ForwardDynamicsContactsKokkevis (emulated_model, emuQ, emuQDot, emuTau, constraint_set_emu, emuQDDot);
  VectorNd emu_force_kokkevis = constraint_set_emu.force;
  ForwardDynamicsContactsKokkevis (multdof3_model, sphQ, sphQDot, sphTau, constraint_set_sph, sphQDDot);
  VectorNd sph_force_kokkevis = constraint_set_sph.force;

  EXPECT_TRUE(EIGEN_MATRIX_NEAR (emu_force_kokkevis, sph_force_kokkevis, TEST_PREC));
}

TEST_F (SphericalJoint, TestEulerZYXvsEmulatedLagrangian ) {
  emuQ[0] = 1.1;
  emuQ[1] = 1.2;
  emuQ[2] = 1.3;
  emuQ[3] = 1.4;
  emuQ[4] = 1.5;

  emuQDot[0] = 1.;
  emuQDot[1] = 2.;
  emuQDot[2] = 3.;
  emuQDot[3] = 4.;
  emuQDot[4] = 5.;

  emuTau[0] = 5.;
  emuTau[1] = 4.;
  emuTau[2] = 7.;
  emuTau[3] = 3.;
  emuTau[4] = 2.;

  VectorNd QDDot_emu = VectorNd::Zero (emulated_model.qdot_size);
  VectorNd QDDot_eulerzyx = VectorNd::Zero (eulerzyx_model.qdot_size);

  ForwardDynamicsLagrangian (emulated_model, emuQ, emuQDot, emuTau, QDDot_emu);
  ForwardDynamicsLagrangian (eulerzyx_model, emuQ, emuQDot, emuTau, QDDot_eulerzyx);

  EXPECT_TRUE(EIGEN_MATRIX_NEAR (QDDot_emu, QDDot_eulerzyx,  TEST_PREC));
}

TEST_F (SphericalJoint, TestEulerZYXvsEmulatedArticulatedBodyAlgorithm ) {
  emuQ[0] = 1.1;
  emuQ[1] = 1.2;
  emuQ[2] = 1.3;
  emuQ[3] = 1.4;
  emuQ[4] = 1.5;

  emuQDot[0] = 1.;
  emuQDot[1] = 2.;
  emuQDot[2] = 3.;
  emuQDot[3] = 4.;
  emuQDot[4] = 5.;

  emuTau[0] = 5.;
  emuTau[1] = 4.;
  emuTau[2] = 7.;
  emuTau[3] = 3.;
  emuTau[4] = 2.;

  VectorNd QDDot_emu = VectorNd::Zero (emulated_model.qdot_size);
  VectorNd QDDot_eulerzyx = VectorNd::Zero (eulerzyx_model.qdot_size);

  ForwardDynamics (emulated_model, emuQ, emuQDot, emuTau, QDDot_emu);
  ForwardDynamics (eulerzyx_model, emuQ, emuQDot, emuTau, QDDot_eulerzyx);

  EXPECT_TRUE(EIGEN_MATRIX_NEAR (QDDot_emu, QDDot_eulerzyx,  TEST_PREC));
}

TEST_F (SphericalJoint, TestEulerZYXvsEmulatedContacts ) {
  emuQ[0] = 1.1;
  emuQ[1] = 1.2;
  emuQ[2] = 1.3;
  emuQ[3] = 1.4;
  emuQ[4] = 1.5;

  emuQDot[0] = 1.;
  emuQDot[1] = 2.;
  emuQDot[2] = 3.;
  emuQDot[3] = 4.;
  emuQDot[4] = 5.;

  emuTau[0] = 5.;
  emuTau[1] = 4.;
  emuTau[2] = 7.;
  emuTau[3] = 3.;
  emuTau[4] = 2.;

  VectorNd QDDot_emu = VectorNd::Zero (emulated_model.qdot_size);
  VectorNd QDDot_eulerzyx = VectorNd::Zero (eulerzyx_model.qdot_size);

  ConstraintSet CS_euler;
  CS_euler.AddConstraint (eulerzyx_child_id, Vector3d (1., 1., 1.), Vector3d (1., 0., 0.));
  CS_euler.AddConstraint (eulerzyx_child_id, Vector3d (0., 0., 0.), Vector3d (0., 1., 0.));
  CS_euler.AddConstraint (eulerzyx_child_id, Vector3d (0., 0., 0.), Vector3d (0., 0., 1.));
  CS_euler.Bind (eulerzyx_model);

  ConstraintSet CS_emulated;
  CS_emulated.AddConstraint (emu_child_id, Vector3d (1., 1., 1.), Vector3d (1., 0., 0.));
  CS_emulated.AddConstraint (emu_child_id, Vector3d (0., 0., 0.), Vector3d (0., 1., 0.));
  CS_emulated.AddConstraint (emu_child_id, Vector3d (0., 0., 0.), Vector3d (0., 0., 1.));
  CS_emulated.Bind (emulated_model);

  ForwardDynamicsContactsDirect (emulated_model, emuQ, emuQDot, emuTau, CS_emulated, QDDot_emu);
  ForwardDynamicsContactsDirect (eulerzyx_model, emuQ, emuQDot, emuTau, CS_euler, QDDot_eulerzyx);

  EXPECT_TRUE(EIGEN_MATRIX_NEAR (QDDot_emu, QDDot_eulerzyx,  TEST_PREC));

  ClearLogOutput();

  ForwardDynamicsContactsKokkevis (emulated_model, emuQ, emuQDot, emuTau, CS_emulated, QDDot_emu);
  ForwardDynamicsContactsKokkevis (eulerzyx_model, emuQ, emuQDot, emuTau, CS_euler, QDDot_eulerzyx);

  EXPECT_TRUE(EIGEN_MATRIX_NEAR (QDDot_emu, QDDot_eulerzyx, TEST_PREC * QDDot_emu.norm()));

  ForwardDynamicsContactsKokkevis (emulated_model, emuQ, emuQDot, emuTau, CS_emulated, QDDot_emu);
  ForwardDynamicsContactsDirect (eulerzyx_model, emuQ, emuQDot, emuTau, CS_euler, QDDot_eulerzyx);

  EXPECT_TRUE(EIGEN_MATRIX_NEAR (QDDot_emu, QDDot_eulerzyx, TEST_PREC * QDDot_emu.norm()));
}

TEST_F (SphericalJoint, TestEulerZYXvsEmulatedCRBA ) {
  emuQ[0] = 1.1;
  emuQ[1] = 1.2;
  emuQ[2] = 1.3;
  emuQ[3] = 1.4;
  emuQ[4] = 1.5;

  MatrixNd H_emulated (MatrixNd::Zero (emulated_model.q_size, emulated_model.q_size));
  MatrixNd H_eulerzyx (MatrixNd::Zero (eulerzyx_model.q_size, eulerzyx_model.q_size));

  CompositeRigidBodyAlgorithm (emulated_model, emuQ, H_emulated);
  CompositeRigidBodyAlgorithm (eulerzyx_model, emuQ, H_eulerzyx);

  EXPECT_TRUE(EIGEN_MATRIX_NEAR (H_emulated, H_eulerzyx, TEST_PREC));
}

TEST (MultiDofTests, TestJointTypeTranslationZYX ) {
  Model model_emulated;
  Model model_3dof;

  Body body (1., Vector3d (1., 2., 1.), Matrix3d (1., 0., 0, 0., 1., 0., 0., 0., 1.));
  Joint joint_emulated (
      SpatialVectord (0., 0., 0., 1., 0., 0.),
      SpatialVectord (0., 0., 0., 0., 1., 0.),
      SpatialVectord (0., 0., 0., 0., 0., 1.)
      );
  Joint joint_3dof (JointTypeTranslationXYZ);

  model_emulated.AppendBody (SpatialTransformd (), joint_emulated, body);
  model_3dof.AppendBody (SpatialTransformd (), joint_3dof, body);

  VectorNd q (VectorNd::Zero (model_emulated.q_size));
  VectorNd qdot (VectorNd::Zero (model_emulated.qdot_size));
  VectorNd qddot_emulated (VectorNd::Zero (model_emulated.qdot_size));
  VectorNd qddot_3dof (VectorNd::Zero (model_emulated.qdot_size));
  VectorNd tau (VectorNd::Zero (model_emulated.qdot_size));

  for (int i = 0; i < q.size(); i++) {
    q[i] = 1.1 * (static_cast<double>(i + 1));
    qdot[i] = 0.1 * (static_cast<double>(i + 1));
    qddot_3dof[i] = 0.21 * (static_cast<double>(i + 1));
    tau[i] = 2.1 * (static_cast<double>(i + 1));
  }

  qddot_emulated = qddot_3dof;
  VectorNd id_emulated (tau);
  VectorNd id_3dof(tau);
  InverseDynamics (model_emulated, q, qdot, qddot_emulated, id_emulated);
  InverseDynamics (model_3dof, q, qdot, qddot_emulated, id_3dof);
  EXPECT_TRUE(EIGEN_MATRIX_NEAR (id_emulated, id_3dof, TEST_PREC * id_emulated.norm()));

  ForwardDynamicsLagrangian (model_emulated, q, qdot, tau, qddot_emulated);
  ForwardDynamicsLagrangian (model_3dof, q, qdot, tau, qddot_3dof);

  EXPECT_TRUE(EIGEN_MATRIX_EQUAL_DOUBLE (qddot_emulated, qddot_3dof));

  MatrixNd H_emulated (MatrixNd::Zero (q.size(), q.size()));
  MatrixNd H_3dof (MatrixNd::Zero (q.size(), q.size()));

  CompositeRigidBodyAlgorithm (model_emulated, q, H_emulated);
  CompositeRigidBodyAlgorithm (model_3dof, q, H_3dof);

  EXPECT_TRUE(EIGEN_MATRIX_NEAR (H_emulated, H_3dof, TEST_PREC));
}

TEST (MultiDofTests, TestJointTypeEulerXYZ ) {
  Model model_emulated;
  Model model_3dof;

  Body body (1., Vector3d (1., 2., 1.), Matrix3d (1., 0., 0, 0., 1., 0., 0., 0., 1.));
  Joint joint_emulated (
      SpatialVectord (1., 0., 0., 0., 0., 0.),
      SpatialVectord (0., 1., 0., 0., 0., 0.),
      SpatialVectord (0., 0., 1., 0., 0., 0.)
      );
  Joint joint_3dof (JointTypeEulerXYZ);

  model_emulated.AppendBody (SpatialTransformd (), joint_emulated, body);
  model_3dof.AppendBody (SpatialTransformd (), joint_3dof, body);

  VectorNd q (VectorNd::Zero (model_emulated.q_size));
  VectorNd qdot (VectorNd::Zero (model_emulated.qdot_size));
  VectorNd qddot_emulated (VectorNd::Zero (model_emulated.qdot_size));
  VectorNd qddot_3dof (VectorNd::Zero (model_emulated.qdot_size));
  VectorNd tau (VectorNd::Zero (model_emulated.qdot_size));

  for (int i = 0; i < q.size(); i++) {
    q[i] = 1.1 * (static_cast<double>(i + 1));
    qdot[i] = 0.55* (static_cast<double>(i + 1));
    qddot_emulated[i] = 0.23 * (static_cast<double>(i + 1));
    qddot_3dof[i] = 0.22 * (static_cast<double>(i + 1));
    tau[i] = 2.1 * (static_cast<double>(i + 1));
  }

  UpdateKinematicsCustom (model_emulated, &q, &qdot, &qddot_emulated);
  UpdateKinematicsCustom (model_3dof, &q, &qdot, &qddot_emulated);

  EXPECT_TRUE(EIGEN_MATRIX_EQUAL_DOUBLE (model_emulated.X_base[3].E, model_3dof.X_base[1].E));
  EXPECT_TRUE(EIGEN_MATRIX_EQUAL_DOUBLE (model_emulated.model_data.v[3], model_3dof.model_data.v[1]));

  ForwardDynamicsLagrangian (model_emulated, q, qdot, tau, qddot_emulated);
  ForwardDynamicsLagrangian (model_3dof, q, qdot, tau, qddot_3dof);

  EXPECT_TRUE(EIGEN_MATRIX_NEAR (qddot_emulated, qddot_3dof,  TEST_PREC));

  MatrixNd H_emulated (MatrixNd::Zero (q.size(), q.size()));
  MatrixNd H_3dof (MatrixNd::Zero (q.size(), q.size()));

  CompositeRigidBodyAlgorithm (model_emulated, q, H_emulated);
  CompositeRigidBodyAlgorithm (model_3dof, q, H_3dof);

  EXPECT_TRUE(EIGEN_MATRIX_NEAR (H_emulated, H_3dof, TEST_PREC));
}

TEST (MultiDofTests, TestJointTypeEulerYXZ ) {
  Model model_emulated;
  Model model_3dof;

  Body body (1., Vector3d (1., 2., 1.), Matrix3d (1., 0., 0, 0., 1., 0., 0., 0., 1.));
  Joint joint_emulated (
      SpatialVectord (0., 1., 0., 0., 0., 0.),
      SpatialVectord (1., 0., 0., 0., 0., 0.),
      SpatialVectord (0., 0., 1., 0., 0., 0.)
      );
  Joint joint_3dof (JointTypeEulerYXZ);

  model_emulated.AppendBody (SpatialTransformd (), joint_emulated, body);
  model_3dof.AppendBody (SpatialTransformd (), joint_3dof, body);

  VectorNd q (VectorNd::Zero (model_emulated.q_size));
  VectorNd qdot (VectorNd::Zero (model_emulated.qdot_size));
  VectorNd qddot_emulated (VectorNd::Zero (model_emulated.qdot_size));
  VectorNd qddot_3dof (VectorNd::Zero (model_emulated.qdot_size));
  VectorNd tau (VectorNd::Zero (model_emulated.qdot_size));

  for (int i = 0; i < q.size(); i++) {
    q[i] = 0.5 * M_PI * static_cast<double>(rand()) / static_cast<double>(RAND_MAX);
    qdot[i] = 0.5 * M_PI * static_cast<double>(rand()) / static_cast<double>(RAND_MAX);
    tau[i] = 0.5 * M_PI * static_cast<double>(rand()) / static_cast<double>(RAND_MAX);
    qddot_3dof[i] = 0.5 * M_PI * static_cast<double>(rand()) / static_cast<double>(RAND_MAX);
  }

  qddot_emulated = qddot_3dof;

  UpdateKinematicsCustom (model_emulated, &q, &qdot, &qddot_emulated);
  UpdateKinematicsCustom (model_3dof, &q, &qdot, &qddot_emulated);

  EXPECT_TRUE(EIGEN_MATRIX_NEAR (model_emulated.X_base[3].E, model_3dof.X_base[1].E,  TEST_PREC));
  EXPECT_TRUE(EIGEN_MATRIX_NEAR (model_emulated.model_data.v[3], model_3dof.model_data.v[1], TEST_PREC));

  VectorNd id_emulated (tau);
  VectorNd id_3dof(tau);
  InverseDynamics (model_emulated, q, qdot, qddot_emulated, id_emulated);
  InverseDynamics (model_3dof, q, qdot, qddot_emulated, id_3dof);

  EXPECT_TRUE(EIGEN_MATRIX_NEAR (id_emulated, id_3dof, TEST_PREC));

  ForwardDynamicsLagrangian (model_emulated, q, qdot, tau, qddot_emulated);
  ForwardDynamicsLagrangian (model_3dof, q, qdot, tau, qddot_3dof);

  EXPECT_TRUE(EIGEN_MATRIX_NEAR (qddot_emulated, qddot_3dof, TEST_PREC));

  MatrixNd H_emulated (MatrixNd::Zero (q.size(), q.size()));
  MatrixNd H_3dof (MatrixNd::Zero (q.size(), q.size()));

  CompositeRigidBodyAlgorithm (model_emulated, q, H_emulated);
  CompositeRigidBodyAlgorithm (model_3dof, q, H_3dof);

  EXPECT_TRUE(EIGEN_MATRIX_NEAR (H_emulated, H_3dof,TEST_PREC));
}

TEST_F  (Human36, TestUpdateKinematics) {
  for (unsigned int i = 0; i < q.size(); i++) {
    q[i] = 0.5 * M_PI * static_cast<double>(rand()) / static_cast<double>(RAND_MAX);
    qdot[i] = 0.5 * M_PI * static_cast<double>(rand()) / static_cast<double>(RAND_MAX);
    qddot[i] = 0.5 * M_PI * static_cast<double>(rand()) / static_cast<double>(RAND_MAX);
    tau[i] = 0.5 * M_PI * static_cast<double>(rand()) / static_cast<double>(RAND_MAX);
  }

  VectorNd id_emulated (tau);
  VectorNd id_3dof (tau);

  UpdateKinematics (*model_emulated, q, qdot, qddot);
  UpdateKinematics (*model_3dof, q, qdot, qddot);

  EXPECT_TRUE(EIGEN_MATRIX_NEAR (model_emulated->model_data.v[body_id_emulated[BodyPelvis]], model_3dof->model_data.v[body_id_3dof[BodyPelvis]], TEST_PREC));
  EXPECT_TRUE(EIGEN_MATRIX_NEAR (model_emulated->model_data.v[body_id_emulated[BodyThighRight]], model_3dof->model_data.v[body_id_3dof[BodyThighRight]],  TEST_PREC));
  EXPECT_TRUE(EIGEN_MATRIX_NEAR (model_emulated->model_data.v[body_id_emulated[BodyShankRight]], model_3dof->model_data.v[body_id_3dof[BodyShankRight]], TEST_PREC));
}

TEST_F  (Human36, TestInverseDynamics) {
  for (unsigned int i = 0; i < q.size(); i++) {
    q[i] = 0.5 * M_PI * static_cast<double>(rand()) / static_cast<double>(RAND_MAX);
    qdot[i] = 0.5 * M_PI * static_cast<double>(rand()) / static_cast<double>(RAND_MAX);
    qddot[i] = 0.5 * M_PI * static_cast<double>(rand()) / static_cast<double>(RAND_MAX);
    tau[i] = 0.5 * M_PI * static_cast<double>(rand()) / static_cast<double>(RAND_MAX);
  }

  VectorNd id_emulated (tau);
  VectorNd id_3dof (tau);

  ClearLogOutput();
  InverseDynamics (*model_emulated, q, qdot, qddot, id_emulated);
  InverseDynamics (*model_3dof, q, qdot, qddot, id_3dof);

  EXPECT_TRUE(EIGEN_MATRIX_NEAR (id_emulated, id_3dof, TEST_PREC * id_emulated.norm()));
}

TEST_F  (Human36, TestNonlinearEffects) {
  for (unsigned int i = 0; i < q.size(); i++) {
    q[i] = 0.5 * M_PI * static_cast<double>(rand()) / static_cast<double>(RAND_MAX);
    qdot[i] = 0.5 * M_PI * static_cast<double>(rand()) / static_cast<double>(RAND_MAX);
    qddot[i] = 0.5 * M_PI * static_cast<double>(rand()) / static_cast<double>(RAND_MAX);
    tau[i] = 0.5 * M_PI * static_cast<double>(rand()) / static_cast<double>(RAND_MAX);
  }

  VectorNd nle_emulated (tau);
  VectorNd nle_3dof (tau);

  ClearLogOutput();
  NonlinearEffects (*model_emulated, q, qdot, nle_emulated);
  NonlinearEffects (*model_3dof, q, qdot, nle_3dof);

  EXPECT_TRUE(EIGEN_MATRIX_NEAR (nle_emulated, nle_3dof, TEST_PREC * nle_emulated.norm()));
}

TEST_F  (Human36, TestContactsEmulatedLagrangianKokkevis) {
  for (unsigned int i = 0; i < q.size(); i++) {
    q[i] = 0.5 * M_PI * static_cast<double>(rand()) / static_cast<double>(RAND_MAX);
    qdot[i] = 0.5 * M_PI * static_cast<double>(rand()) / static_cast<double>(RAND_MAX);
    tau[i] = 0.5 * M_PI * static_cast<double>(rand()) / static_cast<double>(RAND_MAX);
  }

  VectorNd qddot_lagrangian (qddot_emulated);
  VectorNd qddot_kokkevis (qddot_emulated);

  ForwardDynamicsContactsDirect (*model_emulated, q, qdot, tau, constraints_1B1C_emulated, qddot_lagrangian);
  ForwardDynamicsContactsKokkevis (*model_emulated, q, qdot, tau, constraints_1B1C_emulated, qddot_kokkevis);
  EXPECT_TRUE(EIGEN_MATRIX_NEAR (qddot_lagrangian, qddot_kokkevis, TEST_PREC * qddot_lagrangian.norm()));

  ForwardDynamicsContactsDirect (*model_emulated, q, qdot, tau, constraints_1B4C_emulated, qddot_lagrangian);
  ForwardDynamicsContactsKokkevis (*model_emulated, q, qdot, tau, constraints_1B4C_emulated, qddot_kokkevis);
  EXPECT_TRUE(EIGEN_MATRIX_NEAR (qddot_lagrangian, qddot_kokkevis, TEST_PREC * qddot_lagrangian.norm()));

  ForwardDynamicsContactsDirect (*model_emulated, q, qdot, tau, constraints_4B4C_emulated, qddot_lagrangian);
  ForwardDynamicsContactsKokkevis (*model_emulated, q, qdot, tau, constraints_4B4C_emulated, qddot_kokkevis);

  EXPECT_TRUE(EIGEN_MATRIX_NEAR (qddot_lagrangian, qddot_kokkevis, TEST_PREC * qddot_lagrangian.norm()));
}

TEST_F  (Human36, TestContactsEmulatedLagrangianSparse) {
  for (unsigned int i = 0; i < q.size(); i++) {
    q[i] = 0.5 * M_PI * static_cast<double>(rand()) / static_cast<double>(RAND_MAX);
    qdot[i] = 0.5 * M_PI * static_cast<double>(rand()) / static_cast<double>(RAND_MAX);
    tau[i] = -0.5 * M_PI * static_cast<double>(rand()) / static_cast<double>(RAND_MAX);
  }

  VectorNd qddot_lagrangian (qddot_emulated);
  VectorNd qddot_sparse (qddot_emulated);

  ForwardDynamicsContactsDirect (*model_emulated, q, qdot, tau, constraints_1B1C_emulated, qddot_lagrangian);
  ForwardDynamicsContactsRangeSpaceSparse (*model_emulated, q, qdot, tau, constraints_1B1C_emulated, qddot_sparse);
  EXPECT_TRUE(EIGEN_MATRIX_NEAR (qddot_lagrangian, qddot_sparse, TEST_PREC * qddot_lagrangian.norm()));

  ForwardDynamicsContactsDirect (*model_emulated, q, qdot, tau, constraints_1B4C_emulated, qddot_lagrangian);
  ForwardDynamicsContactsRangeSpaceSparse (*model_emulated, q, qdot, tau, constraints_1B4C_emulated, qddot_sparse);
  EXPECT_TRUE(EIGEN_MATRIX_NEAR (qddot_lagrangian, qddot_sparse, TEST_PREC * qddot_lagrangian.norm()));

  ForwardDynamicsContactsDirect (*model_emulated, q, qdot, tau, constraints_4B4C_emulated, qddot_lagrangian);
  ForwardDynamicsContactsRangeSpaceSparse (*model_emulated, q, qdot, tau, constraints_4B4C_emulated, qddot_sparse);
  EXPECT_TRUE(EIGEN_MATRIX_NEAR (qddot_lagrangian, qddot_sparse, TEST_PREC * qddot_lagrangian.norm()));
}

TEST_F  (Human36, TestContactsEmulatedLagrangianNullSpace) {
  for (unsigned int i = 0; i < q.size(); i++) {
    q[i] = 0.5 * M_PI * static_cast<double>(rand()) / static_cast<double>(RAND_MAX);
    qdot[i] = 0.5 * M_PI * static_cast<double>(rand()) / static_cast<double>(RAND_MAX);
    tau[i] = -0.5 * M_PI * static_cast<double>(rand()) / static_cast<double>(RAND_MAX);
  }

  VectorNd qddot_lagrangian (qddot_emulated);
  VectorNd qddot_nullspace (qddot_emulated);

  ForwardDynamicsContactsDirect (*model_emulated, q, qdot, tau, constraints_1B1C_emulated, qddot_lagrangian);
  ForwardDynamicsContactsNullSpace (*model_emulated, q, qdot, tau, constraints_1B1C_emulated, qddot_nullspace);
  EXPECT_TRUE(EIGEN_MATRIX_NEAR (qddot_lagrangian, qddot_nullspace, TEST_PREC * qddot_lagrangian.norm()));

  ForwardDynamicsContactsDirect (*model_emulated, q, qdot, tau, constraints_1B4C_emulated, qddot_lagrangian);
  ForwardDynamicsContactsNullSpace (*model_emulated, q, qdot, tau, constraints_1B4C_emulated, qddot_nullspace);
  EXPECT_TRUE(EIGEN_MATRIX_NEAR (qddot_lagrangian, qddot_nullspace,  TEST_PREC * qddot_lagrangian.norm()));

  ForwardDynamicsContactsDirect (*model_emulated, q, qdot, tau, constraints_4B4C_emulated, qddot_lagrangian);
  ForwardDynamicsContactsNullSpace (*model_emulated, q, qdot, tau, constraints_4B4C_emulated, qddot_nullspace);
  EXPECT_TRUE(EIGEN_MATRIX_NEAR (qddot_lagrangian, qddot_nullspace, TEST_PREC * qddot_lagrangian.norm()));
}

TEST_F  (Human36, TestContactsEmulatedMultdofLagrangian) {
  for (unsigned int i = 0; i < q.size(); i++) {
    q[i] = 0.5 * M_PI * static_cast<double>(rand()) / static_cast<double>(RAND_MAX);
    qdot[i] = 0.5 * M_PI * static_cast<double>(rand()) / static_cast<double>(RAND_MAX);
    qddot[i] = 0.5 * M_PI * static_cast<double>(rand()) / static_cast<double>(RAND_MAX);
    tau[i] = 0.5 * M_PI * static_cast<double>(rand()) / static_cast<double>(RAND_MAX);
  }

  ForwardDynamicsContactsDirect (*model_emulated, q, qdot, tau, constraints_1B1C_emulated, qddot_emulated);
  ForwardDynamicsContactsDirect (*model_3dof, q, qdot, tau, constraints_1B1C_3dof, qddot_3dof);
  EXPECT_TRUE(EIGEN_MATRIX_NEAR (qddot_emulated, qddot_3dof,  TEST_PREC * qddot_emulated.norm()));

  ForwardDynamicsContactsDirect (*model_emulated, q, qdot, tau, constraints_1B4C_emulated, qddot_emulated);
  ForwardDynamicsContactsDirect (*model_3dof, q, qdot, tau, constraints_1B4C_3dof, qddot_3dof);
  EXPECT_TRUE(EIGEN_MATRIX_NEAR (qddot_emulated, qddot_3dof, TEST_PREC * qddot_emulated.norm()));

  ForwardDynamicsContactsDirect (*model_emulated, q, qdot, tau, constraints_4B4C_emulated, qddot_emulated);
  ForwardDynamicsContactsDirect (*model_3dof, q, qdot, tau, constraints_4B4C_3dof, qddot_3dof);
  EXPECT_TRUE(EIGEN_MATRIX_NEAR (qddot_emulated, qddot_3dof, TEST_PREC * qddot_emulated.norm()));
}

TEST_F  (Human36, TestContactsEmulatedMultdofSparse) {
  for (unsigned int i = 0; i < q.size(); i++) {
    q[i] = 0.5 * M_PI * static_cast<double>(rand()) / static_cast<double>(RAND_MAX);
    qdot[i] = 0.5 * M_PI * static_cast<double>(rand()) / static_cast<double>(RAND_MAX);
    tau[i] = 0.5 * M_PI * static_cast<double>(rand()) / static_cast<double>(RAND_MAX);
  }

  ForwardDynamicsContactsRangeSpaceSparse (*model_emulated, q, qdot, tau, constraints_1B1C_emulated, qddot_emulated);

  for (unsigned int i = 0; i < q.size(); i++) {
    EXPECT_EQ (model_emulated->lambda_q[i], model_3dof->lambda_q[i]);
  }

  ForwardDynamicsContactsRangeSpaceSparse (*model_3dof, q, qdot, tau, constraints_1B1C_3dof, qddot_3dof);
  EXPECT_TRUE(EIGEN_MATRIX_NEAR (qddot_emulated, qddot_3dof, TEST_PREC * qddot_emulated.norm()));

  ForwardDynamicsContactsRangeSpaceSparse (*model_emulated, q, qdot, tau, constraints_1B4C_emulated, qddot_emulated);
  ForwardDynamicsContactsRangeSpaceSparse (*model_3dof, q, qdot, tau, constraints_1B4C_3dof, qddot_3dof);
  EXPECT_TRUE(EIGEN_MATRIX_NEAR (qddot_emulated, qddot_3dof, TEST_PREC * qddot_emulated.norm()));

  ForwardDynamicsContactsRangeSpaceSparse (*model_emulated, q, qdot, tau, constraints_4B4C_emulated, qddot_emulated);
  ForwardDynamicsContactsRangeSpaceSparse (*model_3dof, q, qdot, tau, constraints_4B4C_3dof, qddot_3dof);
  EXPECT_TRUE(EIGEN_MATRIX_NEAR (qddot_emulated, qddot_3dof, TEST_PREC * qddot_emulated.norm()));
}

TEST_F  (Human36, TestContactsEmulatedMultdofKokkevisSparse) {
  for (unsigned int i = 0; i < q.size(); i++) {
    q[i] = 0.5 * M_PI * static_cast<double>(rand()) / static_cast<double>(RAND_MAX);
    qdot[i] = 0.5 * M_PI * static_cast<double>(rand()) / static_cast<double>(RAND_MAX);
    tau[i] = 0.5 * M_PI * static_cast<double>(rand()) / static_cast<double>(RAND_MAX);
  }

  ForwardDynamicsContactsRangeSpaceSparse (*model_emulated, q, qdot, tau, constraints_1B1C_emulated, qddot_emulated);

  for (unsigned int i = 0; i < q.size(); i++) {
    EXPECT_EQ (model_emulated->lambda_q[i], model_3dof->lambda_q[i]);
  }

  VectorNd qddot_sparse (qddot_emulated);
  VectorNd qddot_kokkevis (qddot_emulated);

  ForwardDynamicsContactsRangeSpaceSparse (*model_3dof, q, qdot, tau, constraints_1B1C_3dof, qddot_sparse);
  ForwardDynamicsContactsKokkevis (*model_3dof, q, qdot, tau, constraints_1B1C_3dof, qddot_kokkevis);
  EXPECT_TRUE(EIGEN_MATRIX_NEAR (qddot_sparse, qddot_kokkevis, TEST_PREC * qddot_sparse.norm()));

  ForwardDynamicsContactsRangeSpaceSparse (*model_3dof, q, qdot, tau, constraints_1B4C_3dof, qddot_sparse);
  ForwardDynamicsContactsKokkevis (*model_3dof, q, qdot, tau, constraints_1B4C_3dof, qddot_kokkevis);
  EXPECT_TRUE(EIGEN_MATRIX_NEAR (qddot_sparse, qddot_kokkevis, TEST_PREC * qddot_sparse.norm()));

  ForwardDynamicsContactsRangeSpaceSparse (*model_3dof, q, qdot, tau, constraints_4B4C_3dof, qddot_sparse);
  ForwardDynamicsContactsKokkevis (*model_3dof, q, qdot, tau, constraints_4B4C_3dof, qddot_kokkevis);
  EXPECT_TRUE(EIGEN_MATRIX_NEAR (qddot_sparse, qddot_kokkevis, TEST_PREC * qddot_sparse.norm()));
}

TEST_F  (Human36, TestContactsEmulatedMultdofKokkevisMultiple ) {
  for (unsigned int i = 0; i < q.size(); i++) {
    q[i] = 0.5 * M_PI * static_cast<double>(rand()) / static_cast<double>(RAND_MAX);
    qdot[i] = 0.5 * M_PI * static_cast<double>(rand()) / static_cast<double>(RAND_MAX);
    tau[i] = 0.5 * M_PI * static_cast<double>(rand()) / static_cast<double>(RAND_MAX);
  }

  VectorNd qddot_kokkevis (qddot_emulated);
  VectorNd qddot_kokkevis_2 (qddot_emulated);

  ForwardDynamicsContactsKokkevis (*model_3dof, q, qdot, tau, constraints_1B1C_3dof, qddot_kokkevis);
  ForwardDynamicsContactsKokkevis (*model_3dof, q, qdot, tau, constraints_1B1C_3dof, qddot_kokkevis_2);
  EXPECT_TRUE(EIGEN_MATRIX_NEAR (qddot_kokkevis, qddot_kokkevis_2, TEST_PREC * qddot_kokkevis.norm()));

  ForwardDynamicsContactsKokkevis (*model_3dof, q, qdot, tau, constraints_1B4C_3dof, qddot_kokkevis);
  ForwardDynamicsContactsKokkevis (*model_3dof, q, qdot, tau, constraints_1B4C_3dof, qddot_kokkevis_2);
  EXPECT_TRUE(EIGEN_MATRIX_NEAR (qddot_kokkevis, qddot_kokkevis_2, TEST_PREC * qddot_kokkevis.norm()));

  ForwardDynamicsContactsKokkevis (*model_3dof, q, qdot, tau, constraints_4B4C_3dof, qddot_kokkevis);
  ForwardDynamicsContactsKokkevis (*model_3dof, q, qdot, tau, constraints_4B4C_3dof, qddot_kokkevis_2);
  EXPECT_TRUE(EIGEN_MATRIX_NEAR (qddot_kokkevis, qddot_kokkevis_2,TEST_PREC * qddot_kokkevis.norm()));
}

TEST_F (SphericalJoint, TestEulerZYXvsEmulated ) {
  emuQ[0] = 1.1;
  emuQ[1] = 1.2;
  emuQ[2] = 1.3;
  emuQ[3] = 1.4;
  emuQ[4] = 1.5;

  emuQDot[0] = 1.;
  emuQDot[1] = 2.;
  emuQDot[2] = 3.;
  emuQDot[3] = 4.;
  emuQDot[4] = 5.;

  emuTau[0] = 5.;
  emuTau[1] = 4.;
  emuTau[2] = 7.;
  emuTau[3] = 3.;
  emuTau[4] = 2.;

  VectorNd QDDot_emu = VectorNd::Zero (emulated_model.qdot_size);
  VectorNd QDDot_eulerzyx = VectorNd::Zero (eulerzyx_model.qdot_size);

  ForwardDynamicsLagrangian (emulated_model, emuQ, emuQDot, emuTau, QDDot_emu);
  ForwardDynamicsLagrangian (eulerzyx_model, emuQ, emuQDot, emuTau, QDDot_eulerzyx);

  EXPECT_TRUE(EIGEN_MATRIX_NEAR (QDDot_emu, QDDot_eulerzyx, TEST_PREC));
}

TEST_F  ( Human36, SolveMInvTimesTau) {
  for (unsigned int i = 0; i < model->dof_count; i++) {
    q[i] = rand() / static_cast<double>(RAND_MAX);
    tau[i] = rand() / static_cast<double>(RAND_MAX);
  }

  MatrixNd M (MatrixNd::Zero(model->dof_count, model->dof_count));
  CompositeRigidBodyAlgorithm (*model, q, M);

  VectorNd qddot_solve_llt = M.llt().solve(tau);

  VectorNd qddot_minv (q);
  CalcMInvTimesTau (*model, q, tau, qddot_minv);

  EXPECT_TRUE(EIGEN_MATRIX_NEAR (qddot_solve_llt, qddot_minv, TEST_PREC * qddot_minv.norm()));
}

TEST_F  ( Human36, SolveMInvTimesTauReuse) {
  for (unsigned int i = 0; i < model->dof_count; i++) {
    q[i] = rand() / static_cast<double>(RAND_MAX);
    tau[i] = rand() / static_cast<double>(RAND_MAX);
  }

  MatrixNd M (MatrixNd::Zero(model->dof_count, model->dof_count));
  CompositeRigidBodyAlgorithm (*model, q, M);

  VectorNd qddot_solve_llt = M.llt().solve(tau);

  VectorNd qddot_minv (q);
  CalcMInvTimesTau (*model, q, tau, qddot_minv);

  for (unsigned int j = 0; j < 1; j++) {
    for (unsigned int i = 0; i < model->dof_count; i++) {
      tau[i] = rand() / static_cast<double>(RAND_MAX);
    }

    CompositeRigidBodyAlgorithm (*model, q, M);
    qddot_solve_llt = M.llt().solve(tau);

    CalcMInvTimesTau (*model, q, tau, qddot_minv);

    EXPECT_TRUE(EIGEN_MATRIX_NEAR (qddot_solve_llt, qddot_minv, TEST_PREC * qddot_solve_llt.norm()));
  }
}

int main( int argc, char** argv)
{
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
