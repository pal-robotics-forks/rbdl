#include <gtest/gtest.h>
#include <eigen_checks/gtest.h>

#include <iostream>

#include "FixturesTests.h"
#include "rbdl/rbdl_mathutils.h"
#include "rbdl/rbdl_utils.h"
#include "rbdl/Logging.h"

#include "rbdl/Model.h"
#include "rbdl/Kinematics.h"
#include "rbdl/Dynamics.h"

using namespace std;
using namespace RigidBodyDynamics;
using namespace RigidBodyDynamics::Math;

TEST_F(FloatingBase12DoF, TestKineticEnergy) {
  VectorNd q = VectorNd::Zero(model->q_size);
  VectorNd qdot = VectorNd::Zero(model->q_size);

  for (unsigned int i = 0; i < q.size(); i++) {
    q[i] = 0.1 * i;
    qdot[i] = 0.3 * i;
  }

  MatrixNd H = MatrixNd::Zero (model->q_size, model->q_size);
  CompositeRigidBodyAlgorithm (*model, q, H, true);

  double kinetic_energy_ref = 0.5 * qdot.transpose() * H * qdot;
  double kinetic_energy = Utils::CalcKineticEnergy (*model, q, qdot);

  EXPECT_EQ (kinetic_energy_ref, kinetic_energy);
}

TEST(UtilsTests, TestPotentialEnergy) {
  Model model;
  Matrix3d inertia = Matrix3d::Zero(3,3);
  Body body (0.5, Vector3d (0., 0., 0.), inertia);
  Joint joint (
      SpatialVectord (0., 0., 0., 1., 0., 0.),
      SpatialVectord (0., 0., 0., 0., 1., 0.),
      SpatialVectord (0., 0., 0., 0., 0., 1.)
      );

  model.AppendBody (Xtrans<double> (Vector3d::Zero()), joint, body);

  VectorNd q = VectorNd::Zero(model.q_size);
  double potential_energy_zero = Utils::CalcPotentialEnergy (model, q);
  EXPECT_EQ (0., potential_energy_zero);

  q[1] = 1.;
  double potential_energy_lifted = Utils::CalcPotentialEnergy (model, q);
  EXPECT_EQ (4.905, potential_energy_lifted);
}

TEST(UtilsTests, TestCOMSimple) {
  Model model;
  Matrix3d inertia = Matrix3d::Zero(3,3);
  Body body (123., Vector3d (0., 0., 0.), inertia);
  Joint joint (
      SpatialVectord (0., 0., 0., 1., 0., 0.),
      SpatialVectord (0., 0., 0., 0., 1., 0.),
      SpatialVectord (0., 0., 0., 0., 0., 1.)
      );

  model.AppendBody (Xtrans<double> (Vector3d::Zero()), joint, body);

  VectorNd q = VectorNd::Zero(model.q_size);
  VectorNd qdot = VectorNd::Zero(model.qdot_size);

  double mass;
  Vector3d com;
  Vector3d com_velocity;
  Utils::CalcCenterOfMass (model, q, qdot, mass, com, &com_velocity);

  EXPECT_EQ (123., mass);
  EXPECT_EQ (Vector3d (0., 0., 0.), com);
  EXPECT_EQ (Vector3d (0., 0., 0.), com_velocity);

  q[1] = 1.;
  Utils::CalcCenterOfMass (model, q, qdot, mass, com, &com_velocity);
  EXPECT_EQ (Vector3d (0., 1., 0.), com);
  EXPECT_EQ (Vector3d (0., 0., 0.), com_velocity);

  qdot[1] = 1.;
  Utils::CalcCenterOfMass (model, q, qdot, mass, com, &com_velocity);
  EXPECT_EQ (Vector3d (0., 1., 0.), com);
  EXPECT_EQ (Vector3d (0., 1., 0.), com_velocity);
}

TEST(UtilsTests, TestAngularMomentumSimple) {
  Model model;
  Matrix3d inertia = Matrix3d::Zero(3,3);
  inertia(0,0) = 1.1;
  inertia(1,1) = 2.2;
  inertia(2,2) = 3.3;

  Body body (0.5, Vector3d (1., 0., 0.), inertia);
  Joint joint (
      SpatialVectord (1., 0., 0., 0., 0., 0.),
      SpatialVectord (0., 1., 0., 0., 0., 0.),
      SpatialVectord (0., 0., 1., 0., 0., 0.)
      );

  model.AppendBody (Xtrans (Vector3d(0., 0., 0.)), joint, body);

  VectorNd q = VectorNd::Zero(model.q_size);
  VectorNd qdot = VectorNd::Zero(model.qdot_size);

  double mass;
  Vector3d com;
  Vector3d angular_momentum;

  qdot << 1., 0., 0.;
  Utils::CalcCenterOfMass (model, q, qdot, mass, com, NULL, &angular_momentum);
  EXPECT_EQ (Vector3d (1.1, 0., 0.), angular_momentum);

  qdot << 0., 1., 0.;
  Utils::CalcCenterOfMass (model, q, qdot, mass, com, NULL, &angular_momentum);
  EXPECT_EQ (Vector3d (0., 2.2, 0.), angular_momentum);

  qdot << 0., 0., 1.;
  Utils::CalcCenterOfMass (model, q, qdot, mass, com, NULL, &angular_momentum);
  EXPECT_EQ (Vector3d (0., 0., 3.3), angular_momentum);
}

TEST_F (TwoArms12DoF, TestAngularMomentumSimple) {
  double mass;
  Vector3d com;
  Vector3d angular_momentum;

  Utils::CalcCenterOfMass (*model, q, qdot, mass, com, NULL, &angular_momentum);

  EXPECT_EQ (Vector3d (0., 0., 0.), angular_momentum);

  qdot[0] = 1.;
  qdot[1] = 2.;
  qdot[2] = 3.;

  Utils::CalcCenterOfMass (*model, q, qdot, mass, com, NULL, &angular_momentum);

  // only a rough guess from test calculation
  EXPECT_TRUE(EIGEN_MATRIX_NEAR(Vector3d (3.3, 2.54, 1.5), angular_momentum, 1.0e-1));

  qdot[3] = -qdot[0];
  qdot[4] = -qdot[1];
  qdot[5] = -qdot[2];

  ClearLogOutput();
  Utils::CalcCenterOfMass (*model, q, qdot, mass, com, NULL, &angular_momentum);

  EXPECT_TRUE (angular_momentum[0] == 0);
  EXPECT_TRUE (angular_momentum[1] < 0);
  EXPECT_TRUE (angular_momentum[2] == 0.);
}

int main( int argc, char** argv)
{
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
