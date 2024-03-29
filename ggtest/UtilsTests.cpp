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
#include "rbdl/Energy.h"

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
  CompositeRigidBodyAlgorithm (*model, *model_data, q, H, true);

  double kinetic_energy_ref = 0.5 * qdot.transpose() * H * qdot;
  double kinetic_energy = Utils::CalcKineticEnergy (*model, *model_data, q, qdot);

  EXPECT_EQ (kinetic_energy_ref, kinetic_energy);
}

TEST(UtilsTests, TestPotentialEnergy) {
  ModelDatad model_data;
  Model model(model_data);
  Matrix3d inertia = Matrix3d::Zero(3,3);
  Body body (0.5, Vector3d (0., 0., 0.), inertia);
  Joint joint (
      SpatialVectord (0., 0., 0., 1., 0., 0.),
      SpatialVectord (0., 0., 0., 0., 1., 0.),
      SpatialVectord (0., 0., 0., 0., 0., 1.)
      );

  model.AppendBody (model_data, Xtrans<double> (Vector3d::Zero()), joint, body);

  VectorNd q = VectorNd::Zero(model.q_size);
  double potential_energy_zero = Utils::CalcPotentialEnergy (model, model_data, q);
  EXPECT_EQ (0., potential_energy_zero);

  q[1] = 1.;
  double potential_energy_lifted = Utils::CalcPotentialEnergy (model, model_data, q);
  EXPECT_EQ (4.905, potential_energy_lifted);
}

TEST(UtilsTests, TestCOMSimple) {
  ModelDatad model_data;
  Model model(model_data);
  Matrix3d inertia = Matrix3d::Zero(3,3);
  Body body (123., Vector3d (0., 0., 0.), inertia);
  Joint joint (
      SpatialVectord (0., 0., 0., 1., 0., 0.),
      SpatialVectord (0., 0., 0., 0., 1., 0.),
      SpatialVectord (0., 0., 0., 0., 0., 1.)
      );

  model.AppendBody (model_data, Xtrans<double> (Vector3d::Zero()), joint, body);

  VectorNd q = VectorNd::Zero(model.q_size);
  VectorNd qdot = VectorNd::Zero(model.qdot_size);

  double mass;
  Vector3d com;
  Vector3d com_velocity;
  Utils::CalcCenterOfMass (model, model_data, q, qdot, mass, com, &com_velocity);

  EXPECT_EQ (123., mass);
  EXPECT_EQ (Vector3d (0., 0., 0.), com);
  EXPECT_EQ (Vector3d (0., 0., 0.), com_velocity);

  q[1] = 1.;
  Utils::CalcCenterOfMass (model, model_data, q, qdot, mass, com, &com_velocity);
  EXPECT_EQ (Vector3d (0., 1., 0.), com);
  EXPECT_EQ (Vector3d (0., 0., 0.), com_velocity);

  qdot[1] = 1.;
  Utils::CalcCenterOfMass (model, model_data, q, qdot, mass, com, &com_velocity);
  EXPECT_EQ (Vector3d (0., 1., 0.), com);
  EXPECT_EQ (Vector3d (0., 1., 0.), com_velocity);
}

TEST(UtilsTests, TestAngularMomentumSimple) {
  ModelDatad model_data;
  Model model(model_data);
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

  model.AppendBody (model_data, Xtrans (Vector3d(0., 0., 0.)), joint, body);

  VectorNd q = VectorNd::Zero(model.q_size);
  VectorNd qdot = VectorNd::Zero(model.qdot_size);

  double mass;
  Vector3d com;
  Vector3d angular_momentum;

  qdot << 1., 0., 0.;
  Utils::CalcCenterOfMass (model, model_data, q, qdot, mass, com, NULL, &angular_momentum);
  EXPECT_EQ (Vector3d (1.1, 0., 0.), angular_momentum);

  qdot << 0., 1., 0.;
  Utils::CalcCenterOfMass (model, model_data, q, qdot, mass, com, NULL, &angular_momentum);
  EXPECT_EQ (Vector3d (0., 2.2, 0.), angular_momentum);

  qdot << 0., 0., 1.;
  Utils::CalcCenterOfMass (model, model_data, q, qdot, mass, com, NULL, &angular_momentum);
  EXPECT_EQ (Vector3d (0., 0., 3.3), angular_momentum);
}

TEST_F (TwoArms12DoF, TestAngularMomentumSimple) {
  double mass;
  Vector3d com;
  Vector3d angular_momentum;

  Utils::CalcCenterOfMass (*model, *model_data, q, qdot, mass, com, NULL, &angular_momentum);

  EXPECT_EQ (Vector3d (0., 0., 0.), angular_momentum);

  qdot[0] = 1.;
  qdot[1] = 2.;
  qdot[2] = 3.;

  Utils::CalcCenterOfMass (*model, *model_data, q, qdot, mass, com, NULL, &angular_momentum);

  // only a rough guess from test calculation
  EXPECT_TRUE(EIGEN_MATRIX_NEAR(Vector3d (3.3, 2.54, 1.5), angular_momentum, 1.0e-1));

  qdot[3] = -qdot[0];
  qdot[4] = -qdot[1];
  qdot[5] = -qdot[2];

  ClearLogOutput();
  Utils::CalcCenterOfMass (*model, *model_data, q, qdot, mass, com, NULL, &angular_momentum);

  EXPECT_TRUE (angular_momentum[0] == 0);
  EXPECT_TRUE (angular_momentum[1] < 0);
  EXPECT_TRUE (angular_momentum[2] == 0.);
}

TEST(UtilsTests, TestSpatialInteria) {
  Model model;
  Matrix3d inertia = Matrix3d::Zero(3, 3);
  inertia(0, 0) = 1.1;
  inertia(1, 1) = 2.2;
  inertia(2, 2) = 3.3;

  Body body(0.5, Vector3d(1., 0., 0.), inertia);
  Joint joint(SpatialVectord(1., 0., 0., 0., 0., 0.), SpatialVectord(0., 1., 0., 0., 0., 0.),
              SpatialVectord(0., 0., 1., 0., 0., 0.));

  model.AppendBody(*model.getModelData(), Xtrans(Vector3d(0., 0., 0.)), joint, body);

  VectorNd q = VectorNd::Zero(model.q_size);
  VectorNd qdot_zero = VectorNd::Zero(model.qdot_size);
  VectorNd qdot_nonzero_1 = VectorNd::Zero(model.qdot_size);
  VectorNd qdot_nonzero_2 = VectorNd::Zero(model.qdot_size);
  VectorNd qdot_nonzero_3 = VectorNd::Zero(model.qdot_size);
  qdot_nonzero_1.setRandom();
  qdot_nonzero_2.setRandom();
  qdot_nonzero_3.setRandom();

  // Computed with the info inside CalcCenterOfMass method
  Matrix3d world_inertia = Matrix3d::Zero(3, 3);
  world_inertia(0, 0) = 1.1;
  world_inertia(1, 1) = 2.7;
  world_inertia(2, 2) = 3.8;
  double mass;
  Vector3d com;
  Vector3d angular_momentum;
  Vector3d computed_angular_momentum;

  for (const VectorNd& qdot : { qdot_zero, qdot_nonzero_1, qdot_nonzero_2, qdot_nonzero_3 })
  {
    // Test to compute Inertia and Angular Momentum along CoM and a different point

    // Compute Angular Momentum with standard method
    Utils::CalcCenterOfMass(model, q, qdot, mass, com, NULL, &angular_momentum);

    // Main method for the compuration of inertia along CoM
    Eigen::Matrix3d com_inertia = calcGlobalInertiaTensorFromCOM(model, q, true);
    EXPECT_TRUE(EIGEN_MATRIX_NEAR(inertia, com_inertia, 1e-8));

    // Computation at Center of Mass with a specific method
    SpatialRigidBodyInertiad computed_com_inertia =
        SpatialRigidBodyInertiad(0., Vector3d(0., 0., 0.), Matrix3d::Zero(3, 3));
    Utils::CalcCentroidalInertiaMatrix(model, q, qdot, computed_com_inertia,
                                       &computed_angular_momentum, true);
    EXPECT_TRUE(EIGEN_MATRIX_NEAR(inertia,
                                  computed_com_inertia.toMatrix().block(0, 0, 3, 3), 1e-8));
    EXPECT_TRUE(EIGEN_MATRIX_NEAR(Eigen::Vector3d::Zero(), computed_com_inertia.h, 1e-8));
    EXPECT_TRUE(EIGEN_MATRIX_NEAR(angular_momentum, computed_angular_momentum, 1e-8));
    EXPECT_NEAR(body.mMass, computed_com_inertia.m, 1e-8);

    // Computation at a point same as the center of mass
    computed_com_inertia =
        SpatialRigidBodyInertiad(0., Vector3d(0., 0., 0.), Matrix3d::Zero(3, 3));
    Utils::CalcPointSpatialInertiaMatrix(model, *model.getModelData(), q, qdot,
                                         body.mCenterOfMass, computed_com_inertia,
                                         &computed_angular_momentum, true);
    EXPECT_TRUE(EIGEN_MATRIX_NEAR(inertia,
                                  computed_com_inertia.toMatrix().block(0, 0, 3, 3), 1e-8));
    Math::Vector3d computed_com_pos = computed_com_inertia.h / computed_com_inertia.m;
    EXPECT_TRUE(EIGEN_MATRIX_NEAR(Eigen::Vector3d::Zero(), computed_com_pos, 1e-8));
    EXPECT_TRUE(EIGEN_MATRIX_NEAR(angular_momentum, computed_angular_momentum, 1e-8));
    EXPECT_NEAR(body.mMass, computed_com_inertia.m, 1e-8);

    // Computation with respect to world
    computed_com_inertia =
        SpatialRigidBodyInertiad(0., Vector3d(0., 0., 0.), Matrix3d::Zero(3, 3));
    Utils::CalcPointSpatialInertiaMatrix(model, *model.getModelData(), q, qdot,
                                         Math::Vector3d(0, 0, 0), computed_com_inertia,
                                         &computed_angular_momentum, true);
    EXPECT_TRUE(EIGEN_MATRIX_NEAR(world_inertia,
                                  computed_com_inertia.toMatrix().block(0, 0, 3, 3), 1e-8));
    computed_com_pos = computed_com_inertia.h / computed_com_inertia.m;
    // Check with standard formulation of Angular Momentum = Inertia * Angular Velocity
    EXPECT_TRUE(EIGEN_MATRIX_NEAR(computed_com_inertia.toMatrix().block(0, 0, 3, 3) * qdot,
                                  computed_angular_momentum, 1e-8));
    EXPECT_TRUE(EIGEN_MATRIX_NEAR(body.mCenterOfMass,
                                  computed_com_inertia.h / computed_com_inertia.m, 1e-8));
    EXPECT_NEAR(body.mMass, computed_com_inertia.m, 1e-8);
  }

  // Addition of new body
  Matrix3d inertia_1 = Matrix3d::Zero(3, 3);
  inertia_1(0, 0) = 2.2;
  inertia_1(1, 1) = 4.4;
  inertia_1(2, 2) = 6.6;
  Body body_1(2.0, Vector3d(1., 0., 0.), inertia_1);
  Joint joint_1(JointTypeRevoluteX);
  model.AppendBody(*model.getModelData(), Xtrans(Vector3d(1., 0., 0.)), joint_1, body_1);
  q = VectorNd::Zero(model.q_size);
  qdot_zero = VectorNd::Zero(model.qdot_size);
  qdot_nonzero_1 = VectorNd::Zero(model.qdot_size);
  qdot_nonzero_2 = VectorNd::Zero(model.qdot_size);
  qdot_nonzero_3 = VectorNd::Zero(model.qdot_size);
  qdot_nonzero_1.setRandom();
  qdot_nonzero_2.setRandom();
  qdot_nonzero_3.setRandom();

  // Computed with the info inside CalcCenterOfMass method
  world_inertia(0, 0) = 3.3;
  world_inertia(1, 1) = 15.1;
  world_inertia(2, 2) = 18.4;

  Matrix3d combined_inertia;
  combined_inertia.setZero();
  combined_inertia(0,0) = 3.3;
  combined_inertia(1,1) = 7;
  combined_inertia(2,2) = 10.3;

  for (const VectorNd& qdot : { qdot_zero, qdot_nonzero_1, qdot_nonzero_2, qdot_nonzero_3 })
  {
    // Test to compute Inertia and Angular Momentum along CoM and a different point

    // Compute Angular Momentum with standard method
    Utils::CalcCenterOfMass(model, q, qdot, mass, com, NULL, &angular_momentum);

    // Main method for the compuration of inertia along CoM
    Eigen::Matrix3d com_inertia = calcGlobalInertiaTensorFromCOM(model, q, true);
    EXPECT_TRUE(EIGEN_MATRIX_NEAR(combined_inertia, com_inertia, 1e-8));

    // Computation at Center of Mass with a specific method
    SpatialRigidBodyInertiad computed_com_inertia =
        SpatialRigidBodyInertiad(0., Vector3d(0., 0., 0.), Matrix3d::Zero(3, 3));
    Utils::CalcCentroidalInertiaMatrix(model, q, qdot, computed_com_inertia,
                                       &computed_angular_momentum, true);
    EXPECT_TRUE(EIGEN_MATRIX_NEAR(combined_inertia,
                                  computed_com_inertia.toMatrix().block(0, 0, 3, 3), 1e-8));
    EXPECT_TRUE(EIGEN_MATRIX_NEAR(Eigen::Vector3d::Zero(), computed_com_inertia.h, 1e-8));
    EXPECT_TRUE(EIGEN_MATRIX_NEAR(angular_momentum, computed_angular_momentum, 1e-8));
    EXPECT_NEAR(mass, computed_com_inertia.m, 1e-8);

    // Computation at a point same as the center of mass
    computed_com_inertia =
        SpatialRigidBodyInertiad(0., Vector3d(0., 0., 0.), Matrix3d::Zero(3, 3));
    Utils::CalcPointSpatialInertiaMatrix(model, *model.getModelData(), q, qdot,
                                         com, computed_com_inertia,
                                         &computed_angular_momentum, true);
    EXPECT_TRUE(EIGEN_MATRIX_NEAR(combined_inertia,
                                  computed_com_inertia.toMatrix().block(0, 0, 3, 3), 1e-8));
    Math::Vector3d computed_com_pos = computed_com_inertia.h / computed_com_inertia.m;
    EXPECT_TRUE(EIGEN_MATRIX_NEAR(Eigen::Vector3d::Zero(), computed_com_pos, 1e-8));
    EXPECT_TRUE(EIGEN_MATRIX_NEAR(angular_momentum, computed_angular_momentum, 1e-8));
    EXPECT_NEAR(mass, computed_com_inertia.m, 1e-8);

    // Computation with respect to world
    computed_com_inertia =
        SpatialRigidBodyInertiad(0., Vector3d(0., 0., 0.), Matrix3d::Zero(3, 3));
    Utils::CalcPointSpatialInertiaMatrix(model, *model.getModelData(), q, qdot,
                                         Math::Vector3d(0, 0, 0), computed_com_inertia,
                                         &computed_angular_momentum, true);
    EXPECT_TRUE(EIGEN_MATRIX_NEAR(world_inertia,
                                  computed_com_inertia.toMatrix().block(0, 0, 3, 3), 1e-8));
    EXPECT_TRUE(EIGEN_MATRIX_NEAR(com,
                                  computed_com_inertia.h / computed_com_inertia.m, 1e-8));
    EXPECT_NEAR(mass, computed_com_inertia.m, 1e-8);
  }
}

int main( int argc, char** argv)
{
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
