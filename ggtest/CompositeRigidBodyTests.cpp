
#include <gtest/gtest.h>
#include <eigen_checks/gtest.h>

#include <iostream>

#include "rbdl/Logging.h"
#include "rbdl/Model.h"
#include "rbdl/Dynamics.h"

#include "FixturesTests.h"

using namespace std;
using namespace RigidBodyDynamics;
using namespace RigidBodyDynamics::Math;

const double TEST_PREC = 1.0e-12;

class CompositeRigidBodyFixture : public ::testing::Test {
protected:
  virtual void SetUp () {
    ClearLogOutput();
    model_data = new ModelData;
    model = new Model(*model_data);
    model->gravity = Vector3d (0., -9.81, 0.);
  }
  virtual void TearDown () {
    delete model;
  }
  ModelDatad *model_data;
  Model *model;
};

TEST_F(CompositeRigidBodyFixture, TestCompositeRigidBodyForwardDynamicsFloatingBase) {
  Body base_body(1., Vector3d (1., 0., 0.), Vector3d (1., 1., 1.));

  model->AddBody (*model_data, 0, SpatialTransformd(),
      Joint (
        SpatialVectord  (0., 0., 0., 1., 0., 0.),
        SpatialVectord  (0., 0., 0., 0., 1., 0.),
        SpatialVectord  (0., 0., 0., 0., 0., 1.),
        SpatialVectord  (0., 0., 1., 0., 0., 0.),
        SpatialVectord  (0., 1., 0., 0., 0., 0.),
        SpatialVectord  (1., 0., 0., 0., 0., 0.)
        ),
      base_body);

  // Initialization of the input vectors
  VectorNd Q = VectorNd::Constant ((size_t) model->dof_count, 0.);
  VectorNd QDot = VectorNd::Constant ((size_t) model->dof_count, 0.);
  VectorNd QDDot = VectorNd::Constant ((size_t) model->dof_count, 0.);
  VectorNd Tau = VectorNd::Constant ((size_t) model->dof_count, 0.);
  VectorNd TauInv = VectorNd::Constant ((size_t) model->dof_count, 0.);

  MatrixNd H = MatrixNd::Constant ((size_t) model->dof_count, (size_t) model->dof_count, 0.);
  VectorNd C = VectorNd::Constant ((size_t) model->dof_count, 0.);
  VectorNd QDDot_zero = VectorNd::Constant ((size_t) model->dof_count, 0.);
  VectorNd QDDot_crba = VectorNd::Constant ((size_t) model->dof_count, 0.);

  Q[0] = 1.1;
  Q[1] = 1.2;
  Q[2] = 1.3;
  Q[3] = 0.1;
  Q[4] = 0.2;
  Q[5] = 0.3;

  QDot[0] = 1.1;
  QDot[1] = -1.2;
  QDot[2] = 1.3;
  QDot[3] = -0.1;
  QDot[4] = 0.2;
  QDot[5] = -0.3;

  Tau[0] = 2.1;
  Tau[1] = 2.2;
  Tau[2] = 2.3;
  Tau[3] = 1.1;
  Tau[4] = 1.2;
  Tau[5] = 1.3;

  ForwardDynamics(*model, *model_data, Q, QDot, Tau, QDDot);

  ClearLogOutput();
  CompositeRigidBodyAlgorithm (*model, *model_data, Q, H);
  // cout << LogOutput.str() << endl;

  InverseDynamics (*model, *model_data, Q, QDot, QDDot_zero, C);

  EXPECT_TRUE (LinSolveGaussElimPivot (H, C * -1. + Tau, QDDot_crba));

  EXPECT_TRUE(EIGEN_MATRIX_NEAR (QDDot, QDDot_crba, TEST_PREC));
}

TEST_F(FloatingBase12DoF, TestCRBAFloatingBase12DoF) {
  MatrixNd H = MatrixNd::Zero ((size_t) model->dof_count, (size_t) model->dof_count);

  VectorNd C = VectorNd::Constant ((size_t) model->dof_count, 0.);
  VectorNd QDDot_zero = VectorNd::Constant ((size_t) model->dof_count, 0.);
  VectorNd QDDot_crba = VectorNd::Constant ((size_t) model->dof_count, 0.);


  Q[ 0] = 1.1;
  Q[ 1] = 1.2;
  Q[ 2] = 1.3;
  Q[ 3] = 0.1;
  Q[ 4] = 0.2;
  Q[ 5] = 0.3;
  Q[ 6] = -1.3;
  Q[ 7] = -1.4;
  Q[ 8] = -1.5;
  Q[ 9] = -0.3;
  Q[10] = -0.4;
  Q[11] = -0.5;

  QDot[ 0] =  1.1;
  QDot[ 1] = -1.2;
  QDot[ 2] =  1.3;
  QDot[ 3] = -0.1;
  QDot[ 4] =  0.2;
  QDot[ 5] = -0.3;
  QDot[ 6] = -1.1;
  QDot[ 7] =  1.2;
  QDot[ 8] = -1.3;
  QDot[ 9] =  0.1;
  QDot[10] = -0.2;
  QDot[11] =  0.3;

  Tau[ 0] = -1.1;
  Tau[ 1] =  1.2;
  Tau[ 2] = -1.3;
  Tau[ 3] =  1.1;
  Tau[ 4] = -1.2;
  Tau[ 5] =  1.3;
  Tau[ 6] =  0.1;
  Tau[ 7] = -0.2;
  Tau[ 8] =  0.3;
  Tau[ 9] = -0.1;
  Tau[10] =  0.2;
  Tau[11] = -0.3;

  ForwardDynamics(*model, *model_data, Q, QDot, Tau, QDDot);
  ClearLogOutput();
  CompositeRigidBodyAlgorithm (*model, *model_data, Q, H);
  // cout << LogOutput.str() << endl;
  InverseDynamics (*model, *model_data, Q, QDot, QDDot_zero, C);

  EXPECT_TRUE (LinSolveGaussElimPivot (H, C * -1. + Tau, QDDot_crba));

  EXPECT_TRUE(EIGEN_MATRIX_NEAR (QDDot, QDDot_crba,TEST_PREC));
}

TEST_F(FloatingBase12DoF, TestCRBAFloatingBase12DoFInverseDynamics) {
  MatrixNd H_crba = MatrixNd::Zero ((size_t) model->dof_count, (size_t) model->dof_count);
  MatrixNd H_id = MatrixNd::Zero ((size_t) model->dof_count, (size_t) model->dof_count);

  Q[ 0] = 1.1;
  Q[ 1] = 1.2;
  Q[ 2] = 1.3;
  Q[ 3] = 0.1;
  Q[ 4] = 0.2;
  Q[ 5] = 0.3;
  Q[ 6] = -1.3;
  Q[ 7] = -1.4;
  Q[ 8] = -1.5;
  Q[ 9] = -0.3;
  Q[10] = -0.4;
  Q[11] = -0.5;

  QDot.setZero();

  assert (model->dof_count == 12);

  UpdateKinematicsCustom (*model, *model_data, &Q, NULL, NULL);
  CompositeRigidBodyAlgorithm (*model, *model_data, Q, H_crba, false);

  VectorNd H_col = VectorNd::Zero (model->dof_count);
  VectorNd QDDot_zero = VectorNd::Zero (model->dof_count);

  unsigned int i;
  for (i = 0; i < model->dof_count; i++) {
    // compute each column
    VectorNd delta_a = VectorNd::Zero (model->dof_count);
    delta_a[i] = 1.;
    // cout << delta_a << endl;

    // compute ID (model, q, qdot, delta_a)
    VectorNd id_delta = VectorNd::Zero (model->dof_count);
    InverseDynamics (*model, *model_data, Q, QDot, delta_a, id_delta);

    // compute ID (model, q, qdot, zero)
    VectorNd id_zero = VectorNd::Zero (model->dof_count);
    InverseDynamics (*model, *model_data, Q, QDot, QDDot_zero, id_zero);

    H_col = id_delta - id_zero;
    // cout << "H_col = " << H_col << endl;
    H_id.block<12, 1>(0, i) = H_col;
  }

  //	cout << "H (crba) = " << endl << H_crba << endl;
  //	cout << "H (id) = " << endl << H_id << endl;

  EXPECT_TRUE(EIGEN_MATRIX_NEAR (H_crba, H_id, TEST_PREC));
}

TEST_F(FixedBase6DoF, TestCRBAFloatingBase12DoFInverseDynamics) {
  MatrixNd H_crba = MatrixNd::Zero ((size_t) model->dof_count, (size_t) model->dof_count);
  MatrixNd H_id = MatrixNd::Zero ((size_t) model->dof_count, (size_t) model->dof_count);

  Q[ 0] = 1.1;
  Q[ 1] = 1.2;
  Q[ 2] = 1.3;
  Q[ 3] = 0.1;
  Q[ 4] = 0.2;
  Q[ 5] = 0.3;

  QDot.setZero();

  assert (model->dof_count == 6);

  UpdateKinematicsCustom (*model, *model_data, &Q, NULL, NULL);
  CompositeRigidBodyAlgorithm (*model, *model_data, Q, H_crba, false);

  VectorNd H_col = VectorNd::Zero (model->dof_count);
  VectorNd QDDot_zero = VectorNd::Zero (model->dof_count);

  unsigned int i;
  for (i = 0; i < 6; i++) {
    // compute each column
    VectorNd delta_a = VectorNd::Zero (model->dof_count);
    delta_a[i] = 1.;

    ClearLogOutput();
    // compute ID (model, q, qdot, delta_a)
    VectorNd id_delta = VectorNd::Zero (model->dof_count);
    InverseDynamics (*model, *model_data, Q, QDot, delta_a, id_delta);

    // compute ID (model, q, qdot, zero)
    VectorNd id_zero = VectorNd::Zero (model->dof_count);
    InverseDynamics (*model, *model_data, Q, QDot, QDDot_zero, id_zero);

    H_col.setZero();
    H_col = id_delta - id_zero;

    H_id.block<6, 1>(0, i) = H_col;
  }

  EXPECT_TRUE(EIGEN_MATRIX_NEAR (H_crba, H_id, TEST_PREC));
}

TEST_F(CompositeRigidBodyFixture, TestCompositeRigidBodyForwardDynamicsSpherical) {
  Body base_body(1., Vector3d (0., 0., 0.), Vector3d (1., 2., 3.));

  model->AddBody(*model_data, 0, SpatialTransformd(), Joint(JointTypeSpherical), base_body);
  VectorNd Q = VectorNd::Constant ((size_t) model->q_size, 0.);
  model->SetQuaternion (1, Quaterniond(), Q);
  MatrixNd H = MatrixNd::Constant ((size_t) model->qdot_size, (size_t) model->qdot_size, 0.);
  CompositeRigidBodyAlgorithm (*model, *model_data, Q, H, true);

  Matrix3d H_ref (
      1., 0., 0.,
      0., 2., 0.,
      0., 0., 3.
      );

  EXPECT_TRUE(EIGEN_MATRIX_NEAR (H_ref, H, TEST_PREC));
}

int main(int argc, char **argv)
{
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
