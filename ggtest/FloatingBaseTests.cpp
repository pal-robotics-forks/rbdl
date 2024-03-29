
#include <gtest/gtest.h>
#include <eigen_checks/gtest.h>

#include <iostream>

#include "rbdl/rbdl_mathutils.h"
#include "rbdl/Logging.h"

#include "rbdl/Model.h"
#include "rbdl/Kinematics.h"
#include "rbdl/Dynamics.h"

using namespace std;
using namespace RigidBodyDynamics;
using namespace RigidBodyDynamics::Math;

const double TEST_PREC = 1.0e-14;

class FloatingBaseFixture : public ::testing::Test {
protected:
  virtual void SetUp () {
    ClearLogOutput();
    model_data = new ModelDatad;
    model = new Model(*model_data);
    model->gravity = Vector3d (0., -9.81, 0.);

    base = Body (1., Vector3d (1., 0., 0.), Vector3d (1., 1., 1.));

  }
  virtual void TearDown () {
    delete model;
  }

  ModelDatad *model_data;
  Model *model;
  Body base;
  unsigned int base_body_id;

  VectorNd q, qdot, qddot, tau;
};

TEST_F ( FloatingBaseFixture, TestCalcPointTransformation ) {
  base_body_id = model->AddBody (*model_data, 0, SpatialTransformd(),
      Joint (
        SpatialVectord (0., 0., 0., 1., 0., 0.),
        SpatialVectord (0., 0., 0., 0., 1., 0.),
        SpatialVectord (0., 0., 0., 0., 0., 1.),
        SpatialVectord (0., 0., 1., 0., 0., 0.),
        SpatialVectord (0., 1., 0., 0., 0., 0.),
        SpatialVectord (1., 0., 0., 0., 0., 0.)
        ),
      base);

  q = VectorNd::Constant(model->dof_count, 0.);
  qdot = VectorNd::Constant(model->dof_count, 0.);
  qddot = VectorNd::Constant(model->dof_count, 0.);
  tau = VectorNd::Constant(model->dof_count, 0.);

  q[1] = 1.;
  ForwardDynamics (*model, *model_data, q, qdot, tau, qddot);

  Vector3d test_point;

  test_point = CalcBaseToBodyCoordinates (*model, *model_data, q, base_body_id, Vector3d (0., 0., 0.), false);
  EXPECT_TRUE(EIGEN_MATRIX_NEAR (Vector3d (0., -1., 0.), test_point,TEST_PREC));
}

TEST_F(FloatingBaseFixture, TestCalcDynamicFloatingBaseDoubleImplicit) {
  // floating base
  base_body_id = model->AddBody (*model_data, 0, SpatialTransformd(),
      Joint (
        SpatialVectord (0., 0., 0., 1., 0., 0.),
        SpatialVectord (0., 0., 0., 0., 1., 0.),
        SpatialVectord (0., 0., 0., 0., 0., 1.),
        SpatialVectord (0., 0., 1., 0., 0., 0.),
        SpatialVectord (0., 1., 0., 0., 0., 0.),
        SpatialVectord (1., 0., 0., 0., 0., 0.)
        ),
      base);

  // body_a
  Body body_a (1., Vector3d (1., 0., 0), Vector3d (1., 1., 1.));
  Joint joint_a ( SpatialVectord (0., 0., 1., 0., 0., 0.));

  model->AddBody(*model_data, base_body_id, Xtrans(Vector3d(2., 0., 0.)), joint_a, body_a);

  // Initialization of the input vectors
  VectorNd Q = VectorNd::Zero ((size_t) model->dof_count);
  VectorNd QDot = VectorNd::Zero  ((size_t) model->dof_count);
  VectorNd QDDot = VectorNd::Zero  ((size_t) model->dof_count);
  VectorNd Tau = VectorNd::Zero  ((size_t) model->dof_count);

  ForwardDynamics(*model, *model_data, Q, QDot, Tau, QDDot);

  unsigned int i;
  for (i = 0; i < QDDot.size(); i++) {
    LOG << "QDDot[" << i << "] = " << QDDot[i] << endl;
  }

  for (i = 0; i < model_data->a.size(); i++) {
    LOG << "a[" << i << "]     = " << model_data->a.at(i) << endl;
  }

  //	std::cout << LogOutput.str() << std::endl;

  EXPECT_NEAR ( 0.0000, QDDot[0], TEST_PREC);
  EXPECT_NEAR (-9.8100, QDDot[1], TEST_PREC);
  EXPECT_NEAR ( 0.0000, QDDot[2], TEST_PREC);
  EXPECT_NEAR ( 0.0000, QDDot[3], TEST_PREC);
  EXPECT_NEAR ( 0.0000, QDDot[4], TEST_PREC);
  EXPECT_NEAR ( 0.0000, QDDot[5], TEST_PREC);
  EXPECT_NEAR ( 0.0000, QDDot[6], TEST_PREC);

  // We rotate the base... let's see what happens...
  Q[3] = 0.8;
  ForwardDynamics(*model, *model_data, Q, QDot, Tau, QDDot);

  for (i = 0; i < QDDot.size(); i++) {
    LOG << "QDDot[" << i << "] = " << QDDot[i] << endl;
  }

  for (i = 0; i < model_data->a.size(); i++) {
    LOG << "a[" << i << "]     = " << model_data->a.at(i) << endl;
  }

  //	std::cout << LogOutput.str() << std::endl;

  EXPECT_NEAR ( 0.0000, QDDot[0], TEST_PREC);
  EXPECT_NEAR (-9.8100, QDDot[1], TEST_PREC);
  EXPECT_NEAR ( 0.0000, QDDot[2], TEST_PREC);
  EXPECT_NEAR ( 0.0000, QDDot[3], TEST_PREC);
  EXPECT_NEAR ( 0.0000, QDDot[4], TEST_PREC);
  EXPECT_NEAR ( 0.0000, QDDot[5], TEST_PREC);
  EXPECT_NEAR ( 0.0000, QDDot[6], TEST_PREC);

  // We apply a torqe let's see what happens...
  Q[3] = 0.;
  /*
     rot_B[0] = 0.0;
     X_B = XtransRotZYXEuler(pos_B, rot_B);
     */

  Tau[6] = 1.;

  ForwardDynamics(*model, *model_data, Q, QDot, Tau, QDDot);

  for (i = 0; i < QDDot.size(); i++) {
    LOG << "QDDot[" << i << "] = " << QDDot[i] << endl;
  }

  for (i = 0; i < model_data->a.size(); i++) {
    LOG << "a[" << i << "]     = " << model_data->a.at(i) << endl;
  }

  //	std::cout << LogOutput.str() << std::endl;

  EXPECT_NEAR ( 0.0000, QDDot[0], TEST_PREC);
  EXPECT_NEAR (-8.8100, QDDot[1], TEST_PREC);
  EXPECT_NEAR ( 0.0000, QDDot[2], TEST_PREC);
  EXPECT_NEAR (-1.0000, QDDot[3], TEST_PREC);
  EXPECT_NEAR ( 0.0000, QDDot[4], TEST_PREC);
  EXPECT_NEAR ( 0.0000, QDDot[5], TEST_PREC);
  EXPECT_NEAR ( 2.0000, QDDot[6], TEST_PREC);
}

TEST_F(FloatingBaseFixture, TestCalcPointVelocityFloatingBaseSimple) {
  // floating base
  base_body_id = model->AddBody (*model_data, 0, SpatialTransformd(),
      Joint (
        SpatialVectord (0., 0., 0., 1., 0., 0.),
        SpatialVectord (0., 0., 0., 0., 1., 0.),
        SpatialVectord (0., 0., 0., 0., 0., 1.),
        SpatialVectord (0., 0., 1., 0., 0., 0.),
        SpatialVectord (0., 1., 0., 0., 0., 0.),
        SpatialVectord (1., 0., 0., 0., 0., 0.)
        ),
      base);

  VectorNd Q = VectorNd::Zero (model->dof_count);
  VectorNd QDot = VectorNd::Zero (model->dof_count);
  VectorNd QDDot = VectorNd::Zero (model->dof_count);
  VectorNd Tau = VectorNd::Zero (model->dof_count);

  unsigned int ref_body_id = base_body_id;

  // first we calculate the velocity when moving along the X axis
  QDot[0] = 1.;
  Vector3d point_position(1., 0., 0.);
  Vector3d point_velocity;

  point_velocity = CalcPointVelocity(*model, *model_data, Q, QDot, ref_body_id, point_position);

  EXPECT_NEAR(1., point_velocity[0], TEST_PREC);
  EXPECT_NEAR(0., point_velocity[1], TEST_PREC);
  EXPECT_NEAR(0., point_velocity[2], TEST_PREC);

  LOG << "Point velocity = " << point_velocity << endl;
  //	cout << LogOutput.str() << endl;

  ClearLogOutput();

  // Now we calculate the velocity when rotating around the Z axis
  QDot[0] = 0.;
  QDot[3] = 1.;

  point_velocity = CalcPointVelocity(*model, *model_data, Q, QDot, ref_body_id, point_position);

  EXPECT_NEAR(0., point_velocity[0], TEST_PREC);
  EXPECT_NEAR(1., point_velocity[1], TEST_PREC);
  EXPECT_NEAR(0., point_velocity[2], TEST_PREC);

  LOG << "Point velocity = " << point_velocity << endl;
  //	cout << LogOutput.str() << endl;

  // Now we calculate the velocity when rotating around the Z axis and the
  // base is rotated around the z axis by 90 degrees
  ClearLogOutput();
  Q[3] = M_PI * 0.5;
  QDot[3] = 1.;

  point_velocity = CalcPointVelocity(*model, *model_data, Q, QDot, ref_body_id, point_position);

  EXPECT_NEAR(-1., point_velocity[0], TEST_PREC);
  EXPECT_NEAR(0., point_velocity[1], TEST_PREC);
  EXPECT_NEAR(0., point_velocity[2], TEST_PREC);

  LOG << "Point velocity = " << point_velocity << endl;
  //	cout << LogOutput.str() << endl;
}

TEST_F(FloatingBaseFixture, TestCalcPointVelocityCustom) {
  // floating base
  base = Body (1., Vector3d (0., 1., 0.), Vector3d (1., 1., 1.));
  base_body_id = model->AddBody (*model_data, 0, SpatialTransformd(),
      Joint (
        SpatialVectord (0., 0., 0., 1., 0., 0.),
        SpatialVectord (0., 0., 0., 0., 1., 0.),
        SpatialVectord (0., 0., 0., 0., 0., 1.),
        SpatialVectord (0., 0., 1., 0., 0., 0.),
        SpatialVectord (0., 1., 0., 0., 0., 0.),
        SpatialVectord (1., 0., 0., 0., 0., 0.)
        ),
      base);

  VectorNd Q = VectorNd::Zero (model->dof_count);
  VectorNd QDot = VectorNd::Zero (model->dof_count);
  VectorNd QDDot = VectorNd::Zero (model->dof_count);
  VectorNd Tau = VectorNd::Zero (model->dof_count);

  unsigned int ref_body_id = base_body_id;

  Q[0] = 0.1;
  Q[1] = 1.1;
  Q[2] = 1.2;
  Q[3] = 1.3;
  Q[4] = 1.5;
  Q[5] = 1.7;

  QDot[0] = 0.1;
  QDot[1] = 1.1;
  QDot[2] = 1.2;
  QDot[3] = 1.3;
  QDot[4] = 1.5;
  QDot[5] = 1.7;

  // first we calculate the velocity when rotating around the Z axis
  Vector3d point_body_position (1., 0., 0.);
  Vector3d point_base_position;
  Vector3d point_base_velocity;
  Vector3d point_base_velocity_reference;

  ForwardDynamics(*model, *model_data, Q, QDot, Tau, QDDot);

  point_base_velocity = CalcPointVelocity (*model, *model_data, Q, QDot, ref_body_id, point_body_position);

  point_base_velocity_reference = Vector3d (
      -3.888503432977729e-01,
      -3.171179347202455e-01,
      1.093894197498446e+00
      );

  EXPECT_TRUE(EIGEN_MATRIX_NEAR (point_base_velocity_reference, point_base_velocity, TEST_PREC));
}

/** \brief Compares computation of acceleration values for zero qddot
 *
 * Ensures that computation of position, velocity, and acceleration of a
 * point produce the same values as in an equivalent model that was
 * created with the HuMAnS toolbox
 *    http://www.inrialpes.fr/bipop/software/humans/ .
 * Here we omit the term of the generalized acceleration by setting qddot
 * to zero.
 */
TEST_F(FloatingBaseFixture, TestCalcPointAccelerationNoQDDot) {
  // floating base
  base = Body (1., Vector3d (0., 1., 0.), Vector3d (1., 1., 1.));
  base_body_id = model->AddBody (*model_data, 0, SpatialTransformd(),
      Joint (
        SpatialVectord (0., 0., 0., 1., 0., 0.),
        SpatialVectord (0., 0., 0., 0., 1., 0.),
        SpatialVectord (0., 0., 0., 0., 0., 1.),
        SpatialVectord (0., 0., 1., 0., 0., 0.),
        SpatialVectord (0., 1., 0., 0., 0., 0.),
        SpatialVectord (1., 0., 0., 0., 0., 0.)
        ),
      base);

  VectorNd Q = VectorNd::Zero (model->dof_count);
  VectorNd QDot = VectorNd::Zero (model->dof_count);
  VectorNd QDDot = VectorNd::Zero (model->dof_count);
  VectorNd Tau = VectorNd::Zero (model->dof_count);

  unsigned int ref_body_id = base_body_id;

  Q[0] = 0.1;
  Q[1] = 1.1;
  Q[2] = 1.2;
  Q[3] = 1.3;
  Q[4] = 1.5;
  Q[5] = 1.7;

  QDot[0] = 0.1;
  QDot[1] = 1.1;
  QDot[2] = 1.2;
  QDot[3] = 1.3;
  QDot[4] = 1.5;
  QDot[5] = 1.7;

  // first we calculate the velocity when rotating around the Z axis
  Vector3d point_body_position (-1.9, -1.8, 0.);
  Vector3d point_world_position;
  Vector3d point_world_velocity;
  Vector3d point_world_acceleration;

  // call ForwardDynamics to update the model
  ForwardDynamics(*model, *model_data, Q, QDot, Tau, QDDot);
  QDDot = VectorNd::Zero (QDDot.size());

  QDot = QDot;

  point_world_position = CalcBodyToBaseCoordinates (*model, *model_data, Q, ref_body_id, point_body_position, false);
  point_world_velocity = CalcPointVelocity (*model, *model_data, Q, QDot, ref_body_id, point_body_position);

  // we set the generalized acceleration to zero

  ClearLogOutput();

  point_world_acceleration = CalcPointAcceleration (*model, *model_data, Q, QDot, QDDot, ref_body_id, point_body_position);

  Vector3d humans_point_position (
      -6.357089363622626e-01, -6.831041744630977e-01, 2.968974805916970e+00
      );
  Vector3d humans_point_velocity (
      3.091226260907569e-01, 3.891012095550828e+00, 4.100277995030419e+00
      );
  Vector3d humans_point_acceleration (
      -5.302760158847160e+00, 6.541369639625232e+00, -4.795115077652286e+00
      );

  //	cout << LogOutput.str() << endl;
  //
  //	cout << "q     = " << q << endl;
  //	cout << "qdot  = " << qdot << endl;
  //	cout << "qddot = " << qddot << endl;
  //
  //	cout << "body_coords = " << point_body_position << endl;
  //	cout << "world_pos   = " << point_world_position << endl;
  //	cout << "world_vel   = " << point_world_velocity << endl;
  //	cout << "world_accel = " << point_world_acceleration << endl;


  EXPECT_TRUE(EIGEN_MATRIX_NEAR (humans_point_position, point_world_position, TEST_PREC));
  EXPECT_TRUE(EIGEN_MATRIX_NEAR (humans_point_velocity, point_world_velocity, TEST_PREC));
  EXPECT_TRUE(EIGEN_MATRIX_NEAR (humans_point_acceleration, point_world_acceleration, TEST_PREC));
}

/** \brief Compares computation of acceleration values for zero q and qdot
 *
 * Ensures that computation of position, velocity, and acceleration of a
 * point produce the same values as in an equivalent model that was
 * created with the HuMAnS toolbox
 *    http://www.inrialpes.fr/bipop/software/humans/ .
 *
 * Here we set q and qdot to zero and only take into account values that
 * are dependent on qddot.
 */
TEST_F(FloatingBaseFixture, TestCalcPointAccelerationOnlyQDDot) {
  // floating base
  base = Body (1., Vector3d (0., 1., 0.), Vector3d (1., 1., 1.));
  base_body_id = model->AddBody (*model_data, 0, SpatialTransformd(),
      Joint (
        SpatialVectord (0., 0., 0., 1., 0., 0.),
        SpatialVectord (0., 0., 0., 0., 1., 0.),
        SpatialVectord (0., 0., 0., 0., 0., 1.),
        SpatialVectord (0., 0., 1., 0., 0., 0.),
        SpatialVectord (0., 1., 0., 0., 0., 0.),
        SpatialVectord (1., 0., 0., 0., 0., 0.)
        ),
      base);

  VectorNd Q = VectorNd::Zero (model->dof_count);
  VectorNd QDot = VectorNd::Zero (model->dof_count);
  VectorNd QDDot = VectorNd::Zero (model->dof_count);
  VectorNd Tau = VectorNd::Zero (model->dof_count);

  unsigned int ref_body_id = base_body_id;

  // first we calculate the velocity when rotating around the Z axis
  Vector3d point_body_position (-1.9, -1.8, 0.);
  Vector3d point_world_position;
  Vector3d point_world_velocity;
  Vector3d point_world_acceleration;

  ForwardDynamics(*model, *model_data, Q, QDot, Tau, QDDot);

  QDDot = VectorNd::Zero (QDDot.size());

  QDDot[0] = 0.1;
  QDDot[1] = 1.1;
  QDDot[2] = 1.2;
  QDDot[3] = 1.3;
  QDDot[4] = 1.5;
  QDDot[5] = 1.7;

  //	cout << "ref_body_id = " << ref_body_id << endl;
  //	cout << "point_body_position = " << point_body_position << endl;
  point_world_position = CalcBodyToBaseCoordinates (*model, *model_data, Q, ref_body_id, point_body_position, false);
  point_world_velocity = CalcPointVelocity (*model, *model_data, Q, QDot, ref_body_id, point_body_position);

  ClearLogOutput();

  point_world_acceleration = CalcPointAcceleration (*model, *model_data, Q, QDot, QDDot, ref_body_id, point_body_position);

  Vector3d humans_point_position (
      -1.900000000000000e+00, -1.800000000000000e+00, 0.000000000000000e+00
      );
  Vector3d humans_point_velocity (
      0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00
      );
  Vector3d humans_point_acceleration (
      2.440000000000000e+00, -1.370000000000000e+00, 9.899999999999999e-01
      );

  //	cout << LogOutput.str() << endl;
  //
  //	cout << "q     = " << q << endl;
  //	cout << "qdot  = " << qdot << endl;
  //	cout << "qddot = " << qddot << endl;
  //
  //	cout << "body_coords = " << point_body_position << endl;
  //	cout << "world_pos   = " << point_world_position << endl;
  //	cout << "world_vel   = " << point_world_velocity << endl;
  //	cout << "world_accel = " << point_world_acceleration << endl;

  EXPECT_TRUE(EIGEN_MATRIX_NEAR (humans_point_position, point_world_position,  TEST_PREC));
  EXPECT_TRUE(EIGEN_MATRIX_NEAR (humans_point_velocity, point_world_velocity,  TEST_PREC));
  EXPECT_TRUE(EIGEN_MATRIX_NEAR (humans_point_acceleration, point_world_acceleration,  TEST_PREC));
}

/** \brief Compares computation of acceleration values for zero q and qdot
 *
 * Ensures that computation of position, velocity, and acceleration of a
 * point produce the same values as in an equivalent model that was
 * created with the HuMAnS toolbox
 *    http://www.inrialpes.fr/bipop/software/humans/ .
 *
 * Here we set q and qdot to zero and only take into account values that
 * are dependent on qddot.
 */
TEST_F(FloatingBaseFixture, TestCalcPointAccelerationFull) {
  // floating base
  base = Body (1., Vector3d (0., 1., 0.), Vector3d (1., 1., 1.));
  base_body_id = model->AddBody (*model_data, 0, SpatialTransformd(),
      Joint (
        SpatialVectord (0., 0., 0., 1., 0., 0.),
        SpatialVectord (0., 0., 0., 0., 1., 0.),
        SpatialVectord (0., 0., 0., 0., 0., 1.),
        SpatialVectord (0., 0., 1., 0., 0., 0.),
        SpatialVectord (0., 1., 0., 0., 0., 0.),
        SpatialVectord (1., 0., 0., 0., 0., 0.)
        ),
      base);

  VectorNd Q = VectorNd::Zero (model->dof_count);
  VectorNd QDot = VectorNd::Zero (model->dof_count);
  VectorNd QDDot = VectorNd::Zero (model->dof_count);
  VectorNd Tau = VectorNd::Zero (model->dof_count);

  unsigned int ref_body_id = base_body_id;

  // first we calculate the velocity when rotating around the Z axis
  Vector3d point_body_position (-1.9, -1.8, 0.);
  Vector3d point_world_position;
  Vector3d point_world_velocity;
  Vector3d point_world_acceleration;

  Q[0] = 0.1;
  Q[1] = 1.1;
  Q[2] = 1.2;
  Q[3] = 1.3;
  Q[4] = 1.5;
  Q[5] = 1.7;

  QDot[0] = 0.1;
  QDot[1] = 1.1;
  QDot[2] = 1.2;
  QDot[3] = 1.3;
  QDot[4] = 1.5;
  QDot[5] = 1.7;

  ForwardDynamics(*model, *model_data, Q, QDot, Tau, QDDot);

  QDDot[0] = 0.1;
  QDDot[1] = 1.1;
  QDDot[2] = 1.2;
  QDDot[3] = 1.3;
  QDDot[4] = 1.5;
  QDDot[5] = 1.7;

  //	cout << "ref_body_id = " << ref_body_id << endl;
  //	cout << "point_body_position = " << point_body_position << endl;
  point_world_position = CalcBodyToBaseCoordinates (*model, *model_data, Q, ref_body_id, point_body_position, false);
  point_world_velocity = CalcPointVelocity (*model, *model_data, Q, QDot, ref_body_id, point_body_position);

  ClearLogOutput();

  point_world_acceleration = CalcPointAcceleration (*model, *model_data, Q, QDot, QDDot, ref_body_id, point_body_position);

  Vector3d humans_point_position (
      -6.357089363622626e-01, -6.831041744630977e-01, 2.968974805916970e+00
      );
  Vector3d humans_point_velocity (
      3.091226260907569e-01, 3.891012095550828e+00, 4.100277995030419e+00
      );
  Vector3d humans_point_acceleration (
      -4.993637532756404e+00, 1.043238173517606e+01, -6.948370826218673e-01
      );

  //	cout << LogOutput.str() << endl;
  //
  //	cout << "q     = " << q << endl;
  //	cout << "qdot  = " << qdot << endl;
  //	cout << "qddot = " << qddot << endl;
  //
  //	cout << "body_coords = " << point_body_position << endl;
  //	cout << "world_pos   = " << point_world_position << endl;
  //	cout << "world_vel   = " << point_world_velocity << endl;
  //	cout << "world_accel = " << point_world_acceleration << endl;

  EXPECT_TRUE(EIGEN_MATRIX_NEAR (humans_point_position, point_world_position,  TEST_PREC));
  EXPECT_TRUE(EIGEN_MATRIX_NEAR (humans_point_velocity, point_world_velocity,  TEST_PREC));
  EXPECT_TRUE(EIGEN_MATRIX_NEAR (humans_point_acceleration, point_world_acceleration, TEST_PREC));
}


int main( int argc, char** argv)
{
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
