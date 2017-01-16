#include <gtest/gtest.h>
#include <eigen_checks/gtest.h>

#include <iostream>

#include "rbdl/Logging.h"

#include "rbdl/Model.h"
#include "rbdl/Kinematics.h"

using namespace std;
using namespace RigidBodyDynamics;
using namespace RigidBodyDynamics::Math;

const double TEST_PREC = 1.0e-14;

class ModelVelocitiesFixture : public ::testing::Test {

protected:

  virtual void SetUp () {
    ClearLogOutput();
    model_data = new ModelData;
    model = new Model(*model_data);

    body_a = Body (1., Vector3d (1., 0., 0.), Vector3d (1., 1., 1.));
    Joint joint_a ( SpatialVectord (0., 0., 1., 0., 0., 0.));

    body_a_id = model->AddBody(*model_data, 0, Xtrans(Vector3d(0., 0., 0.)), joint_a, body_a);

    body_b = Body (1., Vector3d (0., 1., 0.), Vector3d (1., 1., 1.));
    Joint joint_b ( SpatialVectord (0., 1., 0., 0., 0., 0.));

    body_b_id = model->AddBody(*model_data, 1, Xtrans(Vector3d(1., 0., 0.)), joint_b, body_b);

    body_c = Body (1., Vector3d (1., 0., 0.), Vector3d (1., 1., 1.));
    Joint joint_c ( SpatialVectord (1., 0., 0., 0., 0., 0.));

    body_c_id = model->AddBody(*model_data, 2, Xtrans(Vector3d(0., 1., 0.)), joint_c, body_c);

    Q = VectorNd::Constant ((size_t) model->dof_count, 0.);
    QDot = VectorNd::Constant ((size_t) model->dof_count, 0.);

    point_position = Vector3d::Zero(3);
    point_velocity = Vector3d::Zero(3);

    ref_body_id = 0;

    ClearLogOutput();
  }
  virtual void TearDown () {
    delete model;
  }

  ModelData *model_data;
  Model *model;

  unsigned int body_a_id, body_b_id, body_c_id, ref_body_id;
  Body body_a, body_b, body_c;
  Joint joint_a, joint_b, joint_c;

  VectorNd Q;
  VectorNd QDot;

  Vector3d point_position, point_velocity;
};

TEST_F(ModelVelocitiesFixture, TestCalcPointSimple) {
  ref_body_id = 1;
  QDot[0] = 1.;
  point_position = Vector3d (1., 0., 0.);
  point_velocity = CalcPointVelocity(*model, *model_data, Q, QDot, ref_body_id, point_position);

  EXPECT_NEAR(0., point_velocity[0], TEST_PREC);
  EXPECT_NEAR(1., point_velocity[1], TEST_PREC);
  EXPECT_NEAR(0., point_velocity[2], TEST_PREC);

  LOG << "Point velocity = " << point_velocity << endl;
  //	cout << LogOutput.str() << endl;
}

TEST_F(ModelVelocitiesFixture, TestCalcPointRotatedBaseSimple) {
  // rotated first joint

  ref_body_id = 1;
  Q[0] = M_PI * 0.5;
  QDot[0] = 1.;
  point_position = Vector3d (1., 0., 0.);
  point_velocity = CalcPointVelocity(*model, *model_data, Q, QDot, ref_body_id, point_position);

  EXPECT_NEAR(-1., point_velocity[0], TEST_PREC);
  EXPECT_NEAR( 0., point_velocity[1], TEST_PREC);
  EXPECT_NEAR( 0., point_velocity[2], TEST_PREC);

  //	cout << LogOutput.str() << endl;
}

TEST_F(ModelVelocitiesFixture, TestCalcPointRotatingBodyB) {
  // rotating second joint, point at third body

  ref_body_id = 3;
  QDot[1] = 1.;
  point_position = Vector3d (1., 0., 0.);
  point_velocity = CalcPointVelocity(*model, *model_data, Q, QDot, ref_body_id, point_position);

  //	cout << LogOutput.str() << endl;

  EXPECT_NEAR( 0., point_velocity[0], TEST_PREC);
  EXPECT_NEAR( 0., point_velocity[1], TEST_PREC);
  EXPECT_NEAR(-1., point_velocity[2], TEST_PREC);
}

TEST_F(ModelVelocitiesFixture, TestCalcPointRotatingBaseXAxis) {
  // also rotate the first joint and take a point that is
  // on the X direction

  ref_body_id = 3;
  QDot[0] = 1.;
  QDot[1] = 1.;
  point_position = Vector3d (1., -1., 0.);
  point_velocity = CalcPointVelocity(*model, *model_data, Q, QDot, ref_body_id, point_position);

  //	cout << LogOutput.str() << endl;

  EXPECT_NEAR( 0., point_velocity[0], TEST_PREC);
  EXPECT_NEAR( 2., point_velocity[1], TEST_PREC);
  EXPECT_NEAR(-1., point_velocity[2], TEST_PREC);
}

TEST_F(ModelVelocitiesFixture, TestCalcPointRotatedBaseXAxis) {
  // perform the previous test with the first joint rotated by pi/2
  // upwards
  ClearLogOutput();

  ref_body_id = 3;
  point_position = Vector3d (1., -1., 0.);

  Q[0] = M_PI * 0.5;
  QDot[0] = 1.;
  QDot[1] = 1.;
  point_velocity = CalcPointVelocity(*model, *model_data, Q, QDot, ref_body_id, point_position);

  //	cout << LogOutput.str() << endl;

  EXPECT_NEAR(-2., point_velocity[0], TEST_PREC);
  EXPECT_NEAR( 0., point_velocity[1], TEST_PREC);
  EXPECT_NEAR(-1., point_velocity[2], TEST_PREC);
}

TEST_F(ModelVelocitiesFixture, TestCalcPointBodyOrigin) {
  // Checks whether the computation is also correct for points at the origin
  // of a body

  ref_body_id = body_b_id;
  point_position = Vector3d (0., 0., 0.);

  Q[0] = 0.;
  QDot[0] = 1.;

  point_velocity = CalcPointVelocity(*model, *model_data, Q, QDot, ref_body_id, point_position);

  // cout << LogOutput.str() << endl;

  EXPECT_NEAR( 0., point_velocity[0], TEST_PREC);
  EXPECT_NEAR( 1., point_velocity[1], TEST_PREC);
  EXPECT_NEAR( 0., point_velocity[2], TEST_PREC);
}

TEST(CalcVelocitiesTest, FixedJointCalcPointVelocity ) {
  // the standard modeling using a null body
  Body body(1., Vector3d (1., 0.4, 0.4), Vector3d (1., 1., 1.));
  Body fixed_body(1., Vector3d (1., 0.4, 0.4), Vector3d (1., 1., 1.));

  ModelData model_data;
  Model model(model_data);

  Joint joint_rot_z ( SpatialVectord (0., 0., 1., 0., 0., 0.));
  model.AddBody (model_data, 0, Xtrans(Vector3d(0., 0., 0.)), joint_rot_z, body);

  SpatialTransformd  transform = Xtrans (Vector3d (1., 0., 0.));
  unsigned int fixed_body_id = model.AppendBody (model_data, transform, Joint(JointTypeFixed), fixed_body, "fixed_body");

  VectorNd Q = VectorNd::Zero (model.dof_count);
  VectorNd QDot = VectorNd::Zero (model.dof_count);

  QDot[0] = 1.;

  ClearLogOutput();
  Vector3d point0_velocity = CalcPointVelocity (model, model_data, Q, QDot, fixed_body_id, Vector3d (0., 0., 0.));
  // cout << LogOutput.str() << endl;
  Vector3d point1_velocity = CalcPointVelocity (model, model_data, Q, QDot, fixed_body_id, Vector3d (1., 0., 0.));

   EXPECT_TRUE(EIGEN_MATRIX_NEAR (Vector3d (0., 1., 0.), point0_velocity, TEST_PREC));
   EXPECT_TRUE(EIGEN_MATRIX_NEAR (Vector3d (0., 2., 0.), point1_velocity, TEST_PREC));
}

TEST (CalcVelocitiesTest, FixedJointCalcPointVelocityRotated ) {
  // the standard modeling using a null body
  Body body(1., Vector3d (1., 0.4, 0.4), Vector3d (1., 1., 1.));
  Body fixed_body(1., Vector3d (1., 0.4, 0.4), Vector3d (1., 1., 1.));

  ModelData model_data;
  Model model(model_data);

  Joint joint_rot_z ( SpatialVectord (0., 0., 1., 0., 0., 0.));
  model.AddBody (model_data, 0, Xtrans(Vector3d(0., 0., 0.)), joint_rot_z, body);

  SpatialTransformd  transform = Xtrans (Vector3d (1., 0., 0.));
  unsigned int fixed_body_id = model.AppendBody (model_data, transform, Joint(JointTypeFixed), fixed_body, "fixed_body");

  VectorNd Q = VectorNd::Zero (model.dof_count);
  VectorNd QDot = VectorNd::Zero (model.dof_count);

  Q[0] = M_PI * 0.5;
  QDot[0] = 1.;

  ClearLogOutput();
  Vector3d point0_velocity = CalcPointVelocity (model, model_data, Q, QDot, fixed_body_id, Vector3d (0., 0., 0.));
  // cout << LogOutput.str() << endl;
  Vector3d point1_velocity = CalcPointVelocity (model, model_data, Q, QDot, fixed_body_id, Vector3d (1., 0., 0.));

   EXPECT_TRUE(EIGEN_MATRIX_NEAR (Vector3d (-1., 0., 0.), point0_velocity, TEST_PREC));
   EXPECT_TRUE(EIGEN_MATRIX_NEAR (Vector3d (-2., 0., 0.), point1_velocity, TEST_PREC));
}

int main( int argc, char** argv)
{
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
