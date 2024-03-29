#include <eigen_checks/gtest.h>
#include <gtest/gtest.h>

#include <iostream>

#include "rbdl/rbdl_mathutils.h"
#include "rbdl/Logging.h"

#include "rbdl/Model.h"
#include "rbdl/Kinematics.h"
#include "rbdl/Dynamics.h"

#include "Human36FixtureGTest.h"

using namespace std;
using namespace RigidBodyDynamics;
using namespace RigidBodyDynamics::Math;

const double TEST_PREC = 1.0e-12;

class KinematicsFixture : public ::testing::Test{
protected:
  virtual void SetUp () {
    ClearLogOutput();
    model_data = new ModelDatad;
    model = new Model(*model_data);

    /* Basically a model like this, where X are the Center of Masses
     * and the CoM of the last (3rd) body comes out of the Y=X=0 plane.
     *
     *                X
     *                *
     *              _/
     *            _/  (-Z)
     *      Z    /
     *      *---*
     *      |
     *      |
     *  Z   |
     *  O---*
     *      Y
     */

    body_a = Body (1., Vector3d (1., 0., 0.), Vector3d (1., 1., 1.));
    joint_a = Joint( SpatialVectord (0., 0., 1., 0., 0., 0.));

    body_a_id = model->AddBody(*model_data, 0, Xtrans(Vector3d(0., 0., 0.)), joint_a, body_a);

    body_b = Body (1., Vector3d (0., 1., 0.), Vector3d (1., 1., 1.));
    joint_b = Joint ( SpatialVectord (0., 1., 0., 0., 0., 0.));

    body_b_id = model->AddBody(*model_data, body_a_id, Xtrans(Vector3d(1., 0., 0.)), joint_b, body_b);

    body_c = Body (1., Vector3d (0., 0., 1.), Vector3d (1., 1., 1.));
    joint_c = Joint ( SpatialVectord (0., 0., 1., 0., 0., 0.));

    body_c_id = model->AddBody(*model_data, body_b_id, Xtrans(Vector3d(0., 1., 0.)), joint_c, body_c);

    body_d = Body (1., Vector3d (1., 0., 0.), Vector3d (1., 1., 1.));
    joint_c = Joint ( SpatialVectord (1., 0., 0., 0., 0., 0.));

    body_d_id = model->AddBody(*model_data, body_c_id, Xtrans(Vector3d(0., 0., -1.)), joint_c, body_d);

    Q = VectorNd::Constant ((size_t) model->dof_count, 0.);
    QDot = VectorNd::Constant ((size_t) model->dof_count, 0.);
    QDDot = VectorNd::Constant ((size_t) model->dof_count, 0.);
    Tau = VectorNd::Constant ((size_t) model->dof_count, 0.);

    ClearLogOutput();
  }

  virtual void TearDown () {
    delete model;
  }

  ModelDatad *model_data;
  Model *model;

  unsigned int body_a_id, body_b_id, body_c_id, body_d_id;
  Body body_a, body_b, body_c, body_d;
  Joint joint_a, joint_b, joint_c, joint_d;

  VectorNd Q;
  VectorNd QDot;
  VectorNd QDDot;
  VectorNd Tau;
};

class KinematicsFixture6DoF : public ::testing::Test{
protected:
  virtual void SetUp () {
    ClearLogOutput();

    model_data = new ModelDatad;
    model = new Model(*model_data);

    model->gravity = Vector3d  (0., -9.81, 0.);

    /*
     *
     *          X Contact point (ref child)
     *          |
     *    Base  |
     *   / body |
     *  O-------*
     *           \
     *             Child body
     */

    // base body (3 DoF)
    base = Body (
        1.,
        Vector3d (0.5, 0., 0.),
        Vector3d (1., 1., 1.)
        );
    joint_rotzyx = Joint (
        SpatialVectord (0., 0., 1., 0., 0., 0.),
        SpatialVectord (0., 1., 0., 0., 0., 0.),
        SpatialVectord (1., 0., 0., 0., 0., 0.)
        );
    base_id = model->AddBody (*model_data, 0, Xtrans (Vector3d (0., 0., 0.)), joint_rotzyx, base);

    // child body (3 DoF)
    child = Body (
        1.,
        Vector3d (0., 0.5, 0.),
        Vector3d (1., 1., 1.)
        );
    child_id = model->AddBody (*model_data, base_id, Xtrans (Vector3d (1., 0., 0.)), joint_rotzyx, child);

    Q = VectorNd::Constant (model->mBodies.size() - 1, 0.);
    QDot = VectorNd::Constant (model->mBodies.size() - 1, 0.);
    QDDot = VectorNd::Constant (model->mBodies.size() - 1, 0.);
    Tau = VectorNd::Constant (model->mBodies.size() - 1, 0.);

    ClearLogOutput();
  }

  virtual void TearDown() {
    delete model;
  }

  ModelDatad *model_data;
  Model *model;

  unsigned int base_id, child_id;
  Body base, child;
  Joint joint_rotzyx;

  VectorNd Q;
  VectorNd QDot;
  VectorNd QDDot;
  VectorNd Tau;
};



TEST_F(KinematicsFixture, TestPositionNeutral) {
  // We call ForwardDynamics() as it updates the spatial transformation
  // matrices
  ForwardDynamics(*model, *model_data, Q, QDot, Tau, QDDot);

  Vector3d body_position;

  EXPECT_TRUE(EIGEN_MATRIX_NEAR (Vector3d (0., 0., 0.), CalcBodyToBaseCoordinates(*model, *model_data,  Q, body_a_id, Vector3d (0., 0., 0.), true), TEST_PREC) );
  EXPECT_TRUE(EIGEN_MATRIX_NEAR (Vector3d (1., 0., 0.), CalcBodyToBaseCoordinates(*model, *model_data,  Q, body_b_id, Vector3d (0., 0., 0.), true),  TEST_PREC) );
  EXPECT_TRUE(EIGEN_MATRIX_NEAR (Vector3d (1., 1., 0.), CalcBodyToBaseCoordinates(*model, *model_data,  Q, body_c_id, Vector3d (0., 0., 0.), true),  TEST_PREC) );
  EXPECT_TRUE(EIGEN_MATRIX_NEAR (Vector3d (1., 1., -1.), CalcBodyToBaseCoordinates(*model, *model_data,  Q, body_d_id, Vector3d (0., 0., 0.), true), TEST_PREC) );
}

TEST_F(KinematicsFixture, TestPositionBaseRotated90Deg) {
  // We call ForwardDynamics() as it updates the spatial transformation
  // matrices

  Q[0] = 0.5 * M_PI;
  ForwardDynamics(*model, *model_data, Q, QDot, Tau, QDDot);

  Vector3d body_position;

  //	cout << LogOutput.str() << endl;
  EXPECT_TRUE(EIGEN_MATRIX_NEAR (Vector3d (0., 0., 0.), CalcBodyToBaseCoordinates(*model, *model_data,  Q, body_a_id, Vector3d (0., 0., 0.), true), TEST_PREC) );
  EXPECT_TRUE(EIGEN_MATRIX_NEAR (Vector3d (0., 1., 0.), CalcBodyToBaseCoordinates(*model, *model_data,  Q, body_b_id, Vector3d (0., 0., 0.), true),  TEST_PREC) );
  EXPECT_TRUE(EIGEN_MATRIX_NEAR (Vector3d (-1., 1., 0.),CalcBodyToBaseCoordinates(*model, *model_data,  Q, body_c_id, Vector3d (0., 0., 0.), true),  TEST_PREC) );
  EXPECT_TRUE(EIGEN_MATRIX_NEAR (Vector3d (-1., 1., -1.), CalcBodyToBaseCoordinates(*model, *model_data,  Q, body_d_id, Vector3d (0., 0., 0.), true),  TEST_PREC) );
}

TEST_F(KinematicsFixture, TestPositionBaseRotatedNeg45Deg) {
  // We call ForwardDynamics() as it updates the spatial transformation
  // matrices

  Q[0] = -0.25 * M_PI;
  ForwardDynamics(*model, *model_data, Q, QDot, Tau, QDDot);

  Vector3d body_position;

  //	cout << LogOutput.str() << endl;
  EXPECT_TRUE(EIGEN_MATRIX_NEAR (Vector3d (0., 0., 0.), CalcBodyToBaseCoordinates(*model, *model_data,  Q, body_a_id, Vector3d (0., 0., 0.), true), TEST_PREC) );
  EXPECT_TRUE(EIGEN_MATRIX_NEAR (Vector3d (0.707106781186547, -0.707106781186547, 0.), CalcBodyToBaseCoordinates(*model, *model_data,  Q, body_b_id, Vector3d (0., 0., 0.), true),  TEST_PREC) );
  EXPECT_TRUE(EIGEN_MATRIX_NEAR (Vector3d (sqrt(2.0), 0., 0.),CalcBodyToBaseCoordinates(*model, *model_data,  Q, body_c_id, Vector3d (0., 0., 0.), true), TEST_PREC) );
  EXPECT_TRUE(EIGEN_MATRIX_NEAR (Vector3d (sqrt(2.0), 0., -1.), CalcBodyToBaseCoordinates(*model, *model_data,  Q, body_d_id, Vector3d (0., 0., 0.), true),  TEST_PREC) );
}

TEST_F(KinematicsFixture, TestPositionBodyBRotated90Deg) {
  // We call ForwardDynamics() as it updates the spatial transformation
  // matrices
  Q[1] = 0.5 * M_PI;
  ForwardDynamics(*model, *model_data, Q, QDot, Tau, QDDot);

  Vector3d body_position;

  EXPECT_TRUE(EIGEN_MATRIX_NEAR (Vector3d (0., 0., 0.), CalcBodyToBaseCoordinates(*model, *model_data,  Q, body_a_id, Vector3d (0., 0., 0.), true),  TEST_PREC) );
  EXPECT_TRUE(EIGEN_MATRIX_NEAR (Vector3d (1., 0., 0.), CalcBodyToBaseCoordinates(*model, *model_data,  Q, body_b_id, Vector3d (0., 0., 0.), true),  TEST_PREC) );
  EXPECT_TRUE(EIGEN_MATRIX_NEAR (Vector3d (1., 1., 0.),CalcBodyToBaseCoordinates(*model, *model_data,  Q, body_c_id, Vector3d (0., 0., 0.), true),  TEST_PREC) );
  EXPECT_TRUE(EIGEN_MATRIX_NEAR (Vector3d (0., 1., 0.),CalcBodyToBaseCoordinates(*model, *model_data,  Q, body_d_id, Vector3d (0., 0., 0.), true),  TEST_PREC) );
}

TEST_F(KinematicsFixture, TestPositionBodyBRotatedNeg45Deg) {
  // We call ForwardDynamics() as it updates the spatial transformation
  // matrices
  Q[1] = -0.25 * M_PI;
  ForwardDynamics(*model, *model_data, Q, QDot, Tau, QDDot);

  Vector3d body_position;

  EXPECT_TRUE(EIGEN_MATRIX_NEAR (Vector3d (0., 0., 0.), CalcBodyToBaseCoordinates(*model, *model_data,  Q, body_a_id, Vector3d (0., 0., 0.), true),  TEST_PREC) );
  EXPECT_TRUE(EIGEN_MATRIX_NEAR (Vector3d (1., 0., 0.), CalcBodyToBaseCoordinates(*model, *model_data,  Q, body_b_id, Vector3d (0., 0., 0.), true),  TEST_PREC) );
  EXPECT_TRUE(EIGEN_MATRIX_NEAR (Vector3d (1., 1., 0.),CalcBodyToBaseCoordinates(*model, *model_data,  Q, body_c_id, Vector3d (0., 0., 0.), true),  TEST_PREC) );
  EXPECT_TRUE(EIGEN_MATRIX_NEAR (Vector3d (1 + 0.707106781186547, 1., -0.707106781186547), CalcBodyToBaseCoordinates(*model, *model_data,  Q, body_d_id, Vector3d (0., 0., 0.), true),  TEST_PREC) );
}

TEST_F(KinematicsFixture, TestCalcBodyToBaseCoordinates) {
  // We call ForwardDynamics() as it updates the spatial transformation
  // matrices
  ForwardDynamics(*model, *model_data, Q, QDot, Tau, QDDot);

  EXPECT_TRUE(EIGEN_MATRIX_NEAR (
      Vector3d (1., 2., 0.),
      CalcBodyToBaseCoordinates(*model, *model_data,  Q, body_c_id, Vector3d (0., 1., 0.)),
      TEST_PREC)
      );
}

TEST_F(KinematicsFixture, TestCalcBodyToBaseCoordinatesRotated) {
  Q[2] = 0.5 * M_PI;

  // We call ForwardDynamics() as it updates the spatial transformation
  // matrices
  ForwardDynamics(*model, *model_data, Q, QDot, Tau, QDDot);

  EXPECT_TRUE(EIGEN_MATRIX_NEAR (
      Vector3d (1., 1., 0.),
      CalcBodyToBaseCoordinates(*model, *model_data,  Q, body_c_id, Vector3d (0., 0., 0.), false),
      TEST_PREC)
      );

  EXPECT_TRUE(EIGEN_MATRIX_NEAR (
      Vector3d (0., 1., 0.),
      CalcBodyToBaseCoordinates(*model, *model_data,  Q, body_c_id, Vector3d (0., 1., 0.), false),
       TEST_PREC)
      );

  // Rotate the other way round
  Q[2] = -0.5 * M_PI;

  // We call ForwardDynamics() as it updates the spatial transformation
  // matrices
  ForwardDynamics(*model, *model_data, Q, QDot, Tau, QDDot);

  EXPECT_TRUE(EIGEN_MATRIX_NEAR (
      Vector3d (1., 1., 0.),
      CalcBodyToBaseCoordinates(*model, *model_data,  Q, body_c_id, Vector3d (0., 0., 0.), false),
      TEST_PREC)
      );

  EXPECT_TRUE(EIGEN_MATRIX_NEAR (
      Vector3d (2., 1., 0.),
      CalcBodyToBaseCoordinates(*model, *model_data,  Q, body_c_id, Vector3d (0., 1., 0.), false),
      TEST_PREC)
      );

  // Rotate around the base
  Q[0] = 0.5 * M_PI;
  Q[2] = 0.;

  // We call ForwardDynamics() as it updates the spatial transformation
  // matrices
  ForwardDynamics(*model, *model_data, Q, QDot, Tau, QDDot);

  EXPECT_TRUE(EIGEN_MATRIX_NEAR (
      Vector3d (-1., 1., 0.),
      CalcBodyToBaseCoordinates(*model, *model_data,  Q, body_c_id, Vector3d (0., 0., 0.), false),
      TEST_PREC)
      );

  EXPECT_TRUE(EIGEN_MATRIX_NEAR (
      Vector3d (-2., 1., 0.),
      CalcBodyToBaseCoordinates(*model, *model_data,  Q, body_c_id, Vector3d (0., 1., 0.), false),
       TEST_PREC)
      );

  //	cout << LogOutput.str() << endl;
}

TEST(KinematcisTests, TestCalcPointJacobian) {

  ModelDatad model_data;
  Model model(model_data);
  Body base_body (1., Vector3d (0., 0., 0.), Vector3d (1., 1., 1.));

  unsigned int base_body_id = model.AddBody (model_data, 0, SpatialTransformd(),
      Joint (
        SpatialVectord (0., 0., 0., 1., 0., 0.),
        SpatialVectord (0., 0., 0., 0., 1., 0.),
        SpatialVectord (0., 0., 0., 0., 0., 1.),
        SpatialVectord (0., 0., 1., 0., 0., 0.),
        SpatialVectord (0., 1., 0., 0., 0., 0.),
        SpatialVectord (1., 0., 0., 0., 0., 0.)
        ),
      base_body);

  VectorNd Q = VectorNd::Constant ((size_t) model.dof_count, 0.);
  VectorNd QDot = VectorNd::Constant ((size_t) model.dof_count, 0.);
  MatrixNd G = MatrixNd::Constant (3, model.dof_count, 0.);
  Vector3d point_position (1.1, 1.2, 2.1);
  Vector3d point_velocity_ref;
  Vector3d point_velocity;

  Q[0] = 1.1;
  Q[1] = 1.2;
  Q[2] = 1.3;
  Q[3] = 0.7;
  Q[4] = 0.8;
  Q[5] = 0.9;

  QDot[0] = -1.1;
  QDot[1] = 2.2;
  QDot[2] = 1.3;
  QDot[3] = -2.7;
  QDot[4] = 1.8;
  QDot[5] = -2.9;

  // Compute the reference velocity
  point_velocity_ref = CalcPointVelocity (model, model_data,  Q, QDot, base_body_id, point_position);

  G.setZero();
  CalcPointJacobian (model, model_data, Q, base_body_id, point_position, G);

  point_velocity = G * QDot;

  EXPECT_TRUE(EIGEN_MATRIX_NEAR (
      point_velocity_ref,
      point_velocity,
       TEST_PREC)
      );
}

TEST_F(KinematicsFixture, TestInverseKinematicSimple) {
  std::vector<unsigned int> body_ids;
  std::vector<Vector3d> body_points;
  std::vector<Vector3d> target_pos;

  Q[0] = 0.2;
  Q[1] = 0.1;
  Q[2] = 0.1;

  VectorNd Qres = VectorNd::Zero ((size_t) model->dof_count);

  unsigned int body_id = body_d_id;
  Vector3d body_point = Vector3d (1., 0., 0.);
  Vector3d target (1.3, 0., 0.);

  body_ids.push_back (body_d_id);
  body_points.push_back (body_point);
  target_pos.push_back (target);

  ClearLogOutput();
  bool res = InverseKinematics (*model, *model_data,  Q, body_ids, body_points, target_pos, Qres);
  //	cout << LogOutput.str() << endl;
  EXPECT_TRUE(res);

  UpdateKinematicsCustom<double> (*model, *model_data, &Qres, nullptr, nullptr);

  Vector3d effector;
  effector = CalcBodyToBaseCoordinates(*model, *model_data,  Qres, body_id, body_point, false);

  EXPECT_TRUE(EIGEN_MATRIX_NEAR (target, effector,  TEST_PREC));
}

TEST_F(KinematicsFixture6DoF, TestInverseKinematicUnreachable) {
  std::vector<unsigned int> body_ids;
  std::vector<Vector3d> body_points;
  std::vector<Vector3d> target_pos;

  Q[0] = 0.2;
  Q[1] = 0.1;
  Q[2] = 0.1;

  VectorNd Qres = VectorNd::Zero ((size_t) model->dof_count);

  unsigned int body_id = child_id;
  Vector3d body_point = Vector3d (1., 0., 0.);
  Vector3d target (2.2, 0., 0.);

  body_ids.push_back (body_id);
  body_points.push_back (body_point);
  target_pos.push_back (target);

  ClearLogOutput();
  bool res = InverseKinematics (*model, *model_data, Q, body_ids, body_points, target_pos, Qres, 1.0e-8, 0.9, 1000);
  //	cout << LogOutput.str() << endl;
  EXPECT_TRUE(res);

  UpdateKinematicsCustom<double>(*model, *model_data, &Qres, nullptr, nullptr);

  Vector3d effector;
  effector = CalcBodyToBaseCoordinates(*model, *model_data,  Qres, body_id, body_point, false);

  EXPECT_TRUE(EIGEN_MATRIX_NEAR (Vector3d (2.0, 0., 0.), effector,  1.0e-7));
}

TEST_F(KinematicsFixture6DoF, TestInverseKinematicTwoPoints) {
  std::vector<unsigned int> body_ids;
  std::vector<Vector3d> body_points;
  std::vector<Vector3d> target_pos;

  Q[0] = 0.2;
  Q[1] = 0.1;
  Q[2] = 0.1;

  VectorNd Qres = VectorNd::Zero ((size_t) model->dof_count);

  unsigned int body_id = child_id;
  Vector3d body_point = Vector3d (1., 0., 0.);
  Vector3d target (2., 0., 0.);

  body_ids.push_back (body_id);
  body_points.push_back (body_point);
  target_pos.push_back (target);

  body_ids.push_back (base_id);
  body_points.push_back (Vector3d (0.6, 1.0, 0.));
  target_pos.push_back (Vector3d (0.5, 1.1, 0.));

  ClearLogOutput();
  bool res = InverseKinematics (*model, *model_data, Q, body_ids, body_points, target_pos, Qres, 1.0e-3, 0.9, 200);
  EXPECT_TRUE( res);

  //	cout << LogOutput.str() << endl;
  UpdateKinematicsCustom<double> (*model, *model_data, &Qres, nullptr, nullptr);

  Vector3d effector;

  // testing with very low precision
  effector = CalcBodyToBaseCoordinates(*model, *model_data,  Qres, body_ids[0], body_points[0], false);
  EXPECT_TRUE(EIGEN_MATRIX_NEAR (target_pos[0], effector, 1.0e-1));

  effector = CalcBodyToBaseCoordinates(*model, *model_data,  Qres, body_ids[1], body_points[1], false);
  EXPECT_TRUE(EIGEN_MATRIX_NEAR (target_pos[1], effector,  1.0e-1));
}

TEST ( KinematicsTests, FixedJointBodyCalcBodyToBase ) {
  // the standard modeling using a null body
  Body null_body;
  Body body(1., Vector3d (1., 0.4, 0.4), Vector3d (1., 1., 1.));
  Body fixed_body(1., Vector3d (1., 0.4, 0.4), Vector3d (1., 1., 1.));

  ModelDatad model_data;
  Model model(model_data);

  Joint joint_rot_z ( SpatialVectord (0., 0., 1., 0., 0., 0.));
  model.AddBody (model_data, 0, Xtrans(Vector3d(0., 0., 0.)), joint_rot_z, body);
  unsigned int fixed_body_id = model.AppendBody (model_data, Xtrans(Vector3d(0., 1., 0.)), Joint(JointTypeFixed), fixed_body);

  VectorNd Q_zero = VectorNd::Zero (model.dof_count);
  Vector3d base_coords = CalcBodyToBaseCoordinates (model, model_data,  Q_zero, fixed_body_id, Vector3d (1., 1., 0.1));

  EXPECT_TRUE(EIGEN_MATRIX_NEAR (Vector3d (1., 2., 0.1), base_coords, TEST_PREC));
}

TEST ( KinematicsTests, FixedJointBodyCalcBodyToBaseRotated ) {
  // the standard modeling using a null body
  Body null_body;
  Body body(1., Vector3d (1., 0.4, 0.4), Vector3d (1., 1., 1.));
  Body fixed_body(1., Vector3d (1., 0.4, 0.4), Vector3d (1., 1., 1.));

  ModelDatad model_data;
  Model model(model_data);

  Joint joint_rot_z ( SpatialVectord(0., 0., 1., 0., 0., 0.));
  model.AddBody (model_data, 0, Xtrans(Vector3d(0., 0., 0.)), joint_rot_z, body);
  unsigned int fixed_body_id = model.AppendBody (model_data, Xtrans(Vector3d(1., 0., 0.)), Joint(JointTypeFixed), fixed_body);

  VectorNd Q = VectorNd::Zero (model.dof_count);

  ClearLogOutput();
  Q[0] = M_PI * 0.5;
  Vector3d base_coords = CalcBodyToBaseCoordinates (model, model_data,  Q, fixed_body_id, Vector3d (1., 0., 0.));
  //	cout << LogOutput.str() << endl;

  EXPECT_TRUE(EIGEN_MATRIX_NEAR (Vector3d (0., 2., 0.), base_coords,  TEST_PREC));
}

TEST ( KinematicsTests, FixedJointBodyCalcBaseToBody ) {
  // the standard modeling using a null body
  Body null_body;
  Body body(1., Vector3d (1., 0.4, 0.4), Vector3d (1., 1., 1.));
  Body fixed_body(1., Vector3d (1., 0.4, 0.4), Vector3d (1., 1., 1.));

  ModelDatad model_data;
  Model model(model_data);

  Joint joint_rot_z ( SpatialVectord (0., 0., 1., 0., 0., 0.));
  model.AddBody (model_data, 0, Xtrans(Vector3d(0., 0., 0.)), joint_rot_z, body);
  unsigned int fixed_body_id = model.AppendBody (model_data, Xtrans(Vector3d(0., 1., 0.)), Joint(JointTypeFixed), fixed_body);

  VectorNd Q_zero = VectorNd::Zero (model.dof_count);
  Vector3d base_coords = CalcBaseToBodyCoordinates (model, model_data, Q_zero, fixed_body_id, Vector3d (1., 2., 0.1));

  EXPECT_TRUE(EIGEN_MATRIX_NEAR (Vector3d (1., 1., 0.1), base_coords, TEST_PREC));
}

TEST ( KinematicsTests, FixedJointBodyCalcBaseToBodyRotated ) {
  // the standard modeling using a null body
  Body null_body;
  Body body(1., Vector3d (1., 0.4, 0.4), Vector3d (1., 1., 1.));
  Body fixed_body(1., Vector3d (1., 0.4, 0.4), Vector3d (1., 1., 1.));

  ModelDatad model_data;
  Model model(model_data);

  Joint joint_rot_z ( SpatialVectord (0., 0., 1., 0., 0., 0.));
  model.AddBody (model_data, 0, Xtrans(Vector3d(0., 0., 0.)), joint_rot_z, body);
  unsigned int fixed_body_id = model.AppendBody (model_data, Xtrans(Vector3d(1., 0., 0.)), Joint(JointTypeFixed), fixed_body);

  VectorNd Q = VectorNd::Zero (model.dof_count);

  ClearLogOutput();
  Q[0] = M_PI * 0.5;
  Vector3d base_coords = CalcBaseToBodyCoordinates (model, model_data, Q, fixed_body_id, Vector3d (0., 2., 0.));
  // cout << LogOutput.str() << endl;

  EXPECT_TRUE(EIGEN_MATRIX_NEAR (Vector3d (1., 0., 0.), base_coords,  TEST_PREC));
}

TEST (KinematicsTests, FixedJointBodyWorldOrientation ) {
  // the standard modeling using a null body
  Body null_body;
  Body body(1., Vector3d (1., 0.4, 0.4), Vector3d (1., 1., 1.));
  Body fixed_body(1., Vector3d (1., 0.4, 0.4), Vector3d (1., 1., 1.));

  ModelDatad model_data;
  Model model(model_data);

  Joint joint_rot_z ( SpatialVectord (0., 0., 1., 0., 0., 0.));
  model.AddBody (model_data, 0, Xtrans(Vector3d(0., 0., 0.)), joint_rot_z, body);

  SpatialTransformd transform = Xrotz(0.25) * Xtrans (Vector3d (1., 2., 3.));
  unsigned int fixed_body_id = model.AppendBody (model_data, transform, Joint(JointTypeFixed), fixed_body);

  VectorNd Q_zero = VectorNd::Zero (model.dof_count);
  Matrix3d orientation = CalcBodyWorldOrientation (model, model_data, Q_zero, fixed_body_id);

  Matrix3d reference = transform.E;

  EXPECT_TRUE(EIGEN_MATRIX_NEAR (reference, orientation,  TEST_PREC));
}

TEST ( KinematicsTests, FixedJointCalcPointJacobian ) {
  // the standard modeling using a null body
  Body null_body;
  Body body(1., Vector3d (1., 0.4, 0.4), Vector3d (1., 1., 1.));
  Body fixed_body(1., Vector3d (1., 0.4, 0.4), Vector3d (1., 1., 1.));

  ModelDatad model_data;
  Model model(model_data);

  Joint joint_rot_z ( SpatialVectord (0., 0., 1., 0., 0., 0.));
  model.AddBody (model_data, 0, Xtrans(Vector3d(0., 0., 0.)), joint_rot_z, body);

  SpatialTransformd transform = Xrotz(0.25) * Xtrans (Vector3d (1., 2., 3.));
  unsigned int fixed_body_id = model.AppendBody (model_data, transform, Joint(JointTypeFixed), fixed_body);

  VectorNd Q = VectorNd::Zero (model.dof_count);
  VectorNd QDot = VectorNd::Zero (model.dof_count);

  Q[0] = 1.1;
  QDot[0] = 1.2;

  Vector3d point_position (1., 0., 0.);

  MatrixNd G = MatrixNd::Zero (3, model.dof_count);
  CalcPointJacobian (model, model_data, Q, fixed_body_id, point_position, G);
  Vector3d point_velocity_jacobian = G * QDot;
  Vector3d point_velocity_reference = CalcPointVelocity (model, model_data,  Q, QDot, fixed_body_id, point_position);

  EXPECT_TRUE(EIGEN_MATRIX_NEAR (point_velocity_reference, point_velocity_jacobian, TEST_PREC));
}

TEST_F ( Human36, SpatialJacobianSimple ) {
  randomizeStates();

  unsigned int foot_r_id = model->GetBodyId ("foot_r");
  MatrixNd G (MatrixNd::Zero (6, model->dof_count));

  CalcBodySpatialJacobian (*model, *model_data, q, foot_r_id, G);

  UpdateKinematicsCustom<double> (*model, *model_data, &q, &qdot, nullptr);
  SpatialVectord v_body = SpatialVectord(G * qdot);

  EXPECT_TRUE(EIGEN_MATRIX_NEAR (model_data->v[foot_r_id], v_body,  TEST_PREC));
}

TEST_F ( Human36, SpatialJacobianFixedBody ) {
  randomizeStates();

  unsigned int uppertrunk_id = model->GetBodyId ("uppertrunk");
  MatrixNd G (MatrixNd::Zero (6, model->dof_count));

  CalcBodySpatialJacobian (*model, *model_data, q, uppertrunk_id, G);

  unsigned int fixed_body_id = uppertrunk_id - model->fixed_body_discriminator;
  unsigned int movable_parent = model->mFixedBodies[fixed_body_id].mMovableParent;

  UpdateKinematicsCustom<double>(*model, *model_data, &q, &qdot, nullptr);
  SpatialVectord v_body = SpatialVectord(G * qdot);

  SpatialVectord v_fixed_body = model->mFixedBodies[fixed_body_id].mParentTransform.apply (model_data->v[movable_parent]);

  EXPECT_TRUE(EIGEN_MATRIX_NEAR (v_fixed_body, v_body,  TEST_PREC));
}

TEST_F ( Human36, CalcPointJacobian6D ) {
  randomizeStates();

  unsigned int foot_r_id = model->GetBodyId ("foot_r");
  Vector3d point_local (1.1, 2.2, 3.3);

  // Compute the 6-D velocity using the 6-D Jacobian
  MatrixNd G (MatrixNd::Zero (6, model->dof_count));
  CalcPointJacobian6D (*model, *model_data, q, foot_r_id, point_local, G);
  SpatialVectord v_foot_0_jac = SpatialVectord (G * qdot);

  // Compute the 6-D velocity by transforming the body velocity to the
  // reference point and aligning it with the base coordinate system
  Vector3d r_point = CalcBodyToBaseCoordinates (*model, *model_data, q, foot_r_id, point_local);
  SpatialTransformd X_foot (Matrix3d::Identity(), r_point);
  UpdateKinematicsCustom<double>(*model, *model_data, &q, &qdot, nullptr);
  SpatialVectord v_foot_0_ref = X_foot.apply(model_data->X_base[foot_r_id].inverse().apply(model_data->v[foot_r_id]));

  EXPECT_TRUE(EIGEN_MATRIX_NEAR (v_foot_0_ref, v_foot_0_jac,  TEST_PREC));
}

TEST_F ( Human36, CalcPointVelocity6D ) {
  randomizeStates();

  unsigned int foot_r_id = model->GetBodyId ("foot_r");
  Vector3d point_local (1.1, 2.2, 3.3);

  // Compute the 6-D velocity
  SpatialVectord v_foot_0 = CalcPointVelocity6D (*model, *model_data, q, qdot, foot_r_id, point_local);

  // Compute the 6-D velocity by transforming the body velocity to the
  // reference point and aligning it with the base coordinate system
  Vector3d r_point = CalcBodyToBaseCoordinates (*model, *model_data, q, foot_r_id, point_local);
  SpatialTransformd X_foot (Matrix3d::Identity(), r_point);
  UpdateKinematicsCustom<double>(*model, *model_data, &q, &qdot, nullptr);
  SpatialVectord v_foot_0_ref = X_foot.apply(model_data->X_base[foot_r_id].inverse().apply(model_data->v[foot_r_id]));

  EXPECT_TRUE(EIGEN_MATRIX_NEAR (v_foot_0_ref, v_foot_0,  TEST_PREC));
}

TEST_F ( Human36, CalcPointAcceleration6D ) {
  randomizeStates();

  unsigned int foot_r_id = model->GetBodyId ("foot_r");
  Vector3d point_local (1.1, 2.2, 3.3);

  // Compute the 6-D acceleration
  SpatialVectord a_foot_0 = CalcPointAcceleration6D (*model, *model_data, q, qdot, qddot, foot_r_id, point_local);

  // Compute the 6-D acceleration by adding the coriolis term to the
  // acceleration of the body and transforming the result to the
  // point and align it with the base coordinate system.
  Vector3d r_point = CalcBodyToBaseCoordinates (*model, *model_data, q, foot_r_id, point_local);
  Vector3d v_foot_0 = CalcPointVelocity (*model, *model_data, q, qdot, foot_r_id, point_local);
  SpatialVectord rdot (0., 0., 0., v_foot_0[0], v_foot_0[1], v_foot_0[2]);

  SpatialTransformd X_foot (Matrix3d::Identity(), r_point);
  UpdateKinematicsCustom<double> (*model, *model_data, &q, &qdot, nullptr);
  SpatialVectord a_foot_0_ref = X_foot.apply(
      model_data->X_base[foot_r_id].inverse().apply(model_data->a[foot_r_id])
      - crossm(rdot,
        model_data->X_base[foot_r_id].inverse().apply(model_data->v[foot_r_id])
        )
      );

  EXPECT_TRUE(EIGEN_MATRIX_NEAR (a_foot_0_ref, a_foot_0,  TEST_PREC));
}

int main( int argc, char** argv)
{
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
