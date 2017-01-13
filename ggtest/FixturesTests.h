#ifndef FIXTURESTESTS_H
#define FIXTURESTESTS_H

#include "rbdl/rbdl.h"
#include <gtest/gtest.h>

class FixedBase3DoF : public ::testing::Test {
protected:
  virtual void SetUp() {
    using namespace RigidBodyDynamics;
    using namespace RigidBodyDynamics::Math;

    ClearLogOutput();
    model_data = new ModelData;
    model = new Model(*model_data);

    /* Basically a model like this, where X are the Center of Masses
     * and the CoM of the last (3rd) body comes out of the Y=X=0 plane.
     *
     *      Z
     *      *---*
     *      |
     *      |
     *  Z   |
     *  O---*
     *      Y
     */

    body_a = Body (1., Vector3d (1., 0., 0.), Vector3d (1., 1., 1.));
    joint_a = Joint( SpatialVectord(0., 0., 1., 0., 0., 0.));

    body_a_id = model->AddBody(*model_data, 0, Xtrans(Vector3d(0., 0., 0.)), joint_a, body_a);

    body_b = Body (1., Vector3d (0., 1., 0.), Vector3d (1., 1., 1.));
    joint_b = Joint ( SpatialVectord (0., 1., 0., 0., 0., 0.));

    body_b_id = model->AddBody(*model_data, 1, Xtrans(Vector3d(1., 0., 0.)), joint_b, body_b);

    body_c = Body (1., Vector3d (0., 0., 1.), Vector3d (1., 1., 1.));
    joint_c = Joint ( SpatialVectord (0., 0., 1., 0., 0., 0.));

    body_c_id = model->AddBody(*model_data, 2, Xtrans(Vector3d(0., 1., 0.)), joint_c, body_c);

    Q = VectorNd::Constant ((size_t) model->dof_count, 0.);
    QDot = VectorNd::Constant ((size_t) model->dof_count, 0.);
    QDDot = VectorNd::Constant ((size_t) model->dof_count, 0.);
    Tau = VectorNd::Constant ((size_t) model->dof_count, 0.);

    point_position = Vector3d::Zero (3);
    point_acceleration = Vector3d::Zero (3);

    ref_body_id = 0;

    ClearLogOutput();
  }
  virtual void TearDown () {
    delete model;
  }

  RigidBodyDynamics::ModelData *model_data;
  RigidBodyDynamics::Model *model;

  unsigned int body_a_id, body_b_id, body_c_id, ref_body_id;
  RigidBodyDynamics::Body body_a, body_b, body_c;
  RigidBodyDynamics::Joint joint_a, joint_b, joint_c;

  RigidBodyDynamics::Math::VectorNd Q;
  RigidBodyDynamics::Math::VectorNd QDot;
  RigidBodyDynamics::Math::VectorNd QDDot;
  RigidBodyDynamics::Math::VectorNd Tau;

  RigidBodyDynamics::Math::Vector3d point_position, point_acceleration;
};

class FixedBase6DoF : public ::testing::Test {
protected:
  virtual void SetUp() {
    using namespace RigidBodyDynamics;
    using namespace RigidBodyDynamics::Math;

    ClearLogOutput();
    model_data = new ModelData;
    model = new Model(*model_data);

    model->gravity = Vector3d  (0., -9.81, 0.);

    /* 3 DoF (rot.) joint at base
     * 3 DoF (rot.) joint child origin
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
    base_rot_z = Body (
        0.,
        Vector3d (0., 0., 0.),
        Vector3d (0., 0., 0.)
        );
    joint_base_rot_z = Joint ( SpatialVectord (0., 0., 1., 0., 0., 0.));
    base_rot_z_id = model->AddBody (*model_data, 0, Xtrans (Vector3d (0., 0., 0.)), joint_base_rot_z, base_rot_z);

    base_rot_y = Body (
        0.,
        Vector3d (0., 0., 0.),
        Vector3d (0., 0., 0.)
        );
    joint_base_rot_y = Joint ( SpatialVectord (0., 1., 0., 0., 0., 0.));
    base_rot_y_id = model->AppendBody (*model_data, Xtrans (Vector3d (0., 0., 0.)), joint_base_rot_y, base_rot_y);

    base_rot_x = Body (
        1.,
        Vector3d (0.5, 0., 0.),
        Vector3d (1., 1., 1.)
        );
    joint_base_rot_x = Joint ( SpatialVectord (1., 0., 0., 0., 0., 0.));
    base_rot_x_id = model->AddBody (*model_data, base_rot_y_id, Xtrans (Vector3d (0., 0., 0.)), joint_base_rot_x, base_rot_x);

    // child body (3 DoF)
    child_rot_z = Body (
        0.,
        Vector3d (0., 0., 0.),
        Vector3d (0., 0., 0.)
        );
    joint_child_rot_z = Joint ( SpatialVectord (0., 0., 1., 0., 0., 0.));
    child_rot_z_id = model->AddBody (*model_data, base_rot_x_id, Xtrans (Vector3d (1., 0., 0.)), joint_child_rot_z, child_rot_z);

    child_rot_y = Body (
        0.,
        Vector3d (0., 0., 0.),
        Vector3d (0., 0., 0.)
        );
    joint_child_rot_y = Joint ( SpatialVectord (0., 1., 0., 0., 0., 0.));
    child_rot_y_id = model->AddBody (*model_data, child_rot_z_id, Xtrans (Vector3d (0., 0., 0.)), joint_child_rot_y, child_rot_y);

    child_rot_x = Body (
        1.,
        Vector3d (0., 0.5, 0.),
        Vector3d (1., 1., 1.)
        );
    joint_child_rot_x = Joint ( SpatialVectord (1., 0., 0., 0., 0., 0.));
    child_rot_x_id = model->AddBody (*model_data, child_rot_y_id, Xtrans (Vector3d (0., 0., 0.)), joint_child_rot_x, child_rot_x);

    Q = VectorNd::Constant (model->mBodies.size() - 1, 0.);
    QDot = VectorNd::Constant (model->mBodies.size() - 1, 0.);
    QDDot = VectorNd::Constant (model->mBodies.size() - 1, 0.);
    Tau = VectorNd::Constant (model->mBodies.size() - 1, 0.);

    contact_body_id = child_rot_x_id;
    contact_point = Vector3d  (0.5, 0.5, 0.);
    contact_normal = Vector3d  (0., 1., 0.);

    ClearLogOutput();
  }

  virtual void TearDown() {
    delete model;
  }

  RigidBodyDynamics::ModelData *model_data;
  RigidBodyDynamics::Model *model;

  unsigned int base_rot_z_id, base_rot_y_id, base_rot_x_id,
               child_rot_z_id, child_rot_y_id, child_rot_x_id,
               base_body_id;

  RigidBodyDynamics::Body base_rot_z, base_rot_y, base_rot_x,
    child_rot_z, child_rot_y, child_rot_x;

  RigidBodyDynamics::Joint joint_base_rot_z, joint_base_rot_y, joint_base_rot_x,
    joint_child_rot_z, joint_child_rot_y, joint_child_rot_x;

  RigidBodyDynamics::Math::VectorNd Q;
  RigidBodyDynamics::Math::VectorNd QDot;
  RigidBodyDynamics::Math::VectorNd QDDot;
  RigidBodyDynamics::Math::VectorNd Tau;

  unsigned int contact_body_id;
  RigidBodyDynamics::Math::Vector3d contact_point;
  RigidBodyDynamics::Math::Vector3d contact_normal;
  RigidBodyDynamics::ConstraintSet constraint_set;
};

class FloatingBase12DoF : public ::testing::Test {
protected:
  virtual void SetUp () {
    using namespace RigidBodyDynamics;
    using namespace RigidBodyDynamics::Math;

    ClearLogOutput();
    model_data = new ModelData;
    model = new Model(*model_data);

    model->gravity = Vector3d  (0., -9.81, 0.);

    /* 3 DoF (rot.) joint at base
     * 3 DoF (rot.) joint child origin
     *
     *          X Contact point (ref child)
     *          |
     *    Base  |
     *   / body |
     *  O-------*
     *           \
     *             Child body
     */

    base_rot_x = Body (
        1.,
        Vector3d (0.5, 0., 0.),
        Vector3d (1., 1., 1.)
        );
    base_rot_x_id = model->AddBody (*model_data, 0, SpatialTransformd(),
        Joint (
          SpatialVectord (0., 0., 0., 1., 0., 0.),
          SpatialVectord (0., 0., 0., 0., 1., 0.),
          SpatialVectord (0., 0., 0., 0., 0., 1.),
          SpatialVectord (0., 0., 1., 0., 0., 0.),
          SpatialVectord (0., 1., 0., 0., 0., 0.),
          SpatialVectord (1., 0., 0., 0., 0., 0.)
          ),
        base_rot_x);

    // child body 1 (3 DoF)
    child_rot_z = Body (
        0.,
        Vector3d (0., 0., 0.),
        Vector3d (0., 0., 0.)
        );
    joint_child_rot_z = Joint ( SpatialVectord (0., 0., 1., 0., 0., 0.));
    child_rot_z_id = model->AddBody (*model_data, base_rot_x_id, Xtrans (Vector3d (1., 0., 0.)), joint_child_rot_z, child_rot_z);

    child_rot_y = Body (
        0.,
        Vector3d (0., 0., 0.),
        Vector3d (0., 0., 0.)
        );
    joint_child_rot_y = Joint ( SpatialVectord (0., 1., 0., 0., 0., 0.));
    child_rot_y_id = model->AddBody (*model_data, child_rot_z_id, Xtrans (Vector3d (0., 0., 0.)), joint_child_rot_y, child_rot_y);

    child_rot_x = Body (
        1.,
        Vector3d (0., 0.5, 0.),
        Vector3d (1., 1., 1.)
        );
    joint_child_rot_x = Joint ( SpatialVectord (1., 0., 0., 0., 0., 0.));
    child_rot_x_id = model->AddBody (*model_data, child_rot_y_id, Xtrans (Vector3d (0., 0., 0.)), joint_child_rot_x, child_rot_x);

    // child body (3 DoF)
    child_2_rot_z = Body (
        0.,
        Vector3d (0., 0., 0.),
        Vector3d (0., 0., 0.)
        );
    joint_child_2_rot_z = Joint ( SpatialVectord (0., 0., 1., 0., 0., 0.));
    child_2_rot_z_id = model->AddBody (*model_data, child_rot_x_id, Xtrans (Vector3d (1., 0., 0.)), joint_child_2_rot_z, child_2_rot_z);

    child_2_rot_y = Body (
        0.,
        Vector3d (0., 0., 0.),
        Vector3d (0., 0., 0.)
        );
    joint_child_2_rot_y = Joint ( SpatialVectord (0., 1., 0., 0., 0., 0.));
    child_2_rot_y_id = model->AddBody (*model_data, child_2_rot_z_id, Xtrans (Vector3d (0., 0., 0.)), joint_child_2_rot_y, child_2_rot_y);

    child_2_rot_x = Body (
        1.,
        Vector3d (0., 0.5, 0.),
        Vector3d (1., 1., 1.)
        );
    joint_child_2_rot_x = Joint ( SpatialVectord (1., 0., 0., 0., 0., 0.));
    child_2_rot_x_id = model->AddBody (*model_data, child_2_rot_y_id, Xtrans (Vector3d (0., 0., 0.)), joint_child_2_rot_x, child_2_rot_x);

    Q = VectorNd::Constant (model->dof_count, 0.);
    QDot = VectorNd::Constant (model->dof_count, 0.);
    QDDot = VectorNd::Constant (model->dof_count, 0.);
    Tau = VectorNd::Constant (model->dof_count, 0.);

    ClearLogOutput();
  }

  virtual void TearDown () {
    delete model;
  }

  RigidBodyDynamics::ModelData *model_data;
  RigidBodyDynamics::Model *model;

  unsigned int base_rot_z_id, base_rot_y_id, base_rot_x_id,
               child_rot_z_id, child_rot_y_id, child_rot_x_id,
               child_2_rot_z_id, child_2_rot_y_id,child_2_rot_x_id,
               base_body_id;

  RigidBodyDynamics::Body base_rot_z, base_rot_y, base_rot_x,
    child_rot_z, child_rot_y, child_rot_x,
    child_2_rot_z, child_2_rot_y, child_2_rot_x;

  RigidBodyDynamics::Joint joint_base_rot_z, joint_base_rot_y, joint_base_rot_x,
    joint_child_rot_z, joint_child_rot_y, joint_child_rot_x,
    joint_child_2_rot_z, joint_child_2_rot_y, joint_child_2_rot_x;

  RigidBodyDynamics::Math::VectorNd Q;
  RigidBodyDynamics::Math::VectorNd QDot;
  RigidBodyDynamics::Math::VectorNd QDDot;
  RigidBodyDynamics::Math::VectorNd Tau;
};

class SimpleFixture : public ::testing::Test {
protected:
  virtual void SetUp () {
    ClearLogOutput();
    model_data = new RigidBodyDynamics::ModelData;
    model = new RigidBodyDynamics::Model(*model_data);
    model->gravity = RigidBodyDynamics::Math::Vector3d (0., -9.81, 0.);
  }
  virtual void TearDown() {
    delete model;
  }
  void ResizeVectors () {
    Q = RigidBodyDynamics::Math::VectorNd::Zero (model->dof_count);
    QDot = RigidBodyDynamics::Math::VectorNd::Zero (model->dof_count);
    QDDot = RigidBodyDynamics::Math::VectorNd::Zero (model->dof_count);
    Tau = RigidBodyDynamics::Math::VectorNd::Zero (model->dof_count);
  }

  RigidBodyDynamics::ModelData *model_data;
  RigidBodyDynamics::Model *model;

  RigidBodyDynamics::Math::VectorNd Q;
  RigidBodyDynamics::Math::VectorNd QDot;
  RigidBodyDynamics::Math::VectorNd QDDot;
  RigidBodyDynamics::Math::VectorNd Tau;
};

class FixedJoint2DoF : public ::testing::Test {
protected:
  virtual void SetUp () {
    using namespace RigidBodyDynamics;
    using namespace RigidBodyDynamics::Math;

    ClearLogOutput();
    model_data = new ModelData;
    model = new Model(*model_data);

    /* Basically a model like this, where X are the Center of Masses
     * and the CoM of the last (3rd) body comes out of the Y=X=0 plane.
     *
     *      Z
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
    joint_b = Joint (JointTypeFixed);

    body_b_id = model->AddBody(*model_data, 1, Xtrans(Vector3d(1., 0., 0.)), joint_b, body_b);

    body_c = Body (1., Vector3d (0., 0., 1.), Vector3d (1., 1., 1.));
    joint_c = Joint ( SpatialVectord (0., 0., 1., 0., 0., 0.));

    body_c_id = model->AddBody(*model_data, 2, Xtrans(Vector3d(0., 1., 0.)), joint_c, body_c);

    Q = VectorNd::Constant ((size_t) model->dof_count, 0.);
    QDot = VectorNd::Constant ((size_t) model->dof_count, 0.);
    QDDot = VectorNd::Constant ((size_t) model->dof_count, 0.);

    point_position = Vector3d::Zero (3);
    point_acceleration = Vector3d::Zero (3);

    ref_body_id = 0;

    ClearLogOutput();
  }
  virtual void TearDown () {
    delete model;
  }

  RigidBodyDynamics::ModelData *model_data;
  RigidBodyDynamics::Model *model;

  unsigned int body_a_id, body_b_id, body_c_id, ref_body_id;
  RigidBodyDynamics::Body body_a, body_b, body_c;
  RigidBodyDynamics::Joint joint_a, joint_b, joint_c;

  RigidBodyDynamics::Math::VectorNd Q;
  RigidBodyDynamics::Math::VectorNd QDot;
  RigidBodyDynamics::Math::VectorNd QDDot;

  RigidBodyDynamics::Math::Vector3d point_position, point_acceleration;
};

/** \brief Fixture that contains two models of which one has one joint fixed.
*/
class FixedAndMovableJoint : public ::testing::Test{
protected:
  virtual void SetUp () {
    using namespace RigidBodyDynamics;
    using namespace RigidBodyDynamics::Math;

    ClearLogOutput();

    model_movable_data = new ModelData;
    model_movable = new Model(*model_movable_data);

    /* Basically a model like this, where X are the Center of Masses
     * and the CoM of the last (3rd) body comes out of the Y=X=0 plane.
     *
     *      Z
     *      *---*
     *      |
     *      |
     *  Z   |
     *  O---*
     *      Y
     */

    body_a = Body (1., Vector3d (1., 0., 0.), Vector3d (1., 1., 1.));
    joint_a = Joint( SpatialVectord (0., 0., 1., 0., 0., 0.));

    body_a_id = model_movable->AddBody(*model_movable_data, 0, Xtrans(Vector3d(0., 0., 0.)), joint_a, body_a);

    body_b = Body (1., Vector3d (0., 1., 0.), Vector3d (1., 1., 1.));
    joint_b = Joint ( SpatialVectord (0., 1., 0., 0., 0., 0.));

    body_b_id = model_movable->AddBody(*model_movable_data, body_a_id, Xtrans(Vector3d(1., 0., 0.)), joint_b, body_b);

    body_c = Body (1., Vector3d (0., 0., 1.), Vector3d (1., 1., 1.));
    joint_c = Joint ( SpatialVectord (0., 0., 1., 0., 0., 0.));

    body_c_id = model_movable->AddBody(*model_movable_data, body_b_id, Xtrans(Vector3d(0., 1., 0.)), joint_c, body_c);

    Q = VectorNd::Constant ((size_t) model_movable->dof_count, 0.);
    QDot = VectorNd::Constant ((size_t) model_movable->dof_count, 0.);
    QDDot = VectorNd::Constant ((size_t) model_movable->dof_count, 0.);
    Tau = VectorNd::Constant ((size_t) model_movable->dof_count, 0.);
    C_movable = VectorNd::Zero ((size_t) model_movable->dof_count);
    H_movable = MatrixNd::Zero ((size_t) model_movable->dof_count, (size_t) model_movable->dof_count);

    // Assemble the fixed joint model
    model_fixed_data = new ModelData;
    model_fixed = new Model(*model_fixed_data);

    body_a_fixed_id = model_fixed->AddBody(0, Xtrans(Vector3d(0., 0., 0.)), joint_a, body_a);
    Joint joint_fixed (JointTypeFixed);
    body_b_fixed_id = model_fixed->AddBody(body_a_fixed_id, Xtrans(Vector3d(1., 0., 0.)), joint_fixed, body_b);
    body_c_fixed_id = model_fixed->AddBody(body_b_fixed_id, Xtrans(Vector3d(0., 1., 0.)), joint_c, body_c);

    Q_fixed = VectorNd::Constant ((size_t) model_fixed->dof_count, 0.);
    QDot_fixed = VectorNd::Constant ((size_t) model_fixed->dof_count, 0.);
    QDDot_fixed = VectorNd::Constant ((size_t) model_fixed->dof_count, 0.);
    Tau_fixed = VectorNd::Constant ((size_t) model_fixed->dof_count, 0.);
    C_fixed = VectorNd::Zero ((size_t) model_fixed->dof_count);
    H_fixed = MatrixNd::Zero ((size_t) model_fixed->dof_count, (size_t) model_fixed->dof_count);

    point_position = Vector3d::Zero (3);
    point_acceleration = Vector3d::Zero (3);

    ref_body_id = 0;

    ClearLogOutput();
  }

  virtual void TearDown () {
    delete model_movable;
    delete model_fixed;
  }
  RigidBodyDynamics::Math::VectorNd CreateDofVectorFromReducedVector (const RigidBodyDynamics::Math::VectorNd &q_fixed) {
    assert (q_fixed.size() == model_fixed->dof_count);

    RigidBodyDynamics::Math::VectorNd q_movable (model_movable->dof_count);

    q_movable[0] = q_fixed[0];
    q_movable[1] = 0.;
    q_movable[2] = q_fixed[1];

    return q_movable;
  }

  RigidBodyDynamics::Math::MatrixNd CreateReducedInertiaMatrix(const RigidBodyDynamics::Math::MatrixNd &H_movable) {
    assert (H_movable.rows() == model_movable->dof_count);
    assert (H_movable.cols() == model_movable->dof_count);
    RigidBodyDynamics::Math::MatrixNd H (model_fixed->dof_count, model_fixed->dof_count);

    H (0,0) = H_movable(0,0); H (0,1) = H_movable(0,2);
    H (1,0) = H_movable(2,0); H (1,1) = H_movable(2,2);

    return H;
  }

  RigidBodyDynamics::ModelData *model_fixed_data;
  RigidBodyDynamics::Model *model_fixed;

  RigidBodyDynamics::ModelData *model_movable_data;
  RigidBodyDynamics::Model *model_movable;

  unsigned int body_a_id, body_b_id, body_c_id, ref_body_id;
  unsigned int body_a_fixed_id, body_b_fixed_id, body_c_fixed_id;

  RigidBodyDynamics::Body body_a, body_b, body_c;
  RigidBodyDynamics::Joint joint_a, joint_b, joint_c;

  RigidBodyDynamics::Math::VectorNd Q;
  RigidBodyDynamics::Math::VectorNd QDot;
  RigidBodyDynamics::Math::VectorNd QDDot;
  RigidBodyDynamics::Math::VectorNd Tau;
  RigidBodyDynamics::Math::VectorNd C_movable;
  RigidBodyDynamics::Math::MatrixNd H_movable;

  RigidBodyDynamics::Math::VectorNd Q_fixed;
  RigidBodyDynamics::Math::VectorNd QDot_fixed;
  RigidBodyDynamics::Math::VectorNd QDDot_fixed;
  RigidBodyDynamics::Math::VectorNd Tau_fixed;
  RigidBodyDynamics::Math::VectorNd C_fixed;
  RigidBodyDynamics::Math::MatrixNd H_fixed;

  RigidBodyDynamics::Math::Vector3d point_position, point_acceleration;
};

/** Model with two moving bodies and one fixed body
*/
class RotZRotZYXFixed : public ::testing::Test{
protected:
  virtual void SetUp() {
    using namespace RigidBodyDynamics;
    using namespace RigidBodyDynamics::Math;

    ClearLogOutput();
    model = new Model;

    Joint joint_rot_z ( SpatialVectord (0., 0., 1., 0., 0., 0.));

    Joint joint_rot_zyx (
        SpatialVectord (0., 0., 1., 0., 0., 0.),
        SpatialVectord (0., 1., 0., 0., 0., 0.),
        SpatialVectord (1., 0., 0., 0., 0., 0.)
        );

    Body body_a(1., RigidBodyDynamics::Math::Vector3d (1., 0.4, 0.4), RigidBodyDynamics::Math::Vector3d (1., 1., 1.));
    Body body_b(2., RigidBodyDynamics::Math::Vector3d (1., 0.4, 0.4), RigidBodyDynamics::Math::Vector3d (1., 1., 1.));
    Body body_fixed(10., RigidBodyDynamics::Math::Vector3d (1., 0.4, 0.4), RigidBodyDynamics::Math::Vector3d (1., 1., 1.));

    fixture_transform_a = Xtrans (RigidBodyDynamics::Math::Vector3d(1., 2., 3.));
    fixture_transform_b = Xtrans (RigidBodyDynamics::Math::Vector3d(4., 5., 6.));
    fixture_transform_fixed = Xtrans (RigidBodyDynamics::Math::Vector3d(-1., -2., -3.));

    body_a_id = model->AddBody (*model_data, 0, fixture_transform_a, joint_rot_z, body_a);
    body_b_id = model->AppendBody (fixture_transform_b, joint_rot_zyx, body_b);
    body_fixed_id = model->AppendBody (fixture_transform_fixed, Joint(JointTypeFixed), body_fixed);

    ClearLogOutput();
  }
  virtual void TearDown() {
    delete model;
  }

  RigidBodyDynamics::Model *model;

  unsigned int body_a_id, body_b_id, body_fixed_id;

  RigidBodyDynamics::Math::SpatialTransformd fixture_transform_a;
  RigidBodyDynamics::Math::SpatialTransformd fixture_transform_b;
  RigidBodyDynamics::Math::SpatialTransformd fixture_transform_fixed;
};

class TwoArms12DoF : public ::testing::Test {
protected:
  virtual void SetUp() {
    using namespace RigidBodyDynamics;
    using namespace RigidBodyDynamics::Math;

    ClearLogOutput();
    model = new Model;

    /* Basically a model like this, where X are the Center of Masses
     * and the CoM of the last (3rd) body comes out of the Y=X=0 plane.
     *
     * *----O----*
     * |         |
     * |         |
     * *         *
     * |         |
     * |         |
     *
     */

    Body body_upper = Body (1., Vector3d (0., -0.2, 0.), Vector3d (1.1, 1.3, 1.5));
    Body body_lower = Body (0.5, Vector3d(0., -0.15, 0.), Vector3d (0.3, 0.5, 0.2));

    Joint joint_zyx = Joint (
        SpatialVectord (0., 0., 1., 0., 0., 0.),
        SpatialVectord (0., 1., 0., 0., 0., 0.),
        SpatialVectord (1., 0., 0., 0., 0., 0.)
        );

    right_upper_arm = model->AppendBody (Xtrans (Vector3d (0., 0., -0.3)), joint_zyx, body_upper, "RightUpper");
    //		model->AppendBody (Xtrans (Vector3d (0., -0.4, 0.)), joint_zyx, body_lower, "RightLower");
    left_upper_arm = model->AddBody (*model_data, 0, Xtrans (Vector3d (0., 0., 0.3)), joint_zyx, body_upper, "LeftUpper");
    //		model->AppendBody (Xtrans (Vector3d (0., -0.4, 0.)), joint_zyx, body_lower, "LeftLower");

    q = VectorNd::Constant ((size_t) model->dof_count, 0.);
    qdot = VectorNd::Constant ((size_t) model->dof_count, 0.);
    qddot = VectorNd::Constant ((size_t) model->dof_count, 0.);
    tau = VectorNd::Constant ((size_t) model->dof_count, 0.);

    ClearLogOutput();
  }
  virtual void TearDown() {
    delete model;
  }

  RigidBodyDynamics::Model *model;

  RigidBodyDynamics::Math::VectorNd q;
  RigidBodyDynamics::Math::VectorNd qdot;
  RigidBodyDynamics::Math::VectorNd qddot;
  RigidBodyDynamics::Math::VectorNd tau;

  unsigned int right_upper_arm, left_upper_arm;

};

class FixedBase6DoF12DoFFloatingBase : public ::testing::Test {
protected:
  virtual void SetUp () {
    using namespace RigidBodyDynamics;
    using namespace RigidBodyDynamics::Math;

    ClearLogOutput();
    model_data = new ModelData;
    model = new Model(*model_data);

    model->gravity = Vector3d  (0., -9.81, 0.);

    /* 3 DoF (rot.) joint at base
     * 3 DoF (rot.) joint child origin
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
    base_id = model->AddBody (*model_data, 0, SpatialTransformd(),
        Joint (
          SpatialVectord (0., 0., 0., 1., 0., 0.),
          SpatialVectord (0., 0., 0., 0., 1., 0.),
          SpatialVectord (0., 0., 0., 0., 0., 1.),
          SpatialVectord (0., 0., 1., 0., 0., 0.),
          SpatialVectord (0., 1., 0., 0., 0., 0.),
          SpatialVectord (1., 0., 0., 0., 0., 0.)
          ),
        base);

    // child body 1 (3 DoF)
    child = Body (
        1.,
        Vector3d (0., 0.5, 0.),
        Vector3d (1., 1., 1.)
        );
    joint_rotzyx = Joint (
        SpatialVectord (0., 0., 1., 0., 0., 0.),
        SpatialVectord (0., 1., 0., 0., 0., 0.),
        SpatialVectord (1., 0., 0., 0., 0., 0.)
        );
    child_id = model->AddBody (*model_data, base_id, Xtrans (Vector3d (1., 0., 0.)), joint_rotzyx, child);

    // child body (3 DoF)
    child_2 = Body (
        1.,
        Vector3d (0., 0.5, 0.),
        Vector3d (1., 1., 1.)
        );
    child_2_id = model->AddBody (*model_data, child_id, Xtrans (Vector3d (1., 0., 0.)), joint_rotzyx, child_2);

    Q = VectorNd::Constant (model->dof_count, 0.);
    QDot = VectorNd::Constant (model->dof_count, 0.);
    QDDot = VectorNd::Constant (model->dof_count, 0.);
    Tau = VectorNd::Constant (model->dof_count, 0.);

    contact_body_id = child_id;
    contact_point = Vector3d  (0.5, 0.5, 0.);
    contact_normal = Vector3d  (0., 1., 0.);

    ClearLogOutput();
  }

  virtual void TearDown () {
    delete model;
  }

  RigidBodyDynamics::ModelData *model_data;
  RigidBodyDynamics::Model *model;

  unsigned int base_id, child_id, child_2_id;

  RigidBodyDynamics::Body base, child, child_2;

  RigidBodyDynamics::Joint joint_rotzyx;

  RigidBodyDynamics::Math::VectorNd Q;
  RigidBodyDynamics::Math::VectorNd QDot;
  RigidBodyDynamics::Math::VectorNd QDDot;
  RigidBodyDynamics::Math::VectorNd Tau;

  unsigned int contact_body_id;
  RigidBodyDynamics::Math::Vector3d contact_point;
  RigidBodyDynamics::Math::Vector3d contact_normal;
  RigidBodyDynamics::ConstraintSet constraint_set;
};

#endif // FIXTURESTESTS_H
