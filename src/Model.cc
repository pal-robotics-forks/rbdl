/*
 * RBDL - Rigid Body Dynamics Library
 * Copyright (c) 2011-2016 Martin Felis <martin.felis@iwr.uni-heidelberg.de>
 *
 * Licensed under the zlib license. See LICENSE for more details.
 */

#include <iostream>
#include <limits>
#include <assert.h>

#include "rbdl/rbdl_mathutils.h"

#include "rbdl/Logging.h"

#include "rbdl/Model.h"
#include "rbdl/Body.h"
#include "rbdl/Joint.h"

#include "rbdl/rbdl_utils.h"

using namespace RigidBodyDynamics;
using namespace RigidBodyDynamics::Math;

Model::Model() : model_data_(new ModelDatad())
{
  ModelDatad &model_data = *model_data_.get();

  Body root_body;
  Joint root_joint;

  Vector3d zero_position(0., 0., 0.);
  SpatialVectord zero_spatial(0., 0., 0., 0., 0., 0.);

  // structural information
  lambda.push_back(0);
  lambda_q.push_back(0);
  mu.push_back(std::vector<unsigned int>());
  dof_count = 0;
  q_size = 0;
  qdot_size = 0;
  previously_added_body_id = 0;

  gravity = Vector3d(0., -9.81, 0.);

  // state information
  model_data.v.push_back(zero_spatial);
  model_data.a.push_back(zero_spatial);
  model_data.a_bias.push_back(zero_spatial);
  SpatialMatrixd zero_spatial_matrix;
  zero_spatial_matrix.setZero();
  model_data.acumulated_mass.push_back(zero_spatial_matrix);

  // Joints
  mJoints.push_back(root_joint);
  model_data.S.push_back(zero_spatial);
  X_T.push_back(SpatialTransformd());

  model_data.X_J.push_back(SpatialTransformd());
  model_data.v_J.push_back(zero_spatial);
  model_data.c_J.push_back(zero_spatial);

  // Spherical joints
  model_data.multdof3_S.push_back(Matrix63d::Zero());
  model_data.multdof3_U.push_back(Matrix63d::Zero());
  model_data.multdof3_Dinv.push_back(Matrix3d::Zero());
  model_data.multdof3_u.push_back(Vector3d::Zero());
  multdof3_w_index.push_back(0);

  // Dynamic variables
  model_data.c.push_back(zero_spatial);
  model_data.IA.push_back(SpatialMatrixd::Identity());
  model_data.pA.push_back(zero_spatial);
  model_data.U.push_back(zero_spatial);

  model_data.u = VectorNd::Zero(1);
  model_data.d = VectorNd::Zero(1);

  model_data.f.push_back(zero_spatial);
  SpatialRigidBodyInertiad rbi(0., Vector3d(0., 0., 0.), Matrix3d::Zero(3, 3));
  model_data.Ic.push_back(rbi);
  model_data.I.push_back(rbi);
  model_data.hc.push_back(zero_spatial);

  // Bodies
  model_data.X_lambda.push_back(SpatialTransformd());
  model_data.X_base.push_back(SpatialTransformd());

  mBodies.push_back(root_body);
  mBodyNameMap["ROOT"] = 0;

  fixed_body_discriminator = std::numeric_limits<unsigned int>::max() / 2;
}

Model::Model(ModelDatad &model_data)
{
  Body root_body;
  Joint root_joint;

  Vector3d zero_position(0., 0., 0.);
  SpatialVectord zero_spatial(0., 0., 0., 0., 0., 0.);

  // structural information
  lambda.push_back(0);
  lambda_q.push_back(0);
  mu.push_back(std::vector<unsigned int>());
  dof_count = 0;
  q_size = 0;
  qdot_size = 0;
  previously_added_body_id = 0;

  gravity = Vector3d(0., -9.81, 0.);

  // state information
  model_data.v.push_back(zero_spatial);
  model_data.a.push_back(zero_spatial);
  model_data.a_bias.push_back(zero_spatial);
  SpatialMatrixd zero_spatial_matrix;
  zero_spatial_matrix.setZero();
  model_data.acumulated_mass.push_back(zero_spatial_matrix);

  // Joints
  mJoints.push_back(root_joint);
  model_data.S.push_back(zero_spatial);
  X_T.push_back(SpatialTransformd());

  model_data.X_J.push_back(SpatialTransformd());
  model_data.v_J.push_back(zero_spatial);
  model_data.c_J.push_back(zero_spatial);

  // Spherical joints
  model_data.multdof3_S.push_back(Matrix63d::Zero());
  model_data.multdof3_U.push_back(Matrix63d::Zero());
  model_data.multdof3_Dinv.push_back(Matrix3d::Zero());
  model_data.multdof3_u.push_back(Vector3d::Zero());
  multdof3_w_index.push_back(0);

  // Dynamic variables
  model_data.c.push_back(zero_spatial);
  model_data.IA.push_back(SpatialMatrixd::Identity());
  model_data.pA.push_back(zero_spatial);
  model_data.U.push_back(zero_spatial);

  model_data.u = VectorNd::Zero(1);
  model_data.d = VectorNd::Zero(1);

  model_data.f.push_back(zero_spatial);
  SpatialRigidBodyInertiad rbi(0., Vector3d(0., 0., 0.), Matrix3d::Zero(3, 3));
  model_data.Ic.push_back(rbi);
  model_data.I.push_back(rbi);
  model_data.hc.push_back(zero_spatial);

  // Bodies
  model_data.X_lambda.push_back(SpatialTransformd());
  model_data.X_base.push_back(SpatialTransformd());

  mBodies.push_back(root_body);
  mBodyNameMap["ROOT"] = 0;

  fixed_body_discriminator = std::numeric_limits<unsigned int>::max() / 2;
}

unsigned int AddBodyFixedJoint(Model &model, ModelDatad &model_data, const unsigned int parent_id,
                               const SpatialTransformd &joint_frame, const Joint &joint,
                               const Body &body, std::string body_name)
{
  FixedBody fbody = FixedBody::CreateFromBody(body);
  fbody.mMovableParent = parent_id;
  fbody.mParentTransform = joint_frame;

  if (model.IsFixedBodyId(parent_id))
  {
    FixedBody fixed_parent = model.mFixedBodies[parent_id - model.fixed_body_discriminator];

    fbody.mMovableParent = fixed_parent.mMovableParent;
    fbody.mParentTransform = joint_frame * fixed_parent.mParentTransform;
  }

  // merge the two bodies
  Body parent_body = model.mBodies[fbody.mMovableParent];
  parent_body.Join(fbody.mParentTransform, body);
  model.mBodies[fbody.mMovableParent] = parent_body;
  model_data.I[fbody.mMovableParent] = SpatialRigidBodyInertiad::createFromMassComInertiaC(
      parent_body.mMass, parent_body.mCenterOfMass, parent_body.mInertia);

  model.mFixedBodies.push_back(fbody);

  if (model.mFixedBodies.size() >
      std::numeric_limits<unsigned int>::max() - model.fixed_body_discriminator)
  {
    std::cerr << "Error: cannot add more than "
              << std::numeric_limits<unsigned int>::max() - model.mFixedBodies.size()
              << " fixed bodies. You need to modify "
              << "Model::fixed_body_discriminator for this." << std::endl;
    assert(0);
    abort();
  }

  if (body_name.size() != 0)
  {
    if (model.mBodyNameMap.find(body_name) != model.mBodyNameMap.end())
    {
      std::cerr << "Error: Body with name '" << body_name << "' already exists!" << std::endl;
      assert(0);
      abort();
    }
    model.mBodyNameMap[body_name] =
        model.mFixedBodies.size() + model.fixed_body_discriminator - 1;
  }

  return model.mFixedBodies.size() + model.fixed_body_discriminator - 1;
}

unsigned int AddBodyMultiDofJoint(Model &model, ModelDatad &model_data,
                                  const unsigned int parent_id,
                                  const SpatialTransformd &joint_frame, const Joint &joint,
                                  const Body &body, std::string body_name)
{
  // Here we emulate multi DoF joints by simply adding nullbodies. This
  // allows us to use fixed size elements for S,v,a, etc. which is very
  // fast in Eigen.
  unsigned int joint_count = 0;
  if (joint.mJointType == JointType1DoF)
    joint_count = 1;
  else if (joint.mJointType == JointType2DoF)
    joint_count = 2;
  else if (joint.mJointType == JointType3DoF)
    joint_count = 3;
  else if (joint.mJointType == JointType4DoF)
    joint_count = 4;
  else if (joint.mJointType == JointType5DoF)
    joint_count = 5;
  else if (joint.mJointType == JointType6DoF)
    joint_count = 6;
  else if (joint.mJointType == JointTypeFloatingBase)
  // no action required
  {
  }
  else
  {
    std::cerr << "Error: Invalid joint type: " << joint.mJointType << std::endl;

    assert(0 && !"Invalid joint type!");
  }

  Body null_body(0., Vector3d(0., 0., 0.), Vector3d(0., 0., 0.));
  null_body.mIsVirtual = true;

  unsigned int null_parent = parent_id;
  SpatialTransformd joint_frame_transform;

  if (joint.mJointType == JointTypeFloatingBase)
  {
    null_parent =
        model.AddBody(model_data, parent_id, joint_frame, JointTypeTranslationXYZ, null_body);

    return model.AddBody(model_data, null_parent, SpatialTransformd(), JointTypeSpherical,
                         body, body_name);
  }

  Joint single_dof_joint;
  unsigned int j;

  // Here we add multiple virtual bodies that have no mass or inertia for
  // which each is attached to the model with a single degree of freedom
  // joint.
  for (j = 0; j < joint_count; j++)
  {
    single_dof_joint = Joint(joint.mJointAxes[j]);

    if (single_dof_joint.mJointType == JointType1DoF)
    {
      Vector3d rotation(joint.mJointAxes[j][0], joint.mJointAxes[j][1], joint.mJointAxes[j][2]);
      Vector3d translation(joint.mJointAxes[j][3], joint.mJointAxes[j][4],
                           joint.mJointAxes[j][5]);

      if (rotation == Vector3d(0., 0., 0.))
      {
        single_dof_joint = Joint(JointTypePrismatic, translation);
      }
      else if (translation == Vector3d(0., 0., 0.))
      {
        single_dof_joint = Joint(JointTypeRevolute, rotation);
      }
      else
      {
        std::cerr << "Invalid joint axis: " << joint.mJointAxes[0].transpose()
                  << ". Helical joints not (yet) supported." << std::endl;
        abort();
      }
    }

    // the first joint has to be transformed by joint_frame, all the
    // others must have a null transformation
    if (j == 0)
      joint_frame_transform = joint_frame;
    else
      joint_frame_transform = SpatialTransformd();

    if (j == joint_count - 1)
      // if we are at the last we must add the real body
      break;
    else
    {
      // otherwise we just add an intermediate body
      null_parent = model.AddBody(model_data, null_parent, joint_frame_transform,
                                  single_dof_joint, null_body);
    }
  }

  return model.AddBody(model_data, null_parent, joint_frame_transform, single_dof_joint,
                       body, body_name);
}

unsigned int Model::AddBody(ModelDatad &model_data, const unsigned int parent_id,
                            const SpatialTransformd &joint_frame, const Joint &joint,
                            const Body &body, std::string body_name)
{
  assert(lambda.size() > 0);
  assert(joint.mJointType != JointTypeUndefined);

  if (joint.mJointType == JointTypeFixed)
  {
    previously_added_body_id =
        AddBodyFixedJoint(*this, model_data, parent_id, joint_frame, joint, body, body_name);

    return previously_added_body_id;
  }
  else if ((joint.mJointType == JointTypeSpherical) || (joint.mJointType == JointTypeEulerZYX) ||
           (joint.mJointType == JointTypeEulerXYZ) || (joint.mJointType == JointTypeEulerYXZ) ||
           (joint.mJointType == JointTypeTranslationXYZ)
           //|| (joint.mJointType == JointTypeCustom)
           )
  {
    // no action required
  }
  else if (joint.mJointType != JointTypePrismatic &&
           joint.mJointType != JointTypeRevolute && joint.mJointType != JointTypeRevoluteX &&
           joint.mJointType != JointTypeRevoluteY && joint.mJointType != JointTypeRevoluteZ)
  {
    previously_added_body_id = AddBodyMultiDofJoint(*this, model_data, parent_id,
                                                    joint_frame, joint, body, body_name);
    return previously_added_body_id;
  }

  // If we add the body to a fixed body we have to make sure that we
  // actually add it to its movable parent.
  unsigned int movable_parent_id = parent_id;
  SpatialTransformd movable_parent_transform;

  if (IsFixedBodyId(parent_id))
  {
    unsigned int fbody_id = parent_id - fixed_body_discriminator;
    movable_parent_id = mFixedBodies[fbody_id].mMovableParent;
    movable_parent_transform = mFixedBodies[fbody_id].mParentTransform;
  }

  // structural information
  lambda.push_back(movable_parent_id);
  unsigned int lambda_q_last = mJoints[mJoints.size() - 1].q_index;

  if (mJoints[mJoints.size() - 1].mDoFCount > 0)
  {
    //&& mJoints[mJoints.size() - 1].mJointType != JointTypeCustom) {
    lambda_q_last = lambda_q_last + mJoints[mJoints.size() - 1].mDoFCount;
  } /* else if (mJoints[mJoints.size() - 1].mJointType == JointTypeCustom) {
     unsigned int custom_index = mJoints[mJoints.size() - 1].custom_joint_index;
     lambda_q_last = lambda_q_last
       + mCustomJoints[mCustomJoints.size() - 1]->mDoFCount;
   }*/

  for (unsigned int i = 0; i < joint.mDoFCount; i++)
  {
    lambda_q.push_back(lambda_q_last + i);
  }
  mu.push_back(std::vector<unsigned int>());
  mu.at(movable_parent_id).push_back(mBodies.size());

  // Bodies
  model_data.X_lambda.push_back(SpatialTransformd());
  model_data.X_base.push_back(SpatialTransformd());
  mBodies.push_back(body);

  if (body_name.size() != 0)
  {
    if (mBodyNameMap.find(body_name) != mBodyNameMap.end())
    {
      std::cerr << "Error: Body with name '" << body_name << "' already exists!" << std::endl;
      assert(0);
      abort();
    }
    mBodyNameMap[body_name] = mBodies.size() - 1;
  }

  // state information
  model_data.v.push_back(SpatialVectord(0., 0., 0., 0., 0., 0.));
  model_data.a.push_back(SpatialVectord(0., 0., 0., 0., 0., 0.));
  model_data.a_bias.push_back(SpatialVectord(0., 0., 0., 0., 0., 0.));
  SpatialMatrixd zero_spatial_matrix;
  zero_spatial_matrix.setZero();
  model_data.acumulated_mass.push_back(zero_spatial_matrix);

  // Joints
  unsigned int prev_joint_index = mJoints.size() - 1;
  mJoints.push_back(joint);

  // if (mJoints[prev_joint_index].mJointType != JointTypeCustom) {
  mJoints[mJoints.size() - 1].q_index =
      mJoints[prev_joint_index].q_index + mJoints[prev_joint_index].mDoFCount;
  //} else {
  //  mJoints[mJoints.size() - 1].q_index =
  //      mJoints[prev_joint_index].q_index + mJoints[prev_joint_index].mDoFCount;
  //}

  model_data.S.push_back(joint.mJointAxes[0]);

  // Joint state variables
  model_data.X_J.push_back(SpatialTransformd());
  model_data.v_J.push_back(joint.mJointAxes[0]);
  model_data.c_J.push_back(SpatialVectord(0., 0., 0., 0., 0., 0.));

  // workspace for joints with 3 dof
  model_data.multdof3_S.push_back(Matrix63d::Zero(6, 3));
  model_data.multdof3_U.push_back(Matrix63d::Zero());
  model_data.multdof3_Dinv.push_back(Matrix3d::Zero());
  model_data.multdof3_u.push_back(Vector3d::Zero());
  multdof3_w_index.push_back(0);

  dof_count = dof_count + joint.mDoFCount;

  // update the w components of the Quaternions. They are stored at the end
  // of the q vector
  int multdof3_joint_counter = 0;
  int mCustomJoint_joint_counter = 0;
  for (unsigned int i = 1; i < mJoints.size(); i++)
  {
    if (mJoints[i].mJointType == JointTypeSpherical)
    {
      multdof3_w_index[i] = dof_count + multdof3_joint_counter;
      multdof3_joint_counter++;
    }
  }

  q_size = dof_count + multdof3_joint_counter;

  qdot_size = qdot_size + joint.mDoFCount;

  // we have to invert the transformation as it is later always used from the
  // child bodies perspective.
  X_T.push_back(joint_frame * movable_parent_transform);

  // Dynamic variables
  model_data.c.push_back(SpatialVectord(0., 0., 0., 0., 0., 0.));
  model_data.IA.push_back(SpatialMatrixd::Zero());
  model_data.pA.push_back(SpatialVectord(0., 0., 0., 0., 0., 0.));
  model_data.U.push_back(SpatialVectord(0., 0., 0., 0., 0., 0.));

  model_data.d = VectorNd::Zero(mBodies.size());
  model_data.u = VectorNd::Zero(mBodies.size());

  model_data.f.push_back(SpatialVectord(0., 0., 0., 0., 0., 0.));

  SpatialRigidBodyInertiad rbi = SpatialRigidBodyInertiad::createFromMassComInertiaC(
      body.mMass, body.mCenterOfMass, body.mInertia);

  model_data.Ic.push_back(rbi);
  model_data.I.push_back(rbi);
  model_data.hc.push_back(SpatialVectord(0., 0., 0., 0., 0., 0.));

  if (mBodies.size() == fixed_body_discriminator)
  {
    std::cerr << "Error: cannot add more than " << fixed_body_discriminator
              << " movable bodies. You need to modify "
                 "Model::fixed_body_discriminator for this."
              << std::endl;
    assert(0);
    abort();
  }

  previously_added_body_id = mBodies.size() - 1;

  mJointUpdateOrder.clear();

  // update the joint order computation
  std::vector<std::pair<JointType, unsigned int> > joint_types;
  for (unsigned int i = 0; i < mJoints.size(); i++)
  {
    joint_types.push_back(std::pair<JointType, unsigned int>(mJoints[i].mJointType, i));
    mJointUpdateOrder.push_back(mJoints[i].mJointType);
  }

  mJointUpdateOrder.clear();
  JointType current_joint_type = JointTypeUndefined;
  while (joint_types.size() != 0)
  {
    current_joint_type = joint_types[0].first;

    std::vector<std::pair<JointType, unsigned int> >::iterator type_iter = joint_types.begin();

    while (type_iter != joint_types.end())
    {
      if (type_iter->first == current_joint_type)
      {
        mJointUpdateOrder.push_back(type_iter->second);
        type_iter = joint_types.erase(type_iter);
      }
      else
      {
        ++type_iter;
      }
    }
  }

  //  for (unsigned int i = 0; i < mJointUpdateOrder.size(); i++) {
  //    std::cout << "i = " << i << ": joint_id = " << mJointUpdateOrder[i]
  // << " joint_type = " << mJoints[mJointUpdateOrder[i]].mJointType << std::endl;
  //  }

  return previously_added_body_id;
}

unsigned int Model::AppendBody(ModelDatad &model_data, const Math::SpatialTransformd &joint_frame,
                               const Joint &joint, const Body &body, std::string body_name)
{
  return Model::AddBody(model_data, previously_added_body_id, joint_frame, joint, body, body_name);
}
/*
unsigned int Model::AddBodyCustomJoint (
    ModelDatad &model_data,
    const unsigned int parent_id,
    const Math::SpatialTransformd &joint_frame,
    CustomJoint *custom_joint,
    const Body &body,
    std::string body_name) {
  Joint proxy_joint (JointTypeCustom, custom_joint->mDoFCount);
  proxy_joint.custom_joint_index = mCustomJoints.size();
  //proxy_joint.mDoFCount = custom_joint->mDoFCount; //MM added. Otherwise
  //model.q_size = 0, which is not good.

  mCustomJoints.push_back (custom_joint);

  unsigned int body_id = AddBody (model_data, parent_id,
      joint_frame,
      proxy_joint,
      body,
      body_name);

  return body_id;
}
*/

// Math::Quaterniond Model::GetQuaternion (unsigned int i,
//    const Math::VectorNd &Q) const {
////  assert (mJoints[i].mJointType == JointTypeSpherical);
////  unsigned int q_index = mJoints[i].q_index;
////  return Math::Quaterniond ( Q[q_index],
////      Q[q_index + 1],
////      Q[q_index + 2],
////      Q[multdof3_w_index[i]]);
//}


void Model::SetQuaternion(unsigned int i, const Math::Quaterniond &quat, Math::VectorNd &Q) const
{
  assert(mJoints[i].mJointType == JointTypeSpherical);
  unsigned int q_index = mJoints[i].q_index;

  Q[q_index] = quat[0];
  Q[q_index + 1] = quat[1];
  Q[q_index + 2] = quat[2];
  Q[multdof3_w_index[i]] = quat[3];
}


unsigned int Model::GetBodyId(const char *body_name) const
{
  return GetBodyId(body_name, *model_data_.get());
}
unsigned int Model::GetBodyId(const char *body_name, const ModelDatad &model_data) const
{
  auto it = mBodyNameMap.find(body_name);
  if (it == mBodyNameMap.end())
  {
    std::stringstream ss;
    ss << "GET BODY ID: ID does not exist: " << body_name << std::endl;
    ss << Utils::GetModelHierarchy(*this, model_data);
    throw std::runtime_error(ss.str());
    assert(mBodyNameMap.count(body_name) == 0);
    return std::numeric_limits<unsigned int>::max();
  }
  return it->second;
}
