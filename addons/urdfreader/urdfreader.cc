#include <rbdl/rbdl.h>

#include <rbdl/addons/urdfreader/urdfreader.h>

#include <assert.h>
#include <iostream>
#include <fstream>
#include <map>
#include <stack>

#include <urdf_model/model.h>
#include <urdf_parser/urdf_parser.h>
#include "ros/ros.h"

typedef boost::shared_ptr<urdf::Link> LinkPtr;
typedef const boost::shared_ptr<const urdf::Link> ConstLinkPtr;
typedef boost::shared_ptr<urdf::Joint> JointPtr;

using namespace std;


/* Check if the specified link is inside this urdf submodel */
bool isLinkInUrdfModel(LinkPtr link, std::string tip)
{
  if (tip == link->name)
  {
    return true;
  }
  else
  {
    if (link->child_links.empty())
    {
      return false;
    }
    else
    {
      bool found = false;
      for (unsigned int i = 0; i < link->child_links.size() && !found; ++i)
      {
        found = found || isLinkInUrdfModel(link->child_links[i], tip);
      }
      return found;
    }
  }
}

namespace RigidBodyDynamics
{
namespace Addons
{
using namespace Math;

typedef vector<LinkPtr> URDFLinkVector;
typedef vector<JointPtr> URDFJointVector;
typedef map<string, LinkPtr> URDFLinkMap;
typedef map<string, JointPtr> URDFJointMap;



int initialize_root(Model &rbdl_model, ModelDatad &model_data, ConstLinkPtr urdf_link,
                    const FloatingBaseType floatingBaseType)
{
  int new_id = 0;

  // If a floating base is desired we need to compute the inertia data of the floating
  // base link,
  // if the floating base has no inertial details, then it makes no sense to make it
  // floating
  // base (Like kdl). (The root is the only link allowed to be floating base

  if (floatingBaseType != +FloatingBaseType::FixedBase)
  {
    double link_inertial_mass = 0;
    Vector3d link_inertial_position;
    Vector3d link_inertial_rpy;
    Matrix3d link_inertial_inertia;

    link_inertial_position.setZero();
    link_inertial_rpy.setZero();
    link_inertial_inertia.setZero();

    if (urdf_link->inertial)
    {
      link_inertial_mass = urdf_link->inertial->mass;

      link_inertial_position.set(urdf_link->inertial->origin.position.x,
                                 urdf_link->inertial->origin.position.y,
                                 urdf_link->inertial->origin.position.z);
      urdf_link->inertial->origin.rotation.getRPY(link_inertial_rpy[0], link_inertial_rpy[1],
                                                  link_inertial_rpy[2]);

      link_inertial_inertia(0, 0) = urdf_link->inertial->ixx;
      link_inertial_inertia(0, 1) = urdf_link->inertial->ixy;
      link_inertial_inertia(0, 2) = urdf_link->inertial->ixz;

      link_inertial_inertia(1, 0) = urdf_link->inertial->ixy;
      link_inertial_inertia(1, 1) = urdf_link->inertial->iyy;
      link_inertial_inertia(1, 2) = urdf_link->inertial->iyz;

      link_inertial_inertia(2, 0) = urdf_link->inertial->ixz;
      link_inertial_inertia(2, 1) = urdf_link->inertial->iyz;
      link_inertial_inertia(2, 2) = urdf_link->inertial->izz;

      if (link_inertial_rpy != Vector3d(0., 0., 0.))
      {
        cerr << "Error while processing body '" << urdf_link->name
             << "': rotation of body frames not yet supported. Please rotate the joint frame instead."
             << endl;
      }
    }

    Body base = Body(link_inertial_mass, link_inertial_position, link_inertial_inertia);
    // make floating base
    if (floatingBaseType == +FloatingBaseType::XY_Yaw)
    {
      ROS_DEBUG_STREAM("CREATING PLANAR FLOATING BASE");
      Joint floating_base_joint(SpatialVectord(0., 0., 0., 1., 0., 0.),
                                SpatialVectord(0., 0., 0., 0., 1., 0.),
                                SpatialVectord(0., 0., 1., 0., 0., 0.));

      new_id = rbdl_model.AddBody(model_data, 0, Xtrans(Vector3d(0., 0., 0.)),
                                  floating_base_joint, base, urdf_link->name);
    }
    else if (floatingBaseType == +FloatingBaseType::XYZ_RollPitchYaw)
    {
      ROS_DEBUG_STREAM("Creating 6D FLOATING BASE");
      Joint floating_base_joint(
          SpatialVectord(0., 0., 0., 1., 0., 0.), SpatialVectord(0., 0., 0., 0., 1., 0.),
          SpatialVectord(0., 0., 0., 0., 0., 1.), SpatialVectord(0., 0., 1., 0., 0., 0.),
          SpatialVectord(0., 1., 0., 0., 0., 0.), SpatialVectord(1., 0., 0., 0., 0., 0.));

      new_id = rbdl_model.AddBody(model_data, 0, Xtrans(Vector3d(0., 0., 0.)),
                                  floating_base_joint, base, urdf_link->name);
    }
    else if (floatingBaseType == +FloatingBaseType::XYZ_Quaternion)
    {
      ROS_DEBUG_STREAM("Creating 6D QUATERNION FLOATING BASE");
      Joint root_joint = JointTypeFloatingBase;

      SpatialTransformd root_joint_frame = SpatialTransformd();

      new_id = rbdl_model.AddBody(model_data, 0, root_joint_frame, root_joint, base,
                                  urdf_link->name);
    }
    else
    {
      throw std::runtime_error("Floating base type not supported");
    }
  }
  else
  {
    // Add alias for the root node
    // AS it is necessary to delede the "ROOT" key since the same
    // map is used for reverse lookup, duh
    rbdl_model.mBodyNameMap.erase("ROOT");
    rbdl_model.mBodyNameMap[urdf_link->name] = 0;
    // AS it would be useful to keep root inertial parameters as
    // well, but who knows what kind of side effects this may have.
  }

  return (new_id);
}


int initialize_link(Model &rbdl_model, ModelDatad &model_data, ConstLinkPtr urdf_link,
                    const int parent_id)
{
  int new_id = 0;

  JointPtr urdf_joint = urdf_link->parent_joint;
  // create the joint
  Joint rbdl_joint;
  if ((urdf_joint->type == urdf::Joint::REVOLUTE || urdf_joint->type == urdf::Joint::CONTINUOUS) &&
      !urdf_joint->mimic)
  {
    rbdl_joint = Joint(SpatialVectord(urdf_joint->axis.x, urdf_joint->axis.y,
                                      urdf_joint->axis.z, 0., 0., 0.));
  }
  else if (urdf_joint->type == urdf::Joint::PRISMATIC && !urdf_joint->mimic)
  {
    rbdl_joint = Joint(SpatialVectord(0., 0., 0., urdf_joint->axis.x, urdf_joint->axis.y,
                                      urdf_joint->axis.z));
  }
  else if (urdf_joint->type == urdf::Joint::FIXED || urdf_joint->mimic)
  {
    rbdl_joint = Joint(JointTypeFixed);
  }

  // compute the joint transformation
  Vector3d joint_rpy;
  Vector3d joint_translation;
  urdf_joint->parent_to_joint_origin_transform.rotation.getRPY(joint_rpy[0], joint_rpy[1],
                                                               joint_rpy[2]);
  joint_translation.set(urdf_joint->parent_to_joint_origin_transform.position.x,
                        urdf_joint->parent_to_joint_origin_transform.position.y,
                        urdf_joint->parent_to_joint_origin_transform.position.z);

  SpatialTransformd rbdl_joint_frame =
      Xrot(joint_rpy[0], Vector3d(1., 0., 0.)) * Xrot(joint_rpy[1], Vector3d(0., 1., 0.)) *
      Xrot(joint_rpy[2], Vector3d(0., 0., 1.)) * Xtrans(Vector3d(joint_translation));

  // assemble the body
  Vector3d link_inertial_position;
  Vector3d link_inertial_rpy;
  Matrix3d link_inertial_inertia = Matrix3d::Zero();
  double link_inertial_mass = 0;

  // but only if we actually have inertial data
  if (urdf_link->inertial)
  {
    link_inertial_mass = urdf_link->inertial->mass;

    link_inertial_position.set(urdf_link->inertial->origin.position.x,
                               urdf_link->inertial->origin.position.y,
                               urdf_link->inertial->origin.position.z);
    urdf_link->inertial->origin.rotation.getRPY(link_inertial_rpy[0], link_inertial_rpy[1],
                                                link_inertial_rpy[2]);

    link_inertial_inertia(0, 0) = urdf_link->inertial->ixx;
    link_inertial_inertia(0, 1) = urdf_link->inertial->ixy;
    link_inertial_inertia(0, 2) = urdf_link->inertial->ixz;

    link_inertial_inertia(1, 0) = urdf_link->inertial->ixy;
    link_inertial_inertia(1, 1) = urdf_link->inertial->iyy;
    link_inertial_inertia(1, 2) = urdf_link->inertial->iyz;

    link_inertial_inertia(2, 0) = urdf_link->inertial->ixz;
    link_inertial_inertia(2, 1) = urdf_link->inertial->iyz;
    link_inertial_inertia(2, 2) = urdf_link->inertial->izz;


    if (link_inertial_rpy != Vector3d(0., 0., 0.))
    {
      cerr << "Error while processing body '" << urdf_link->name
           << "': rotation of body frames not yet supported. Please rotate the joint frame instead."
           << endl;
    }
  }
  Body rbdl_body = Body(link_inertial_mass, link_inertial_position, link_inertial_inertia);

  new_id = rbdl_model.AddBody(model_data, parent_id, rbdl_joint_frame, rbdl_joint,
                              rbdl_body, urdf_link->name);

  return (new_id);
}


bool construct_model(Model &rbdl_model, ModelDatad &model_data, ConstLinkPtr urdf_link, int parent_id)
{
  for (unsigned int i = 0; i < urdf_link->child_links.size(); ++i)
  {
    int new_id = initialize_link(rbdl_model, model_data, urdf_link->child_links[i], parent_id);
    construct_model(rbdl_model, model_data, urdf_link->child_links[i], new_id);
  }

  return true;
}


bool construct_model_cuttips(Model &rbdl_model, ModelDatad &model_data, ConstLinkPtr urdf_link,
                             const int parent_id, const std::vector<std::string> &tipLinks)
{
  /* The current link is a desired tip */
  bool tip_found = false;
  for (unsigned int i = 0; i < tipLinks.size() && !tip_found; ++i)
  {
    if (tipLinks[i] == urdf_link->name)
    {
      tip_found = true;
    }
  }

  if (!tip_found)
  {
    for (unsigned int i = 0; i < urdf_link->child_links.size(); ++i)
    {
      /* Check that the branch contains atleast on of the specified links */
      for (unsigned int j = 0; j < tipLinks.size(); ++j)
      {
        if (true == isLinkInUrdfModel(urdf_link->child_links[i], tipLinks[j]))
        {
          int new_id =
              initialize_link(rbdl_model, model_data, urdf_link->child_links[i], parent_id);
          construct_model_cuttips(rbdl_model, model_data, urdf_link->child_links[i],
                                  new_id, tipLinks);
          break;
        }
      }
    }
  }

  return true;
}


// ============================================================
// from URDF
// ============================================================

// basic version
bool URDFReadFromURDF(urdf::Model &urdf_model, Model *model, ModelDatad &model_data,
                      FloatingBaseType floatingBaseType, bool verbose)
{
  boost::shared_ptr<urdf::Link> root(boost::const_pointer_cast<urdf::Link>(urdf_model.getRoot()));

  int new_id = initialize_root(*model, model_data, root, floatingBaseType);

  if (false == construct_model(*model, model_data, root, new_id))
  {
    cerr << "Error constructing model from urdf file." << endl;
    return false;
  }

  model->gravity.set(0., 0., -9.81);
  return true;
}


// cut tips
bool URDFReadFromURDF(urdf::Model &urdf_model, Model *model, ModelDatad &model_data,
                      const std::vector<std::string> &tip_links,
                      FloatingBaseType floatingBaseType, bool verbose)
{
  boost::shared_ptr<const urdf::Link> root(urdf_model.getRoot());

  int new_id = initialize_root(*model, model_data, root, floatingBaseType);

  if (false == construct_model_cuttips(*model, model_data, root, new_id, tip_links))
  {
    cerr << "Error constructing model from urdf file." << endl;
    return false;
  }

  model->gravity.set(0., 0., -9.81);
  return true;
}


// basic version with extra info
bool URDFReadFromURDF(urdf::Model &urdf_model, Model *model, ModelDatad &model_data,
                      FloatingBaseType floatingBaseType, std::vector<std::string> &joint_names,
                      std::vector<double> &position_min, std::vector<double> &position_max,
                      std::vector<double> &vel_min, std::vector<double> &vel_max,
                      std::vector<double> &damping, std::vector<double> &friction,
                      std::vector<double> &max_effort, bool verbose)
{
  bool urdfOK = URDFReadFromURDF(urdf_model, model, model_data, floatingBaseType, verbose);
  bool extraOK = parseExtraInformation(urdf_model, model, joint_names, position_min,
                                       position_max, vel_min, vel_max, damping, friction,
                                       max_effort, floatingBaseType);

  return urdfOK && extraOK;
}


// cut tips + extra info
bool URDFReadFromURDF(urdf::Model &urdf_model, Model *model, ModelDatad &model_data,
                      FloatingBaseType floatingBaseType, const std::vector<std::string> &cut_tips,
                      std::vector<std::string> &joint_names,
                      std::vector<double> &position_min, std::vector<double> &position_max,
                      std::vector<double> &vel_min, std::vector<double> &vel_max,
                      std::vector<double> &damping, std::vector<double> &friction,
                      std::vector<double> &max_effort, bool verbose)
{
  bool urdfOK =
      URDFReadFromURDF(urdf_model, model, model_data, cut_tips, floatingBaseType, verbose);
  bool extraOK = parseExtraInformation(urdf_model, model, joint_names, position_min,
                                       position_max, vel_min, vel_max, damping, friction,
                                       max_effort, floatingBaseType);

  return urdfOK && extraOK;
}


// compatibility
bool URDFReadFromURDF(urdf::Model &urdf_model, Model *model,
                      FloatingBaseType floatingBaseType, std::vector<std::string> &joint_names,
                      std::vector<double> &position_min, std::vector<double> &position_max,
                      std::vector<double> &vel_min, std::vector<double> &vel_max,
                      std::vector<double> &damping, std::vector<double> &friction,
                      std::vector<double> &max_effort, bool verbose)
{
  bool urdfOK =
      URDFReadFromURDF(urdf_model, model, *model->getModelData(), floatingBaseType, verbose);
  bool extraOK = parseExtraInformation(urdf_model, model, joint_names, position_min,
                                       position_max, vel_min, vel_max, damping, friction,
                                       max_effort, floatingBaseType);

  return urdfOK && extraOK;
}



// ============================================================
// from a string
// ============================================================

bool initializeURDFModelFromString(urdf::Model &urdf_model, const char *model_xml_string)
{
  if (!urdf_model.initString(model_xml_string))
  {
    ROS_ERROR("Failed to parse urdf file");
    return false;
  }
  return (true);
}


bool URDFReadFromString(const char *model_xml_string, Model *model, ModelDatad &model_data,
                        FloatingBaseType floatingBaseType, bool verbose)
{
  urdf::Model urdf_model;
  bool result = initializeURDFModelFromString(urdf_model, model_xml_string);

  if (result)
  {
    result = URDFReadFromURDF(urdf_model, model, model_data, floatingBaseType, verbose);
  }

  return result;
}


bool URDFReadFromString(const char *model_xml_string, Model *model, ModelDatad &model_data,
                        FloatingBaseType floatingBaseType, std::vector<std::string> &joint_names,
                        std::vector<double> &position_min, std::vector<double> &position_max,
                        std::vector<double> &vel_min, std::vector<double> &vel_max,
                        std::vector<double> &damping, std::vector<double> &friction,
                        std::vector<double> &max_effort, bool verbose)
{
  urdf::Model urdf_model;
  bool result = initializeURDFModelFromString(urdf_model, model_xml_string);

  if (result)
  {
    result = URDFReadFromURDF(urdf_model, model, model_data, floatingBaseType,
                              joint_names, position_min, position_max, vel_min, vel_max,
                              damping, friction, max_effort, verbose);
  }
  return (result);
}


// compatibility
bool URDFReadFromString(const char *model_xml_string, Model *model,
                        FloatingBaseType floatingBaseType, std::vector<std::string> &joint_names,
                        std::vector<double> &position_min, std::vector<double> &position_max,
                        std::vector<double> &vel_min, std::vector<double> &vel_max,
                        std::vector<double> &damping, std::vector<double> &friction,
                        std::vector<double> &max_effort, bool verbose)
{
  return (URDFReadFromString(model_xml_string, model, *model->getModelData(),
                             floatingBaseType, joint_names, position_min, position_max,
                             vel_min, vel_max, damping, friction, max_effort, verbose));
}


// compatibility
bool URDFReadFromString(const char *model_xml_string, Model *model,
                        FloatingBaseType floatingBaseType, bool verbose)
{
  return URDFReadFromString(model_xml_string, model, *model->getModelData(),
                            floatingBaseType, verbose);
}



// ============================================================
// from the parameter server
// ============================================================

bool initializeURDFModelFromParamServer(urdf::Model &urdf_model)
{
  if (!urdf_model.initParam("robot_description"))
  {
    ROS_ERROR("Failed to parse urdf file");
    return false;
  }
  return (true);
}


bool URDFReadFromParamServer(Model *model, ModelDatad &model_data,
                             FloatingBaseType floatingBaseType, bool verbose)
{
  urdf::Model urdf_model;
  bool result = initializeURDFModelFromParamServer(urdf_model);
  if (result)
  {
    result = URDFReadFromURDF(urdf_model, model, model_data, floatingBaseType, verbose);
  }
  return (result);
}


bool URDFReadFromParamServer(Model *model, ModelDatad &model_data,
                             FloatingBaseType floatingBaseType,
                             std::vector<std::string> &joint_names,
                             std::vector<double> &position_min, std::vector<double> &position_max,
                             std::vector<double> &vel_min, std::vector<double> &vel_max,
                             std::vector<double> &damping, std::vector<double> &friction,
                             std::vector<double> &max_effort, bool verbose)
{
  urdf::Model urdf_model;
  bool result = initializeURDFModelFromParamServer(urdf_model);
  if (result)
  {
    result = URDFReadFromURDF(urdf_model, model, model_data, floatingBaseType,
                              joint_names, position_min, position_max, vel_min, vel_max,
                              damping, friction, max_effort, verbose);
  }
  return (result);
}


// cut tips
bool URDFReadFromParamServer(Model *model, ModelDatad &model_data,
                             const std::vector<std::string> &tipLinks,
                             FloatingBaseType floatingBaseType, bool verbose)
{
  urdf::Model urdf_model;
  bool result = initializeURDFModelFromParamServer(urdf_model);
  if (result)
  {
    result = URDFReadFromURDF(urdf_model, model, model_data, tipLinks, floatingBaseType, verbose);
  }

  return result;
}


bool URDFReadFromParamServer(Model *model, ModelDatad &model_data,
                             FloatingBaseType floatingBaseType,
                             const std::vector<std::string> &tipLinks,
                             std::vector<std::string> &joint_names,
                             std::vector<double> &position_min, std::vector<double> &position_max,
                             std::vector<double> &vel_min, std::vector<double> &vel_max,
                             std::vector<double> &damping, std::vector<double> &friction,
                             std::vector<double> &max_effort, bool verbose)
{
  urdf::Model urdf_model;
  bool result = initializeURDFModelFromParamServer(urdf_model);
  if (result)
  {
    result = URDFReadFromURDF(urdf_model, model, model_data, floatingBaseType, tipLinks,
                              joint_names, position_min, position_max, vel_min, vel_max,
                              damping, friction, max_effort, verbose);
  }
  return (result);
}



// compatibility
bool URDFReadFromParamServer(Model *model, FloatingBaseType floatingBaseType, bool verbose)
{
  return URDFReadFromParamServer(model, *model->getModelData(), floatingBaseType, verbose);
}

// compatibility
bool URDFReadFromParamServer(Model *model, const std::vector<std::string> &tipLinks,
                             FloatingBaseType floatingBaseType, bool verbose)
{
  return URDFReadFromParamServer(model, *model->getModelData(), tipLinks,
                                 floatingBaseType, verbose);
}

// compatibility
bool URDFReadFromParamServer(Model *model, FloatingBaseType floatingBaseType,
                             std::vector<std::string> &joint_names,
                             std::vector<double> &position_min, std::vector<double> &position_max,
                             std::vector<double> &vel_min, std::vector<double> &vel_max,
                             std::vector<double> &damping, std::vector<double> &friction,
                             std::vector<double> &max_effort, bool verbose)
{
  return URDFReadFromParamServer(model, *model->getModelData(), floatingBaseType,
                                 joint_names, position_min, position_max, vel_min,
                                 vel_max, damping, friction, max_effort, verbose);
}

// compatibility
bool URDFReadFromParamServer(Model *model, FloatingBaseType floatingBaseType,
                             const std::vector<std::string> &tipLinks,
                             std::vector<std::string> &joint_names,
                             std::vector<double> &position_min, std::vector<double> &position_max,
                             std::vector<double> &vel_min, std::vector<double> &vel_max,
                             std::vector<double> &damping, std::vector<double> &friction,
                             std::vector<double> &max_effort, bool verbose)
{
  return (URDFReadFromParamServer(model, *model->getModelData(), floatingBaseType, tipLinks,
                                  joint_names, position_min, position_max, vel_min,
                                  vel_max, damping, friction, max_effort, verbose));
}


// ============================================================
// from a file
// ============================================================

bool initializeURDFModelFromFile(urdf::Model &urdf_model, const char *filename)
{
  ifstream model_file(filename);
  if (!model_file)
  {
    cerr << "Error opening file '" << filename << "'." << endl;
    return (false);
  }

  // reserve memory for the contents of the file
  string model_xml_string;
  model_file.seekg(0, std::ios::end);
  model_xml_string.reserve(model_file.tellg());
  model_file.seekg(0, std::ios::beg);
  model_xml_string.assign((std::istreambuf_iterator<char>(model_file)),
                          std::istreambuf_iterator<char>());

  model_file.close();

  return (initializeURDFModelFromString(urdf_model, model_xml_string.c_str()));
}


bool URDFReadFromFile(const char *filename, Model *model, ModelDatad &model_data,
                      FloatingBaseType floatingBaseType, bool verbose)
{
  urdf::Model urdf_model;

  bool result = initializeURDFModelFromFile(urdf_model, filename);

  if (result)
  {
    result = URDFReadFromURDF(urdf_model, model, model_data, floatingBaseType, verbose);
  }

  return (result);
}


bool URDFReadFromFile(const char *filename, Model *model, ModelDatad &model_data,
                      FloatingBaseType floatingBaseType, std::vector<std::string> &joint_names,
                      std::vector<double> &position_min, std::vector<double> &position_max,
                      std::vector<double> &vel_min, std::vector<double> &vel_max,
                      std::vector<double> &damping, std::vector<double> &friction,
                      std::vector<double> &max_effort, bool verbose)
{
  urdf::Model urdf_model;

  bool result = initializeURDFModelFromFile(urdf_model, filename);

  if (result)
  {
    result = URDFReadFromURDF(urdf_model, model, model_data, floatingBaseType,
                              joint_names, position_min, position_max, vel_min, vel_max,
                              damping, friction, max_effort, verbose);
  }
  return (result);
}


// compatibility
bool URDFReadFromFile(const char *filename, Model *model,
                      FloatingBaseType floatingBaseType, bool verbose)
{
  return URDFReadFromFile(filename, model, *model->getModelData(), floatingBaseType, verbose);
}


// compatibility
bool URDFReadFromFile(const char *filename, Model *model, FloatingBaseType floatingBaseType,
                      std::vector<std::string> &joint_names,
                      std::vector<double> &position_min, std::vector<double> &position_max,
                      std::vector<double> &vel_min, std::vector<double> &vel_max,
                      std::vector<double> &damping, std::vector<double> &friction,
                      std::vector<double> &max_effort, bool verbose)
{
  return URDFReadFromFile(filename, model, floatingBaseType, joint_names, position_min,
                          position_max, vel_min, vel_max, damping, friction, max_effort,
                          verbose);
}


///////////////////////////

bool parseExtraInformation(urdf::Model &urdf_model, Model *rbdl_model,
                           std::vector<std::string> &joint_names,
                           std::vector<double> &position_min, std::vector<double> &position_max,
                           std::vector<double> &vel_min, std::vector<double> &vel_max,
                           std::vector<double> &damping, std::vector<double> &friction,
                           std::vector<double> &max_effort, FloatingBaseType floatingBaseType)
{
  int nDof;
  if (floatingBaseType == +FloatingBaseType::XYZ_RollPitchYaw)
  {
    nDof = rbdl_model->qdot_size - 6;
  }
  else if (floatingBaseType == +FloatingBaseType::XYZ_Quaternion)
  {
    nDof = rbdl_model->qdot_size - 6;
  }
  else if (floatingBaseType == +FloatingBaseType::XY_Yaw)
  {
    nDof = rbdl_model->qdot_size - 3;
  }
  else if (floatingBaseType == +FloatingBaseType::FixedBase)
  {
    nDof = rbdl_model->qdot_size;
  }
  else
  {
    throw std::runtime_error("Floating base type not supported");
    return false;
  }
  joint_names.resize(nDof);
  position_min.resize(nDof);
  position_max.resize(nDof);
  vel_min.resize(nDof);
  vel_max.resize(nDof);
  damping.resize(nDof);
  friction.resize(nDof);
  max_effort.resize(nDof);


  std::string rootName = urdf_model.getRoot()->name;

  for (auto it = rbdl_model->mBodyNameMap.begin(); it != rbdl_model->mBodyNameMap.end(); ++it)
  {
    if (it->first != "ROOT" && it->first != rootName)
    {
      boost::shared_ptr<urdf::Link> urdf_link;
      // As the joints, the bodyes with movable joints start ad +1
      // Since the movable bodies have a small numbering and the fixed bodyes have
      // a high numbering will get n joint names in the lower numbers

      std::string body_name = it->first;
      urdf_model.getLink(body_name, urdf_link);

      if ((urdf_link->parent_joint->type != urdf::Joint::FIXED) &&
          !urdf_link->parent_joint->mimic)
      {
        unsigned id;
        if (floatingBaseType == +FloatingBaseType::XYZ_RollPitchYaw)
        {
          id = rbdl_model->GetBodyId(body_name.c_str()) - 7;
        }
        else if (floatingBaseType == +FloatingBaseType::XYZ_Quaternion)
        {
          id = rbdl_model->GetBodyId(body_name.c_str()) - 3;
        }
        else if (floatingBaseType == +FloatingBaseType::XY_Yaw)
        {
          id = rbdl_model->GetBodyId(body_name.c_str()) - 4;
        }
        else if (floatingBaseType == +FloatingBaseType::FixedBase)
        {
          id = rbdl_model->GetBodyId(body_name.c_str()) - 1;
        }
        else
        {
          throw std::runtime_error("Floating base type not supported");
          return false;
        }

        ROS_DEBUG_STREAM("rbdl body name: " << body_name);

        joint_names[id] = (urdf_link->parent_joint->name);
        ROS_DEBUG_STREAM(" names " << urdf_link->parent_joint->name);
        // Store the joint limits position
        if (urdf_link->parent_joint->type != urdf::Joint::CONTINUOUS)
        {
          position_min[id] = (urdf_link->parent_joint->limits->lower);
          position_max[id] = (urdf_link->parent_joint->limits->upper);
        }
        else
        {
          position_min[id] = (-3.14);
          position_max[id] = (3.14);
        }
        // Store the joint limits velocity
        if (urdf_link->parent_joint->type != urdf::Joint::CONTINUOUS)
        {
          vel_min[id] = (-urdf_link->parent_joint->limits->velocity);
          vel_max[id] = (urdf_link->parent_joint->limits->velocity);
        }
        else
        {
          /// Random high value
          vel_min[id] = (-100.0);
          vel_max[id] = (100.0);
        }
        // Store joint damping
        if (urdf_link->parent_joint->dynamics)
        {
          damping[id] = (urdf_link->parent_joint->dynamics->damping);
        }
        else
        {
          damping[id] = (0.0);
        }
        // Store joint friction
        if (urdf_link->parent_joint->dynamics)
        {
          friction[id] = (urdf_link->parent_joint->dynamics->friction);
        }
        else
        {
          friction[id] = (0.0);
        }
        // Store joint max_effor
        if (urdf_link->parent_joint->limits)
        {
          max_effort[id] = (urdf_link->parent_joint->limits->effort);
        }
        else
        {
          max_effort[id] = (0.0);
        }
      }
    }
  }
  ROS_DEBUG_STREAM("parsed the rbdl model succesfully");
  return true;
}
}
}
