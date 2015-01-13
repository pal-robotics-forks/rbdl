#include <string.h>
#include <vector>
#include <iostream>

#include <boost/algorithm/string.hpp>
#include <boost/shared_ptr.hpp>

#include <ros/ros.h>
#include <urdf_model/link.h>
#include <urdf/model.h>

/*
#include <tf/transform_listener.h>
#include <tf/transform_broadcaster.h>
#include <tf/tf.h>
#include <tf_conversions/tf_kdl.h>
*/

#include <sensor_msgs/JointState.h>
#include <rbdl/addons/rbdlUrdfParser.h>

using namespace std;
using namespace RigidBodyDynamics;
using namespace RigidBodyDynamics::Math;

typedef boost::shared_ptr<urdf::Link> LinkPtr;
typedef boost::shared_ptr<urdf::Joint> JointPtr;

/* Check if the specified link is inside this urdf submodel */
bool isLinkInUrdfModel(LinkPtr link, std::string tip){
  if(tip == link->name){
    return true;
  }
  else{
    if(link->child_links.empty()){
      return false;
    }
    else{
      bool found = false;
      for(unsigned int i=0; i<link->child_links.size() && !found; ++i){
        found = found || isLinkInUrdfModel(link->child_links[i], tip);
      }
      return found;
    }
  }
}

/*This function will be called recursively adding the Links from the urdf to rbdl */
void constructRBDLfromURDF(Model &rbdl_model, LinkPtr urdf_link, int parent_id, bool floating_base){

  int new_id = 0;
  if(urdf_link->parent_joint.get()){
    JointPtr urdf_joint = urdf_link->parent_joint;
    // create the joint
    Joint rbdl_joint;
    if (urdf_joint->type == urdf::Joint::REVOLUTE || urdf_joint->type == urdf::Joint::CONTINUOUS) {
      rbdl_joint = Joint (SpatialVector (urdf_joint->axis.x, urdf_joint->axis.y, urdf_joint->axis.z, 0., 0., 0.));

    }
    else if (urdf_joint->type == urdf::Joint::PRISMATIC) {
      rbdl_joint = Joint (SpatialVector (0., 0., 0., urdf_joint->axis.x, urdf_joint->axis.y, urdf_joint->axis.z));
      //rbdl_joint = Joint (JointTypeFixed);
    }
    else if (urdf_joint->type == urdf::Joint::FIXED) {
      rbdl_joint = Joint (JointTypeFixed);
    }

    // compute the joint transformation
    Vector3d joint_rpy;
    Vector3d joint_translation;
    urdf_joint->parent_to_joint_origin_transform.rotation.getRPY (joint_rpy[0], joint_rpy[1], joint_rpy[2]);
    joint_translation.set (
          urdf_joint->parent_to_joint_origin_transform.position.x,
          urdf_joint->parent_to_joint_origin_transform.position.y,
          urdf_joint->parent_to_joint_origin_transform.position.z
          );
    SpatialTransform rbdl_joint_frame =
        Xrot (joint_rpy[0], Vector3d (1., 0., 0.))
        * Xrot (joint_rpy[1], Vector3d (0., 1., 0.))
        * Xrot (joint_rpy[2], Vector3d (0., 0., 1.))
        * Xtrans (Vector3d (
                    joint_translation
                    ));

    // assemble the body
    Vector3d link_inertial_position;
    Vector3d link_inertial_rpy;
    Matrix3d link_inertial_inertia = Matrix3d::Zero();
    double link_inertial_mass;

    // but only if we actually have inertial data
    if (urdf_link->inertial) {
      link_inertial_mass = urdf_link->inertial->mass;

      link_inertial_position.set (
            urdf_link->inertial->origin.position.x,
            urdf_link->inertial->origin.position.y,
            urdf_link->inertial->origin.position.z
            );
      urdf_link->inertial->origin.rotation.getRPY (link_inertial_rpy[0], link_inertial_rpy[1], link_inertial_rpy[2]);

      link_inertial_inertia(0,0) = urdf_link->inertial->ixx;
      link_inertial_inertia(0,1) = urdf_link->inertial->ixy;
      link_inertial_inertia(0,2) = urdf_link->inertial->ixz;

      link_inertial_inertia(1,0) = urdf_link->inertial->ixy;
      link_inertial_inertia(1,1) = urdf_link->inertial->iyy;
      link_inertial_inertia(1,2) = urdf_link->inertial->iyz;

      link_inertial_inertia(2,0) = urdf_link->inertial->ixz;
      link_inertial_inertia(2,1) = urdf_link->inertial->iyz;
      link_inertial_inertia(2,2) = urdf_link->inertial->izz;


      if (link_inertial_rpy != Vector3d (0., 0., 0.)) {
        cerr << "Error while processing body '" << urdf_link->name <<
                "': rotation of body frames not yet supported. Please rotate the joint frame instead." << endl;
      }
    }
    Body rbdl_body = Body (link_inertial_mass, link_inertial_position, link_inertial_inertia);

    new_id = rbdl_model.AddBody (parent_id, rbdl_joint_frame, rbdl_joint, rbdl_body, urdf_link->name);
  }
  else
    //If a floating base is desired we need to compute the inertia data of the floating base link,
    //if the floating base has no inertial details, then it makes no sense to make it floating
    //base (Like kdl). (The root is the only link allowed to be floating base

    if(floating_base){
      double link_inertial_mass;
      Vector3d link_inertial_position;
      Vector3d link_inertial_rpy;
      Matrix3d link_inertial_inertia = Matrix3d::Zero();

      if (urdf_link->inertial) {
        link_inertial_mass = urdf_link->inertial->mass;

        link_inertial_position.set (
              urdf_link->inertial->origin.position.x,
              urdf_link->inertial->origin.position.y,
              urdf_link->inertial->origin.position.z
              );
        urdf_link->inertial->origin.rotation.getRPY (link_inertial_rpy[0], link_inertial_rpy[1], link_inertial_rpy[2]);

        link_inertial_inertia(0,0) = urdf_link->inertial->ixx;
        link_inertial_inertia(0,1) = urdf_link->inertial->ixy;
        link_inertial_inertia(0,2) = urdf_link->inertial->ixz;

        link_inertial_inertia(1,0) = urdf_link->inertial->ixy;
        link_inertial_inertia(1,1) = urdf_link->inertial->iyy;
        link_inertial_inertia(1,2) = urdf_link->inertial->iyz;

        link_inertial_inertia(2,0) = urdf_link->inertial->ixz;
        link_inertial_inertia(2,1) = urdf_link->inertial->iyz;
        link_inertial_inertia(2,2) = urdf_link->inertial->izz;


        if (link_inertial_rpy != Vector3d (0., 0., 0.)) {
          cerr << "Error while processing body '" << urdf_link->name <<
                  "': rotation of body frames not yet supported. Please rotate the joint frame instead." << endl;
        }
      }

      Body base = Body (link_inertial_mass, link_inertial_position, link_inertial_inertia);
      //make floating base
      new_id = rbdl_model.SetFloatingBaseBody(base, urdf_link->name);
    }

  for(unsigned int i=0; i< urdf_link->child_links.size(); ++i){
    constructRBDLfromURDF(rbdl_model, urdf_link->child_links[i], new_id, false);
  }
}

/*This function will be called recursively adding the Links from the urdf to rbdl, But only with the tips that are specified inside
  the tips vectors*/

void constructRBDLfromURDF(Model &rbdl_model, LinkPtr urdf_link, int parent_id, std::vector<std::string> tips,
                           bool floating_base){

  int new_id = 0;
  if(urdf_link->parent_joint.get()){
    JointPtr urdf_joint = urdf_link->parent_joint;
    // create the joint
    Joint rbdl_joint;
    if (urdf_joint->type == urdf::Joint::REVOLUTE || urdf_joint->type == urdf::Joint::CONTINUOUS) {
      rbdl_joint = Joint (SpatialVector (urdf_joint->axis.x, urdf_joint->axis.y, urdf_joint->axis.z, 0., 0., 0.));

    }
    else if (urdf_joint->type == urdf::Joint::PRISMATIC) {
      //	rbdl_joint = Joint (SpatialVector (0., 0., 0., urdf_joint->axis.x, urdf_joint->axis.y, urdf_joint->axis.z));
      rbdl_joint = Joint (JointTypeFixed);
    }
    else if (urdf_joint->type == urdf::Joint::FIXED) {
      rbdl_joint = Joint (JointTypeFixed);
    }

    // compute the joint transformation
    Vector3d joint_rpy;
    Vector3d joint_translation;
    urdf_joint->parent_to_joint_origin_transform.rotation.getRPY (joint_rpy[0], joint_rpy[1], joint_rpy[2]);
    joint_translation.set (
          urdf_joint->parent_to_joint_origin_transform.position.x,
          urdf_joint->parent_to_joint_origin_transform.position.y,
          urdf_joint->parent_to_joint_origin_transform.position.z
          );
    SpatialTransform rbdl_joint_frame =
        Xrot (joint_rpy[0], Vector3d (1., 0., 0.))
        * Xrot (joint_rpy[1], Vector3d (0., 1., 0.))
        * Xrot (joint_rpy[2], Vector3d (0., 0., 1.))
        * Xtrans (Vector3d (
                    joint_translation
                    ));

    // assemble the body
    Vector3d link_inertial_position;
    Vector3d link_inertial_rpy;
    Matrix3d link_inertial_inertia = Matrix3d::Zero();
    double link_inertial_mass;

    // but only if we actually have inertial data
    if (urdf_link->inertial) {
      link_inertial_mass = urdf_link->inertial->mass;

      link_inertial_position.set (
            urdf_link->inertial->origin.position.x,
            urdf_link->inertial->origin.position.y,
            urdf_link->inertial->origin.position.z
            );
      urdf_link->inertial->origin.rotation.getRPY (link_inertial_rpy[0], link_inertial_rpy[1], link_inertial_rpy[2]);

      link_inertial_inertia(0,0) = urdf_link->inertial->ixx;
      link_inertial_inertia(0,1) = urdf_link->inertial->ixy;
      link_inertial_inertia(0,2) = urdf_link->inertial->ixz;

      link_inertial_inertia(1,0) = urdf_link->inertial->ixy;
      link_inertial_inertia(1,1) = urdf_link->inertial->iyy;
      link_inertial_inertia(1,2) = urdf_link->inertial->iyz;

      link_inertial_inertia(2,0) = urdf_link->inertial->ixz;
      link_inertial_inertia(2,1) = urdf_link->inertial->iyz;
      link_inertial_inertia(2,2) = urdf_link->inertial->izz;


      if (link_inertial_rpy != Vector3d (0., 0., 0.)) {
        cerr << "Error while processing body '" << urdf_link->name <<
                "': rotation of body frames not yet supported. Please rotate the joint frame instead." << endl;
      }
    }
    Body rbdl_body = Body (link_inertial_mass, link_inertial_position, link_inertial_inertia);

    new_id = rbdl_model.AddBody (parent_id, rbdl_joint_frame, rbdl_joint, rbdl_body, urdf_link->name);
  }
  else
    //If a floating base is desired we need to compute the inertia data of the floating base link,
    //if the floating base has no inertial details, then it makes no sense to make it floating
    //base (Like kdl). (The root is the only link allowed to be floating base

    if(floating_base){
      double link_inertial_mass;
      Vector3d link_inertial_position;
      Vector3d link_inertial_rpy;
      Matrix3d link_inertial_inertia = Matrix3d::Zero();

      if (urdf_link->inertial) {
        link_inertial_mass = urdf_link->inertial->mass;

        link_inertial_position.set (
              urdf_link->inertial->origin.position.x,
              urdf_link->inertial->origin.position.y,
              urdf_link->inertial->origin.position.z
              );
        urdf_link->inertial->origin.rotation.getRPY (link_inertial_rpy[0], link_inertial_rpy[1], link_inertial_rpy[2]);

        link_inertial_inertia(0,0) = urdf_link->inertial->ixx;
        link_inertial_inertia(0,1) = urdf_link->inertial->ixy;
        link_inertial_inertia(0,2) = urdf_link->inertial->ixz;

        link_inertial_inertia(1,0) = urdf_link->inertial->ixy;
        link_inertial_inertia(1,1) = urdf_link->inertial->iyy;
        link_inertial_inertia(1,2) = urdf_link->inertial->iyz;

        link_inertial_inertia(2,0) = urdf_link->inertial->ixz;
        link_inertial_inertia(2,1) = urdf_link->inertial->iyz;
        link_inertial_inertia(2,2) = urdf_link->inertial->izz;


        if (link_inertial_rpy != Vector3d (0., 0., 0.)) {
          cerr << "Error while processing body '" << urdf_link->name <<
                  "': rotation of body frames not yet supported. Please rotate the joint frame instead." << endl;
        }
      }

      Body base = Body (link_inertial_mass, link_inertial_position, link_inertial_inertia);
      //make floating base
      new_id = rbdl_model.SetFloatingBaseBody(base, urdf_link->name);
    }

  /* The current link is a desired tip */
  bool tip_found = false;
  for(unsigned int i=0; i<tips.size() && !tip_found; ++i){
    if(tips[i] == urdf_link->name){
      tip_found = true;
    }
  }

  if(!tip_found){
    for(unsigned int i=0; i< urdf_link->child_links.size(); ++i){
      /* Check that the branch contains atleast on of the specified links */
      bool found = false;
      for(unsigned int j=0; j<tips.size() && !found; ++j){
        found = isLinkInUrdfModel(urdf_link->child_links[i], tips[j]);
      }
      if(found){
        constructRBDLfromURDF(rbdl_model, urdf_link->child_links[i], new_id, tips, false);
      }
    }
  }
}


// WTF Happened here, that in was not found in the .so when linked?
bool parseUrdf(urdf::Model &urdf_model, Model &rbdl_model, std::vector<std::string> tips, bool floating_base) {

  boost::shared_ptr<urdf::Link> root(boost::const_pointer_cast<urdf::Link>(urdf_model.getRoot()));
  constructRBDLfromURDF(rbdl_model, root, 0, tips, floating_base);
  return true;
}

bool parseUrdf(urdf::Model &urdf_model, Model &rbdl_model, bool floating_base) {

  boost::shared_ptr<urdf::Link> root(boost::const_pointer_cast<urdf::Link>(urdf_model.getRoot()));
  constructRBDLfromURDF(rbdl_model, root, 0, floating_base);
  return true;
}

bool parseUrdfParamServer(RigidBodyDynamics::Model &rbdl_model, std::vector<std::string> &joint_names, bool floating_base){

  ros::NodeHandle n;
  std::string urdf_name, full_urdf_xml;

  n.param("urdf_xml_model",urdf_name,std::string("robot_description"));
  n.searchParam(urdf_name,full_urdf_xml);

  urdf::Model urdf_model;
  if (!urdf_model.initParam(full_urdf_xml)){
    ROS_ERROR("Failed to parse urdf file");
    return false;
  }
  ROS_INFO("Successfully parsed urdf file");
  bool res = parseUrdf(urdf_model, rbdl_model, floating_base);
  //We need a link->name to joint->name map;
  for(unsigned int i=0; i< rbdl_model.dof_count; ++i){
    boost::shared_ptr<urdf::Link> urdf_link;
    //As the joints, the bodyes with movable joints start ad +1
    //Since the movable bodies have a small numbering and the fixed bodyes have
    //a high numbering will get n joint names in the lower numbers
    std::string body_name = rbdl_model.GetBodyName(i+1);
    urdf_model.getLink(body_name, urdf_link);
    joint_names.push_back(urdf_link->parent_joint->name);
  }

  ROS_INFO_STREAM("parsed the urdf into rbdl succesfully");
  return res;
}

bool parseUrdfParamServerParameters(RigidBodyDynamics::Model &rbdl_model, bool floating_base){
  ros::NodeHandle n;
  std::string urdf_name, full_urdf_xml;

  ROS_INFO_STREAM("parsing urdf into rbdl from param server");

  n.param("urdf_xml_model",urdf_name,std::string("robot_description"));
  n.searchParam(urdf_name,full_urdf_xml);
  urdf::Model urdf_model;

  if (!urdf_model.initParam(full_urdf_xml)){
    ROS_ERROR("Failed to parse urdf file");
    return false;
  }

  ROS_INFO("Successfully parsed urdf file");
  bool res = parseUrdf(urdf_model, rbdl_model, floating_base);
  return res;
}

bool parseUrdfParamServerParameters(RigidBodyDynamics::Model &rbdl_model, std::vector<std::string> tips, bool floating_base){
  ros::NodeHandle n;
  std::string urdf_name, full_urdf_xml;

  ROS_INFO_STREAM("parsing urdf into rbdl from param server");

  n.param("urdf_xml_model",urdf_name,std::string("robot_description"));
  n.searchParam(urdf_name,full_urdf_xml);
  urdf::Model urdf_model;

  if (!urdf_model.initParam(full_urdf_xml)){
    ROS_ERROR("Failed to parse urdf file");
    return false;
  }

  ROS_INFO("Successfully parsed urdf file");
  bool res = parseUrdf(urdf_model, rbdl_model, tips, floating_base);
  return res;
}


bool parseUrdfParamServerParameters(RigidBodyDynamics::Model &rbdl_model, std::vector<std::string> &joint_names,
                                    std::vector<double> &position_min,  std::vector<double> &position_max,
                                    std::vector<double> &vel_min,  std::vector<double> &vel_max,
                                    std::vector<double> &damping, std::vector<double> &friction,
                                    std::vector<double> &max_effort,
                                    bool floating_base){

  ros::NodeHandle n;
  std::string urdf_name, full_urdf_xml;

  ROS_INFO_STREAM("parsing urdf into rbdl from param server");

  n.param("urdf_xml_model",urdf_name,std::string("robot_description"));
  n.searchParam(urdf_name,full_urdf_xml);
  urdf::Model urdf_model;

  if (!urdf_model.initParam(full_urdf_xml)){
    ROS_ERROR("Failed to parse urdf file");
    return false;
  }

  ROS_INFO("Successfully parsed urdf file");
  bool res = parseUrdf(urdf_model, rbdl_model, floating_base);

  //We need a link->name to joint->name map;

  unsigned int i=0;
  if(floating_base) i=6; //Floating base has been added

  for(; i< rbdl_model.dof_count; ++i){
    boost::shared_ptr<urdf::Link> urdf_link;
    //As the joints, the bodyes with movable joints start ad +1
    //Since the movable bodies have a small numbering and the fixed bodyes have
    //a high numbering will get n joint names in the lower numbers
    std::string body_name = rbdl_model.GetBodyName(i+1);

    ROS_DEBUG_STREAM("rbdl body name: "<<body_name);

    urdf_model.getLink(body_name, urdf_link);
    joint_names.push_back(urdf_link->parent_joint->name);
    ROS_DEBUG_STREAM(" names "<<urdf_link->parent_joint->name);
    //Store the joint limits position
    if ( urdf_link->parent_joint->type != urdf::Joint::CONTINUOUS ) {
      position_min.push_back(urdf_link->parent_joint->limits->lower);
      position_max.push_back(urdf_link->parent_joint->limits->upper);
    }
    else{
      position_min.push_back(-3.14);
      position_max.push_back(3.14);
    }
    //Store the joint limits velocity
    if ( urdf_link->parent_joint->type != urdf::Joint::CONTINUOUS ) {
      vel_min.push_back(-urdf_link->parent_joint->limits->velocity);
      vel_max.push_back(urdf_link->parent_joint->limits->velocity);
    }
    else{
      /// Random high value
      vel_min.push_back(-100.0);
      vel_max.push_back(100.0);
    }
    //Store joint damping
    if (urdf_link->parent_joint->dynamics) {
      damping.push_back(urdf_link->parent_joint->dynamics->damping);
    }
    else{
      damping.push_back(0.0);
    }
    //Store joint friction
    if (urdf_link->parent_joint->dynamics) {
      friction.push_back(urdf_link->parent_joint->dynamics->friction);
    }
    else{
      friction.push_back(0.0);
    }
    //Store joint max_effor
    if (urdf_link->parent_joint->limits) {
      max_effort.push_back(urdf_link->parent_joint->limits->effort);
    }
    else{
      max_effort.push_back(0.0);
    }

  }
  ROS_INFO_STREAM("parsed the urdf into rbdl succesfully");
  return res;

}


bool parseUrdfParamServerParameters(RigidBodyDynamics::Model &rbdl_model, std::vector<std::string> &joint_names,
                                    std::vector<double> &position_min,  std::vector<double> &position_max,
                                    std::vector<double> &vel_min,  std::vector<double> &vel_max,
                                    std::vector<double> &damping, std::vector<double> &friction,
                                    std::vector<double> &max_effort,
                                    bool floating_base,
                                    std::vector<std::string> tip_links){

  ros::NodeHandle n;
  std::string urdf_name, full_urdf_xml;

  ROS_INFO_STREAM("parsing urdf into rbdl from param server");

  n.param("urdf_xml_model",urdf_name,std::string("robot_description"));
  n.searchParam(urdf_name,full_urdf_xml);
  urdf::Model urdf_model;

  if (!urdf_model.initParam(full_urdf_xml)){
    ROS_ERROR("Failed to parse urdf file");
    return false;
  }

  ROS_INFO("Successfully parsed urdf file");
  bool res = parseUrdf(urdf_model, rbdl_model, tip_links, floating_base);

  //We need a link->name to joint->name map;

  unsigned int i=0;
  if(floating_base) i=6; //Floating base has been added

  for(; i< rbdl_model.dof_count; ++i){
    boost::shared_ptr<urdf::Link> urdf_link;
    //As the joints, the bodyes with movable joints start ad +1
    //Since the movable bodies have a small numbering and the fixed bodyes have
    //a high numbering will get n joint names in the lower numbers
    std::string body_name = rbdl_model.GetBodyName(i+1);

    ROS_DEBUG_STREAM("rbdl body name: "<<body_name);

    urdf_model.getLink(body_name, urdf_link);
    joint_names.push_back(urdf_link->parent_joint->name);
    ROS_DEBUG_STREAM(" names "<<urdf_link->parent_joint->name);
    //Store the joint limits position
    if ( urdf_link->parent_joint->type != urdf::Joint::CONTINUOUS ) {
      position_min.push_back(urdf_link->parent_joint->limits->lower);
      position_max.push_back(urdf_link->parent_joint->limits->upper);
    }
    else{
      position_min.push_back(-3.14);
      position_max.push_back(3.14);
    }
    //Store the joint limits velocity
    if ( urdf_link->parent_joint->type != urdf::Joint::CONTINUOUS ) {
      vel_min.push_back(-urdf_link->parent_joint->limits->velocity);
      vel_max.push_back(urdf_link->parent_joint->limits->velocity);
    }
    else{
      /// Random high value
      vel_min.push_back(-100.0);
      vel_max.push_back(100.0);
    }
    //Store joint damping
    if (urdf_link->parent_joint->dynamics) {
      damping.push_back(urdf_link->parent_joint->dynamics->damping);
    }
    else{
      damping.push_back(0.0);
    }
    //Store joint friction
    if (urdf_link->parent_joint->dynamics) {
      friction.push_back(urdf_link->parent_joint->dynamics->friction);
    }
    else{
      friction.push_back(0.0);
    }
    //Store joint max_effor
    if (urdf_link->parent_joint->limits) {
      max_effort.push_back(urdf_link->parent_joint->limits->effort);
    }
    else{
      max_effort.push_back(0.0);
    }

  }
  ROS_INFO_STREAM("parsed the urdf into rbdl succesfully");
  return res;

}

bool parseUrdfFromFile(RigidBodyDynamics::Model &rbdl_model, std::string file_path, bool floating_base){
  ROS_INFO_STREAM("parsing urdf into rbdl from file");

  urdf::Model urdf_model;
  if (!urdf_model.initFile(file_path)){
    ROS_ERROR("Failed to parse urdf file");
    return false;
  }
  ROS_INFO("Successfully parsed urdf file");

  return parseUrdf(urdf_model, rbdl_model, floating_base);

}

bool parseUrdfFromFile(RigidBodyDynamics::Model &rbdl_model, std::vector<std::string> &joint_names,
                       std::vector<double> &position_min,  std::vector<double> &position_max,
                       std::vector<double> &vel_min, std::vector<double> &vel_max,
                       std::vector<double> &damping, std::vector<double> &friction,
                       std::vector<double> &max_effort, bool floating_base,
                       std::string file_path){

  ROS_INFO_STREAM("parsing urdf into rbdl from file");

  urdf::Model urdf_model;

  if (!urdf_model.initFile(file_path)){
    ROS_ERROR("Failed to parse urdf file");
    return false;
  }

  ROS_INFO("Successfully parsed urdf file");
  bool res = parseUrdf(urdf_model, rbdl_model, floating_base);

  unsigned int i=0;
  if(floating_base) i=6; //Floating base has been added

  for(; i< rbdl_model.dof_count; ++i){
    boost::shared_ptr<urdf::Link> urdf_link;
    //As the joints, the bodyes with movable joints start ad +1
    //Since the movable bodies have a small numbering and the fixed bodyes have
    //a high numbering will get n joint names in the lower numbers
    std::string body_name = rbdl_model.GetBodyName(i+1);

    ROS_DEBUG_STREAM("rbdl body name: "<<body_name);

    urdf_model.getLink(body_name, urdf_link);
    joint_names.push_back(urdf_link->parent_joint->name);
    ROS_DEBUG_STREAM(" names "<<urdf_link->parent_joint->name);
    //Store the joint limits position
    if ( urdf_link->parent_joint->type != urdf::Joint::CONTINUOUS ) {
      position_min.push_back(urdf_link->parent_joint->limits->lower);
      position_max.push_back(urdf_link->parent_joint->limits->upper);
    }
    else{
      position_min.push_back(-3.14);
      position_max.push_back(3.14);
    }
    //Store the joint limits velocity
    if ( urdf_link->parent_joint->type != urdf::Joint::CONTINUOUS ) {
      vel_min.push_back(-urdf_link->parent_joint->limits->velocity);
      vel_max.push_back(urdf_link->parent_joint->limits->velocity);
    }
    else{
      /// Random high value
      vel_min.push_back(-100.0);
      vel_max.push_back(100.0);
    }
    //Store joint damping
    if (urdf_link->parent_joint->dynamics) {
      damping.push_back(urdf_link->parent_joint->dynamics->damping);
    }
    else{
      damping.push_back(0.0);
    }
    //Store joint friction
    if (urdf_link->parent_joint->dynamics) {
      friction.push_back(urdf_link->parent_joint->dynamics->friction);
    }
    else{
      friction.push_back(0.0);
    }
    //Store joint max_effor
    if (urdf_link->parent_joint->limits) {
      max_effort.push_back(urdf_link->parent_joint->limits->effort);
    }
    else{
      max_effort.push_back(0.0);
    }

  }
  ROS_INFO_STREAM("parsed the urdf into rbdl succesfully");
  return res;

}


RigidBodyDynamics::Model getSubTree(RigidBodyDynamics::Model &rbdl_model, std::vector<std::string> tips, std::string root){
  RigidBodyDynamics::Model subtree;

  ROS_ERROR_STREAM("not implemented");

  return subtree;
}

/* Function to publish the poses of all the links inside the model */
/*
void publish_link_poses(RigidBodyDynamics::Model &model){
   ros::NodeHandle n;
  static tf::TransformBroadcaster br;

  VectorNd Q (VectorNd::Zero(model.dof_count));
 // UpdateKinematicsCustom (model, &Q, NULL, NULL);

  for (unsigned int body_id = 0; body_id < model.mBodies.size(); body_id++) {
    std::string body_name = model.GetBodyName (body_id);
    if (body_name.size() == 0)
      continue;

    Vector3d position = CalcBodyToBaseCoordinates (model, Q, body_id, Vector3d (0., 0., 0.), false);
    Eigen::Matrix3d orientation = RigidBodyDynamics::CalcBodyWorldOrientation(model, Q, body_id,  false).transpose();
    Eigen::Quaterniond quaternion (orientation);


    tf::Transform transform;
    transform.setOrigin( tf::Vector3(position[0],
                                     position[1],
                                     position[2]) );

    transform.setRotation( tf::Quaternion(quaternion.x(),
                                          quaternion.y(),
                                          quaternion.z(),
                                          quaternion.w()) );
    br.sendTransform(tf::StampedTransform(transform, ros::Time::now(), "world", body_name + "rbdl"));

  }
}


void publish_link_com_poses(RigidBodyDynamics::Model &model){
   ros::NodeHandle n;
  static tf::TransformBroadcaster br;

  VectorNd Q (VectorNd::Zero(model.dof_count));
 // UpdateKinematicsCustom (model, &Q, NULL, NULL);

  for (unsigned int body_id = 0; body_id < model.mBodies.size(); body_id++) {
    std::string body_name = model.GetBodyName (body_id);
    if (body_name.size() == 0)
      continue;

    Vector3d position = RigidBodyDynamics::CalcBodyToBaseCoordinates(model, Q, body_id, model.mBodies[body_id].mCenterOfMass, false);

    //Vector3d position = CalcBodyToBaseCoordinates (model, Q, body_id, Vector3d (0., 0., 0.), false);
    Eigen::Matrix3d orientation = RigidBodyDynamics::CalcBodyWorldOrientation(model, Q, body_id,  false).transpose();
    Eigen::Quaterniond quaternion (orientation);


    tf::Transform transform;
    transform.setOrigin( tf::Vector3(position[0],
                                     position[1],
                                     position[2]) );

    transform.setRotation( tf::Quaternion(quaternion.x(),
                                          quaternion.y(),
                                          quaternion.z(),
                                          quaternion.w()) );
    br.sendTransform(tf::StampedTransform(transform, ros::Time::now(), "world", body_name + "_COM_rbdl"));

  }
}
*/
