//Rbdl stuff
#include <rbdl/rbdl.h>
#include <vector>

//THIS FUCKS UP THE LINKING
//bool parseUrdf(urdf::Model &urdf_model, Model &rbdl_model);

/**
  Parse rbdl model from param server
*/
bool parseUrdfParamServerParameters(RigidBodyDynamics::Model &rbdl_model, bool floating_base);

/**
  Parse rbdl model from param server, only getting the branches of the tree that have their tips listed in the vector tips
*/
bool parseUrdfParamServerParameters(RigidBodyDynamics::Model &rbdl_model, std::vector<std::string> tips, bool floating_base);

/**
  Parse rbdl model from param server with aditional parameters
*/
bool parseUrdfParamServerParameters(RigidBodyDynamics::Model &rbdl_model, std::vector<std::string> &joint_names,
                                    std::vector<double> &position_min,  std::vector<double> &position_max,
                                    std::vector<double> &vel_min, std::vector<double> &vel_max,
                                    std::vector<double> &damping, std::vector<double> &friction,
                                    std::vector<double> &max_effort,
                                    bool floating_base);

/**
  Parse rbdl model from param server with aditional parameters, only the subtree specified in the tip links
*/
bool parseUrdfParamServerParameters(RigidBodyDynamics::Model &rbdl_model, std::vector<std::string> &joint_names,
                                    std::vector<double> &position_min,  std::vector<double> &position_max,
                                    std::vector<double> &vel_min,  std::vector<double> &vel_max,
                                    std::vector<double> &damping, std::vector<double> &friction,
                                    std::vector<double> &max_effort,
                                    bool floating_base,
                                    std::vector<std::string> tip_links);

/**
  Parse rbdl model from file
*/
bool parseUrdfFromFile(RigidBodyDynamics::Model &rbdl_model, std::string file_path, bool floating_base);

/**
  Parse rbdl model from file with adition parameters
*/
bool parseUrdfFromFile(RigidBodyDynamics::Model &rbdl_model, std::vector<std::string> &joint_names,
                                    std::vector<double> &position_min,  std::vector<double> &position_max,
                                    std::vector<double> &vel_min, std::vector<double> &vel_max,
                                    std::vector<double> &damping, std::vector<double> &friction,
                                    std::vector<double> &max_effort, bool floating_base,
                                    std::string file_path);

/* Return a subtree that includes specified root, and the specified tips */
RigidBodyDynamics::Model getSubTree(RigidBodyDynamics::Model &rbdl_model, std::vector<std::string> tips, std::string root);

/* Publish the internal kinematic state of the robot */
void publish_link_poses(RigidBodyDynamics::Model &model);

/* Publish the internal COM kinematic state of the robot */
void publish_link_com_poses(RigidBodyDynamics::Model &model);

