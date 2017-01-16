#ifndef RBDL_URDFREADER_H
#define RBDL_URDFREADER_H

#include <ros/ros.h>
#include <urdf/model.h>
#include <rbdl/Model.h>

namespace RigidBodyDynamics {

  struct Model;

  namespace Addons {

    /// @todo urdf should be const

    bool URDFReadFromFile (const char* filename, Model* model, FloatingBaseType floatingBaseType, bool verbose = false);
    bool URDFReadFromString (const char* model_xml_string, Model* model, FloatingBaseType floatingBaseType, bool verbose = false);
    bool URDFReadFromURDF(urdf::Model &urdf_model, Model* model, FloatingBaseType floatingBaseType, bool verbose = false);
    bool URDFReadFromParamServer(Model* model, FloatingBaseType floatingBaseType, bool verbose = false);

    bool URDFReadFromFile (const char* filename, Model* model, FloatingBaseType floatingBaseType,
                           std::vector<std::string> &joint_names,
                           std::vector<double> &position_min,  std::vector<double> &position_max,
                           std::vector<double> &vel_min,  std::vector<double> &vel_max,
                           std::vector<double> &damping, std::vector<double> &friction,
                           std::vector<double> &max_effort,
                           bool verbose = false);

    bool URDFReadFromString (const char* model_xml_string, Model* model, FloatingBaseType floatingBaseType,
                             std::vector<std::string> &joint_names,
                             std::vector<double> &position_min,  std::vector<double> &position_max,
                             std::vector<double> &vel_min,  std::vector<double> &vel_max,
                             std::vector<double> &damping, std::vector<double> &friction,
                             std::vector<double> &max_effort,
                             bool verbose = false);

    bool URDFReadFromURDF(urdf::Model &urdf_model, Model* model, FloatingBaseType floatingBaseType,
                          std::vector<std::string> &joint_names,
                          std::vector<double> &position_min,  std::vector<double> &position_max,
                          std::vector<double> &vel_min,  std::vector<double> &vel_max,
                          std::vector<double> &damping, std::vector<double> &friction,
                          std::vector<double> &max_effort,
                          bool verbose = false);

    bool URDFReadFromParamServer(Model* model, ModelDatad &model_data, FloatingBaseType floatingBaseType,
                                 std::vector<std::string> &joint_names,
                                 std::vector<double> &position_min,  std::vector<double> &position_max,
                                 std::vector<double> &vel_min,  std::vector<double> &vel_max,
                                 std::vector<double> &damping, std::vector<double> &friction,
                                 std::vector<double> &max_effort,
                                 bool verbose = false);

    bool parseExtraInformation(urdf::Model &urdf_model, Model *rbdl_model, std::vector<std::string> &joint_names,
                               std::vector<double> &position_min,  std::vector<double> &position_max,
                               std::vector<double> &vel_min,  std::vector<double> &vel_max,
                               std::vector<double> &damping, std::vector<double> &friction,
                               std::vector<double> &max_effort,
                               FloatingBaseType floatingBaseType);



    // With subtree parsing

    bool URDFReadFromParamServer(Model* model, ModelDatad &model_data, FloatingBaseType floatingBaseType, const std::vector<std::string> &tipLinks,
                                 std::vector<std::string> &joint_names,
                                 std::vector<double> &position_min,  std::vector<double> &position_max,
                                 std::vector<double> &vel_min,  std::vector<double> &vel_max,
                                 std::vector<double> &damping, std::vector<double> &friction,
                                 std::vector<double> &max_effort,
                                 bool verbose = false);
  }
}

/* _RBDL_URDFREADER_H */
#endif
