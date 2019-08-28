#ifndef RBDL_URDFREADER_H
#define RBDL_URDFREADER_H

#include <ros/ros.h>
#include <urdf/model.h>
#include <rbdl/Model.h>

namespace RigidBodyDynamics
{
struct Model;

namespace Addons
{

const double INF = 1e10;

/// @todo urdf should be const
bool initializeURDFModelFromFile(urdf::Model &urdf_model, const char *filename);
bool initializeURDFModelFromString(urdf::Model &urdf_model, const char *model_xml_string);
bool initializeURDFModelFromParamServer(urdf::Model &urdf_model);

bool URDFReadFromFile(const char *filename, Model *model, ModelDatad &model_data,
                      const FloatingBaseType floatingBaseType, const bool verbose = false);
bool URDFReadFromFile(const char *filename, Model *model,
                      const FloatingBaseType floatingBaseType, const bool verbose = false);
bool URDFReadFromString(const char *model_xml_string, Model *model, ModelDatad &model_data,
                        const FloatingBaseType floatingBaseType, const bool verbose = false);
bool URDFReadFromString(const char *model_xml_string, Model *model,
                        const FloatingBaseType floatingBaseType, const bool verbose = false);
bool URDFReadFromURDF(urdf::Model &urdf_model, Model *model, ModelDatad &model_data,
                      const FloatingBaseType floatingBaseType, const bool verbose = false);
bool URDFReadFromURDF(Model *model, ModelDatad &model_data,
                      urdf::LinkConstSharedPtr root,
                      const FloatingBaseType floatingBaseType, const bool verbose = false);
bool URDFReadFromURDF(urdf::Model &urdf_model, Model *model, ModelDatad &model_data,
                      const std::vector<std::string> &tip_links,
                      const FloatingBaseType floatingBaseType, const bool verbose = false);
bool URDFReadFromParamServer(Model *model, ModelDatad &model_data,
                             const FloatingBaseType floatingBaseType,
                             const bool verbose = false);

bool URDFReadFromFile(const char *filename, Model *model, ModelDatad &model_data,
                      const FloatingBaseType floatingBaseType,
                      std::vector<std::string> &joint_names,
                      std::vector<double> &position_min, std::vector<double> &position_max,
                      std::vector<double> &vel_min, std::vector<double> &vel_max,
                      std::vector<double> &damping, std::vector<double> &friction,
                      std::vector<double> &max_effort, const bool verbose = false);

bool URDFReadFromString(const char *model_xml_string, Model *model,
                        ModelDatad &model_data, const FloatingBaseType floatingBaseType,
                        std::vector<std::string> &joint_names,
                        std::vector<double> &position_min, std::vector<double> &position_max,
                        std::vector<double> &vel_min, std::vector<double> &vel_max,
                        std::vector<double> &damping, std::vector<double> &friction,
                        std::vector<double> &max_effort, const bool verbose = false);

bool URDFReadFromURDF(urdf::Model &urdf_model, Model *model, ModelDatad &model_data,
                      urdf::LinkConstSharedPtr root,
                      const FloatingBaseType floatingBaseType,
                      std::vector<std::string> &joint_names,
                      std::vector<double> &position_min, std::vector<double> &position_max,
                      std::vector<double> &vel_min, std::vector<double> &vel_max,
                      std::vector<double> &damping, std::vector<double> &friction,
                      std::vector<double> &max_effort, const bool verbose = false);

bool URDFReadFromURDF(urdf::Model &urdf_model, Model *model, ModelDatad &model_data,
                      const FloatingBaseType floatingBaseType,
                      std::vector<std::string> &joint_names,
                      std::vector<double> &position_min, std::vector<double> &position_max,
                      std::vector<double> &vel_min, std::vector<double> &vel_max,
                      std::vector<double> &damping, std::vector<double> &friction,
                      std::vector<double> &max_effort, const bool verbose = false);

bool URDFReadFromURDF(urdf::Model &urdf_model, Model *model, ModelDatad &model_data,
                      const FloatingBaseType floatingBaseType,
                      const std::vector<std::string> &tip_links,
                      std::vector<std::string> &joint_names,
                      std::vector<double> &position_min, std::vector<double> &position_max,
                      std::vector<double> &vel_min, std::vector<double> &vel_max,
                      std::vector<double> &damping, std::vector<double> &friction,
                      std::vector<double> &max_effort, const bool verbose = false);

bool URDFReadFromParamServer(Model *model, ModelDatad &model_data,
                             const FloatingBaseType floatingBaseType,
                             std::vector<std::string> &joint_names,
                             std::vector<double> &position_min, std::vector<double> &position_max,
                             std::vector<double> &vel_min, std::vector<double> &vel_max,
                             std::vector<double> &damping, std::vector<double> &friction,
                             std::vector<double> &max_effort, const bool verbose = false);

bool URDFReadFromParamServer(Model *model, const std::vector<std::string> &tipLinks,
                             const FloatingBaseType floatingBaseType,
                             const bool verbose = false);

// With subtree parsing

bool URDFReadFromParamServer(Model *model, ModelDatad &model_data,
                             const std::vector<std::string> &tipLinks,
                             const FloatingBaseType floatingBaseType,
                             const bool verbose = false);

bool URDFReadFromParamServer(Model *model, ModelDatad &model_data,
                             const FloatingBaseType floatingBaseType,
                             const std::vector<std::string> &tipLinks,
                             std::vector<std::string> &joint_names,
                             std::vector<double> &position_min, std::vector<double> &position_max,
                             std::vector<double> &vel_min, std::vector<double> &vel_max,
                             std::vector<double> &damping, std::vector<double> &friction,
                             std::vector<double> &max_effort, const bool verbose = false);

// With defautl model data allocation (compatibility)

bool URDFReadFromURDF(urdf::Model &urdf_model, Model *model,
                      const FloatingBaseType floatingBaseType,
                      std::vector<std::string> &joint_names,
                      std::vector<double> &position_min, std::vector<double> &position_max,
                      std::vector<double> &vel_min, std::vector<double> &vel_max,
                      std::vector<double> &damping, std::vector<double> &friction,
                      std::vector<double> &max_effort, const bool verbose = false);


bool URDFReadFromString(const char *model_xml_string, Model *model,
                        const FloatingBaseType floatingBaseType,
                        std::vector<std::string> &joint_names,
                        std::vector<double> &position_min, std::vector<double> &position_max,
                        std::vector<double> &vel_min, std::vector<double> &vel_max,
                        std::vector<double> &damping, std::vector<double> &friction,
                        std::vector<double> &max_effort, const bool verbose = false);

bool URDFReadFromParamServer(Model *model, const FloatingBaseType floatingBaseType,
                             const bool verbose = false);

bool URDFReadFromParamServer(Model *model, const FloatingBaseType floatingBaseType,
                             std::vector<std::string> &joint_names,
                             std::vector<double> &position_min, std::vector<double> &position_max,
                             std::vector<double> &vel_min, std::vector<double> &vel_max,
                             std::vector<double> &damping, std::vector<double> &friction,
                             std::vector<double> &max_effort, const bool verbose = false);

bool URDFReadFromParamServer(Model *model, const FloatingBaseType floatingBaseType,
                             const std::vector<std::string> &tipLinks,
                             std::vector<std::string> &joint_names,
                             std::vector<double> &position_min, std::vector<double> &position_max,
                             std::vector<double> &vel_min, std::vector<double> &vel_max,
                             std::vector<double> &damping, std::vector<double> &friction,
                             std::vector<double> &max_effort, const bool verbose = false);

bool URDFReadFromFile(const char *filename, Model *model, const FloatingBaseType floatingBaseType,
                      std::vector<std::string> &joint_names,
                      std::vector<double> &position_min, std::vector<double> &position_max,
                      std::vector<double> &vel_min, std::vector<double> &vel_max,
                      std::vector<double> &damping, std::vector<double> &friction,
                      std::vector<double> &max_effort, const bool verbose = false);

bool parseExtraInformation(urdf::Model &urdf_model, Model *rbdl_model,
                           std::vector<std::string> &joint_names, std::vector<double> &position_min,
                           std::vector<double> &position_max, std::vector<double> &vel_min,
                           std::vector<double> &vel_max, std::vector<double> &damping,
                           std::vector<double> &friction, std::vector<double> &max_effort,
                           const FloatingBaseType floatingBaseType);
}
}

/* _RBDL_URDFREADER_H */
#endif
