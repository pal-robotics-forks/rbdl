/*
 * RBDL - Rigid Body Dynamics Library
 * Copyright (c) 2011-2016 Martin Felis <martin.felis@iwr.uni-heidelberg.de>
 *
 * Licensed under the zlib license. See LICENSE for more details.
 */

#ifndef RBDL_UTILS_H
#define RBDL_UTILS_H

#include <string>
#include <rbdl/rbdl_config.h>
#include <rbdl/rbdl_math.h>
#include <rbdl/ModelData.h>

namespace RigidBodyDynamics
{
struct Model;

/** \brief Namespace that contains optional helper functions */
namespace Utils
{
/** \brief Creates a human readable overview of the model. */
RBDL_DLLAPI std::string GetModelHierarchy(const Model &model, const ModelDatad &model_data);
/** \brief Creates a human readable overview of the Degrees of Freedom. */
RBDL_DLLAPI std::string GetModelDOFOverview(const Model &model, const ModelDatad &model_data);
/** \brief Creates a human readable overview of the locations of all bodies that have
 * names. */
RBDL_DLLAPI std::string GetNamedBodyOriginsOverview(const Model &model, ModelDatad &model_data);

/** \brief Computes the Center of Mass (COM) and optionally its linear velocity.
 *
 * When only interested in computing the location of the COM you can use
 * NULL as value for com_velocity.
 *
 * \param model The model for which we want to compute the COM
 * \param q The current joint positions
 * \param qdot The current joint velocities
 * \param mass (output) total mass of the model
 * \param com (output) location of the Center of Mass of the model in base coordinates
 * \param com_velocity (optional output) linear velocity of the COM in base coordinates
 * \param angular_momentum (optional output) angular momentum of the model at the COM in
 * base coordinates
 * \param update_kinematics (optional input) whether the kinematics should be updated
 * (defaults to true)
 */
RBDL_DLLAPI void CalcCenterOfMass(const Model &model, ModelDatad &model_data,
                                  const Math::VectorNd &q, const Math::VectorNd &qdot,
                                  double &mass, Math::Vector3d &com,
                                  Math::Vector3d *com_velocity = NULL,
                                  Math::Vector3d *angular_momentum = NULL,
                                  bool update_kinematics = true);

void CalcCenterOfMass(Model &model, const Math::VectorNd &q, const Math::VectorNd &qdot,
                      double &mass, Math::Vector3d &com, Math::Vector3d *com_velocity = NULL,
                      Math::Vector3d *angular_momentum = NULL, bool update_kinematics = true);

Math::Vector3d CalcZeroMomentPoint(Model &model, const Math::VectorNd &q,
                                   const Math::VectorNd &qdot,
                                   const Math::Vector3d &normal = Math::Vector3d(0., 0., 1.),
                                   const Math::Vector3d &point = Math::Vector3d(0., 0., 0.),
                                   bool update_kinematics = true);

/** \brief Computes the potential energy of the full model. */
RBDL_DLLAPI double CalcPotentialEnergy(Model &model, ModelDatad &model_data,
                                       const Math::VectorNd &q, bool update_kinematics = true);

/** \brief Computes the kinetic energy of the full model. */
RBDL_DLLAPI double CalcKineticEnergy(Model &model, ModelDatad &model_data,
                                     const Math::VectorNd &q, const Math::VectorNd &qdot,
                                     bool update_kinematics = true);

RBDL_DLLAPI void CalcZeroMomentPoint(Model &model, ModelDatad &model_data,
                                     const Math::VectorNd &q, const Math::VectorNd &qdot,
                                     Math::Vector3d *zmp = NULL,
                                     const Math::Vector3d &normal = Math::Vector3d(0., 0., 1.),
                                     const Math::Vector3d &point = Math::Vector3d(0., 0., 0.),
                                     bool update_kinematics = true);

/** \brief Computed the Spatial Inertia Matrix of a given point position (w.r.t base_link)
 */
RBDL_DLLAPI void CalcPointSpatialInertiaMatrix(
    Model &model, ModelDatad &model_data, const Math::VectorNd &q, const Math::VectorNd &qdot,
    const Math::Vector3d &point_position, Math::SpatialRigidBodyInertiad &inertia_matrix,
    Math::Vector3d *angular_momentum, bool update_kinematics = true);

void CalcPointSpatialInertiaMatrix(Model &model, const Math::VectorNd &q,
                                   const Math::VectorNd &qdot, unsigned int body_id,
                                   const Math::Vector3d &point_position,
                                   Math::SpatialRigidBodyInertiad &inertia_matrix,
                                   Math::Vector3d *angular_momentum,
                                   bool update_kinematics = true);

void CalcCentroidalInertiaMatrix(Model &model, const Math::VectorNd &q,
                                 const Math::VectorNd &qdot,
                                 Math::SpatialRigidBodyInertiad &inertia_matrix,
                                 Math::Vector3d *angular_momentum,
                                 bool update_kinematics = true);
}
}

/* RBDL_UTILS_H */
#endif
