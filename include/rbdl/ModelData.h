#ifndef _MODEL_DATA_
#define _MODEL_DATA_

#include <rbdl/rbdl_math.h>

/** \brief Namespace for all structures of the RigidBodyDynamics library
*/
namespace RigidBodyDynamics {

template <class T>
struct FixedBodyData{
  /// \brief Transforms spatial quantities expressed for the parent to the
  // fixed body.
  Math::SpatialTransform<T> mBaseTransform;
};

using FixedBodyDataD = FixedBodyData<double>;

template <typename T>
class ModelData{

public:

  // State information
  /// \brief The spatial velocity of the bodies
  std::vector<Math::SpatialVector<T> > v;
  /// \brief The spatial acceleration of the bodies
  std::vector<Math::SpatialVector<T> > a;
  /// \brief The spatial bias acceleration of the bodies
  std::vector<Math::SpatialVector<T> > a_bias;

  /// \brief For computing COM jacobian efficiently
  std::vector<Math::SpatialMatrix<T> > acumulated_mass;

  /// \brief Transformation from the base to bodies reference frame
  std::vector<Math::SpatialTransform<T> > X_base;

  // Joint state variables
  std::vector<Math::SpatialTransform<T> > X_J;
  std::vector<Math::SpatialVector<T> > v_J;
  std::vector<Math::SpatialVector<T> > c_J;

  /// \brief The joint axis for joint i
  std::vector<Math::SpatialVector<T> > S;

  // Special variables for joints with 3 degrees of freedom
  /// \brief Motion subspace for joints with 3 degrees of freedom
  std::vector<Math::Matrix63<T> > multdof3_S;

  // Bodies

  /** \brief Transformation from the parent body to the current body
       * \f[
       *	X_{\lambda(i)} = {}^{i} X_{\lambda(i)}
       * \f]
       */
  std::vector<Math::SpatialTransform<T> > X_lambda;

  // Dynamics variables
  /// \brief The velocity dependent spatial acceleration
  std::vector<Math::SpatialVector<T> > c;


  std::vector<FixedBodyData<T> > mFixedBodiesData;
};

using ModelDatad = ModelData<double>;

}

#endif
