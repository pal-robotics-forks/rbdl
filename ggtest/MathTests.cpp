#include <gtest/gtest.h>
#include <eigen_checks/gtest.h>
#include "rbdl/Logging.h"
#include "rbdl/rbdl_math.h"
#include "rbdl/rbdl_mathutils.h"
#include <iostream>

const double TEST_PREC = 1.0e-14;

using namespace std;
using namespace RigidBodyDynamics::Math;

class MathFixture : public ::testing::Test{
protected:
  virtual void SetUp(){}
  virtual void TearDown(){}
};

TEST (MathTests, GaussElimPivot) {
  ClearLogOutput();

  MatrixNd A;
  A.resize(3,3);
  VectorNd b(3);
  VectorNd x(3);

  A(0,0) = 0; A(0,1) = 2; A(0,2) = 1;
  A(1,0) = 1; A(1,1) = 1; A(1,2) = 5;
  A(2,0) = 0; A(2,1) = 0; A(2,2) = 1;

  b[0] = 1;
  b[1] = 2;
  b[2] = 3;

  VectorNd test_result (3);

  test_result[0] = -12;
  test_result[1] = -1;
  test_result[2] = 3;

  LinSolveGaussElimPivot (A, b, x);

  EXPECT_TRUE(EIGEN_MATRIX_NEAR (test_result, x, TEST_PREC));

  A(0,0) = 0; A(0,1) = -2; A(0,2) = 1;
  A(1,0) = 1; A(1,1) =  1; A(1,2) = 5;
  A(2,0) = 0; A(2,1) =  0; A(2,2) = 1;

  LinSolveGaussElimPivot (A, b, x);
  test_result[0] = -14;
  test_result[1] = 1;
  test_result[2] = 3;

  EXPECT_TRUE(EIGEN_MATRIX_NEAR (test_result, x,  TEST_PREC));
}

TEST (MathTests, Dynamic_1D_initialize_value) {
  VectorNd myvector_10 = VectorNd::Constant ((size_t) 10, 12.);
  VectorNd test_values(10);

  for (unsigned int i = 0; i < 10; i++)
    test_values[i] = 12.;

  EXPECT_TRUE(EIGEN_MATRIX_EQUAL_DOUBLE(test_values, myvector_10));
}

TEST (MathTests,Dynamic_2D_initialize_value) {
  MatrixNd mymatrix_10x10 = MatrixNd::Constant (10, 10, 12.);
  MatrixNd test_values(10, 10);

  for (unsigned int i = 0; i < 10; i++)
    for (unsigned int j = 0; j < 10; j++)
      test_values(i,j) = 12.;

  EXPECT_TRUE(EIGEN_MATRIX_EQUAL_DOUBLE(test_values, mymatrix_10x10));
}

TEST (MathTests, SpatialMatrix_Multiplication) {
  SpatialMatrixd X_1 (
      1.,  2.,  3.,  4.,  5.,  6.,
      11., 12., 13., 14., 15., 16.,
      21., 22., 23., 24., 25., 26.,
      31., 32., 33., 34., 35., 36.,
      41., 42., 43., 44., 45., 46.,
      51., 52., 53., 54., 55., 56.
      );

  SpatialMatrixd X_2 (X_1);

  X_2 *= 2;

  SpatialMatrixd correct_result (
      1442,    1484,    1526,    1568,    1610,    1652,
      4562,    4724,    4886,    5048,    5210,    5372,
      7682,    7964,    8246,    8528,    8810,    9092,
      10802,   11204,   11606,   12008,   12410,   12812,
      13922,   14444,   14966,   15488,   16010,   16532,
      17042,   17684,   18326,   18968,   19610,   20252
      );

  SpatialMatrixd test_result = X_1 * X_2;

  EXPECT_TRUE(EIGEN_MATRIX_NEAR (correct_result, test_result, TEST_PREC));

  // check the *= operator:
  test_result = X_1;
  test_result *= X_2;

  EXPECT_TRUE(EIGEN_MATRIX_NEAR (correct_result, test_result,  TEST_PREC));
}

int main( int argc, char** argv)
{
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
