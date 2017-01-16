if(CATKIN_ENABLE_TESTING)

# Tests
 catkin_add_gtest(GTestBodyTests BodyTests.cpp)
 target_link_libraries(GTestBodyTests ${catkin_LIBRARIES} rbdl)

 catkin_add_gtest(GTestCalcAccelerations CalcAccelerationsTests.cpp)
 target_link_libraries(GTestCalcAccelerations ${catkin_LIBRARIES} rbdl)

 catkin_add_gtest(GTestCalcVelocitiesTests CalcVelocitiesTests.cpp)
 target_link_libraries(GTestCalcVelocitiesTests ${catkin_LIBRARIES} rbdl)

 catkin_add_gtest(GTestCompositeRigidBodyTests CompositeRigidBodyTests.cpp)
 target_link_libraries(GTestCompositeRigidBodyTests ${catkin_LIBRARIES} rbdl)

 catkin_add_gtest(GTestContactsTests ContactsTests.cpp)
 target_link_libraries(GTestContactsTests ${catkin_LIBRARIES} rbdl)

 #catkin_add_gtest(GTestCustomJointMultiBodyTests CustomJointMultiBodyTests.cpp)
 #target_link_libraries(GTestCustomJointMultiBodyTests ${catkin_LIBRARIES} rbdl)

 #catkin_add_gtest(GTestCustomJointSingleBodyTests CustomJointSingleBodyTests.cpp)
 #target_link_libraries(GTestCustomJointSingleBodyTests ${catkin_LIBRARIES} rbdl)

 #catkin_add_gtest(GTestCustomJointTests CustomJointTests.cpp)
 #target_link_libraries(GTestCustomJointTests ${catkin_LIBRARIES} rbdl)

 #catkin_add_gtest(GTestDynamicsTests DynamicsTests.cpp)
 #target_link_libraries(GTestDynamicsTests ${catkin_LIBRARIES} rbdl)

 #catkin_add_gtest(GTestFloatingBaseTests FloatingBaseTests.cpp)
 #target_link_libraries(GTestFloatingBaseTests ${catkin_LIBRARIES} rbdl)

 #catkin_add_gtest(GTestImpulsesTests ImpulsesTests.cpp)
 #target_link_libraries(GTestImpulsesTests ${catkin_LIBRARIES} rbdl)

 #catkin_add_gtest(GTestInverseDynamicsTests InverseDynamicsTests.cpp)
 #target_link_libraries(GTestInverseDynamicsTests ${catkin_LIBRARIES} rbdl)

 #catkin_add_gtest(GTestKinematicsTests KinematicsTests.cpp)
 #target_link_libraries(GTestKinematicsTests ${catkin_LIBRARIES} rbdl)

 #catkin_add_gtest(GTestMathTests MathTests.cpp)
 #target_link_libraries(GTestMathTests ${catkin_LIBRARIES} rbdl)

 #catkin_add_gtest(GTestModelTests ModelTests.cpp)
 #target_link_libraries(GTestModelTests ${catkin_LIBRARIES} rbdl)

 #catkin_add_gtest(GTestMultiDofTests MultiDofTests.cpp)
 #target_link_libraries(GTestMultiDofTests ${catkin_LIBRARIES} rbdl)

 #catkin_add_gtest(GTestSparseFactorizationTests SparseFactorizationTests.cpp)
 #target_link_libraries(GTestSparseFactorizationTests ${catkin_LIBRARIES} rbdl)

 #catkin_add_gtest(GTestSpatialAlgebraTests SpatialAlgebraTests.cpp)
 #target_link_libraries(GTestSpatialAlgebraTests ${catkin_LIBRARIES} rbdl)

 #catkin_add_gtest(GTestTwolegModelTests TwolegModelTests.cpp)
 #target_link_libraries(GTestTwolegModelTests ${catkin_LIBRARIES} rbdl)

 #catkin_add_gtest(GTestUtilsTests UtilsTests.cpp)
 #target_link_libraries(GTestUtilsTests ${catkin_LIBRARIES} rbdl)

endif()
