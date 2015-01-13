cmake_minimum_required(VERSION 2.8.3)
PROJECT (rbdl)
# Find and use the system's Eigen3 library
find_package(Eigen REQUIRED)

find_package(catkin REQUIRED COMPONENTS urdf)

INCLUDE_DIRECTORIES (
        include
        ${catkin_INCLUDE_DIRS}
        ${Eigen_INCLUDE_DIRS}
)

SET_TARGET_PROPERTIES ( ${PROJECT_EXECUTABLES} PROPERTIES
		LINKER_LANGUAGE CXX
	)

# Perform the proper linking
SET (CMAKE_SKIP_BUILD_RPATH FALSE)
SET (CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
SET (CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
SET (CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# Options
OPTION (BUILD_STATIC "Build the statically linked library" OFF)
OPTION (BUILD_TESTS "Build the test executables" OFF)
OPTION (RBDL_ENABLE_LOGGING "Enable logging (warning: major impact on performance!)" OFF)
OPTION (RBDL_USE_SIMPLE_MATH "Use slow math instead of the fast Eigen3 library (faster compilation)" OFF)
OPTION (RBDL_STORE_VERSION "Enable storing of version information in the library (requires build from valid repository)" OFF)

SET ( RBDL_SOURCES 
	src/rbdl_version.cc
	src/rbdl_mathutils.cc
	src/rbdl_utils.cc
	src/Contacts.cc
        src/Energy.cc
        src/Dynamics.cc
	src/Logging.cc
	src/Joint.cc
	src/Model.cc
	src/Kinematics.cc
	)

ADD_LIBRARY ( rbdl SHARED ${RBDL_SOURCES} )

#IF (BUILD_TESTS)
	ADD_SUBDIRECTORY ( vendor/UnitTest++ )
	INCLUDE_DIRECTORIES ( vendor/UnitTest++/src )
        ADD_SUBDIRECTORY ( tests )
#ENDIF (BUILD_TESTS)

IF (BUILD_STATIC)
  ADD_LIBRARY ( rbdl-static STATIC ${RBDL_SOURCES} )
  SET_TARGET_PROPERTIES ( rbdl-static PROPERTIES PREFIX "lib")
  SET_TARGET_PROPERTIES ( rbdl-static PROPERTIES OUTPUT_NAME "rbdl")

	INSTALL (TARGETS rbdl-static
	  LIBRARY DESTINATION lib
  	ARCHIVE DESTINATION lib
	)
ENDIF (BUILD_STATIC)


IF (RBDL_STORE_VERSION)
	# Set versioning information that can be queried during runtime
	EXEC_PROGRAM("hg" ${CMAKE_CURRENT_SOURCE_DIR} ARGS "id -i"
		OUTPUT_VARIABLE RBDL_BUILD_REVISION)
	EXEC_PROGRAM("hg" ${CMAKE_CURRENT_SOURCE_DIR} ARGS "branch"
		OUTPUT_VARIABLE RBDL_BUILD_BRANCH)
	SET (RBDL_BUILD_TYPE ${CMAKE_BUILD_TYPE})
ELSE (RBDL_STORE_VERSION)
	SET (RBDL_BUILD_REVISION "unknown")
	SET (RBDL_BUILD_BRANCH "unknown")
	SET (RBDL_BUILD_TYPE "unknown")
ENDIF (RBDL_STORE_VERSION)

CONFIGURE_FILE (
	"${CMAKE_CURRENT_SOURCE_DIR}/include/rbdl/rbdl_config.h.cmake" 
	"${CMAKE_CURRENT_SOURCE_DIR}/include/rbdl/rbdl_config.h"
	)

FILE ( GLOB headers 
	${CMAKE_CURRENT_SOURCE_DIR}/include/rbdl/*.h
	${CMAKE_CURRENT_SOURCE_DIR}/include/rbdl/rbdl_config.h
	)

INSTALL ( FILES ${headers} DESTINATION include/rbdl )

# Setup of SimpleMath install settings
IF (RBDL_USE_SIMPLE_MATH)
	INSTALL ( DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/include/rbdl/SimpleMath"
		DESTINATION include/rbdl
	)
ENDIF (RBDL_USE_SIMPLE_MATH)

catkin_package(
  INCLUDE_DIRS include
  LIBRARIES rbdl RBDL_URDFParser
  DEPENDS Eigen
)

add_library(RBDL_URDFParser addons/urdfreader/RBDL_URDFParser.cpp)
target_link_libraries(RBDL_URDFParser ${catkin_LIBRARIES})


# Installing
INSTALL (TARGETS rbdl RBDL_URDFParser
	LIBRARY DESTINATION lib
	ARCHIVE DESTINATION lib
	)

install(DIRECTORY include/
   DESTINATION ${CATKIN_GLOBAL_INCLUDE_DESTINATION}
   FILES_MATCHING PATTERN "*.h"
)

			
