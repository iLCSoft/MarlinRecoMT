########################################################
# cmake file for building MarlinRecoMT
# @author Remi Ete, DESY
CMAKE_MINIMUM_REQUIRED( VERSION 3.5 FATAL_ERROR )
########################################################

####################################################################################################

PROJECT( MarlinRecoMT )

# project version
SET( ${PROJECT_NAME}_VERSION_MAJOR 0 )
SET( ${PROJECT_NAME}_VERSION_MINOR 1 )
SET( ${PROJECT_NAME}_VERSION_PATCH 0 )



### DEPENDENCIES ############################################################

FIND_PACKAGE( ILCUTIL COMPONENTS ILCSOFT_CMAKE_MODULES REQUIRED )

# load default settings from ILCSOFT_CMAKE_MODULES
INCLUDE( ilcsoft_default_settings )

FIND_PACKAGE( Marlin REQUIRED )
FIND_PACKAGE( ROOT 6.02 REQUIRED COMPONENTS MathMore TMVA)
FIND_PACKAGE( DD4hep REQUIRED COMPONENTS DDRec )


FOREACH( pkg Marlin ROOT DD4hep )
  IF( ${pkg}_FOUND )
    INCLUDE_DIRECTORIES( SYSTEM ${${pkg}_INCLUDE_DIRS} )
    LINK_LIBRARIES( ${${pkg}_LIBRARIES} ${${pkg}_COMPONENT_LIBRARIES} )
    ADD_DEFINITIONS( ${${pkg}_DEFINITIONS} )
  ENDIF()
ENDFOREACH()

### LIBRARY ##################################################################

MACRO( ADD_MARLINRECOMT_PKG package_name )
    MESSAGE( STATUS "+=> MarlinRecoMT package: ${package_name}" )
    SET( package_dir ${CMAKE_CURRENT_SOURCE_DIR}/source/${package_name} )
    IF( EXISTS ${package_dir}/include )
      MESSAGE( STATUS "+===> include directory: ${package_dir}/include" )
      INCLUDE_DIRECTORIES( BEFORE ${package_dir}/include )
      INSTALL_DIRECTORY( ${package_dir}/include DESTINATION . FILES_MATCHING PATTERN "*.h" )
    ENDIF()
    IF( EXISTS /${package_dir}/src )
      MESSAGE( STATUS "+===> compile sources: ${package_dir}/src" )
      AUX_SOURCE_DIRECTORY( ${package_dir}/src MarlinRecoMT_SRCS )
    ENDIF()
    IF( EXISTS ${package_dir}/plugins )
      MESSAGE( STATUS "+===> compile plugins: ${package_dir}/plugins" )
      AUX_SOURCE_DIRECTORY( ${package_dir}/plugins MarlinRecoMT_SRCS )
    ENDIF()
ENDMACRO()

ADD_MARLINRECOMT_PKG( Utils )
ADD_MARLINRECOMT_PKG( CaloDigi )
ADD_MARLINRECOMT_PKG( TrackDigi )

# create library
ADD_SHARED_LIBRARY( MarlinRecoMT ${MarlinRecoMT_SRCS} )
INSTALL_SHARED_LIBRARY( MarlinRecoMT DESTINATION lib )


# display some variables and write them to cache
DISPLAY_STD_VARIABLES()

# INCLUDE("cmake/clang-cpp-checks.cmake")
