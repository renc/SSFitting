cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
cmake_policy(VERSION 2.6)

##########
# Add OpenMesh
 
set(OPENMESH_VERSION "OpenMesh-2.0-RC3") 
set(OPENMESH_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/external/${OPENMESH_VERSION}) 

if(CMAKE_BUILD_TYPE STREQUAL "Debug")
	set(OPENMESH_LIB_NAME "OpenMeshCore" "OpenMeshTools")
else(CMAKE_BUILD_TYPE STREQUAL "Debug")
	set(OPENMESH_LIB_NAME "OpenMeshCore" "OpenMeshTools")
endif(CMAKE_BUILD_TYPE STREQUAL "Debug")

include_directories("${OPENMESH_DIRECTORY}/src")
link_directories("${OPENMESH_DIRECTORY}/lib") 


###############
# add taucs 

set(TAUCS_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/external/taucs )

include_directories( ${CMAKE_CURRENT_SOURCE_DIR}/external )
include_directories( ${TAUCS_DIRECTORY}/ ) 
link_directories( ${TAUCS_DIRECTORY}/Debug )

# renc: 20150829, if taucs is built in Release, and used here in Debug, we have error: 
#'_ITERATOR_DEBUG_LEVEL': value '0' doesn't match value '2' 

#2>taucs.lib(mulsymatrix.obj) : warning LNK4075: ignoring '/EDITANDCONTINUE' due to '/SAFESEH' specification
#2>taucs.lib(scnrm2.obj) : error LNK2026: module unsafe for SAFESEH image.
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} /SAFESEH:NO")
set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} /SAFESEH:NO")
set(CMAKE_MODULE_LINKER_FLAGS "${CMAKE_MODULE_LINKER_FLAGS} /SAFESEH:NO")
SET(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} /NODEFAULTLIBS:MSVCRT" )

set(TAUSE_LIB taucs.lib)
  

##############
# Add FLTK

#set(FLTK_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/external/fltk-1.1.10 )
#set(FLTK_LIB ${FLTK_DIRECTORY}/lib/fltkd.lib ${FLTK_DIRECTORY}/lib/fltkgld.lib )

set(FLTK_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/external/fltk-1.3.0 )
set(FLTK_LIB ${FLTK_DIRECTORY}/lib/fltkd.lib ${FLTK_DIRECTORY}/lib/fltkgld.lib )


include_directories( ${FLTK_DIRECTORY} )   
 


###############
set(COMMON_DIR ${CMAKE_CURRENT_SOURCE_DIR}/common/) 
include_directories( ${COMMON_DIR} ) 
 
################################
# Add subdirectory

# comment this out if you want to not compile these projects.

add_subdirectory(src) 

