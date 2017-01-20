# $Id: CMakeLists.txt 68058 2013-03-13 14:47:43Z gcosmo $

#----------------------------------------------------------------------------
# Setup the project
#
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(MCSpin_simple)

set(ROOT_INSTALL /home/camelud/GEANT4/ROOT/install)
if(UNIX)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -std=c++11 -O3")        ## Optimize
    set(CMAKE_EXE_LINKER_FLAGS "-s")  ## Strip binary

endif()

message( STATUS ${PROJECT_SOURCE_DIR})

include_directories(${ROOT_INSTALL}/include)
include_directories(${PROJECT_SOURCE_DIR}/include)
#include_directories(${PROJECT_SOURCE_DIR}/libraries)
include_directories(${PROJECT_SOURCE_DIR}/src)
#----------------------------------------------------------------------------
# Locate sources and headers for this project
# NB: headers are included so they will show up in IDEs
#

file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

#file(GLOB sourcesLibre ${PROJECT_SOURCE_DIR}/libraries/*.cpp)
#file(GLOB headersLibre ${PROJECT_SOURCE_DIR}/libraries/*.h)
#file(GLOB sources ${PROJECT_SOURCE_DIR}/libraries/*.cpp)
#file(GLOB headers ${PROJECT_SOURCE_DIR}/libraries/*.h)

add_executable(MCSpin_simple MCSpin_simple.cc ${sources} ${headers} ${sourcesLibre} ${headersLibre} )
target_link_libraries(MCSpin_simple -L${ROOT_INSTALL}/lib -lCore  -lGpad -lHist -lGraf  -lPostscript -lMathCore -rdynamic)

#-L${ROOT_INSTALL}/lib -lCore -lCint -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lm -ldl -rdynamic

#if(CMAKE_COMPILER_IS_GNUCXX)
#    set(CMAKE_CXX_FLAGS "-O3")        ## Optimize
#    set(CMAKE_EXE_LINKER_FLAGS "-s")  ## Strip binary
#endif()

#----------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we
# build B4a. This is so that we can run the executable directly because it
# relies on these scripts being in the current working directory.
#
#set(EXAMPLEB4A_SCRIPTS
#  INIMaterial
#  INISource
#  Sigmae.txt
#  Sigmap.txt
#  )

#foreach(_script ${EXAMPLEB4A_SCRIPTS})
 # configure_file(
 #   ${PROJECT_SOURCE_DIR}/${_script}
  #  ${PROJECT_BINARY_DIR}/${_script}
   # COPYONLY
   # )
#endforeach()


set(SCRIPTS
  INIMacro
fit.gnu
  )

foreach(_script ${SCRIPTS})
 configure_file(
    ${PROJECT_SOURCE_DIR}/${_script}
   ${PROJECT_BINARY_DIR}/${_script}
  COPYONLY
  )
endforeach()
#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS MCSpin_simple DESTINATION bin)
