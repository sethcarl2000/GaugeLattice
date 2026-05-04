#---------------------------------------------------------
#   I copied most of this CMakeLists.txt from: https://cliutils.gitlab.io/modern-cmake/examples/root-dict/README.html 
cmake_minimum_required(VERSION 3.11)

#---------------------------------------------------------
#   Define the project 
project(GaugeLattice VERSION 0.1
        LANGUAGES CXX)  

#---------------------------------------------------------
#   Define this preprocessor macro to enable debug statements 
option(DEBUG "Enable debug printing" OFF)
if (DEBUG)
    add_compile_definitions(DEBUG)
endif()

#---------------------------------------------------------
#   Which include directories to use. 
set(CMAKE_CXX_STANDARD
    20
    CACHE STRING "C++ standard to use")
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_PLATFORM_INDEPENDENT_CODE ON)

#---------------------------------------------------------
#   Find ROOT (necessary)
find_package(ROOT 6.20 REQUIRED COMPONENTS Core RIO Physics)

#---------------------------------------------------------
#   Define sources and headers
set(sources 
    SU3.cxx
)

set(headers 
    ComplexNum.hxx
    GaugeLattice.hxx
    SU3.hxx
)

#---------------------------------------------------------
#   Add application executable
add_executable(main main.cxx ${sources} ${headers})

message("cmake: ROOT_INCLUDE_DIRS: '${ROOT_INCLUDE_DIRS}'")

target_include_directories(main PRIVATE ${CMAKE_CURRENT_SOURCE_DIR} ${ROOT_INCLUDE_DIRS})

target_link_libraries(main PRIVATE ${ROOT_LIBRARIES})

#---------------------------------------------------------
#   Generate the ROOT dictionary
#root_generate_dictionary(G__GaugeLattice ${headers} LINKDEF LinkDef.h)

#---------------------------------------------------------
#   Add Libraries
#add_library(GaugeLattice SHARED ${headers} ${sources} G__GaugeLattice.cxx)
#target_include_directories(GaugeLattice PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")
#target_link_libraries(GaugeLattice PUBLIC ROOT::Core)


#target_include_directories(pick_sieveholes PRIVATE
#    "${CMAKE_CURRENT_SOURCE_DIR}"
#    "${ROOT_INCLUDE_DIRS}"
#)

#target_link_libraries(pick_sieveholes PRIVATE
#    ${ROOT_LIBRARIES}
#    GaugeLattice
#)
