cmake_minimum_required(VERSION 3.15)
project(libcdgbs)

set(CMAKE_CXX_STANDARD 20)

# ---- Compiler Flags ----
if (CMAKE_CXX_COMPILER_ID MATCHES "GNU")
  # Common warning and PIC flags
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wno-volatile -Wno-reorder -fPIC")

  # Debug build flags
  set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g -DDEBUG")
  if (NOT WIN32 AND NOT APPLE)
    # Enable AddressSanitizer on non-Windows platforms only
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -fsanitize=address -fno-omit-frame-pointer")
    set(CMAKE_EXE_LINKER_FLAGS_DEBUG "${CMAKE_EXE_LINKER_FLAGS_DEBUG} -fsanitize=address -fno-omit-frame-pointer")
  endif()

  # Release build flags
  set(CMAKE_CXX_FLAGS_RELEASE "-O3 -g -DNDEBUG")
endif()

# ---- External Libraries ----

# Eigen (header-only)
include_directories(external/eigen)

set(BUILD_SHARED_LIBS OFF CACHE BOOL "" FORCE)

# Triangle (build subdir with known output location)
set(TRIANGLE_SOURCE_DIR ${CMAKE_SOURCE_DIR}/external/Triangle)
set(TRIANGLE_BINARY_DIR ${CMAKE_BINARY_DIR}/triangle_build)
add_subdirectory(${TRIANGLE_SOURCE_DIR} ${TRIANGLE_BINARY_DIR})
# ---- After pulling in the libigl triangle subdir ----
if (CMAKE_CXX_COMPILER_ID MATCHES "GNU")
  # Make the pointer→unsigned long casts only warnings
  target_compile_options(triangle PRIVATE -fpermissive)
endif()


set(LIBGEOM_SOURCE_DIR ${CMAKE_SOURCE_DIR}/external/libgeom)
set(OPENMESH_SOURCE_DIR ${CMAKE_SOURCE_DIR}/external/OpenMesh/src)

# Avoid building tools with Qt and custom macros

set(OPENMESH_BUILD_UNIT_TESTS OFF CACHE BOOL "" FORCE)
set(OPENMESH_BUILD_DECIMATER OFF CACHE BOOL "" FORCE)
set(OPENMESH_BUILD_NORMAL_DEVIATION OFF CACHE BOOL "" FORCE)
set(OPENMESH_BUILD_PYTHON_BINDINGS OFF CACHE BOOL "" FORCE)
set(OPENMESH_BUILD_APPS OFF CACHE BOOL "" FORCE)
set(OPENMESH_BUILD_DOCS OFF CACHE BOOL "" FORCE)
set(OPENMESH_BUILD_UTILS OFF CACHE BOOL "" FORCE)

# OpenMesh and libgeom
add_subdirectory(external/OpenMesh)
add_subdirectory(external/libgeom)


# ---- Your Library ----
include_directories(include)

file(GLOB_RECURSE LIBCDGBS_SRC "src/*.cpp")

add_library(libcdgbs STATIC ${LIBCDGBS_SRC})

# Precompiled headers for libcdgbs
# Include heavy template headers once for faster builds and lower memory usage
target_precompile_headers(libcdgbs PRIVATE
    <Eigen/Core>
    <OpenMesh/Core/IO/MeshIO.hh>
)

# 🔧 Add include paths for Triangle, including generated headers
target_include_directories(libcdgbs
    PRIVATE
        ${OPENMESH_SOURCE_DIR}
        ${LIBGEOM_SOURCE_DIR}
        ${TRIANGLE_SOURCE_DIR}
)

if(WIN32)
target_compile_definitions(libcdgbs PRIVATE NOMINMAX)
endif()


# 🔗 Link libraries
target_link_libraries(libcdgbs
    triangle
    OpenMeshCore
    OpenMeshTools
    geom
)

# ---- Test Executable ----
add_executable(libcdgbs_test tests/main_test.cpp)

# Precompiled headers for test executable
target_precompile_headers(libcdgbs_test PRIVATE
    <Eigen/Core>
    <OpenMesh/Core/IO/MeshIO.hh>
)

# 🔧 Also add include path for generated headers here
target_include_directories(libcdgbs_test
    PRIVATE
        ${OPENMESH_SOURCE_DIR}
        ${LIBGEOM_SOURCE_DIR}
        ${TRIANGLE_BINARY_DIR}
)

if(WIN32)
target_compile_definitions(libcdgbs_test PRIVATE NOMINMAX)
endif()

target_link_libraries(libcdgbs_test libcdgbs)
