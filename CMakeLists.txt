cmake_minimum_required(VERSION 3.21)
project(REEF3D LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# ── Dependency versions ───────────────────────────────────────────────────────
set(REEF3D_HYPRE_VERSION  "3.2.0")
set(REEF3D_EIGEN_VERSION  "5.0.0")

# Default build type mirrors Makefile's .DEFAULT_GOAL = release
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()
set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS Release Debug Dev All)

# Non-standard build type flag sets (Dev / All)
set(CMAKE_CXX_FLAGS_DEV
  "-O3 -Wall -pedantic -Wpedantic -Wextra -Wshadow -Wcast-align -Wconversion -Wsign-conversion -Wnull-dereference -Wdouble-promotion -Wformat=2"
  CACHE STRING "Flags for Dev build" FORCE)
set(CMAKE_CXX_FLAGS_ALL "-O3 -w -g" CACHE STRING "Flags for All build" FORCE)
mark_as_advanced(CMAKE_CXX_FLAGS_DEV CMAKE_CXX_FLAGS_ALL)

# ── Git version ──────────────────────────────────────────────────────────────

find_package(Git QUIET)
if(Git_FOUND AND EXISTS "${PROJECT_SOURCE_DIR}/.git")
  execute_process(
    COMMAND ${GIT_EXECUTABLE} rev-parse --short=7 HEAD
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
    OUTPUT_VARIABLE GIT_COMMIT OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_QUIET)
  execute_process(
    COMMAND ${GIT_EXECUTABLE} rev-parse --abbrev-ref HEAD
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
    OUTPUT_VARIABLE GIT_BRANCH OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_QUIET)
  execute_process(
    COMMAND ${GIT_EXECUTABLE} diff --quiet --ignore-submodules HEAD --
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
    RESULT_VARIABLE GIT_IS_DIRTY ERROR_QUIET)
  if(NOT GIT_IS_DIRTY EQUAL 0)
    set(GIT_VERSION "${GIT_COMMIT}-dirty")
  else()
    set(GIT_VERSION "${GIT_COMMIT}")
  endif()
else()
  set(GIT_COMMIT "unknown")
  set(GIT_BRANCH "unknown")
  set(GIT_VERSION "unknown")
endif()

# ── MPI ──────────────────────────────────────────────────────────────────────
# MPI is a system dependency; no FetchContent fallback.

find_package(MPI REQUIRED CXX)

# ── HYPRE ────────────────────────────────────────────────────────────────────

set(HYPRE_DIR "/usr/local/hypre" CACHE PATH "HYPRE installation prefix")

find_library(_HYPRE_LIB NAMES HYPRE HINTS "${HYPRE_DIR}/lib" NO_DEFAULT_PATH)
if(_HYPRE_LIB)
  message(STATUS "Found HYPRE: ${_HYPRE_LIB}")
  add_library(HYPRE::HYPRE UNKNOWN IMPORTED)
  set_target_properties(HYPRE::HYPRE PROPERTIES
    IMPORTED_LOCATION             "${_HYPRE_LIB}"
    INTERFACE_INCLUDE_DIRECTORIES "${HYPRE_DIR}/include")
else()
  message(STATUS "HYPRE not found at ${HYPRE_DIR} — downloading v${REEF3D_HYPRE_VERSION} via FetchContent")
  include(FetchContent)
  FetchContent_Declare(hypre_src
    GIT_REPOSITORY https://github.com/hypre-space/hypre.git
    GIT_TAG        v${REEF3D_HYPRE_VERSION}
    GIT_SHALLOW    ON
    SOURCE_SUBDIR  src)
  set(HYPRE_BUILD_SHARED_LIBS  OFF CACHE BOOL "" FORCE)
  set(HYPRE_ENABLE_SHARED      OFF CACHE BOOL "" FORCE)
  FetchContent_MakeAvailable(hypre_src)
  # HYPRE's CMake exposes target "HYPRE"
  add_library(HYPRE::HYPRE ALIAS HYPRE)
endif()

# ── Eigen ─────────────────────────────────────────────────────────────────────

set(EIGEN_LOCAL_DIR "${PROJECT_SOURCE_DIR}/ThirdParty/eigen-${REEF3D_EIGEN_VERSION}")

if(EXISTS "${EIGEN_LOCAL_DIR}/Eigen/Dense")
  message(STATUS "Found Eigen: ${EIGEN_LOCAL_DIR}")
  add_library(Eigen3::Eigen INTERFACE IMPORTED)
  set_target_properties(Eigen3::Eigen PROPERTIES
    INTERFACE_INCLUDE_DIRECTORIES "${EIGEN_LOCAL_DIR}")
else()
  find_package(Eigen3 ${REEF3D_EIGEN_VERSION} QUIET NO_MODULE)
  if(NOT TARGET Eigen3::Eigen)
    message(STATUS "Eigen not found locally — downloading ${REEF3D_EIGEN_VERSION} via FetchContent")
    include(FetchContent)
    FetchContent_Declare(Eigen3
      GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
      GIT_TAG        ${REEF3D_EIGEN_VERSION}
      GIT_SHALLOW    ON)
    set(EIGEN_BUILD_DOC   OFF CACHE BOOL "" FORCE)
    set(BUILD_TESTING     OFF CACHE BOOL "" FORCE)
    FetchContent_MakeAvailable(Eigen3)
  endif()
endif()
# Enforce MPL2-only regardless of how Eigen was found
target_compile_definitions(Eigen3::Eigen INTERFACE EIGEN_MPL2_ONLY)

# ── Executable ───────────────────────────────────────────────────────────────

file(GLOB SOURCES CONFIGURE_DEPENDS src/*.cpp)
add_executable(REEF3D ${SOURCES})

set_target_properties(REEF3D PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/bin")

# ── Compile definitions ───────────────────────────────────────────────────────

target_compile_definitions(REEF3D PRIVATE
  VERSION="${GIT_VERSION}"
  BRANCH="${GIT_BRANCH}"
  $<$<CONFIG:Release>:NDEBUG;EIGEN_NO_DEBUG;BUILD="release">
  $<$<CONFIG:Debug>:BUILD="debug">
  $<$<CONFIG:Dev>:BUILD="dev">
  $<$<CONFIG:All>:BUILD="all">)

# ── Compile flags ─────────────────────────────────────────────────────────────
# CMAKE_CXX_FLAGS_<CONFIG> already covers Dev/All (set above).
# Generator expressions handle Release and Debug overrides.

target_compile_options(REEF3D PRIVATE
  $<$<CONFIG:Release>:-O3;-march=native;-flto;-w>
  $<$<CONFIG:Debug>:-O0;-g;-g3;-Wall>)

target_link_options(REEF3D PRIVATE
  $<$<CONFIG:Release>:-flto=auto>)

# ── Link libraries ────────────────────────────────────────────────────────────

target_link_libraries(REEF3D PRIVATE
  MPI::MPI_CXX
  HYPRE::HYPRE
  Eigen3::Eigen)
