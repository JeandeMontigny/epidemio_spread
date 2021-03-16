# -----------------------------------------------------------------------------
#
# Copyright (C) The BioDynaMo Project.
# All Rights Reserved.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
#
# See the LICENSE file distributed with this work for details.
# See the NOTICE file distributed with this work for additional information
# regarding copyright ownership.
#
# -----------------------------------------------------------------------------
cmake_minimum_required(VERSION 3.2.0)

project(openlb_sim)

include_directories(${CMAKE_SOURCE_DIR}/openlb/olb-1.4r0/src/)
link_directories(${CMAKE_SOURCE_DIR}/openlb/olb-1.4r0/build/precompiled/lib/)
# Define openlb parallelization strategy. Has to be consistent with the openlb build
# If missing, functions defined in openlb headers that are included in BDM will be compiled
# without parallelization, causing a crash.
add_definitions(-DPARALLEL_MODE_MPI)

find_package(BioDynaMo REQUIRED)
find_package(MPI REQUIRED)

include(${BDM_USE_FILE})
include_directories("src")
include_directories(SYSTEM ${MPI_INCLUDE_PATH})

file(GLOB_RECURSE HEADERS src/*.h)
file(GLOB_RECURSE SOURCES src/*.cc)

bdm_add_executable(openlb_sim
                   HEADERS ${HEADERS}
                   SOURCES ${SOURCES}
                   LIBRARIES ${BDM_REQUIRED_LIBRARIES} ${MPI_C_LIBRARIES} ${MPI_CXX_LIBRARIES} libolb.a libz.a)