#!/bin/bash


# set OpenLB
BDM_SCRIPT_DIR=$(readlink -e $(dirname "${BASH_SOURCE[0]}"))
cd $BDM_SCRIPT_DIR
source $BDM_SCRIPT_DIR/../../util/main.sh

cd openlb

wget https://www.openlb.net/wp-content/uploads/2020/11/olb-1.4r0.tgz
tar -xzf olb-1.4r0.tgz
rm olb-1.4r0.tgz
cd olb-1.4r0

# fix for include in openlb/olb-1.4r0/src/particles/particleOperations/particleOperations3D.h (#include "../../functors/functors3D.h")
cp ../particleOperations3D.h src/particles/particleOperations/
# overwrite config file to enable MPI in olb
cp ../config.mk .

# make clean
make -j

# Quick fix for linking error inside docker container:
sudo mv /usr/lib/x86_64-linux-gnu/libz.a /usr/lib/x86_64-linux-gnu/libz.a.bak

# set BioDynaMo
# if 1, do not perform clean bdm build
# export BPE_NO_CLEANBUILD=1

source $BDM_SCRIPT_DIR/../../util/default-compile-script.sh "-Djemalloc=off" "-Djemalloc=off"

export OMP_NUM_THREADS=1
export BPE_REPEAT=1
unset OMP_PROC_BIND

# FIXME Temporary fix for ROOT_INCLUDE_PATH issue. It seems that due to some reason 
# static inititalization of bdm dictionaries occurs after cling inititalization.
export ROOT_INCLUDE_PATH=$PROJECT_ROOT_DIR/bdm-paper-examples/bdm/epidemio_spread/openlb/olb-1.4r0/src:/usr/include:$PROJECT_ROOT_DIR/biodynamo/build/third_party/root/include:$PROJECT_ROOT_DIR/bdm-paper-examples/bdm/epidemio_spread/build/omp:$PROJECT_ROOT_DIR/biodynamo/build/include:$PROJECT_ROOT_DIR/bdm-paper-examples/bdm/epidemio_spread/src:/usr/lib/x86_64-linux-gnu/openmpi/include/openmpi:/usr/lib/x86_64-linux-gnu/openmpi/include:$PROJECT_ROOT_DIR/biodynamo/build/third_party/paraview/include/paraview-5.8:$PROJECT_ROOT_DIR/biodynamo/build/third_party/root/include:$PROJECT_ROOT_DIR/biodynamo/build/include::

mpirun -np 71 ./openlb_sim
