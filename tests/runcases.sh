#!/bin/bash
TOPDIR=${PWD}

#number of processors to use
NPROCS=${1:-16}

#need mpirun command (on eagle it is srun and 
#peregrine it is mpirun)
MPI_RUN_COMMAND=${2:-srun}

MAKEFILE_OPTIONS=${3:-}

declare -a allcases=('HCS' 'hopper_autogen' 'Settle_box'  'Settle_cyl'  'Silo'  'Silo_autogen' 'incline' 'wedge_hopper' 'moving_incline')

cd ../build
#make realclean
make -j ${MAKEFILE_OPTIONS}
cd ${TOPDIR}


#run cases
for case in "${allcases[@]}";
do
        echo ${case}
	cd ${case}
        . ../clean.sh
        $2 -n $1 ../../build/*.ex inputs >& out$case
        cd ${TOPDIR}
done
