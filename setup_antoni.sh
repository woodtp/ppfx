#!/bin/bash


# the default mode is NuMI
export MODE="NUMI"

# set package to be used
source /cvmfs/icarus.opensciencegrid.org/products/icarus/setup_icarus.sh
source /cvmfs/sbnd.opensciencegrid.org/products/sbnd/setup
##ROOT: 
#setup -qe14:prof  -f Linux64bit+3.10-2.17 root v6_08_06g

#BOOST:
setup -q e9:prof -f Linux64bit+3.10-2.17 boost v1_57_0a
export BOOSTROOT=${BOOST_DIR}/source/boost_1_57_0
#DK2NU:
setup -q e9:prof:r5 -f Linux64bit+3.10-2.17 dk2nu v01_03_00c

export DK2NU_INC=${DK2NU}/include/dk2nu/tree
export DK2NU_LIB=${DK2NU}/lib

#GRID
#Note: with jobsub_lite this is no longer needed
#source /cvmfs/fermilab.opensciencegrid.org/products/common/etc/setups.sh
#setup jobsub_client
setup ifdhc #v2_2_3
export IFDH_GRIDFTP_EXTRA="-st 10" #set ifdh cp stall timeout to 10 sec
export IFDH_CP_MAXRETRIES=2

# gets the full path to the location of setup.sh
export PPFX_DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
echo "setting PPFX_DIR=${PPFX_DIR}"

export LD_LIBRARY_PATH=$PPFX_DIR/lib:$LD_LIBRARY_PATH
echo "LD_LIBRARY_PATH=$LD_LIBRARY_PATH"
