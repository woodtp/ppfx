GCC_VERSION=12.2.0
ROOT_VERSION=6.28.06
FIFEUTILS_VERSION=3.7.4

source /cvmfs/larsoft.opensciencegrid.org/spack-packages/setup-env.sh

spack load gcc@$GCC_VERSION root@$ROOT_VERSION fife-utils@$FIFEUTILS_VERSION

export MODE="NUMI"
echo "setting MODE=${MODE}"

export PPFX_DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
echo "setting PPFX_DIR=${PPFX_DIR}"
