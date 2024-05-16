#!/usr/bin/env python3
##################################################
# ppfx grid job submitter
# 2019.07 -- L. Aliaga
# Python3 rework -- A. Wood (2024)
##################################################

import argparse
import os
import random
import shutil
import tarfile
import subprocess
from pathlib import Path
from datetime import datetime

SCRATCH_AREA = Path(f"/pnfs/icarus/scratch/users/{os.getenv('USER')}")
CACHE_PNFS_AREA = SCRATCH_AREA / "grid_cache/"
PWD = Path(os.getcwd())

##################################################
# Job Defaults
##################################################
N_JOBS = 993  # set this equal to the number of input files
# N_JOBS = 495  # set this equal to the number of input files
INPUT_OPTIONS = "scripts/inputs_default.xml"
# OPTION = "QuarterWeight"
# INPUT_OPTIONS      = "scripts/inputs_"+OPTION+".xml"

# define detector positions (x, y, z), underscores needed, they will be removed in the job script

MINERVA = "-56.28_-53.29_103231.9"
NOVA_ND = 3

# IDET               = "326.0_8278.0_79883.0"
POINT1 = "871.8_7904.3_78519.0"  # point1
POINT2 = "92.1_8011.7_80360.1"  # point2
POINT1_1 = "832.81_7909.67_78611.061"  # point1_1 (1 meter towards center from point 1)
POINT2_1 = "131.08_8006.33_80268.05"  # point2_1 (1 meter towards center from point 2)
ICARUS_CENTER = "481.95_7958.00_79439.55"  # point in the middle between 1 and 2

OFF_AXIS_0_5 = (
    "240.97_3979.00_79439.55"  # off axis point twice closer to the beamline than ICARUS
)

ICARUS_22998 = "450.3730_8015.3901_79511.2945"  # SBN-doc-22998

ICARUS_geometrical_center = "450.37_7991.98_79512.66"

ICARUS_Xm358_49 = "120.19_7983.84_79373.28"
ICARUS_Xp358_49 = "780.56_8000.12_79652.05"
ICARUS_Ym158_41 = "450.37_7833.84_79521.90"
ICARUS_Ym79_20 = "450.37_7912.91_79517.28"
ICARUS_Yp79_21 = "450.37_8071.05_79508.04"
ICARUS_Yp158_41 = "450.37_8150.12_79503.42"
ICARUS_Zm894_95 = "798.94_7943.91_78689.78"
ICARUS_Zp894_95 = "101.81_8040.05_80335.54"

TWO_BY_TWO = "0.0_0.0_103648.837"

# IDET=TWO_BY_TWO
IDET = ICARUS_geometrical_center
# IDET=ICARUS_Xm358_49
# IDET=ICARUS_Xp358_49
# IDET=ICARUS_Ym158_41
# IDET=ICARUS_Ym79_20
# IDET=ICARUS_Yp79_21
# IDET=ICARUS_Yp158_41
# IDET=ICARUS_Zm894_95
# IDET=ICARUS_Zp894_95

# name of the dataset
# DATA_TAG = "Nilay_NuMI_1MW_RHC"
DATA_TAG = "userdataset_awood_NuMI_g4_10_4_all"
# DATA_TAG = f"userdataset_awood_numi_g4_10_4_RHC"

TARFILE_NAME = "local_install.tar.gz"

today = datetime.today().strftime("%Y-%m-%d")
OUTDIR = SCRATCH_AREA / f"{today}_ppfx_{DATA_TAG}"

##################################################


def main():
    options = get_options()

    cache_folder = CACHE_PNFS_AREA / str(random.randint(10000, 99999))

    cache_folder.mkdir(parents=True, exist_ok=False)
    options.outdir.mkdir(parents=True,exist_ok=False)

    # os.makedirs(cache_folder, exist_ok=False)

    # os.makedirs(options.outdir, exist_ok=False)

    print("\nTarring up local area...")
    make_tarfile(TARFILE_NAME, ".")

    shutil.move(TARFILE_NAME, cache_folder)
    shutil.copy("ppfx_job_list.sh", cache_folder)

    print("\nTarball of local area:", cache_folder / TARFILE_NAME)

    logfile = options.outdir / f"ppfx_{DATA_TAG}_$PROCESS.log"

    print("\nOutput logfile(s):", logfile)

    submit_command = (
        "jobsub_submit {GRID} {MEMORY} {DISK} -N {NJOBS} -d PPFX {OUTDIR} "
        "-G icarus "
        "-e DATA_TAG={DATA_TAG} "
        "-e INPUT_OPTIONS={INPUT_OPTIONS} "
        "-e IDET={IDET} "
        "-f {TARFILE} "
        "-L {LOGFILE} "
        "--singularity-image=/cvmfs/singularity.opensciencegrid.org/fermilab/fnal-wn-sl7:latest "
        "--mail-always "
        "file://{CACHE}/ppfx_job_list.sh".format(
            GRID=(
                # "--OS=SL7 "
                "--resource-provides=usage_model=DEDICATED,OPPORTUNISTIC "
                "--expected-lifetime 7200 "
                "--role=Analysis "
            ),
            MEMORY="--memory 4000MB ",
            DISK="--disk 2GB ",
            NJOBS=options.n_jobs,
            OUTDIR=options.outdir,
            DATA_TAG=DATA_TAG,
            INPUT_OPTIONS=options.input_options,
            IDET=IDET,
            TARFILE=cache_folder / TARFILE_NAME,
            LOGFILE=logfile,
            CACHE=cache_folder,
        )
    )

    # Ship it
    print("\nSubmitting to grid:\n" + submit_command + "\n")
    status = subprocess.call(submit_command, shell=True)
    return status


def get_options():
    parser = argparse.ArgumentParser(
        prog="ProcessPPFX_list.py", description="PPFX grid job submitter"
    )
    # parser = optparse.OptionParser(usage="usage: %prog [options]")
    # grid_group = optparse.OptionGroup(parser, "Grid Options")

    parser.add_argument(
        "--outdir",
        default=OUTDIR,
        type=Path,
        help="Output flux histograms location. (default: %(default)s)",
    )

    parser.add_argument(
        "--n_jobs", default=N_JOBS, help="Number of g4numi jobs. (default: %(default)s)"
    )

    # beam_group = optparse.OptionGroup(parser, "Beam Options")
    #
    # run_group = optparse.OptionGroup(parser, "Run Options")

    parser.add_argument(
        "--input_options",
        default=INPUT_OPTIONS,
        help="PPFX input: number of universes, MIPP on/off, etc. (default: %(default)s)",
    )

    # parser.add_option_group(grid_group)
    # parser.add_option_group(beam_group)
    # parser.add_option_group(run_group)

    options = parser.parse_args()

    return options


def make_tarfile(output_filename, source_dir):
    tar = tarfile.open(output_filename, "w:gz")

    # We can reduce the tarball size by over 50% by omiting unuseful files, like hidden (.git), html...
    def omit_file(string):
        prefixes = [".git", "html", "src", "latex"]
        for p in prefixes:
            if string[0 : len(p)] == p:
                return True
        return False

    for i in os.listdir(source_dir):
        if not omit_file(i):
            tar.add(i)
        else:
            print("Omitting file", i)
    tar.close()


if __name__ == "__main__":
    main()
