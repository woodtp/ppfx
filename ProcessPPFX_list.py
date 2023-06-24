#!/usr/bin/env python2
##################################################
# ppfx grid job submitter
# 2019.07 -- L. Aliaga
##################################################
import os, optparse, random, shutil, tarfile, sys
import subprocess, string

CACHE_PNFS_AREA = "/pnfs/icarus/scratch/users/{USER}/grid_cache/".format(USER = os.getenv("USER"))
PWD = os.getenv("PWD")

##################################################
# Job Defaults
##################################################
#N_JOBS             = 1000
N_JOBS             = 1
#N_JOBS             = 5000
INPUT_OPTIONS      = "scripts/inputs_default.xml"
#OPTION = "QuarterWeight"
#INPUT_OPTIONS      = "scripts/inputs_"+OPTION+".xml"

#define detector positions (x, y, z), underscores needed, they will be removed in the job script

MINERVA            = "-56.28_-53.29_103231.9"
NOVA_ND            = 3


#IDET               = "326.0_8278.0_79883.0"
POINT1             = "871.8_7904.3_78519.0" # point1
POINT2             = "92.1_8011.7_80360.1" #point2
POINT1_1           = "832.81_7909.67_78611.061" #point1_1 (1 meter towards center from point 1)
POINT2_1           = "131.08_8006.33_80268.05" #point2_1 (1 meter towards center from point 2)
ICARUS_CENTER      = "481.95_7958.00_79439.55" #point in the middle between 1 and 2

OFF_AXIS_0_5       = "240.97_3979.00_79439.55" #off axis point twice closer to the beamline than ICARUS

ICARUS_22998       = "450.3730_8015.3901_79511.2945"   #SBN-doc-22998

ICARUS_geometrical_center = "450.37_7991.98_79512.66"

ICARUS_Xm358_49 = "120.19_7983.84_79373.28"
ICARUS_Xp358_49 = "780.56_8000.12_79652.05"
ICARUS_Ym158_41 = "450.37_7833.84_79521.90"
ICARUS_Ym79_20  = "450.37_7912.91_79517.28"
ICARUS_Yp79_21  = "450.37_8071.05_79508.04"
ICARUS_Yp158_41 = "450.37_8150.12_79503.42"
ICARUS_Zm894_95 = "798.94_7943.91_78689.78"
ICARUS_Zp894_95 = "101.81_8040.05_80335.54"

TWO_BY_TWO      = "0.0_0.0_103648.837"

IDET=TWO_BY_TWO
#IDET=ICARUS_geometrical_center
#IDET=ICARUS_Xm358_49
#IDET=ICARUS_Xp358_49
#IDET=ICARUS_Ym158_41
#IDET=ICARUS_Ym79_20
#IDET=ICARUS_Yp79_21
#IDET=ICARUS_Yp158_41
#IDET=ICARUS_Zm894_95
#IDET=ICARUS_Zp894_95

#name of the dataset
DATA_TAG   = "Nilay_NuMI_1MW_RHC"

TARFILE_NAME       = "local_install.tar.gz"

#OUTDIR             = "/pnfs/icarus/persistent/users/{USER}/".format(USER = os.getenv("USER"))+"ppfx_2023-04_pCQEL10a/"+IDET+"_"+RUN_NUMBER+"_"+BEAMCONFIG+"/"
OUTDIR             = "/pnfs/icarus/scratch/users/{USER}/".format(USER = os.getenv("USER"))+"ppfx_2023-06_" + DATA_TAG +"/"

##################################################

def main():
  options = get_options()
  
  cache_folder = CACHE_PNFS_AREA + str(random.randint(10000,99999)) + "/"
#  os.makedirs(cache_folder, exist_ok=False)   #python3
  os.makedirs(cache_folder)
  
#  os.makedirs(options.outdir, exist_ok=False) #python3
  os.makedirs(options.outdir)

  print("\nTarring up local area...")
  make_tarfile(TARFILE_NAME, ".")

  shutil.move(TARFILE_NAME,    cache_folder) 
  shutil.copy("ppfx_job_list.sh", cache_folder)
  
  print("\nTarball of local area:", cache_folder + TARFILE_NAME)

  logfile = options.outdir + "/ppfx_" + DATA_TAG + "_\$PROCESS.log"
  
  print("\nOutput logfile(s):",logfile)

  submit_command = ("jobsub_submit {GRID} {MEMORY} {DISK} -N {NJOBS} -d PPFX {OUTDIR} "
      "-G icarus "
      "-e DATA_TAG={DATA_TAG} " 
      "-e INPUT_OPTIONS={INPUT_OPTIONS} " 
      "-e IDET={IDET} " 
      "-f {TARFILE} "
      "-L {LOGFILE} "
      "--mail-always "
      "file://{CACHE}/ppfx_job_list.sh".format(
      GRID       = ("--OS=SL7 "
                    "--resource-provides=usage_model=DEDICATED,OPPORTUNISTIC "
#                    "--expected-lifetime=900 "
                    "--expected-lifetime 3600 "
                    "--role=Analysis "),
      MEMORY     = "--memory 500MB ",
      DISK       = "--disk 2GB ",
      NJOBS      = options.n_jobs,
      OUTDIR     = options.outdir,
      DATA_TAG   = DATA_TAG,
      INPUT_OPTIONS = options.input_options,
      IDET       = IDET,
      TARFILE    = cache_folder + TARFILE_NAME,
      LOGFILE    = logfile,
      CACHE      = cache_folder)
  )

  #Ship it
  print("\nSubmitting to grid:\n"+submit_command+"\n")
  status = subprocess.call(submit_command, shell=True)
  return(status)

def get_options():
  parser       = optparse.OptionParser(usage="usage: %prog [options]")
  grid_group   = optparse.OptionGroup(parser, "Grid Options")
  
  grid_group.add_option("--outdir",
                        default = OUTDIR,
                        help    = "Output flux histograms location. Default = %default.")
  
  grid_group.add_option("--n_jobs",
                        default = N_JOBS,
                        help = "Number of g4numi jobs. Default = %default.")
  
  beam_group   = optparse.OptionGroup(parser, "Beam Options")
  
  run_group   = optparse.OptionGroup(parser, "Run Options")
  
  run_group.add_option("--input_options",
                       default = INPUT_OPTIONS,
                       help    = "PPFX input: number of universes, MIPP on/off, etc. Default = %default.")
  
  
  parser.add_option_group(grid_group)
  parser.add_option_group(beam_group)
  parser.add_option_group(run_group)

  options, remainder = parser.parse_args()

  return options

def make_tarfile(output_filename, source_dir):
  tar = tarfile.open(output_filename, "w:gz")

  #We can reduce the tarball size by over 50% by omiting unuseful files, like hidden (.git), html...
  def omit_file(string):
    prefixes = [".git", "html", "src", "latex"]
    for p in prefixes:
      if string[0:len(p)] == p:
        return True
    return False

  for i in os.listdir(source_dir):
    if not omit_file(i):
      tar.add(i)
    else:
      print("Omitting file", i)
  tar.close()

if __name__ == "__main__":
  sys.exit(main())
