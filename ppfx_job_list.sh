#!/bin/bash

echo
echo $(date) "======== cd to CONDOR_DIR_INPUT ========"
cd $CONDOR_DIR_INPUT
pwd

echo
echo $(date) "======== env ========"
env | sort

echo
echo $(date) "======== ls ========"
ls

echo
echo $(date) "======== UNTARRING... ========"
tar xvfz local_install.tar.gz -C ./ > /dev/null

echo
echo $(date) "======== Done untarring. ls ========"
ls

echo
echo $(date) "======== SETUP ROOT, BOOST and DK2NU ========"
echo "source setup_antoni.sh"
source setup_antoni.sh

echo $(date) "======== ups active ========"
ups active

echo
echo $(date) "======== UPDATE g4numi run number to select input ========"
echo PROCESS=$PROCESS
eval process=${PROCESS}
echo process=$process
#INPUT_FILE=$(sed "$((process+1))q;d" filelist_${DATA_TAG})
INPUT_FILE="/pnfs/numix/persistent/users/nbostan/G4NuMI_RHC_new_target/g4numiv6_minervame_me000z-200i_\${PROCESS}_0001.root"
name_file=$(basename $INPUT_FILE)
#eval name_file=$(basename $INPUT_FILE)

echo name_file=$name_file
echo INPUT_FILE=$INPUT_FILE

ifdh ls $INPUT_FILE
ifdh cp "${INPUT_FILE}" "$CONDOR_DIR_INPUT/"
ls -l $INPUT_FILE
if [ $(ifdh ls $INPUT_FILE | wc -l) -eq 0 ]
then
  echo ERROR: File $input_file is empty or does not exist!!! Aborting.
  exit
fi

echo running: ifdh cp "${INPUT_FILE}" "$CONDOR_DIR_INPUT/"

#more than 1 attempt to copy might be needed if the file is not staged
n_attempts=5
while [ $n_attempts -gt 0 ]
do
  echo $n_attempts attempts to copy left
  ifdh cp "${INPUT_FILE}" "$CONDOR_DIR_INPUT/"
  if [ -s $name_file ]
  then
    break
  fi
  n_attempts=$((n_attempts-1))
  sleep 300
done

ls -l

#if [ ! -s $INPUT_FILE ]
#then
#  echo "Can't find input file. Aborting."
#  exit
#fi

OUTPUT_FILE="ppfx_${DATA_TAG}_${name_file%.*}.root"


echo "DATA_TAG=$DATA_TAG"
echo "PROCESS=$PROCESS"
echo "OUTFILE=$OUTPUT_FILE"
eval output_file=$OUTPUT_FILE

echo
echo $(date) "======== EXECUTING ppfx ========"
#deal with \${PROCESS}
#replace underscores in IDET with spaces
idet=$(echo $IDET | sed s/_/\ /g)
cmd="bin/doReweight_dk2nu_original ${name_file} ${output_file} ${INPUT_OPTIONS} ${idet}"
#cmd="bin/doReweight_dk2nu_original ${INPUT_FILE} ${output_file} ${INPUT_OPTIONS} ${idet}"
echo $cmd
eval $cmd

ls -l

echo
echo $(date) "Moving output to CONDOR_DIR_PPFX: "
echo "=> CONDOR_DIR_PPFX: $CONDOR_DIR_PPFX"
rm g4numi*root
mv -v ppfx*.root $CONDOR_DIR_PPFX
