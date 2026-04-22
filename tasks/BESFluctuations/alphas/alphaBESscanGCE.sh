#!/bin/sh

# Number of events from first argument (default 10000)
nevents=${1:-10000}

# input folder prefix
inputFolderPrefix=${2:-"../../../tasks/BESFluctuations/input/"}

energies=(7.7 14.5 19.6 27 39 62.4 200)

outdir="alphasGCE"
mkdir $outdir

for energy in ${energies[@]}; do
  ./BES-alphas-KpLaQ $inputFolderPrefix"input.AuAu."$energy".C0-5" --nevents=$nevents --Bcanonical=0 --ecm=$energy --output_file=$outdir"/alphas.AuAu."$energy".C0-5.dat"  --decays=3
done
