#!/bin/bash
#PBS -l nodes=1:ppn=1
#PBS -l walltime=00:03:00
#PBS -N PMTSummer
#PBS -q debug
#PBS -V

#Call "date" at beginning and end just to get timestamps
date

cd /N/u/frangonz/Carbonate/New_MCA_Analysis/Outputs

# Folders
IN_FOLDER=/N/u/frangonz/BigRed3/runLists/summedPmt
IN_FOLDER_B=/N/u/frangonz/Carbonate/runLists/backgrounds
OUT_FOLDER=/N/u/frangonz/Carbonate/New_MCA_Analysis/data-test/

# This is hardcoded for our gaps
SAVE_LOC=${OUT_FOLDER}
../AnalyzerForeach 06178 50 500 8 3 2 4
../AnalyzerForeach 08284 50 500 8 3 2 4


date
