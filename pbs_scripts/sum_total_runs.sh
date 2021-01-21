#!/bin/bash
#PBS -l nodes=1:ppn=24
#PBS -l walltime=01:00:00
#PBS -N RunSummer
#PBS -q debug
#PBS -V

#Call "date" at beginning and end just to get timestamps
date

cd /N/u/frangonz/Carbonate/New_MCA_Analysis/Outputs

# Folders
IN_FOLDER=/N/u/frangonz/BigRed3/runLists/summedPmt
IN_FOLDER_B=/N/u/frangonz/Carbonate/runLists/backgrounds
OUT_FOLDER=/N/slate/frangonz/summedHists/removeEvts

# This is hardcoded for our gaps
RUN_LISTS=${IN_FOLDER}/4231/
RUN_LISTS_B=${IN_FOLDER_B}/4231/
SAVE_LOC=${OUT_FOLDER}/4231/
mkdir ${SAVE_LOC}
seq 0 23 | parallel ../AnalyzerForeach ${RUN_LISTS}/RunList.txt_ 50 500 8 3 2 4 {}
seq 0 23 | parallel ../AnalyzerForeach ${RUN_LISTS_B}/beamOff.txt_ 50 500 8 3 2 4 {}

RUN_LISTS=${IN_FOLDER}/7332/
RUN_LISTS_B=${IN_FOLDER_B}/7332/
SAVE_LOC=${OUT_FOLDER}/7332/
mkdir ${SAVE_LOC}
seq 0 23 | parallel ../AnalyzerForeach ${RUN_LISTS}/RunList.txt_ 50 500 8 3 2 4 {}
seq 0 23 | parallel ../AnalyzerForeach ${RUN_LISTS_B}/beamOff.txt_ 50 500 8 3 2 4 {}

RUN_LISTS=${IN_FOLDER}/9600/
RUN_LISTS_B=${IN_FOLDER_B}/9600/
SAVE_LOC=${OUT_FOLDER}/9600/
mkdir ${SAVE_LOC}
seq 0 23 | parallel ../AnalyzerForeach ${RUN_LISTS}/RunList.txt_ 50 500 8 3 2 4 {}
seq 0 23 | parallel ../AnalyzerForeach ${RUN_LISTS_B}/beamOff.txt_ 50 500 8 3 2 4 {}

RUN_LISTS=${IN_FOLDER}/13219/
RUN_LISTS_B=${IN_FOLDER_B}/13219/
SAVE_LOC=${OUT_FOLDER}/13219/
mkdir ${SAVE_LOC}
seq 0 23 | parallel ../AnalyzerForeach ${RUN_LISTS}/RunList.txt_ 50 500 8 3 2 4 {}
seq 0 23 | parallel ../AnalyzerForeach ${RUN_LISTS_B}/beamOff.txt_ 50 500 8 3 2 4 {}

date
