#!/bin/bash
#PBS -l nodes=1:ppn=24
#PBS -l walltime=01:00:00
#PBS -N PMTSummer
#PBS -q debug
#PBS -V

#Call "date" at beginning and end just to get timestamps
date

cd /N/u/frangonz/Carbonate/New_MCA_Analysis/Outputs

RUN_LISTS=/N/u/frangonz/BigRed3/runLists/summedPmt/4231/
##SAVE_LOC=/N/dc2/scratch/frangonz/TDepB_Hi/4231/
SAVE_LOC=/N/slate/frangonz/summedHists/rawHits/LowT_Comb/PMT2/4231/
seq 0 23 | parallel ../AnalyzerForeach ${RUN_LISTS}/RunList.txt_ 50 500 8 7 6 5 {}
seq 0 23 | parallel ../AnalyzerForeach ${RUN_LISTS}/RunList.txt_ 50 500 8 7 6 6 {}
seq 0 23 | parallel ../AnalyzerForeach ${RUN_LISTS}/RunList.txt_ 50 500 8 7 6 7 {}
seq 0 23 | parallel ../AnalyzerForeach ${RUN_LISTS}/RunList.txt_ 50 500 8 7 6 8 {}

RUN_LISTS=/N/u/frangonz/BigRed3/runLists/summedPmt/7332/
SAVE_LOC=/N/slate/frangonz/summedHists/rawHits/LowT_Comb/PMT2/7332/
seq 0 23 | parallel ../AnalyzerForeach ${RUN_LISTS}/RunList.txt_ 50 500 8 7 6 5 {}
seq 0 23 | parallel ../AnalyzerForeach ${RUN_LISTS}/RunList.txt_ 50 500 8 7 6 6 {}
seq 0 23 | parallel ../AnalyzerForeach ${RUN_LISTS}/RunList.txt_ 50 500 8 7 6 7 {}
seq 0 23 | parallel ../AnalyzerForeach ${RUN_LISTS}/RunList.txt_ 50 500 8 7 6 8 {}

RUN_LISTS=/N/u/frangonz/BigRed3/runLists/summedPmt/9600/
SAVE_LOC=/N/slate/frangonz/summedHists/rawHits/LowT_Comb/PMT2/9600/
seq 0 23 | parallel ../AnalyzerForeach ${RUN_LISTS}/RunList.txt_ 50 500 8 7 6 5 {}
seq 0 23 | parallel ../AnalyzerForeach ${RUN_LISTS}/RunList.txt_ 50 500 8 7 6 6 {}
seq 0 23 | parallel ../AnalyzerForeach ${RUN_LISTS}/RunList.txt_ 50 500 8 7 6 7 {}
seq 0 23 | parallel ../AnalyzerForeach ${RUN_LISTS}/RunList.txt_ 50 500 8 7 6 8 {}

RUN_LISTS=/N/u/frangonz/BigRed3/runLists/summedPmt/13219/
SAVE_LOC=/N/slate/frangonz/summedHists/rawHits/LowT_Comb/PMT2/13219/
seq 0 23 | parallel ../AnalyzerForeach ${RUN_LISTS}/RunList.txt_ 50 500 8 7 6 5 {}
seq 0 23 | parallel ../AnalyzerForeach ${RUN_LISTS}/RunList.txt_ 50 500 8 7 6 6 {}
seq 0 23 | parallel ../AnalyzerForeach ${RUN_LISTS}/RunList.txt_ 50 500 8 7 6 7 {}
seq 0 23 | parallel ../AnalyzerForeach ${RUN_LISTS}/RunList.txt_ 50 500 8 7 6 8 {}

date
