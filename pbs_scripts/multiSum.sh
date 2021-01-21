#!/bin/bash
#PBS -l nodes=1:ppn=1
#PBS -l walltime=01:00:00
#PBS -N AnalyzerForTau_debug
#PBS -q debug
#PBS -V

#Call "date" at beginning and end just to get timestamps
date

cd /N/u/frangonz/Carbonate/New_MCA_Analysis/data-DeepCleaning
export SAVE_LOC=/N/u/frangonz/Carbonate/New_MCA_Analysis/data-DeepCleaning

#GNU parallel expects a list of items to loop over in a parallel way.
#Here "seq" just provides a list of integers 1 to 24. We then run another script
#Which takes a given one of those integers N and runs an instance of g4simple-naive
#while writing the output to <some file>_N.root. The "{}" is a special instruction
#to GNU parallel to tell it to supply the given integer.
../AnalyzerForeach 05147,05465,05471 50 500 8 4 2 1
mv $SAVE_LOC/summedTraces/summedDipShort.root $SAVE_LOC/summedTraces/40k.root
../AnalyzerForeach 05148,05466 50 500 8 4 2 1
mv $SAVE_LOC/summedTraces/summedDipShort.root $SAVE_LOC/summedTraces/60k.root
../AnalyzerForeach 05149,05468 50 500 8 4 2 1
mv $SAVE_LOC/summedTraces/summedDipShort.root $SAVE_LOC/summedTraces/80k.root 
../AnalyzerForeach 05150,05469 50 500 8 4 2 1
mv $SAVE_LOC/summedTraces/summedDipShort.root $SAVE_LOC/summedTraces/110k.root
../AnalyzerForeach 05151,05472 50 500 8 4 2 1
mv $SAVE_LOC/summedTraces/summedDipShort.root $SAVE_LOC/summedTraces/140k.root
../AnalyzerForeach 05152,05473 50 500 8 4 2 1
mv $SAVE_LOC/summedTraces/summedDipShort.root $SAVE_LOC/summedTraces/250k.root
../AnalyzerForeach 05153,05474 50 500 8 4 2 1
mv $SAVE_LOC/summedTraces/summedDipShort.root $SAVE_LOC/summedTraces/180k.root 
date
