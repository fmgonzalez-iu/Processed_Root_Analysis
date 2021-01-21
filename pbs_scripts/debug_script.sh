#!/bin/bash
#PBS -l nodes=2:ppn=24
#PBS -l walltime=01:00:00
#PBS -N AnalyzerForTau_debug
#PBS -q debug
#PBS -V

#Call "date" at beginning and end just to get timestamps
date

cd /N/u/frangonz/Carbonate/New_MCA_Analysis/data-LowT_bkg

#GNU parallel expects a list of items to loop over in a parallel way.
#Here "seq" just provides a list of integers 1 to 24. We then run another script
#Which takes a given one of those integers N and runs an instance of g4simple-naive
#while writing the output to <some file>_N.root. The "{}" is a special instruction
#to GNU parallel to tell it to supply the given integer.
seq 0 47 | parallel ../AnalyzerForeach /N/u/frangonz/Carbonate/runLists/backgrounds/beamOn.txt_ 50 500 8 3 1 1 {}
#seq 0 47 | parallel ../AnalyzerForeach ../runLists/Carbonate/RunListSplit.txt_ 50 500 8 3 1 1 {}

date
