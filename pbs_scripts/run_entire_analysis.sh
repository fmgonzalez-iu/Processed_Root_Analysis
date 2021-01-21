# Generate Directories
# First load some global variables for analysis
export PBS_DIR=/N/u/frangonz/Carbonate/Processed_Root_Analysis/pbs_scripts # Location of scripts
export ANALYZER=/N/u/frangonz/Carbonate/Processed_Root_Analysis/rootAnalyzer # Location of analyzer
export RUNLISTS=/N/u/frangonz/Carbonate/runLists # Directory for runlists
export OUTPUT_DIR=/N/slate/frangonz/ParsedData/lowT-Full # Where we're saving to
mkdir ${OUTPUT_DIR} # And if it doesn't already exist, save.

# Load the "Coincidence Parameters" we want to use
export INI_WINDOW=50
export TEL_WINDOW=1000
export PROMPT_WIN=1000 # Does nothing right now.
export PE_THRESH=8
export COINC_ALG=9 # Might want to set up this as paired?

# TODO: Move some of these out of the debug queue
# For now assumes we're in
qsub ${PBS_DIR}/long_hold_background.pbs 
qsub ${PBS_DIR}/background_scans.pbs

sleep 15m
# It's faster to run everything in the debug queue than suffering through Carbonate...
# Probably should migrate to BR3 eventually
#sleep 1h # Wait for the background runs to finish before generating data.
# TODO: Run an instance of background_general.py when these are complete
qsub ${PBS_DIR}/generate_data.pbs
