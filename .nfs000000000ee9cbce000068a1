#!/bin/bash

# Deals with the resource or device busy problem
fuser -k ./output
rm -rf output_simulations

# This is the number of times that we want to be able to run our simulations.
SIMULATION_NUMBER=100
FILE_NUMBER=89
make data-cleanup
make reset 
 

# To store all of the resulting files from our simulations, we want to be able 
# to have a directory where we can store everything.
OUTPUT_DIR="./output_simulations" 
mkdir -p "$OUTPUT_DIR"
run_simulation(){
	fuser -k ./output 
	local i=$1 # First argument of the function is taken
	echo "Running simulation ${i}, please wait ..."

	# During each of the simulations, we want to be able to have a 
	# subdirectory to store all of the resulting simulation output files
	# and have a number to associate the simulation number with the file
	# number. 
	SIM_DIR="./simulation_${i}"
	mkdir -p "$SIM_DIR" && echo "Directory created: $SIM_DIR" || { echo "Failed to create $SIM_DIR"; exit 1; }
	
	# Note that now we are inside the OUTPUTDIR/simulation_[iteration]
	# After the directories are made, we want to be able to go through our 
	# our commands that helps us to build our files 	
	make data-cleanup 
	make reset
	make hepatocyte-sample 
	make 
	./hepatocyte-cryopreservation 
	 
	#cp -rf output/*.* "$SIM_DIR/".

	cp -rf output/afterfrozen_Cellnumber_and_stateoftheirneighbours.txt "$SIM_DIR/".
	cp -rf output/beforefrozen_Cellnumber_and_stateoftheirneighbours.txt "$SIM_DIR/".
	cp -rf output/final.svg "$SIM_DIR/".
	cp -rf output/final.xml "$SIM_DIR/".
	cp -rf output/Gillespie_alpha10_rathepatocyte_22ncell_B400.txt "$SIM_DIR/".
	cp -rf output/initial.svg "$SIM_DIR/".
	cp -rf output/initial.xml "$SIM_DIR/".
	cp -rf output/rat_hepatocyte_22cells_B400.txt "$SIM_DIR/".
	cp -rf output/rathepatocyte_freezing_time_B400.txt "$SIM_DIR/".
	cp -rf output/rat_hepatocyte_Tau_22cells_B400.txt "$SIM_DIR/".
	cp -rf output/snapshot*.svg "$SIM_DIR/".
	cp -rf output/temp_points.csv "$SIM_DIR/".
} 
export -f run_simulation
seq 1 $SIMULATION_NUMBER | xargs -P 4 -I {} bash -c 'run_simulation "$@"' _ {}

# Once the simulation is done, one for-loop is ran here to go through all of those 
# simulation folders and add them into one GIANT folder called output_simulations 
for i in $(seq 1 $SIMULATION_NUMBER); do
	cp -rf ./simulation_${i} ./output_simulations/.
	rm -rf ./simulation_${i}
done
