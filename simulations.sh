#!/bin/bash

# This is the number of times that we want to be able to run our simulations.
SIMULATION_NUMBER = 100

# To store all of the resulting files from our simulations, we want to be able 
# to have a directory where we can store everything.
OUTPUT_DIR="outputted_simulations" 
mkdir -p "$OUTPUT_DIR"

for i in $(seq 1 $SIMULATION_NUMBER); do
	# Proof that we got into the for-loop and the simulation is running. 
	echo "Running simulation $i, please wait ..."

	# During each of the simulations, we want to be able to have a 
	# subdirectory to store all of the resulting simulation output files
	# and have a number to associate the simulation number with the file
	# number. 
	SIM_DIR = "$OUTPUT_DIR/simulation_$i"
	mkdir -p "SIM_DIR"

	# Note that now we are inside the OUTPUTDIR/simulation_[iteration]
	# After the directories are made, we want to be able to go through our 
	# our commands that helps us to build our files 	
	make data-cleanup 
	make reset
	make hepatocyte-sample 
	make 
	./hepatocyte-cryopreservation 

	mv afterfrozen_Cellnumber_and_stateoftheirneighbours.txt "$SIM_DIR/"
	mv beforefrozen_Cellnumber_and_stateoftheirneighbours.txt "$SIM_DIR/"
	mv final.svg "$SIM_DIR/"
	mv final.xml "$SIM_DIR/"
	mv Gillespie_alpha10_rathepatocyte_22ncell_B400.txt "$SIM_DIR/"
	mv initial.svg "$SIM_DIR/"
	mv initial.xml "$SIM_DIR/"
	mv rat_hepatocyte_22cells_B400.txt "$SIM_DIR/"
	mv rathepatocyte_freezing_time_B400.txt "$SIM_DIR/"
	mv rat_hepatocyte_Tau_22cells_B400.txt "$SIM_DIR/"
	mv snapshot*.svg "$SIM_DIR/"
	mv temp_points.csv "$SIM_DIR/"
done

echo "Officially completed all the simulations, exiting now..."
