#!/bin/bash

# Deals with the resource or device busy problem
fuser -k ./output

# This is the number of times that we want to be able to run our simulations.
SIMULATION_NUMBER=5
FILE_NUMBER=89
make data-cleanup
make reset 

# To store all of the resulting files from our simulations, we want to be able 
# to have a directory where we can store everything.
OUTPUT_DIR="output_simulations" 
mkdir -p "$OUTPUT_DIR"

for i in $(seq 1 $SIMULATION_NUMBER); do
	fuser -k ./output 
	# Proof that we got into the for-loop and the simulation is running. 
	echo "Running simulation $i, please wait ..."

	# During each of the simulations, we want to be able to have a 
	# subdirectory to store all of the resulting simulation output files
	# and have a number to associate the simulation number with the file
	# number. 
	SIM_DIR="$OUTPUT_DIR/simulation_$i"
	mkdir -p "$SIM_DIR"
	echo "Opening a subdirectory inside output for simulation_$i"

	# Note that now we are inside the OUTPUTDIR/simulation_[iteration]
	# After the directories are made, we want to be able to go through our 
	# our commands that helps us to build our files 	
	make data-cleanup 
	make reset
	make hepatocyte-sample 
	make 
	./hepatocyte-cryopreservation 

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
	
	# This part of the for-loop is to deal with the iteration of the snapshot
	# image files - We want to ensure that we are able to go from 0 to 89th 
	# image file and copy them to their respective simulation numbered folder
	for i in $(seq 0 $FILE_NUMBER); do
		PADDED_NUMBER=$(printf "%02d" $i) # 00, 01, 02 ... 89
		cp -rf output/snapshot000000${PADDED_NUMBER}.svg "$SIM_DIR/".
	done 
done

echo "Officially completed all the simulations, exiting now..."
