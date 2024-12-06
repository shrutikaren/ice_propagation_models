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

	# After the directories are made, we want to be able to go through our 

