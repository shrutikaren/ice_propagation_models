#!/bin/bash 

DIRECTORY=./output_simulations 
TOTAL=100

# This bash files focuses on taking out two of the text files from our simulation
# specific output folder so that we can take these files with us and then perform
# some kind of graphical representation. The first time will be done with Fatima's 
# code then it will be done with the newly added code to see if the graphical 
# representation is the same. 

# All of the Gillespie_alpha10_rathepatocyte_22ncell_B400.txt is stored in 
# ./time_simulated files. All of the rathepatocyte_freezing_time_B400.txt is placed
# in ./freezing_time directory.
NEWDIRECTORY=./time_simulated_files 
NEWDIRECTORY2=./freezing_time

mkdir -p ${NEWDIRECTORY}
mkdir -p ${NEWDIRECTORY2}
for i in $(seq 1 $TOTAL);do
	# mv ./output_simulations/simulation_${i}/rathepatocyte_freezing_time_B400.txt ./output_simulations/simulation_${i}/rathepatocyte_freezing_time_B400_simulation_${i}.txt

	# cp -rf ./output_simulations/simulation_${i}/rathepatocyte_freezing_time_B400_simulation_${i}.txt $NEWDIRECTORY
	# rm /output_simulations/simulation_${i}/rathepatocyte_freezing_time_B400_simulation_${i}.txt

	# mv ./output_simulations/simulation_${i}/Gillespie_alpha10_rathepatocyte_22ncell_B400.txt ./output_simulations/simulation_${i}/Gillespie_alpha10_rathepatocyte_22ncell_B400_simulation_${i}.txt

	cp -rf ./output_simulations/simulation_${i}/Gillespie_alpha10_rathepatocyte_22ncell_B400_simulation_${i}.txt $NEWDIRECTORY2

	rm /output_simulations/simulation_${i}/Gillespie_alpha10_rathepatocyte_22ncell_B400_simulation_${i}.txt

done
