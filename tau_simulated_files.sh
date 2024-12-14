#!/bin/bash 

DIRECTORY=./output_simulations 
TOTAL=100
NEWDIRECTORY=./tau_simulated_files
mkdir -p ${NEWDIRECTORY}
for i in $(seq 1 $TOTAL);do
	rm /output_simulations/simulation_${i}/rat_hepatocyte_Tau_22cells_B400_simulation_${i}.txt
	mv ./output_simulations/simulation_${i}/rat_hepatocyte_Tau_22cells_B400.txt ./output_simulations/simulation_${i}/rat_hepatocyte_Tau_22cells_B400_simulation_${i}.txt

	cp -rf ./output_simulations/simulation_${i}/rat_hepatocyte_Tau_22cells_B400_simulation_${i}.txt $NEWDIRECTORY
done
