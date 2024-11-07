# ice_propagation_models
Utilising a Master's degree student code and further exemplifying it to be able to produce more accurate picture of those frozen cells by making use of the parameter "Tau" in our equations.

There are a few things actually: 
- Could we try and make a python code utilising the two parameter model or Mazur's Temperature Dependent Model? 
- Could we implement the Monte Carlo simulation that the paper discusses? 
- Can we briefly discuss on how can we calculate tau and averaging them as connections break? - what literature would you suggest on this?

# Notes from current undergraduate student working on the code:
This ice_propagation_models directory covers all files from fatima's repository copied here and a directory named Temp_Model to store "Mazur's Temperature Dependent Model". There is one more directory called Notes that will store all my literature reviews, presentation slides with Joseph and important thoughts that I might have along the way of this research project. 

# Notes from the previous Master's student:
This code implementes Gillespie (Monte-Carlo) model of ice formation and propagation in hepatocyte disks monolayer of 22 cells tissue, using  PhysiCell http://physicell.org/.
PhysiCell was structured in several critical subdirectories. Here we use the following subdirectories:

1. BioFVM: This includes a working copy of the BioFVM multi-substrate diffusion code. 

2. config: This directory is used for XML configuration files. Any custom configuration files
should be placed in this directory.

3. core: The core library files for PhysiCell are kept in this directory. Users should not modify functions
in the core library.

5. modules: This is where standard code modules palced for distribution with PhysiCell. 
Currently, it includes pathology, MultiCellDS, and SVG (scalable vector graphic) functions. 

4. custom_modules: This directory includes our main code for this project.

6. sample_projects: This directory includes the main.cpp and Makefile.

7. output: Output files will be here.

PhysiCell should successfully compile and run on any C++11 or later compiler that supports OpenMP.
I use g++ on Linux.

To compile and run the code:
1) make data-cleanup
2) make reset 
3) make hepatocyte-sample
4) make
5) ./hepatocyte-cryopreservation

My 3 main code files are: "hepatocyte-cryopreservation.cpp", "hepatocyte-cryopreservation.h" and "main.cpp".  
The "hepatocyte-cryopreservation.cpp" and "hepatocyte-cryopreservation.h" are located in /Codes/PhysiCell/custom_modules
The "main.cpp" is located in /Codes/PhysiCell/sample_projects/hepatocyte
The results are saved in /Codes/PhysiCell/output as .SVG files. In each .SVG file the red cell indicates the unfrozen cell and the blue cell with a red circle inside indicates the frozen cell in the tissue.
