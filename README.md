# neutralviralpassage
This folder contains all the necessary files to reproduce the neutral simulation results in "Mutation rate of SARS-CoV-2 and emergence of mutators during experimental evolution". Specifically the data for Fig.2a and Fig. S3.
Find the following files:
- "Virus_Evolution_Neutral.R" contains the R script for the neutral simulations without migration. Running such script has it is will generate 2 outcomes: "Virus_Evolution_N*_L*.RData" and "Virus_Evolution_Pool.RData".
The first contains the raw data for the study of the mutation spectrum while the second contains the genome pool that is used for the simulations with migration.
- "Virus_Evolution_Neutral_Migration.R" contains the R script for the simulations with migration. It requires an input files with a pool of genotypes (e.g. the previously generated "Virus_Evolution_Pool.RData").
Running such script as it is will generate "Virus_Evolution_M*_N*_L*.RData" containing the raw data for the study of the mutation spectrum with 10% migration.
- "Virus_Code_Analysis.RMD contains the RStudio script to obtain the tables to be plotted (e.g. Fig.S3). It requires that the raw data files were already generated.

Running these 3 script in the given order will produce the tables necessary for plots in Fig.2a and Fig. S3.

Parameters of the simulations can be changed at the beginning of the *.R files.

For any problem or doubts email: mamicone@igc.gulbenkian.pt    
