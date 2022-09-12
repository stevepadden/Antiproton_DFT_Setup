# Antiproton_DFT_Setup
Setup scripts used in construction force fields required for MDRange simulations of antiprotons using the DFT code ORCA

MAIN handles pretty much all running once the construction has been undertaken, each script runs using the one before it in the line as an input.

Orca_Plot extracts all the required data from a collection of input files, it also shows how to create an ORCA scan file required by the program

Screening_Function takes all datapoints and converts them from a potential energy surface into a screened coulomb potential
  
LeastSqFit uses a Levenberg-Marquadt algorithim to fit a continous function to the screened coulomb potential

Finally Continuous_Potential_MDRange converts back into a potential form as required by MDRange, it also uses some symbolic calculus to find force field calculations, and finally produces a table in the correct format
