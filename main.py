def main():
    import Orca_plot
    import Screening_Function
    import Recaster
    import LeastSqFit
    import Continuous_Potential_MDRange
    import Norlun_Test

'''
How does this work?
step 1 - We take the potentials from Orca outputs
step 2 - Subtract from all potentials the potential at 50 angstroms, ie when very little interaction between the two
step 3 - TO do this adjust the NumPoints in the function Orca_Plot, this determines how many scans were used
step 4 - Plot all the potentials - expect Leonard Jones for matter - matter and strictly "London Dispersion" for antimatter (ie purely repulsive)
step 5 - Construct a screening function from the potentials - phi(r) = V(r)/Coulomb potential -> This is done in Screening_function.py
step 6 - Fit a continuous function to the screenng function - first we pass through Recaster which just maniuplates the dataframe then into Least Square fit
        LeastSqFit uses levenberg marquadt fitting to fit exponential terms to the screening functions, it presents linear/linear scales aswell as log log
        dont get lost in the log log plots!
step 7 - We now use the fitted screening function to construct a continuous potential, previously we had a calculated potential, and fitted a screening function, now we use the fitted screening function to go full circle
        and bring back a continuous potential, this continuous potential is output in to tables suitable for reppot files in MDRange. It also returns the exponential format for the continuous potential.
        It also differentiates the function found wrt distance and returns that, ie from the potential we gather a force profile of the interaction.
        

~notes:
    There are a few implicit conversions in this version of the code - namely the use of joules and meters as a new coulmn, this is to help with the screening function calculation and keep units
for the coulomb repulsion sensible, however using the levenberg-Marquadt algorithim does NOT work well (read at all) using small/large numbers, so although the screening function
is factored in in joules and meters, there is an implicit conversion back to angstroms for the distance, the screening function is simply a scale factor so stays the unit it begins with (ie unitless)
SEE SCREENING_FUNCTION.PY line 35 (by the append to array section) to observe conversion back into Angstroms

The screening function is only ever a way to scale how impactful the coulomb force is! Keep this in mind!

NOTE - There is some useful plotting stuff for formatting ETC in ORCA_PLOT commented out!

'''


if __name__ == '__main__':
    main()

#TODO! Put NEON BACK IN ORCA_PLOT
#TODO! FIX SCREENING FUNCTION UNITS!
#