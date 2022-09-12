#This file constructs the screening functions required. It casts the energy into one which uses a screened coulomb function!

import pandas as pd
from math import pi
import numpy as np
import matplotlib.pyplot as plt
e =1.6e-19 #does this work because we are working in angstroms not meters?
episilon = 8.85e-12# * 1e-10#*5.2918E-11#elec     #Farads per angstrom    - is this is for sure correct?
meters2ang = 1e10
# e=1.6e-19
# episilon = 8.815e-12
fourpi = 4*pi

df = pd.read_pickle("Orca_Values.pkl")


CoulombScreening = []
spacing = []
df2 = pd.DataFrame()


for index,x in df.iterrows(): # For each data set
    #for space,potential in zip(x.Spacing,x.eV):
    for space,potential in zip(x.Meters,x.Joules):  #At each energy and spacing
        Z1 = x.z1 #Charges
        Z2 = x.z2
        var1 = ((1/(fourpi * episilon))*Z1*Z2*((e*e)/(space)))  #Coulomb repulstion at that space
        #print (var1)
        var2 = ((potential)/var1)   #modified by the potential energy surface
        #print (var2)
        if abs(var2)<1e100:
            CoulombScreening.append(float(var2))
            spacing.append(float(space*meters2ang))
    df2 = df2.append([{'Spacing': spacing, 'Potential': CoulombScreening,'Name':index + "_ORCA_MP2"}],ignore_index=True)
    print(df2)
    spacing=[]
    CoulombScreening=[]
#Setting up some pickle arrays for the next section and plotting routines.

df2=df2.set_index('Name')         #Setting the index, makes accessing easier

fig = plt.figure()
ax = plt.gca()
ax.set_yscale('log')
ax.set_xscale('log')
for i,row in df2.iterrows():
    tname=i
    tname=tname.replace("_ORCA_MP2" ,"")
    tname=tname.replace("_","-")
    if "Proton" in tname:
        marker = "o"
        size=35
    else:
        marker = "+"
        size=45
    tname = tname.replace("Pbar", "$\overline{p}$")
    tname = tname.replace("Proton","p")
    ax.scatter(df2.at[i,'Spacing'],df2.at[i,'Potential'],label=tname,marker=marker,s=size)
ax.legend(frameon=False,ncol=2)
ax.set_xlabel("Interatomic Distance ($\AA$)")
ax.set_ylabel("Screening Function $\phi_r$")
#ax.set_ylim([0.0001,1])
#ax.set_yscale("log")
fig.show()
df2.to_pickle("ScreeningFunctions.pkl")