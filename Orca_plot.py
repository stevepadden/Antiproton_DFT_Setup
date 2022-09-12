import sys
from numpy import linspace
import matplotlib.pyplot as plt
import math
import pandas as pd
import os
#Example Orca commands
# /mnt/d/Orca$ ./Orca.exe OrcaScan.inp > ../Orca_Outputs/Pbar_Al.out
# /mnt/d/Orca$ ./Orca.exe OrcaScan.inp > ../Orca_Outputs/Proton_Al.out
# /mnt/d/Orca$ ./Orca.exe OrcaScan.inp > ../Orca_Outputs/Pbar_Si.out
# /mnt/d/Orca$ ./Orca.exe OrcaScan.inp > ../Orca_Outputs/Proton_Si.out
# /mnt/d/Orca$ ./Orca.exe Orca50Scan.inp > ../Orca_Outputs/Proton_Si_50.out

#Example orca scan file - place this into a file called OrcaScan.inp (after adjustment) and run using one of the above commands

# ! MP2RI def2-TZVP def2-TZVP/C TightSCF Opt        #Specifiying the fit types and basis sets (it moller plessent 2 def2-TZVP basis with tightSCF conversion optimised
# %geom Scan        #telling it we want it to scan
# B 0 1 =8,0.007, 50 # Q-AL distance that will be scanned       #Scans from 0.007 bohr to 8 bohr in 50 steps between atoms in order 0 and 1
# end       #need to use 2 ends? dunno why
# end
#
# * xyz 0 1	#Ctype, charge, multiplicity Here we use Internal Co-ords with 0 charge and multiplicty of 1
# #recall multiplicity = (2S + 1) thus for spin half (2*1/2 + 1 = Multiplicity 2) - See pages 37 (68)
# Q -1  0 0 0.000    #point charge (Q) at -1 charge (anti proton) with positions 0x 0y 0z
# Si	  0 0 8      #silicon molecule at 0x 0y 8z but scans DOWN to 0.007
#*

#RUN ORDER
#1) Orca_Plot.py    # Takes the data from the .out files, maniupulates them into ev from hatrees and angstroms from bohr
                    #It is necessary to change the charges here aswell as the names of the files for both the 0-x angstroms (set in orca)
                    #aswell as changing the names required for the 50 angstrom scan which is subtracted from the close spacing to reveal
                    #unscreened nuclear charges     #Creates a pickle file Orca_Values which is used later
#2) ScreeningFunction.py    #This performs the screened coulomb calculations - need to check the values here because this is where the simulation is essentially determined by accuracy
                            #It then returns a dataframe ScreenFunction.pkl (pickle format) which can be used later
#3)Recaster - this takes the pickled dataframe from screening function and makes it into something more suitable for fitting with matlab etc
                #this could really be done in screening function but it is at times useful to be able to just return all the data as a list not a table so i decided to keep seperate
                #returns 2 datasets one called (ScreeningData.Csv) for use with matlab as it can be read in much easier and another called screeningCoulmn.csv for
                #use with the scipy fitting
#4) Least Sq fit - performs fitting through use of a residual function on the defined function, need to read file to understand as too lengthy to put in comment byt
                #essentiall it does Levenberg–Marquardt fitting of least squares to the dataset, returns the fit parameters as a pickle output
                #ScreenFuncParams.pkl and also does some log log plotting of all the input datasets to show how well the fit matches the actual data
from sympy.integrals.meijerint_doc import obj

path= 'Orca_Outputs'        #CHANGE THE PATH TO WHERE THE .OUT FILES ARE LOCATED!
Eh2Ev = 1# 27.2113834      #Hatree to Ev conversion
Bohr2Angstrom = 1.88973#1 #0.5291772083 Data is already in angstroms!!!
# = 0.5291772083
ang2meter =  1e-10
hatree2joules = 4.35974e-18
df = pd.DataFrame()#(columns=["RunElements","Spacing","Distance"])
names = ['Proton_Al','Proton_Si','Pbar_Al','Pbar_Si','Pbar_Ne','Pbar_Ti','Pbar_Cu','Pbar_Ag','Pbar_Au'] #The "file" names that identify each species scan
fifties = ['Proton_Al_50','Proton_Si_50','Pbar_Al_50','Pbar_Si_50','Pbar_Ne_50','Pbar_Ti_50','Pbar_Cu_50','Pbar_Ag_50','Pbar_Au_50'] #The "file" names that contain the potential energy scan at 50 Bohr

z1charges = [1,1,-1,-1,-1,-1,-1,-1,-1]       #Proton vs antiproton charges
z2charges = [13,14,13,14,10,22,29,47,79]      #Nuclear "target" charge
ext='.out'  #The extension on the end of each file
for name,fifty,z1,z2 in zip(names,fifties,z1charges,z2charges): #Extracting the data
    if "Proton" in name:    #Using protons as open circles
        marker = "o"
        size=35
    else:
        marker = "+"    #And antiprotons as + marks.
        size=45
    #Extracting the energy for 50 angstroms
    ang_fifty = open(os.path.join(path,fifty+ext))  #Opening the 50 angstrom scan
    fifty_MP2=11   #data is N lines up we find the data here.
    fifty_counter=-fifty_MP2
    for line in ang_fifty:
        fifty_counter = fifty_counter+1
    ang_fifty = open(os.path.join(path, fifty+ext))
    i=1
    fifty_val = 0
    # print(fifty_counter)
    #Extracting the 50 ang value
    for line in ang_fifty.readlines():
        if i == fifty_counter:
            # print("LINE READING")
            # print(line)
            # print("LINE READ")
            line=line.strip("\n")
            x=line.split('  ')[1]
            fifty_val = float((x.split(' ')[1]))
            # print(fifty_val)

        i=i+1
    # print("The value at fifty is :  {}".format(fifty_val))

    #Extracting the scan data
    Num_points= 49
    counter = 1
    lines=linspace(12,12+Num_points,num=Num_points+1,dtype=int) #Setting the number of lines we are extracting!
    filename = os.path.join(path,name+ext)
    file = open(filename)
    for line in file:
        counter= counter +1
    # print(counter)
    lines = [counter-line for line in lines]        #Selecting the correct lines! going from 12 up from bottom
    # print(lines)
    i=1

    file = open(filename)   #Opening the scan file
    spacing=[]
    hatrees=[]
    meters =[]
    joules = []
    for line in (file.readlines()): #Readint the lines
        if i in lines:  #Selecting only from the lines we want
            # print (line)
            line=line.strip("\n")   #Stripping new line tokens
            x=line.split('   ')[1]  #Splitting into requried data strcuture
            # print(x)
            spacing.append((float(x.split(' ')[0])) * Bohr2Angstrom)          #Be aware of the conversions here! Have to convert to angstrom
            #hatrees.append((float((x.split(' ')[1]))) * Eh2Ev)
            hatrees.append((float( (x.split(' ')[1]) ) - fifty_val)*Eh2Ev)  #Subtracting the 50 val from the scan val, then converting from Hartrees to eV
            meters.append((float(x.split(' ')[0])) * ang2meter)
            joules.append((float( (x.split(' ')[1]) ) - fifty_val)*hatree2joules)
        i=i+1


    file.close()

    spacing= (list(reversed(spacing)))      #largely irrelevant step but i prefer to have the list flipped so we go from smallest seperation to largest
    hatrees= (list(reversed(hatrees)))
    meters = (list(reversed(meters)))
    joules = (list(reversed(joules)))
    fig= plt.figure()   #type: plt.figure()
    ax = plt.gca()
    tname=name
    tname= tname.replace("_","-")
    tname = tname.replace("Pbar", "$\overline{p}$")
    tname = tname.replace("Proton","p")
    plt.scatter(spacing,hatrees,label=tname,marker=marker,color="royalblue",s=size)
    ax.set_xlabel("Interatomic Distance ($\AA$)")
    ax.set_ylabel(r"Potential Energy (eV)")
    ax.axhline(0, color='black')
    #ax.set_xlim(0.03,8)
    ax.set_ylim(-1,0.2)
    ax.set_ylim(-1,1)
    ax.legend(loc="best",frameon=False)
    fig.show()
    savestring=name+".png"
    fig.savefig(savestring)
    dftemp = pd.DataFrame([{"Name":name,"Spacing":spacing,"eV":hatrees,"z1":z1,"z2":z2 , "Joules" : joules , "Meters": meters}])           #Constructing a temporary dataframe
    df=df.append(dftemp,ignore_index=True)          #adding to the main dataframe - notice the need to be explicit to update, ie df = df ...

df=df.set_index('Name')         #Setting the index, makes accessing easier (See ScreeningFunction_


#Use some type hinting here for speed etc.
fig = plt.figure()  # type: plt.figure
ax = plt.gca()      # type: plt.gca
for i,line in df.iterrows():
    tname = i
    tname= tname.replace("_","-")
    if "Proton" in tname:
        marker = "o"
        size = 35
    else:
        marker = "+"
        size=45
    tname = tname.replace("Pbar", "$\overline{p}$")
    tname = tname.replace("Proton","p")
    # print(df.at[i,'Spacing'])
    plt.scatter(df.at[i,'Spacing'],df.at[i,'eV'],label=tname,marker=marker,s=size)
ax.set_xlabel("Interatomic Distance ($\AA$)")
ax.set_ylabel(r"Potential Energy (eV)")
ax.set_xlim(0.04,8)
ax.set_ylim(-1,1)
ax.axhline(0,color='black')
ax.legend(loc="upper right", ncol=2,frameon=False)
fig.savefig("ORCA_Potentials.png")
fig.show()

df.to_pickle("Orca_Values.pkl")

for index,x in df.iterrows():
    print(index)
    print(x.Spacing)
    print(x.eV)


#plt.style.use('dark_background')
# fig = plt.figure()  # type: plt.figure
# ax = plt.gca()      # type: plt.gca
# frameon=False

#plt.grid(color='#666666', linestyle='--', linewidth=0.3)
# plt.grid(color='lightgrey', linestyle='--', linewidth=0.3)
# plt.scatter(df.at['Proton_Si','Spacing'],df.at['Proton_Si','eV'],50,label="Proton - Silicon",marker="1",color="w")
# plt.scatter(df.at['Pbar_Si','Spacing'],df.at['Pbar_Si','eV'],50,label="Antiproton - Silicon",marker="2",color="r")
# plt.scatter(df.at['Pbar_Ne','Spacing'],df.at['Pbar_Ne','eV'],40,label="Antiproton - Neon",marker="^",color="r",alpha=1)
# ax.set_xlabel("Separation ($\AA$)",fontsize=15)
# ax.set_ylabel(r"Potential Energy (eV)",fontsize=15)
# ax.set_xlim(0.04,8)
# ax.set_ylim(-1.5,1.5)
# ax.spines['right'].set_visible(False)
# ax.spines['top'].set_visible(False)
# #ax.axhline(0,color='black')
# ax.axhline(0,color='lightgrey')
# plt.legend(loc="best",frameon=False,fontsize=13)
# fig.savefig("ORCA_Potentials.png")
# fig.show