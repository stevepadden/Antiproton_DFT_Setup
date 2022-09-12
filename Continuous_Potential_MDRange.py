#We now convert back into a constant function, we also use some computed calculus here to find force fields.

from sympy import *
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.pyplot import rcParams
from math import pi as pi
path = r"D:\MDrange_Shared\Reppot_Files" #Path to where the repulsive potential files used within MDRange are saved too
def tablewrite(colname,NNUM,rdist,energyfunc,force):    #Function to write the repulsive potential tables, adjust to required program strucutre as required
    file = open(path+"/reppot_"+colname+".0.1.in", "w+")
    for i in range(len(rdist)):
        file.write("{} {} {}\n".format(rdist[i],energyfunc[i],force[i]))
    file.close()

#Setting up some conversions and numbers required - small lists of constants.
fourpi = 4*pi

epsilon =8.85e-12# 8.85e-12
e = 1.6e-19#1.6e-19

ang2meter = 1e-10

startang = 0.001
endang = 8
NNUM=500

hbar = 6.5821e-16      #hbar ev/s
fine_struc = 1./137.    #Dimensionless
c_ang = 2.998e18    #Speed of light Angstroms per second
hbar_c = 2000
#r_meter= np.linspace(start)



df = pd.read_pickle("ScreenFuncFitParams.pkl")  #reading in screening functions
print("Found Dataframe with columns for : \n {} \n Now proceeding to plot potential energy curves for Pbar events \n \n ".format(df.columns))


#Setting up plots
fig = plt.figure()
ax = plt.gca()
ax.set_xlabel("Interatomic Distance ($\AA$)")
ax.set_ylabel("Potential Energy (eV)")
ax.axhline(0,color='black')

#rcParams['text.usetex'] = True
fig1 = plt.figure()
ax1 = plt.gca()
ax1.set_xlabel("Interatomic Distance ($\AA$)")
ax1.set_ylabel(r"Force   $-\frac{\delta V(r)}{ \delta(r)}$     (N)")
ax1.axhline(0,color='black')
ax1.usetex=True



#Doing the differential - we use sympy but we need to tell it now which we intend to use symbolically
r = Symbol('r')  # Setting the symbol we need to differentiate over!

# def fprime(func,r):
#     return func.diff(r)

#The plan now is to iterate over all screening functions. extracting them as the function they are in the table and appending them
#to a coulomb potential - previously these all saved in the same file but now it goes into its own reppot file
for col in df.columns:  #Each column is named after what it is representing - so we just iterate through these columns

    if "Pbar" in col:       #if its pbar then z = -1
        #Just setting the correct params
        if "Si" in col:
            Z2 = 14
        elif "Al" in col:
            Z2 = 13
        elif "Ne" in col:
            Z2 = 10
        elif "Ti" in col:
            Z2 = 22
        elif "Cu" in col:
            Z2 = 29
        elif "Ag" in col:
            Z2 = 47
        elif "Au" in col:
            Z2 = 79
        else:
            Z2 = input("Please input the target atomic number for : {}".format(col))
        Z1 = -1     #pbar charge=-1

        #Selecting the right dataset and which function it contains
        x =  df.loc[:,col]   #selecting Proton in Aluminium from Dataframe
        if len(x) == 8:
            #exp_func = x[0]*exp(x[1]*r)+x[2]*exp(x[3]*r)+x[4]*exp(x[5]*r)+1-x[0]-x[2]-x[4]*exp(x[7]*r)  #using 4 exponential function
            #exp_func = fourexp_pbar(x,r)
            exp_func = f = x[0] * exp(-x[1] * r) + x[2] * exp(-x[3] * r) + x[4] * exp(-x[5] * r) + ((1-x[0]-x[2]-x[4])* exp(-x[7] * r))
            print("Found 4 exp function")
        elif len(x) == 6:
            exp_func = x[0]*exp(x[1]*r)+x[2]*exp(x[3]*r)+x[4]*exp(x[5]*r)
            print("Found 3 exp function")
        else:
            print("Unkonwn function found! Please see codes in LeastSquaresFit and update VrContinuousWDiff to account!")
            break

        #doing the actual calculation
        func = (hbar_c*fine_struc*Z1*Z2/r)*exp_func
        funcdiff = -func.diff(r)
        print("Function : \n {}".format(func))
        print("Function differentiated : \n {}".format(funcdiff))

        rdist = np.linspace(startang, endang, NNUM)
        print("Plotting between {} and {}: \n ".format(startang,endang))
        # print(rdist)


        f = lambdify (r, func, "numpy")         #telling the symbolic function we can now evaluate r in the symfunc func using numpy
        #print(f(rdist))
        fprime = lambdify(r,funcdiff,"numpy")
        tname = col
        tname = tname.replace("_ORCA_MP2", "")
        tname = tname.replace("_", "-")
        tname = tname.replace("Pbar", "$\overline{p}$")
        tname = tname.replace("Proton", "p")
        ax.plot(rdist,f(rdist),label=tname)      #Using the function f evaluated at rdist array and plotting
        ax1.plot(rdist,fprime(rdist),label=tname)
        df1 = pd.DataFrame(columns=['NNUM', 'ENERGY', 'FORCE'])
        df1.NNUM = rdist
        df1.Energy = func
        df1.FORCE = fprime(rdist)
        df1.to_pickle(col+".pkl")

        y = f(rdist)
        yprime = fprime(rdist)
        #filtered = filter(lambda num: num > 0, y)
        #filtered = [x for x in y if x <= 0]
        index = 0
        """
        for i in y:
            if i >0:
                print(i)
                y = np.delete(y,index,0)
                yprime = np.delete(yprime,index,0)
                rdist = np.delete(rdist,index,0)
            index = index+1
        print("FILTERED")
        print(y)
        """
        #Necessary stops to remove suprious values bigger than 0
        indexbig0 = np.argwhere(y>0)
        print(indexbig0)
        rdist = np.delete(rdist,indexbig0)
        y = np.delete(y,indexbig0)
        yprime = np.delete(yprime,indexbig0)
        print("Filtered")
        print(y)


        tablewrite(col,NNUM,rdist,y,yprime)



    elif "Proton" in col:
        if "Si" in col:
            Z2 = 14
        elif "Al" in col:
            Z2 = 13
        elif "Ne" in col:
            Z2 = 10
        else:
            Z2 = input("Please input the target atomic number for : {}".format(col))


        Z1 = 1
        x = df.loc[:, col]  # selecting Proton in Aluminium from Dataframe
        print(x)
        if len(x) == 8:
            #exp_func = x[0] * exp(x[1] * r) + x[2] * exp(x[3] * r) + x[4] * exp(x[5] * r) + x[6] * exp( x[7] * r)  # using 4 exponential function
            #exp_func = fourexp_ZBL(x,r)
            exp_func = f = x[0] * exp(-x[1] * r) + x[2] * exp(-x[3] * r) + x[4] * exp(-x[5] * r) + x[6] * exp(-x[7] * r)
            print("Found 4 exp function")
        elif len(x) == 6:
            exp_func = x[0] * exp(x[1] * r) + x[2] * exp(x[3] * r) + x[4] * exp(x[5] * r)
            print("Found 3 exp function")
        else:
            print(
                "Unkonwn function found! Please see codes in LeastSquaresFit and update VrContinuousWDiff to account!")
            break
        #func = (1 / fourpi * epsilon) * (Z1 * Z2 * (e ** 2) / (r*ang2meter)) * exp_func
        func = (hbar_c * fine_struc * Z1 * Z2 / r) * exp_func
        # print(func)
        funcdiff = -func.diff(r)
        print("Function differentiated : \n {}".format(funcdiff))

        rdist = np.linspace(startang, endang, NNUM)
        print("Plotting between {} and {}: \n ".format(startang, endang))
        # print(rdist)

        f = lambdify(r, func,
                     "numpy")  # telling the symbolic function we can now evaluate r in the symfunc func using numpy
        # print(f(rdist))
        fprime = lambdify(r,funcdiff,"numpy")
        tname = col
        tname = tname.replace("_ORCA_MP2", "")
        tname = tname.replace("_", "-")
        tname = tname.replace("Pbar", "$\overline{p}$")
        tname = tname.replace("Proton", "p")


        ax.plot(rdist, f(rdist), label=tname)  # Using the function f evaluated at rdist array and plotting
        ax1.plot(rdist,fprime(rdist),label=tname)

        df1 = pd.DataFrame(columns=['NNUM','ENERGY','FORCE'])
        df1.NNUM = rdist
        df1.ENERGY = func
        df1.FORCE = fprime(rdist)
        df1.to_pickle(col+".pkl")

        tablewrite(col,NNUM,rdist,f(rdist),fprime(rdist))


ax.legend(loc='best',frameon=False,ncol=2)
ax.set_xlim([0.1,7])
ax.set_ylim([-6,2])
fig.show()

ax1.legend(loc='best',frameon=False,ncol=2)
ax1.set_xlim([0.2,7])
ax1.set_ylim([-6,2])
fig1.show()


