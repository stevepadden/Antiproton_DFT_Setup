#Actually performing the least square fitting routine, using a LM routine to fit data.
import scipy.optimize as opt
import pandas as pd
from numpy import exp,power,array,linspace,sqrt,longdouble
import numpy as np
#import scipy.special.expit as exp
#from scipy import exp
import matplotlib.pyplot as plt
#import bigfloat



x = [1.,1.,1.,1.,1.,1.,1.,1.]#Initial "guesses"
df = pd.read_pickle("ScreeningColumn.pkl") #reading in the screened functions previously created
print(df)
xdata = df.Spacing
dffit = pd.DataFrame()
#Selection of potential functions that fit the bill!
def threeexp(x,r):
    f = x[0]*exp(x[1]*r)+x[2]*exp(x[3]*r)+x[4]*exp(x[5]*r)
    return f

def fourexp_ZBL(x,r):

    f = x[0] * exp(-x[1] * r) + x[2] * exp(-x[3] * r) + x[4] * exp(-x[5] * r) + x[6] * exp(-x[7] * r)

    return f

def fourexp_pbar(x,r):
    f = x[0] * exp(-x[1] * r) + x[2] * exp(-x[3] * r) + x[4] * exp(-x[5] * r) + ((1-x[0]-x[2]-x[4])* exp(-x[7] * r))
    return f
#Residuals required
def residual_proton(x,xdata,ydata):
    return(ydata-fourexp_ZBL(x,xdata))

def residual_pbar(x,xdata,ydata):
    return(ydata-fourexp_pbar(x,xdata))
#Some basic figure setup
fig = plt.figure()
ax = plt.gca()
ax.set_yscale('log')
ax.set_xscale('log')
#x.set_ylim(10e-3,10)
#ax.set_xlim(10e-1,10e1)
ax.set_xlabel("Interatomic Distance ($\AA$)")
ax.set_ylabel("Screening Function $\phi_r$")


fig1 = plt.figure()
ax1 = plt.gca()
# ax1.set_ylim(10e-4,10)
# ax1.set_xlim(10e-3,10e1)
ax1.set_xlabel("Interatomic Distance ($\AA$)")
ax1.set_ylabel("Screening Function $\phi_r$")

fig2,ax2 = plt.subplots()
ax2.set_yscale('log')
ax2.set_xscale('log')
#x.set_ylim(10e-3,10)
#ax.set_xlim(10e-1,10e1)
ax2.set_xlabel("Interatomic Distance ($\AA$)")
ax2.set_ylabel("Screening Function $\phi_r$")

#For each column in the overall data set
for col in df.columns:
    if col!="Spacing": #If it isnt the spacing
        ydata=df.loc[:,col] #extract the data from that column
        print(ydata)
        if "Proton" in col: #If its protons
            popt = opt.least_squares(residual_proton, x, method='lm',args=(xdata, ydata),max_nfev=30000)  # Fitting to the function residual which calculates the
            # residuals between our function and the data, using initial guess x
            # and arguments xdata and ydata
            xarray = popt.x  # Extracting the x 'guess' array from POPT returned from fit
            print(xarray)
            yn = fourexp_ZBL(xarray,xdata)  # using an the fit params from above (xarray) along with x-coords to plot the function returned n
            print(popt)
            residuals = ydata - fourexp_ZBL(xarray,xdata) #Calculating residuals
            #Below lines calculate the R^2 value
            ss_res = np.sum(residuals**2)
            ss_tot = np.sum((ydata-np.mean(ydata))**2)
            r2 = 1 - (ss_res/ss_tot)
            print("RSQUARED")
            print(r2)
        elif "Pbar" in col: #Same as above in terms of structure, except now uses a different fitting method for the residuals to reflect its antiprotons not protons

            popt = opt.least_squares(residual_pbar,x,method='lm',args=(xdata,ydata),max_nfev=30000)
            xarray = popt.x
            yn = fourexp_pbar(xarray,xdata)
            print(popt)
            residuals = ydata - fourexp_pbar(xarray,xdata)
            ss_res = np.sum(residuals**2)
            ss_tot = np.sum((ydata-np.mean(ydata))**2)
            r2 = 1 - (ss_res/ss_tot)
            print("RSQUARED")
            print(r2)
        else:
            print("Column found without Pbar or Proton in title - unsure of functional form to fit - please adjust tables")
            continue


        #Some plotting routines that show the least squared fits and their ability to fit the data correctly
        tname = col
        if "Proton" in tname:
            marker = "o"
            size = 35
        else:
            marker = "+"
            size = 45
        tname = tname.replace("_ORCA_MP2", "")
        tname = tname.replace("_", "-")
        tname = tname.replace("Pbar", "$\overline{p}$")
        tname = tname.replace("Proton", "p")

        #ax.scatter(xdata,ydata,marker='+',s=20,label=col)
        str2 = tname+" Fit \n$R^2$ : " + str(r2)[:7]
        ax.plot(xdata,yn,label=str2)
        ax.legend(loc='best',frameon=False,ncol=2)

        ax1.scatter(xdata,ydata,marker=marker,s=size)#,label=col)
        ax1.plot(xdata,yn,label=str2)
        ax1.legend(loc='best',frameon=False,ncol=2)
        ax1.set_xlim([0.008,8])

        dffit[col] = xarray     #adding a new column into an existing dataframe
        print("\n \n " +col+" \n \n")
        print(popt.x)

        ax2.plot(xdata,yn,label=str2)
        ax2.scatter(xdata,ydata,marker=marker,s=size)
        ax2.legend(loc='best',frameon=False,ncol=2)

fig.show()
fig1.show()
fig2.show()
print(dffit)
dffit.to_pickle("ScreenFuncFitParams.pkl")
dffit.to_csv("ScreenFuncFitParams.csv")
