# See readme for information about this file.
# The main function is named Main()
#v1.1
#!/usr/bin/python3

import numpy as np
from scipy.optimize import curve_fit
import math
import matplotlib.pyplot as plt

# Semivariogram function
# Function for spherical semivariogram modeling
def spherical(h,Nugget,Range,Sill):
    k=len(h)
    gammas=np.zeros(k)
    for i in range(k):
            if(h[i]<Range):
                gammas[i]=(Nugget+(Sill-Nugget)*(1.5*h[i]/Range-0.5*h[i]*h[i]*h[i]\
                                /(Range*Range*Range)))
            else:
                gammas[i]=Sill
    return gammas

# Function for exponential semivariogram modeling
def exponential(h,Nugget,Range,Sill):
    k=len(h)
    gammae=np.zeros(k)
    for i in range(k):
        gammae[i]=(Nugget+(Sill-Nugget)*(1-math.exp(-3.0*h[i]/Range)))
    return gammae

# Function for gauss semivariogram modeling
def gauss(h,Nugget,Range,Sill):
    k=len(h)
    gammag=np.zeros(k)
    for i in range(k):
        gammag[i]=(Nugget+(Sill-Nugget)*(1-math.exp(-3.0*h[i]*h[i]/(Range*Range))))
    return gammag

#Function to check if the point is inside or outside the region
# True for inside and False for Outside
def ray_tracing_method(x,y,xlim,ylim):
    n = len(xlim)
    inside = False
    p1x=xlim[0]
    p1y=ylim[0]
    for i in range(n+1):
        p2x = xlim[i%n]
        p2y = ylim[i%n]
        if y > min(p1y,p2y):
            if y <= max(p1y,p2y):
                if x <= max(p1x,p2x):
                    if p1y != p2y:
                        xints = (y-p1y)*(p2x-p1x)/(p2y-p1y)+p1x
                    if p1x == p2x or x <= xints:
                        inside = not inside
        p1x=p2x
        p1y=p2y
    return inside

# Function to find the nearest neighbor
def find_PNeigborns(nNeig, xx, yy, pID, x, y, z):
    k=len(x)
    neigDist=np.zeros(nNeig)
    neigID=np.zeros(nNeig)
    neigZ=np.zeros(nNeig)
    neigDist=np.zeros(nNeig)
    for i in range(nNeig):
        neigDist[i]=1E30
    for i in range(k):
        maxDist=neigDist[0]
        maxJ=0
        for j in range(nNeig):
            if (maxDist<neigDist[j]):
                maxDist=neigDist[j]
                maxJ=j
        deltaX=math.fabs(xx-x[i])
        deltaY=math.fabs(yy-y[i])
        distance=math.sqrt(deltaX*deltaX+deltaY*deltaY)
        if(distance<=maxDist):
            neigDist[maxJ]=distance
            neigID[maxJ]=pID[i]
            neigZ[maxJ]=z[i]
    for i in range(nNeig-1):
        minDist=neigDist[i]
        minNum=i
        minID=neigID[i]
        minZ=neigZ[i]
        for i in range(nNeig):
            if(minDist>neigDist[j]):
                minDist=neigDist[j]
                minNum=j
                minID=neigID[j]
                minZ=neigZ[j]
        if(i!=minNum):
            neigDist[minNum]=neigDist[i]
            neigID[minNum]=neigID[i]
            neigZ[minNum]=neigZ[i]
            neigDist[i]=minDist
            neigID[i]=minID
            neigZ[i]=minZ
    return neigID,neigZ,neigDist

# Function to calculate the covariance - Used to elaborate the covariance matrix	
def covar_Funct(Model, Nugget, Range, Sill, Hh):
      output=0.0
      if (Model=="exponential"):
              if(Hh>0.0):
                  output=(Sill-Nugget)*math.exp(-3.0*Hh/Range)
              else:
                  output=Sill
      if (Model=="gauss"):
              if(Hh>0.0):
                  output=(Sill-Nugget)*(1.0-math.exp(-3.0*Hh*Hh/(Range*Range)))
              else:
                  output=Range
      if (Model=="spherical"):
              if(Hh>=Range):
                  output=0.0
              else:
                  if(Hh>0.0 and Hh<Range):
                      output=(Sill-Nugget)*(3*Hh/(2.0*Range)-Hh*Hh*Hh/\
                                          (2.0*Range*Range*Range))
                  else:
                      output=Sill
      return output
    
# Function to elaborate the covariance matrix
def covarianceMatrix(nNeig, Model, Nugget, Range, Sill, neigID, neigZ, neigDist, x, y):
    covM=np.zeros((nNeig+1,nNeig+1))
    xDc=np.zeros(nNeig+1)
    for j in range(0,nNeig):
        for i in range(j+1,nNeig):
            ik=int(neigID[i])
            jk=int(neigID[j])
            xx=math.sqrt((x[ik]-x[jk])*(x[ik]-x[jk])+(y[ik]-y[jk])*(y[ik]-y[jk]))
            gamma=covar_Funct(Model,Nugget,Range,Sill,xx)
            covM[i][j]=gamma
            covM[j][i]=gamma
    for i in range (0,nNeig):
        covM[i][nNeig]=1.0
        covM[nNeig][i]=1.0
        xx=0
        gamma=covar_Funct(Model,Nugget,Range,Sill,xx)
        covM[i][i]=gamma
        xx=neigDist[i]
        gamma=covar_Funct(Model,Nugget,Range,Sill,xx)
        xDc[i]=gamma
    covM[nNeig][nNeig]=0.0
    xDc[nNeig]=1.0
    return covM, xDc

# Function to execute the kriging
def kriging(xx, yy, nNeig, Model, Nugget, Range, Sill, neigID, neigZ, neigDist, x, y):
    negativeW=True
    k=0
    n=nNeig+1-k
    while((negativeW==True) and (nNeig+1-k>=4)):
        n=nNeig+1-k
        dummynNeig=nNeig-k
        covM=np.zeros((n,n))
        xDc=np.zeros(n)
        b=np.zeros(n)
        covM,xDc=covarianceMatrix(dummynNeig,Model,Nugget,Range,Sill,neigID,neigZ,neigDist,x,y)
        for i in range(0,n):
            b[i]=xDc[i]
        xD=np.linalg.solve(covM,b)
        negativeW=False
        for i in range(0,n-1):
            if(xD[i]<0):
                negativeW=True
        if(negativeW==True):
            k=k+4;
    zEstimated=0.0
    eEstimated=Sill
    for i in range(0,n-1):
        zEstimated=zEstimated+xD[i]*neigZ[i]    
        eEstimated=eEstimated-xD[i]*xDc[i]
    eEstimated=eEstimated-xD[n-1]
    if (eEstimated<0):
        eEstimated=0.0
    if(eEstimated>0.0):
        eEstimated=math.sqrt(eEstimated)
    return zEstimated, eEstimated

# Main function to perform kriging
def Main(pixel,nNeig,kk,fator_max_dist,pID,x,y,z,xlim,ylim): # 

    #Calculate the maximum distance between points

    maxdist=0.0
    n = len(z)
    for i in range (0,n):
        for j in range(i,n):
            dist=math.sqrt((x[i]-x[j])*(x[i]-x[j])+(y[i]-y[j])*(y[i]-y[j]))
            if dist>maxdist:
                maxdist=dist

    # Generate the experimental semivariogram
    gamma=np.zeros(kk)
    lag=np.zeros(kk)
    num=np.zeros(kk)
    for i in range (0,n-1):
        for j in range(i,n):
            dist=math.sqrt((x[i]-x[j])*(x[i]-x[j])+(y[i]-y[j])*(y[i]-y[j]))
            k=int(kk*dist/(fator_max_dist*maxdist))   
            if(k<=kk-1):
                gamma[k]=gamma[k]+(z[i]-z[j])*(z[i]-z[j])
                lag[k]=lag[k]+dist
                num[k]=num[k]+1
    for i in range(kk):
        gamma[i]=gamma[i]/(2*num[i])
        lag[i]=lag[i]/num[i]

    #Adjust the theoretical semivariogram model
    semiModel=np.zeros((3,4))
    # Spherical model
    Model="spherical"
    #Pick a random initial value for the nugget
    Nugget=(gamma[1]*lag[0]-gamma[0]*lag[1])/(lag[0]-lag[1])
    if Nugget<0: Nugget=gamma[0]
    #Pick a random initial value for the sill
    Sill=(gamma[kk-4]+gamma[kk-3]+gamma[kk-2]+gamma[kk-1])/4.0             
    #kick the starting value for the range
    Range=lag[int(kk/2)]                                                            
    init_vals = [Nugget, Range, Sill]
    #define the maximum values
    maxlim=[Sill,max(lag),max(gamma)]
    #optios method : ‘lm’, ‘trf’, ‘dogbox’
    best_vals, covar = curve_fit(spherical, lag, gamma,method='trf', p0=init_vals,bounds=(0, maxlim))
    Nugget=best_vals[0]
    Range=best_vals[1]
    Sill=best_vals[2]
    gammaT=np.zeros(kk)
    gammaT=spherical(lag,Nugget,Range,Sill)
    sqr_esf=0.0
    for i in range(kk):
        sqr_esf=sqr_esf+(gamma[i]-gammaT[i])*(gamma[i]-gammaT[i])
    semiModel[0][0]=Nugget
    semiModel[0][1]=Range
    semiModel[0][2]=Sill
    semiModel[0][3]=sqr_esf
    # Gauss model
    Model="gauss"
    best_vals, covar = curve_fit(spherical, lag, gamma,method='trf', p0=init_vals,bounds=(0, maxlim))   
    Nugget=best_vals[0]
    Range=best_vals[1]
    Sill=best_vals[2]
    gammaT=np.zeros(kk)
    gammaT=gauss(lag,Nugget,Range,Sill)
    sqr_gas=0.0
    for i in range(kk):
        sqr_gas=sqr_gas+(gamma[i]-gammaT[i])*(gamma[i]-gammaT[i])
    semiModel[1][0]=Nugget
    semiModel[1][1]=Range
    semiModel[1][2]=Sill
    semiModel[1][3]=sqr_gas
    # Exponential model
    Model="exponential"
    best_vals, covar = curve_fit(spherical, lag, gamma,method='trf', p0=init_vals,bounds=(0, maxlim))   
    Nugget=best_vals[0]
    Range=best_vals[1]
    Sill=best_vals[2]
    gammaT=np.zeros(kk)
    gammaT=exponential(lag,Nugget,Range,Sill)
    sqr_exp=0.0
    for i in range(kk):
        sqr_exp=sqr_exp+(gamma[i]-gammaT[i])*(gamma[i]-gammaT[i])
    semiModel[2][0]=Nugget
    semiModel[2][1]=Range
    semiModel[2][2]=Sill
    semiModel[2][3]=sqr_exp
   # check the best fit model
    Model='spherical' # 
    imodelo=0
    sqr_minimo=sqr_esf
    if(sqr_minimo>sqr_gas):
        sqr_minimo=sqr_gas
        imodelo=1
        Model='gauss'
    if(sqr_minimo>sqr_exp):
        sqr_minimo=sqr_exp
        imodelo=2
        Model='exponential'
    Nugget=semiModel[imodelo][0]
    Range=semiModel[imodelo][1]
    Sill=semiModel[imodelo][2]
    print ("The Semivariogram selected is....")
    print ("Model,Nugget,Range,Sill used in semivariogram")
    print (Model,Nugget,Range,Sill)
    print ("")

    # Generate the grid to start the kriging process
    # array to grid point
    nx=int((max(xlim)-min(xlim))/pixel)+10
    ny=int((max(ylim)-min(ylim))/pixel)+10
    ###Results
    xgrid,ygrid,zkrig,ekrig=[],[],[],[] # z estimated by kriging  error in estimation by kriging
    for i in range(nx):
        for j in range(ny):
            x1=min(xlim)+(i-5)*pixel
            y1=min(ylim)+(j-5)*pixel
            if(ray_tracing_method(x1,y1,xlim,ylim)==True): # if poit grid is inside contour
                # Find the nearest neighbor
                neigID,neigZ,neigDist=find_PNeigborns(nNeig,x1,y1,pID,x,y,z)
                #estimate the value of the variable and its error by kriging
                zEstimated,eEstimated=kriging(x1,y1,nNeig,Model,Nugget,Range,Sill,neigID,neigZ,neigDist,x,y)
                xgrid.append(x1) # e-w coordinate
                ygrid.append(y1) # n-s coordinate
                zkrig.append(zEstimated) # kriging value
                ekrig.append(eEstimated) # # kriging errors
    return lag,gamma,gammaT,xgrid,ygrid,zkrig,ekrig
