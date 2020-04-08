# See readme for information about this file.
'''
Modified version, for reduce processing time

'''
#v1.1
#!/usr/bin/python3

import numpy as np
from scipy.optimize import curve_fit
import math
import time
import pandas as pd
import itertools as it
from scipy import spatial
import matplotlib.path as mplPath
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


# Function to calculate the covariance - Used to elaborate the covariance matrix	
def covar_Funct(Model,Hh):
        
      func=Model[0]
      Nugget=Model[1]
      Range=Model[2]
      Sill=Model[3]
      
      output=0.0
      if (func=="Exponential"):
              if(Hh>0.0):
                  output=(Sill-Nugget)*math.exp(-3.0*Hh/Range)
              else:
                  output=Sill
      if (func=="Gauss"):
              if(Hh>0.0):
                  output=(Sill-Nugget)*(1.0-math.exp(-3.0*Hh*Hh/(Range*Range)))
              else:
                  output=Range
      if (func=="Spherical"):
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
def covarianceMatrix(n,Model, DistP2N, DistN2N):
    
    #n is size of matrix
    #model its type, nugget ,range, sill
    #neigZ it is z value of neigDist
    #DistP2N it is distance of neigDist to grid point
    #DistN2N it is distance of neighboor to neighboor
    #covM it is the covariance matriz between the  neighboor
    C=np.ones((n,n))
    C[n-1][n-1]=0
    #xDc it is the covariance array between and the grid point
    D=np.ones(n)
    #Build the covM matrix
    # The covariance matrix it is estimaded using the fit model, based on distance
    # DistP2N it is n-1 x n-1
    
    
    for i in range(0,n-1):
        for j in range(0,n-1):
            C[i,j]=covar_Funct(Model,DistN2N[i][j])

    for i in range(0,n-1):
        D[i]=covar_Funct(Model,DistP2N[i])

    #Return
    return C, D

# Function to execute the kriging
#Receive the model (type, nugget,range, sill), the distance of neighboor to grid point,
#distance of neighboor to neighboor
# the z value of each neighboor
def exec_kriging(Model, DistP2N, DistN2N, neigZ):
    
    
    #number of neighboor
    nNeig=len(neigZ)
    #This software do not use negative Weight
    negativeW=True
    k=0
    n=nNeig+1-k #matrix n+1 x n+1
    
    while((negativeW==True) and (nNeig+1-k>=4)):
        n=nNeig+1-k
        #covM it is the covariance matriz between the  neighboors
        #xDc it is the covariance array between and the grid point
        #Call function to buid the matrix
        #pass the number size of matriz , the model (type, nugget,range, sill),
        #the distance p2n , n2n and 
        covM,xDc=covarianceMatrix(n,Model, DistP2N, DistN2N)
        # Solve xD=covM-1 *xDC . 
        #xD= weight array and the Lagrange
        xD=np.linalg.solve(covM,xDc)
        # Check if have negative Weight
        negativeW=False
        for i in range(0,n-1):
            if(xD[i]<0): negativeW=True
        if(negativeW==True):
            k=k+4;
            #print ('Was found negative Weight. Reduce in 4 (for) the number of neig neighboors')
    
    # Calculate z and est, acoording GeoStatistics Books
    zEstimated=0.0
    eEstimated=covM[0][0] #sill (COV[i,i]) it is used 
    #
    for i in range(0,n-1):
        zEstimated=zEstimated+xD[i]*neigZ[i]    
        eEstimated=eEstimated-xD[i]*xDc[i]
        
    eEstimated=eEstimated-xD[n-1] #the last element its the Lagrange
    #
    if (eEstimated<0):  eEstimated=0.0
    # This software use standart deviation of estimation
    if(eEstimated>0.0):  eEstimated=math.sqrt(eEstimated)
        
    #Return estimation and standart deviation of estimation
    return zEstimated, eEstimated

# Function to Generate Semivariogram
def SemiVariogram(nlag,fator_max_dist,xyz): # 
    
    #nlag it is the number of lags
    #fator_max_dist it is the active distance used in semivariogram
    #xyz it is x,y and z experimental points
    #Build a dataframe to save experimental semivariogram
    var=pd.DataFrame()
    #distance between all experimental points
    var['lag']=spatial.distance.pdist(xyz.iloc[:,0:2], metric='euclidean')
    #distance between z(i) e z(i+h) ==> (z(i)-z(i+h))**2
    var['gamma']=[(y - x)**2 for x, y in it.combinations(xyz.iloc[:,2], 2)]
    
    #sort in crescent order acoording distance (lag)
    #gamma ij follow the order
    var=var.sort_values(by='lag')
    #remove point if distance > max_dist * factor
    remove_index=var[var['lag'] > fator_max_dist*max(var['lag']) ].index
    var=var.drop(remove_index)
    
    #ranges to split lags
    #max(var['lag'] it is new max distance after remove points
    bins=np.arange(0,max(var['lag']),max(var['lag'])/nlag)
    
    #classify the distances according the range
    ind = np.digitize(var['lag'],bins)
    
    #group lag e and calculate mean
    lag=var['lag'].groupby(ind).mean()
    lag=lag.to_numpy() #convert pandas series to numpy
    #group gamma e calculate mean()/2, acoording matheron estimator
    gamma=var['gamma'].groupby(ind).mean().div(2)
    gamma=gamma.to_numpy() #convert pandas series to numpy
    
    #Adjust the theoretical semivariogram model
   
    #Pick a random initial value for the nugget
    Nugget=(gamma[1]*lag[0]-gamma[0]*lag[1])/(lag[0]-lag[1])
    if Nugget<0: Nugget=gamma[0]
    #Pick a random initial value for the sill
    Sill=(gamma[nlag-4]+gamma[nlag-3]+gamma[nlag-2]+gamma[nlag-1])/4.0             
    #kick the starting value for the range
    Range=lag[int(nlag/2)]                                                            
    init_vals = [Nugget, Range, Sill]
    #define the maximum values
    maxlim=[Sill,max(lag),max(gamma)]
    
    ###############################################################
    #Fit Spherical Models
    #return Nugget, Range , Sill and estimated covariance (not used)
    [Nugget,Range,Sill] , _ = curve_fit(spherical, lag, gamma,method='trf', p0=init_vals,bounds=(0, maxlim))

    
    gammaT=np.zeros(nlag)
    gammaT=spherical(lag,Nugget,Range,Sill)
    #Calculate RMSE of fit
    rmse=0.0
    for i in range(nlag): rmse=rmse+(gamma[i]-gammaT[i])*(gamma[i]-gammaT[i])

    # List of results fit of spherical model
    model=['Spherical',Nugget,Range,Sill,rmse]
    
    print ('Model:',model[0],' Nugget:',model[1],' Range:',model[2],\
           ' Sill:',model[3],' RMSE:',model[4],'\n')
    
    ###############################################################
    #Fit Gauus Models
    [Nugget,Range,Sill] , _  = curve_fit(gauss, lag, gamma,method='trf', p0=init_vals,bounds=(0, maxlim))   

    gammaT=np.zeros(nlag)
    gammaT=gauss(lag,Nugget,Range,Sill)
    #Calculate RMSE of fit
    rmse_2=0.0
    for i in range(nlag):
        rmse_2=rmse_2+(gamma[i]-gammaT[i])*(gamma[i]-gammaT[i])
    
    #If RMSE_2 its smaller than RMSE (model[4] == >  Gauss Fit Model it is better ==> update model list
    # Else Spherical Fit Model it is better
    if rmse_2 < model[4] : 
        model=['Gauss',Nugget,Range,Sill,rmse_2]
        print ('Gauss Fit Model it is better than Spherical Fit Model\n')
        print ('Model:',model[0],' Nugget:',model[1],' Range:',model[2],\
           ' Sill:',model[3],' RMSE:',model[4],'\n')
    
    else:
         print ( model[0], ' Fit Model still it is better \n')
    
    ###############################################################
    #Fit Exponentil Models

    [Nugget,Range,Sill], _ = curve_fit(exponential, lag, gamma,method='trf', p0=init_vals,bounds=(0, maxlim))   
 
    gammaT=np.zeros(nlag)
    gammaT=exponential(lag,Nugget,Range,Sill)
    #Calculate RMSE of fit
    rmse_3=0.0
    for i in range(nlag):
        rmse_3=rmse_3+(gamma[i]-gammaT[i])*(gamma[i]-gammaT[i])
    
    #If RMSE_3 its smaller than RMSE (model[4]) == >  Exponential Fit Model it is better ==> update model list
    # Else Spherical or Gauss Fit Model it is better
    if rmse_3 < model[4] : 
        model=['Exponential',Nugget,Range,Sill,rmse_3]
        print ('Exponential Fit Model it is better than Spherical Fit Model\n')
        print ('Model:',model[0],' Nugget:',model[1],' Range:',model[2],\
           ' Sill:',model[3],' RMSE:',model[4],'\n')
    
    else:
         print ( model[0], ' Fit Model still it is better \n')
    
    return lag,gamma,gammaT,model[0:4] # the RMSE it is not necessary

###############################################################
#Function to generate the grid,
def Grid (grid,xylim):
    
    #grid is the size of grid, in meter
    #xylim it is the pair of points that define the contour
    # Generate the grid to start the kriging process
    #Find min and max of contour
    x_min = xylim.iloc[:,0].min()
    x_max = xylim.iloc[:,0].max()

    y_min = xylim.iloc[:,1].min()
    y_max = xylim.iloc[:,1].max()
    
    #Generate the grid
    gridx = np.arange(x_min, x_max, grid)
    gridy = np.arange(y_min, y_max, grid)
    
   
    # Generate a polygon with contour
    contours = mplPath.Path(np.array(xylim))
    
    # Generate a array with all grid points inside contour
    #gridxy it is a n x 2 , where n its the number o points
    gridxy=[]
    #
    # Run all combination of i and j points of grid
    for i in gridx:
        for j in gridy:
            # if point it is internal of contour
            if contours.contains_point((i,j)):
                gridxy.append([i,j])
    
    #return a nx2 array with pair of point internal of contour
    return np.array(gridxy)
    

#Function to execute kriging
def Kriging (model,nneig,gridxy,xyz):
   
    #model contain type (shperical, exponential ou gauss), nugget, range and sil of
    #fit semivariogram
    #nneig it is the number of neigbor for interpolation
    #gridxy it is the pair of point for interpolation
    #xyz it is x, t and z of experimental points
    #
    #Build the tree to find neigboor
    tree = spatial.KDTree(xyz.iloc[:,0:2])
    
    zkrig,ekrig=[],[] # z estimated by kriging  error in estimation by kriging
    
    # Run all combination of i and j points of grid
    for i in gridxy:
   
        #find distance of neighboor to grid poitn
        #and index of k neigboor
        dist_p2n,index=tree.query((i[0],i[1]),k=nneig)
        #pass the [model,nugget,range,sill,
        #distance of neigbhoor to grid point
        # distance of each neigbhoor with other neigbhoor
        #and z value of neigbhoor
        z_neig=xyz.iloc[index,2]
            #find distance between neighboor
        #find distance between neighboor
        dist_n2n=spatial.distance.pdist(xyz.iloc[index,0:2], metric='euclidean')
        # a square matrix n x n, with distance between the neighboor
        dist_n2n=spatial.distance.squareform(dist_n2n)
        z_est,e_est=exec_kriging(model,dist_p2n,dist_n2n,z_neig.to_numpy())
        zkrig.append(z_est)
        ekrig.append(e_est)
    
    #return it is estimation value and standart deviation of estimation
    return np.array(zkrig),ekrig



