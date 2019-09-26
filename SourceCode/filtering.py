# See readme for information about this file.
#!/usr/bin/python3
#v1.1
import pandas as pd
import math


# Function to find the nearest neighbor
#DIst = Max distance, xx and yy  = point to find neighbor , pId,x,y,z - dataset used
# return neid ID,neigz,neid Dist
def find_PNeigborns(Dist, xx, yy, x, y, z):
    k=len(x)
    neigDist=[]
    neigZ=[]
    for i in range(k):
        deltaX=math.fabs(xx-x[i])
        deltaY=math.fabs(yy-y[i])
        distance=math.sqrt(deltaX*deltaX+deltaY*deltaY)
        if (distance<=Dist):
            neigDist.append(distance)
            neigZ.append(z[i])
    return neigZ,neigDist


# remove data outside of admissible range
def PreFilter(minV,maxV,x,y,z):
    print ("Running pre filtering...")
    n=len(x)
    x_Prefil=[]
    y_Prefil=[]
    z_Prefil=[]
    z_del=[]
    y_del=[]
    x_del=[]
    for i in range(n):
        if ((z[i]>minV) and (z[i]<maxV)): # if is inside of the range defineted
            z_Prefil.append(z[i]) # 
            y_Prefil.append(y[i]) # 
            x_Prefil.append(x[i]) # 
        else:
            z_del.append(z[i])
            y_del.append(y[i])
            x_del.append(x[i])
    print ('Prefilter removed datas : before and after',len(x),len(x_del))
    return x_Prefil,y_Prefil,z_Prefil



# Function to outiler filtering 
def Outlier(Dist,out_factor,local_factor,x,y,z):
    print ("Running outiler filtering...")
    n=len(x)
    # define sd_factor
    # if mean - sd_factor * sd < x[i] < meand + sd_factor * sd => the data is maintainned
    # else , the data is deleted
    # Calculates Statistics for z data... The results are in n[0],mean[0],sd[0],.....
    data=pd.DataFrame(z)
    mean=data.mean()
    sd=data.std()
    # list to filtering data
    z_Outfil=[]
    y_Outfil=[]
    x_Outfil=[]
    z_del=[]
    y_del=[]
    x_del=[]
    #The rule for data removal
    for i in range (n):
         if (z[i] > (mean[0] - out_factor*sd[0])) and (z[i] < (mean[0] + out_factor*sd[0])): # if is inside the range
           # Execute an Oultier local filtering. This is, check if the point is too a local outilier
            xx=x[i]
            yy=y[i]
            neigZ,neidDist=find_PNeigborns(Dist, xx, yy, x, y, z)
            # calc statistic for point and its neighbors
            data=pd.DataFrame(neigZ) # 
            med=data.median()
            ###### The rule for data removal
            if z[i]*(1-local_factor)<med[0] and z[i]*(1+local_factor)>med[0]: # if z*(1-factor) < median of neighobor points  < z*(1+factor) => the point is true. 
                z_Outfil.append(z[i])
                x_Outfil.append(x[i])
                y_Outfil.append(y[i])
            else:
                z_del.append(z[i])
                y_del.append(y[i])
                x_del.append(x[i])
   # return result
    print ('Outiler removed data : before and after',len(x),len(x_del))
    return x_Outfil,y_Outfil,z_Outfil


# Function to inlier filtering (local filtering analysis)
def Inlier(Dist,in_factor,x,y,z):
    print ("Running inlier filtering...")
    n=len(x)
    # number of neigbhor to use
    z_Infil=[]
    y_Infil=[]
    x_Infil=[]
    z_del=[]
    y_del=[]
    x_del=[]
    # For the point i
    for i in range (n):
        xx=x[i]
        yy=y[i]
        # call function that find neighbour
        #nNeig = Number of neighbor, xx and yy  = point to find neighbor , pId,x,y,z - dataset used
        # return neid ID,neigz,neid Dist
        neigZ,neidDist=find_PNeigborns(Dist, xx, yy, x, y, z)
       # calc statistic for point yours neighbor
        data=pd.DataFrame(neigZ) # 
        med=data.median()
       ###### The rule for data removal
        Id=0 # to start in 0
        if z[i]*(1-in_factor)<med[0] and z[i]*(1+in_factor)>med[0]: # if z*(1-factor) < median of neighbor points  < z*(1+factor) => the point is true. 
            Id=Id+1
            z_Infil.append(z[i])
            x_Infil.append(x[i])
            y_Infil.append(y[i])
        else:
            z_del.append(z[i])
            y_del.append(y[i])
            x_del.append(x[i])
    print ('Inlier removed data : before and after',len(x),len(x_del))
    return x_Infil,y_Infil,z_Infil
