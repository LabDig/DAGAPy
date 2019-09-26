# See readme for information about this file.
#!/usr/bin/python3
#v1.1
import numpy as np
from scipy import spatial as sps

def cost(X, R, M, N, D, K, dtype):
    cost = 0
    V=np.zeros((D,D))
    Var=np.zeros(D)
    VI=np.zeros((D,D))
    if ((dtype=='diagonal') | (dtype=='mahalanobis')) :
        V=np.cov(X.transpose())
        for i in range(D):
            Var[i]=V[i][i]
    for i in range(N):
        for j in range(K):
            if(dtype=='euclidean'):
                d1=sps.distance.euclidean(M[j],X[i])
            elif(dtype=='diagonal'):
                d1=sps.distance.seuclidean(M[j],X[i],Var)
            elif(dtype=='mahalanobis'):
                V=np.cov(X.transpose())
                VI=np.linalg.inv(V)
                d1=sps.distance.mahalanobis(M[j],X[i],VI)
            cost=cost+R[i][j]*d1
    return cost

def initmember(scatter,nclass,ndata):
    U=np.zeros((ndata,nclass))
    U= np.random.rand(ndata,nclass)
    cs=np.sum(U, axis=1)
    for i in range(ndata):
        for j in range(nclass):
            U[i][j]=U[i][j]/cs[i]
    return U

def Cluster(X, K, max_iter,beta,dtype): # X is N,D array
    print ("Running Clustering...")
    N, D = X.shape
    M = np.zeros((K, D))
    U = np.zeros((N,K))
    VI = np.identity(D)
    Var = np.zeros(D)
    scatter=0.2
    U=initmember(scatter,K,N)
    UB=np.power(U,beta)
    for j in range(K):
        somaD=0.0
        somaN = np.zeros(D)
        for k in range(N):
            somaD=somaD+UB[k][j]
            for i in range(D):
                somaN[i]=somaN[i]+UB[k][j]*X[k][i]
        for i in range(D):
            M[j][i]=somaN[i]/somaD
    costs = np.zeros(max_iter)
    V=np.cov(X.transpose())
    if D==1:
        Var=V
    elif D>1:
        for i in range(D):
            Var[i]=V[i][i]
    for ii in range(max_iter):
        for i in range(N):
            for j in range(K):
                soma=0.0
                for k in range(K):
                    if(dtype=='euclidean'):
                        d1=sps.distance.euclidean(M[j],X[i])
                        d2=sps.distance.euclidean(M[k],X[i])
                    elif(dtype=='diagonal'):
                        d1=sps.distance.seuclidean(M[j],X[i],Var)
                        d2=sps.distance.seuclidean(M[k],X[i],Var)
                    elif(dtype=='mahalanobis'):
                        V=np.cov(X.transpose())
                        VI=np.linalg.inv(V)
                        d1=sps.distance.mahalanobis(M[j],X[i],VI)
                        d2=sps.distance.mahalanobis(M[k],X[i],VI)
                    soma=soma+(d1/d2)**(2.0/(beta-1.0))
                U[i][j]=1/soma
        UB=np.power(U,beta)
        for j in range(K):
            somaD=0.0
            somaN = np.zeros(D)
            for k in range(N):
                somaD=somaD+UB[k][j]
                for i in range(D):
                    somaN[i]=somaN[i]+UB[k][j]*X[k][i]
            for i in range(D):
                M[j][i]=somaN[i]/somaD
        costs[ii] = cost(X, UB, M,  N, D, K, dtype)
        if ii > 0:
            if np.abs(costs[i] - costs[i-1]) < 10e-5:
                break
        
    # array to generate patterns
    if K==2:
        pt=[[0,0,0],[1,0,0]]
    elif K==3:
        pt=[[0,0,0],[1,0,0],[0,1,0]]
    elif K==4:
        pt=[[0,0,0],[1,0,0],[0,1,0],[0,0,1]]
    elif K==5:
        pt=[[0,0,0],[1,0,0],[0,1,0],[0,0,1],[0,1,1]]

    UB_pt = U.dot(pt)
    UB=np.round_(UB_pt)
    UB.astype(int)
    

   # UB have patterns with 0,1 combination
   # Example : [1,0,0],[0,1,0],[1,1,1] =
   # Each combination is a specific classes
    l=UB.shape[0] # size
    #print (l)
    cluster=np.zeros(l) 
    # store information of each class to be sorted in crescent order
    mean_c=np.zeros((4,K)) 
    # each collum of mean_c : class 1,2,3,4,5
    # row 1 : number of point in each class
    # row 2 : sum of point atribuite
    # row 3 : mean of each class
    # row 4 : original label : 1,2,3,4,5
    # Calculate mean of each class
    #X[i,0] =z atribute use to cluster
    for i in range(l):
        if UB[i,0]==0 and UB[i,1]==0 and UB[i,2]==0:
            mean_c[0,0]=mean_c[0,0]+1
            mean_c[1,0]=mean_c[1,0]+X[i,0]
            mean_c[2,0]=mean_c[1,0]/mean_c[0,0]
            mean_c[3,0]=1.0
            cluster[i]=1.0
        elif UB[i,0]==1 and UB[i,1]==0 and UB[i,2]==0:
            mean_c[0,1]=mean_c[0,1]+1
            mean_c[1,1]=mean_c[1,1]+X[i,0]
            mean_c[2,1]=mean_c[1,1]/mean_c[0,1]
            mean_c[3,1]=2.0
            cluster[i]=2.0
        elif UB[i,0]==0 and UB[i,1]==1 and UB[i,2]==0:
            mean_c[0,2]=mean_c[0,2]+1
            mean_c[1,2]=mean_c[1,2]+X[i,0]
            mean_c[2,2]=mean_c[1,2]/mean_c[0,2]
            mean_c[3,2]=3.0
            cluster[i]=3.0
        elif UB[i,0]==0 and UB[i,1]==0 and UB[i,2]==1:
            mean_c[0,3]=mean_c[0,3]+1
            mean_c[1,3]=mean_c[1,3]+X[i,0]
            mean_c[2,3]=mean_c[1,3]/mean_c[0,3]
            mean_c[3,3]=4.0
            cluster[i]=4.0
        elif UB[i,0]==0 and UB[i,1]==1 and UB[i,2]==1:
            mean_c[0,4]=mean_c[0,4]+1
            mean_c[1,4]=mean_c[1,4]+X[i,0]
            mean_c[2,4]=mean_c[1,4]/mean_c[0,4]
            mean_c[3,4]=5.0
            cluster[i]=5.0

    # return index of sorted array in crescent order of mean classes
    sort_list=np.argsort(mean_c[2,:]) #sort_list[0] is lower mean and sort_list[5] is higher
    # define cluster atributes, with label :1,2,3,4,5,
    cluster_label=np.zeros(K)
    for i in range(K):
        cluster_label[i]=mean_c[3,sort_list[i]] # rename class label (1 is lower and K is higher)

    # Redefine label class to cluster. Label 1 is class that have lower mean and 5 higher mean
    order_cluster=np.zeros(l)  
    for i in range(l):
            if cluster[i]==1.0:
                order_cluster[i]=cluster_label[0]
            elif cluster[i]==2.0:
                order_cluster[i]=cluster_label[1]
            elif cluster[i]==3.0:
                order_cluster[i]=cluster_label[2]
            elif cluster[i]==4.0:
                order_cluster[i]=cluster_label[3]
            elif cluster[i]==5.0:
                order_cluster[i]=cluster_label[4]
    return order_cluster

