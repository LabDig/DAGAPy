# See readme for information about this file.
# import python modules
#!/usr/bin/python3
#V1.2- Modification in layout software and simplifications.

import os
import sys
import numpy as np
from PyQt5 import QtWidgets
from PyQt5.QtWidgets import QFileDialog
from PyQt5 import QtGui
import utm
import matplotlib.pyplot as plt
import matplotlib.cm as cm

# import python files
from guiv12 import Ui_dagapy #graphical user interface script
import kriging as kr #kriging script
import clustering as cl # clustering script
import filtering as fil # filtering script


#
class DagaPy(QtWidgets.QTabWidget,Ui_dagapy): 
    def __init__(self,parent=None):
        super(DagaPy,self).__init__(parent)
        self.setupUi(self)
        self.this_dir=os.path.dirname(os.path.abspath(__file__)) # directory of this file (main.py)

        # create res directory for results 
        if not os.path.exists('res/'):
            os.makedirs('res/')

        ###In tab Main
        #button that start data acquisition
        self.start_read.clicked.connect(self.StartAcquisition)
        #button that select file to data analysis
        self.select_file.clicked.connect(self.SelectFile)
        #button that select file that contain limit points
        self.select_limit_file.clicked.connect(self.LimiteFile)
        #button that read the dataset file
        self.read_file.clicked.connect(self.ReadFile)

        ### In tab Filtering
        # run filter button
        self.fil_run.clicked.connect(self.Filter) 

        ### In tab Kriging        
        # button that run kriging
        self.krig_run.clicked.connect(self.RunKrig)
        # button that plot and show semivariogram
        self.show_semiv.clicked.connect(self.ShowSemi)
        # button that show grid used in kriging
        self.show_grid.clicked.connect(self.Grid)
        # button that show kriging estimated value
        self.show_krig_value.clicked.connect(self.KrigValue)
        # button that show kriging estimated errors
        self.show_krig_error.clicked.connect(self.KrigError)

        ### In Clustering
        self.run_cluster.clicked.connect(self.RunCluster) # Run Cluster

    #Functions
    def StartAcquisition(self): # data acquistion
        #
        pass
        #
        #
        #
        # Can be used by the user to add commands to read the sensor signals and save the data in a dataset file
        # 
        # 
        # 

    def SelectFile(self): # select file to data analysis

        fname = QFileDialog.getOpenFileName(self, 'Select Dataset File', '')

        #read the delimiter char select im combobox
        delimiter=self.delimiter.currentText()
        #
        #self.dataset contain the dataset im matrix data and it is used in ReadFile function
        self.dataset=None # clear variable for security
        # read dataset file
        with open(fname[0], "r",encoding='latin-1') as f:
            self.dataset = f.read().splitlines()
        f.close
        # read header of dataset file   
        ctitle=self.dataset[0].split(delimiter)
        self.nlen=len(self.dataset)
        # clear combox itens
        self.lat.clear()
        self.longi.clear()
        self.atribute.clear()
        # add options of variables in combox
        self.lat.addItems(ctitle) 
        self.longi.addItems(ctitle)
        self.atribute.addItems(ctitle)

    def LimiteFile(self): #
        fname = QFileDialog.getOpenFileName(self, 'Select Limite File', '')
        #read the delimiter char select im combobox
        delimiter=self.delimiter.currentText()
        #self.limite contain the limit points in matrix, and it is used in Read File Function
        self.limite=None # clear variable for security
        # read dataset file
        with open(fname[0], "r",encoding='latin-1') as f:
            self.limite = f.read().splitlines()
        f.close

        
    def ReadFile(self): # read file to data analysis (dataset and limite)

        #Read dataset and organize in vector X , Y and Z
        nlen=int(len(self.dataset)) 
        # read in combox the item selected by user
        cx,cy,cz=self.lat.currentIndex(),self.longi.currentIndex(),self.atribute.currentIndex()
        #read the delimiter char select im combobox
        delimiter=self.delimiter.currentText()
        self.x,self.y,self.z=[],[],[] # temporary x and y
        # Separate the data column
        for i in range(1,nlen-1): #remove header
            Row=self.dataset[i].split(delimiter)
            if self.utm_check.isChecked(): # alread in UTM
                self.x.append(float(Row[cx])) # separe x data in a vector
                self.y.append(float(Row[cy])) # separe y data in a vector
            else: #it is necessary to convert to utm
                utm_conv=utm.from_latlon(float(Row[cx]),float(Row[cy]))
                self.x=utm_conv[0]
                self.y=utm_conv[1]
            self.z.append(float(Row[cz])) # separe z data in a vector
       
        #Read limite and organize in vector Xlim and Ylim
        nlen=int(len(self.limite))
        self.xlim,self.ylim=[],[] #
        for i in range(1,nlen-1): #remove header
            Row=self.limite[i].split(delimiter)
            if self.utm_check.isChecked(): # alread in UTM
                self.xlim.append(float(Row[0])) #First Column
                self.ylim.append(float(Row[1])) #First Column
            else: #it is necessary to convert to utm
                utm_conv=utm.from_latlon(float(Row[0]),float(Row[1]))
                self.xlim=utm_conv[0]
                self.ylim=utm_conv[1]
        
        # show xy experimental point in software
        #Produce graph of the location of the experimental points and countour
        plt2=self.wdgt_main.canvas
        plt2.ax.plot(self.x,self.y, 'b*',label='Collected Point')
        plt2.ax.plot(self.xlim,self.ylim, 'r^',label='Countor Point')
        plt2.ax.set_ylabel("Coordinate UTM N-S, m")
        plt2.ax.set_xlabel("Coordinate UTM E-W, m")
        plt2.ax.legend()
        plt2.draw()

        # show histogram
        plt2=self.wdgt_hist_before.canvas
        plt2.ax.hist(self.z)
        plt2.ax.set_title('Raw Data Histrogram')
        plt2.ax.set_xlabel("Atribute")
        plt2.ax.set_ylabel("Frequency")
        plt2.draw()



        
    def Filter(self): # Function that calculate  filtering => call function  in filering.py script

        
        # read spinbox min and maximum value
        minValue=float(self.minVal.value())
        maxValue=float(self.maxVal.value())
        x_fil,y_fil,z_fil=fil.PreFilter(minValue,maxValue,self.x,self.y,self.z)
        # update the arrays
        self.x,self.y,self.z=x_fil,y_fil,z_fil
        # Outlier filtering
        n_dist=self.neig_dist.value() # maximum distance to search neigbhor
        sd_factor=float(self.sd_ratio.value()) # standat deviation factor
        median_factor=float(self.median_ratio.value()) # median factor
        x_fil,y_fil,z_fil=fil.Outlier(n_dist,sd_factor,median_factor,self.x,self.y,self.z)
        # update variables with filtered data
        self.x,self.y,self.z=x_fil,y_fil,z_fil
        # Inlier Filtering
        x_fil,y_fil,z_fil=fil.Inlier(n_dist,median_factor,self.x,self.y,self.z)
        # update variables with filtered data
        self.x,self.y,self.z=x_fil,y_fil,z_fil
       # Descriptive statistic of dataset
       # show histogram
        plt2=self.wdgt_hist_filter.canvas
        plt2.ax.hist(self.z)
        plt2.ax.set_title('Raw Data Histrogram')
        plt2.ax.set_xlabel("Atribute")
        plt2.ax.set_ylabel("Frequency")
        plt2.draw()
        # Automatically exports dataset filtered to txt file
        #coordinate utm e-w, coordinate utm n-s, estimated value, estimated error
        file_name=os.path.join(self.this_dir,"res/filtering.txt")
        res=open(file_name,"w")
        res.write("coordinate utm e-w, coordinate utm n-s, Variable after filter \n")
        res.close()
        res=open(file_name,"a")
        for i in range(len(self.x)):
            res.write(str(self.x[i])+','+str(self.y[i])+','+str(self.z[i])+'\n')
        res.close()
        
    def RunKrig(self): # Function that run kriging 

        grid=self.gridsize.value() # read in spin box, the grid size
        nneig=self.neig.value() # neig number
        nlag=self.lag.value() # lag number
        max_dist_factor=self.max_dist.value() # max dist factor to semivariogram
        pID=[]
        for i in range (len(self.x)):
            pID.append(i)
        self.lag,self.gamma,self.gammaT,self.xgrid,self.ygrid,self.zkrig,self.ekrig=kr.Main(grid,nneig,nlag,max_dist_factor,pID,self.x,self.y,self.z,\
                                                                                        self.xlim,self.ylim)
        # Call function that plot Semivariogram.
        self.ShowSemi()
        #coordinate utm e-w, coordinate utm n-s, estimated value, estimated error
        this_dir=os.path.dirname(os.path.abspath(__file__)) # directory of this ptyhon file
        file_name=os.path.join(this_dir,"res/kriging_results.txt")
        res=open(file_name,"w")
        res.write("id,coordinate utm e-w, coordinate utm n-s, estimated value, estimated error \n")
        res.close()
        res=open(file_name,"a")
        for i in range(len(self.xgrid)):
            res.write(str(self.xgrid[i])+','+str(self.ygrid[i])+','+str(self.zkrig[i])+','+str(self.ekrig[i])+'\n')
        res.close()
        
    #Function to plot SemiVariogram
    def ShowSemi(self): 

        try:self.cb.remove() #remove colorbar figure
        except:pass
        plt2=self.wdgt_krig.canvas
        plt2.ax.cla()
        plt2.ax.plot(self.lag,self.gamma,'b+',label='Experimental')
        plt2.ax.plot(self.lag,self.gammaT,'b',label='Fit')
        plt2.ax.set_ylabel("Semivariance")
        plt2.ax.set_xlabel("Distance (m)")
        plt2.ax.legend()
        plt2.draw()


    #Function to plot grid point
    def Grid(self):

        try:self.cb.remove() #remove colorbar figure
        except:pass
        plt2=self.wdgt_krig.canvas
        plt2.ax.cla()
        plt2.ax.plot(self.xgrid,self.ygrid, 'b*',label='Grid Point')
        plt2.ax.plot(self.xlim,self.ylim, 'r^',label='Countor Point')
        plt2.ax.set_ylabel("Coordinate UTM N-S, m")
        plt2.ax.set_xlabel("Coordinate UTM E-W, m")
        plt2.draw()

    def KrigValue(self): 

        plt2=self.wdgt_krig.canvas
        plt2.ax.cla()
        try:self.cb.remove()
        except:pass
        ax2=plt2.ax.scatter(self.xgrid,self.ygrid, cmap=cm.jet,c=self.zkrig,marker='s',s=75)
        self.cb=plt2.fig.colorbar(ax2) #global colobar self.cb
        plt2.ax.set_ylabel("Coordinate UTM N-S, m")
        plt2.ax.set_xlabel("Coordinate UTM E-W, m")
        plt2.draw()

    def KrigError(self):

        plt2=self.wdgt_krig.canvas
        plt2.ax.cla()
        try:self.cb.remove()
        except:pass
        ax2=plt2.ax.scatter(self.xgrid,self.ygrid, cmap=cm.jet,c=self.ekrig,marker='s',s=75)
        self.cb=plt2.fig.colorbar(ax2) #global colobar self.cb
        plt2.ax.set_ylabel("Coordinate UTM N-S, m")
        plt2.ax.set_xlabel("Coordinate UTM E-W, m")
        plt2.draw()
        
    def RunCluster(self):

        k=self.n_clust.value()# read number of cluster
        max_iter=self.max_int.value() # max iteration
        beta=float(self.max_int.value()) # fuzzy exponent
        dtype="euclidean"
        n=len(self.xgrid)
        D=1 # only one atribute it is used
        Z=np.zeros((n,D))
        for i in range(n):
            Z[i][0]=self.zkrig[i]  # transforme to N x 1 vector
        # Run cluster    
        self.cluster=cl.Cluster(Z, k, max_iter, beta, dtype) # Z is N,D array. In this case,D=1
        
        #plot cluster graphics
        plt2=self.wdgt_cluster.canvas
        plt2.ax.cla() #clear the graphics
        try:self.cbar.remove() #remove tha actual colorbar, if existent
        except:pass
        cmap = plt.get_cmap('jet', k)
        ax2=plt2.ax.scatter(self.xgrid,self.ygrid, cmap=cm.jet,c=self.cluster,marker='s',s=75)
        
        #add ticks
        tick=[]
        for i in range(k):
            tick.append(i+1)
        cbar=plt2.fig.colorbar(ax2,ticks=tick) #global colobar self.cb
        #add tick labels
        if k==2:
             cbar.ax.set_yticklabels(['Low', 'High'])
        elif k==3:
            cbar.ax.set_yticklabels(['Low', 'Medium', 'High'])
        elif k==4:
            cbar.ax.set_yticklabels(['Very Low','Low', 'High','Very high'])
        elif k==5:
            cbar.ax.set_yticklabels(['Very Low','Low', 'Medium','High','Very high'])
        
        plt2.ax.set_ylabel("Coordinate UTM N-S, m")
        plt2.ax.set_xlabel("Coordinate UTM E-W, m")
        plt2.draw()

        #export cluster results
   
        file_name=os.path.join(self.this_dir,"res/clustering_results.txt")
        res=open(file_name,"w")
        res.write("coordinate utm e-w, coordinate utm n-s, cluster \n")
        res.close()
        res=open(file_name,"a")
        for i in range(len(self.xgrid)):
            res.write(str(self.xgrid[i])+','+str(self.ygrid[i])+','+str(self.cluster[i])+'\n')
        res.close()
       
        

#Run the app:
if __name__ == '__main__':
    if not QtWidgets.QApplication.instance():
        app = QtWidgets.QApplication(sys.argv)
    else:
        app = QtWidgets.QApplication.instance() 
    ex = DagaPy()
    ex.show()
    sys.exit(app.exec_())
    
