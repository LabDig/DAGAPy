# See readme for information about this file.
# import python modules
#!/usr/bin/python3
#V1.3- Modification in layout software and simplifications.

import os
import sys
import numpy as np
from PyQt5 import QtWidgets
from PyQt5.QtWidgets import QFileDialog
import utm
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import pandas as pd
from sklearn.cluster import KMeans

# import python files
from gui import Ui_dagapy #graphical user interface script
import kriging as kr #kriging script
import filtering as fil # filtering script


#
class DagaPy(QtWidgets.QTabWidget,Ui_dagapy): 
   def __init__(self,parent=None):
      super(DagaPy,self).__init__(parent)
      self.setupUi(self)

      # create res directory for results
      #In any cases,Python no make the directory.
      if not os.path.exists('res/'): os.makedirs('res/')
      #
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
      # Can be used by the user to add commands to read the sensor signals 
      #and save the data in a dataset file
      # 
      # 
      # 


       
       
   def SelectFile(self): # select file to data analysis

      fname = QFileDialog.getOpenFileName(self, 'Select Dataset File', '')
      
     
      #read the delimiter char select im combobox
      delimiter=self.delimiter.currentText()
      #
      #Read dataset and transforme to a Pandas Datraframe
      self.dataset=pd.read_csv(fname[0],sep=delimiter)
      #
      # clear combox itens on GUI
      self.lat.clear()
      self.longi.clear()
      self.atribute.clear()
      
      #Read name of columns of datrafame and add the items in combobox
      #list(self.dataset.columns) create a list of name of columns
      self.lat.addItems(list(self.dataset.columns))
      self.longi.addItems(list(self.dataset.columns))
      self.atribute.addItems(list(self.dataset.columns))

   def LimiteFile(self): # select file to limite
      fname = QFileDialog.getOpenFileName(self, 'Select Limite File', '')
      
      
      #read the delimiter char select im combobox
      delimiter=self.delimiter.currentText()
      #Read limites and transforme to a Pandas Datraframe
      self.limit=pd.read_csv(fname[0],sep=delimiter)

   def ToUtm(self,Long,Lat): #
      
      return utm.from_latlon(Lat,Long) #conv contain Xutm,Yutm,Zone Number, Zone Letter

   def ReadFile(self): # read file to data analysis (dataset and limite)

      # read in combox the item selected by user
      cx,cy,cz=self.lat.currentIndex(),self.longi.currentIndex(),self.atribute.currentIndex()
      
      
      #self.atribute_name it is used in plot's and export
      self.atribute_name=self.atribute.currentText() #it is the name of atribute column
      #
      #create x,y and z list from dataset. 
      self.x=self.dataset.iloc[:, cx]  #read colum with index cx
      self.y=self.dataset.iloc[:, cy]
      self.z=self.dataset.iloc[:, cz]
      #creat xlim and ylim list from limite
      self.xlim=self.limit.iloc[:, 0] #x it is in collum 0
      self.ylim=self.limit.iloc[:, 1] #y it is in collum 1
      #
      # if coordinate of dataset and limite is Latitude and Longitude
      if self.utm_check.isChecked() is False:
         #For dataset
         #pass Longitude and Latitude Array for ToUtm Function.  Retur tuple with Xutm,Yutm,Zone Number, Zone Letter
         conv= list(map(self.ToUtm, self.x,self.y)) 
         #conv it's a list of tuple. Transform to dataframe
         conv = pd.DataFrame(conv,columns = ['xutm' , 'ytum', 'Zone N','Zone L'],)
         #get xutm and yutm of dataframe 
         self.x=conv.iloc[:, 0] #xutm it's collum 0
         self.y=conv.iloc[:, 1] #yutm it's collum 0
         #For Limite
         #pass Longitude and Latitude Array for ToUtm Function.  Retur tuple with Xutm,Yutm,Zone Number, Zone Letter
         conv= list(map(self.ToUtm, self.xlim,self.ylim))
         #conv it's a list of tuple. Transform to dataframe
         conv = pd.DataFrame(conv,columns = ['xutm' , 'ytum', 'Zone N','Zone L'],)
         #get xutm and yutm of dataframe 
         self.xlim=conv.iloc[:, 0] #xutm it's collum 0
         self.ylim=conv.iloc[:, 1] #yutm it's collum 0
      #
      # show xy experimental point in software
      #Produce graph of the location of the experimental points and countour
      plt2=self.wdgt_main.canvas
      plt2.ax.plot(self.x,self.y, 'b*',label='Collected Point')
      plt2.ax.plot(self.xlim,self.ylim, 'r^',label='Countor Point')
      plt2.ax.set_ylabel("Coordinate UTM N-S, m")
      plt2.ax.set_xlabel("Coordinate UTM E-W, m")
      plt2.ax.axis('equal') #keep ratio aspect
      plt2.ax.set_yticklabels([]) #hide tick of label
      plt2.ax.set_xticklabels([]) #hide tick of label
      plt2.ax.legend()
      plt2.draw()
      #
      # show histogram in tab filtering
      plt2=self.wdgt_hist_before.canvas
      plt2.ax.hist(self.z)
      plt2.ax.set_title('Raw Data Histrogram')
      plt2.ax.set_xlabel(self.atribute_name) #xlabel it is the name of atribute column
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
      plt2.ax.set_xlabel(self.atribute_name) #xlabel it is the name of atribute column
      plt2.ax.set_ylabel("Frequency")
      plt2.draw()
      #
      #Exports dataset filtered to csv file
      data={'xutm':self.x,'yutm':self.y,self.atribute_name:self.z} #build a dictonaty
      data=pd.DataFrame(data) #dictonary to dataframe
      data.to_csv('res/filtering.csv',index=False) #export


   def RunKrig(self): # Function that run kriging 
       

      # Prepare the dataset to kriging
      xyz=pd.DataFrame({'x':self.x,'y':self.y,'z':self.z})
      
      xylim=pd.DataFrame({'x':self.xlim,'y':self.ylim})
       
      #Read configuration in Tab
      grid=self.gridsize.value() # read in spin box, the grid size
      nneig=self.neig.value() # neig number
      nlag=self.lag.value() # lag number
      fator_max_dist=self.max_dist.value() # max dist factor to semivariogram

      #Call function to generate the semivariogram 
      lag,gamma,gammaT,model=kr.SemiVariogram(nlag,fator_max_dist,xyz)
      #call function to generathe the grid
      #In future, a xygrid file can be upload
      xygrid=kr.Grid(grid,xylim)
      #Call function to kriging
      zkrig,ekrig=kr.Kriging(model,nneig,xygrid,xyz)
         
      # local to global
      self.lag=lag
      self.gamma=gamma
      self.gammaT=gammaT
      self.xgrid=xygrid[:,0]
      self.ygrid=xygrid[:,1]
      self.zkrig=zkrig
      self.ekrig=ekrig
      # Call function that plot Semivariogram.
      self.ShowSemi()
      #Exports kriging points to csv file
      data={'xgrid':self.xgrid,'ygrid':self.ygrid,'estimated value':self.zkrig,\
             'standart deviation of estimation':self.ekrig} #build a dictonaty
      data=pd.DataFrame(data) #dictonary to dataframe
      data.to_csv('res/kriging.csv',index=False) #export

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
      plt2.ax.axis('auto') #keep ratio aspect
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
      plt2.ax.set_yticklabels([]) #hide tick of label
      plt2.ax.set_xticklabels([]) #hide tick of label
      plt2.ax.axis('equal') #keep ratio aspect
      plt2.draw()

   def KrigValue(self): 

      plt2=self.wdgt_krig.canvas
      plt2.ax.cla()
      try:self.cb.remove()
      except:pass
      ax2=plt2.ax.scatter(self.xgrid,self.ygrid, cmap=cm.jet,c=self.zkrig,marker='s',s=75)
      self.cb=plt2.fig.colorbar(ax2) #global colobar self.cb
      self.cb.set_label(label=self.atribute_name)
      plt2.ax.set_ylabel("Coordinate UTM N-S, m")
      plt2.ax.set_xlabel("Coordinate UTM E-W, m")
      plt2.ax.set_yticklabels([]) #hide tick of label
      plt2.ax.set_xticklabels([]) #hide tick of label
      plt2.ax.axis('equal') #keep ratio aspect
      plt2.draw()
     

   def KrigError(self):

      plt2=self.wdgt_krig.canvas
      plt2.ax.cla()
      try:self.cb.remove()
      except:pass
      ax2=plt2.ax.scatter(self.xgrid,self.ygrid, cmap=cm.jet,c=self.ekrig,marker='s',s=75)
      self.cb=plt2.fig.colorbar(ax2) #global colobar self.cb
      self.cb.set_label(label=self.atribute_name)
      plt2.ax.set_ylabel("Coordinate UTM N-S, m")
      plt2.ax.set_xlabel("Coordinate UTM E-W, m")
      plt2.ax.set_yticklabels([]) #hide tick of label
      plt2.ax.set_xticklabels([]) #hide tick of label
      plt2.ax.axis('equal') #keep ratio aspect
      plt2.draw()

   def RunCluster(self):
 
      k=self.n_clust.value()# read number of cluster
      #Running Kmeans
      z=np.asarray(self.zkrig)
      z=z.reshape(-1,1) #necessary to used kmeans for one atribute
      kmeans = KMeans(n_clusters=k)
      z_cluster = kmeans.fit_predict(z)
      #
      # obtain the class of each point
      #plot cluster graphics
      plt2=self.wdgt_cluster.canvas
      plt2.ax.cla() #clear the graphics
      #remove colorbar, if exist
      try:self.cbar.remove()
      except:pass
      #plot the map,with k color
      ax2=plt2.ax.scatter(self.xgrid,self.ygrid, cmap=plt.get_cmap('jet', k),c=z_cluster,marker='s',s=75)
      #ticks for colorvar
      tick=np.arange(k)
      self.cbar=plt2.fig.colorbar(ax2,ticks=tick)
      #define the tick label, acoording number of cluster
      if k==2: self.cbar.ax.set_yticklabels(['Low', 'High'])
      elif k==3: self.cbar.ax.set_yticklabels(['Low', 'Medium', 'High'])
      elif k==4: self.cbar.ax.set_yticklabels(['Very Low','Low', 'High','Very high'])
      elif k==5: self.cbar.ax.set_yticklabels(['Very Low','Low', 'Medium','High','Very high'])
      #plot configurations
      plt2.ax.set_yticklabels([]) #hide tick of label
      plt2.ax.set_xticklabels([]) #hide tick of label
      plt2.ax.set_ylabel("Coordinate UTM N-S, m")
      plt2.ax.set_xlabel("Coordinate UTM E-W, m")
      self.cbar.set_label(label=self.atribute_name)
      plt2.ax.axis('equal') #keep ratio aspect
      plt2.draw()
      #Exports kriging points to csv file
      data={'xgrid':self.xgrid,'ygrid':self.ygrid,'estimated value':self.zkrig,\
             'class':z_cluster} #build a dictonaty
      data=pd.DataFrame(data) #dictonary to dataframe
      data.to_csv('res/clustering.csv',index=False) #export

#Run the app:
if __name__ == '__main__':
   if not QtWidgets.QApplication.instance():
      app = QtWidgets.QApplication(sys.argv)
   else:
      app = QtWidgets.QApplication.instance() 
   ex = DagaPy()
   ex.show()
   sys.exit(app.exec_())

