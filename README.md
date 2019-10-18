# DAGAPy
Data Acquistion and Geostatistic Analysis Python Software

Software developed and published in Journal Computer and Eletronics in Agricultural.

Manuscript title : An open-source spatial analysis system for embedded systems


    DAGAPy is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    DAGAPy is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with DAGAPy.  If not, see <https://www.gnu.org/licenses/>6.
    
**************************************************************
Authors of software

M.Sc Andre Luiz de Freitas Coelho  (Federal University of Vicosa).Email:andre.coelho@ufv.br

PhD Daniel Marcal de Queiroz (Federal University of Vicosa). Email:queiroz@ufv.br

PhD Francisco de Assis de Carvalho Pinto (Federal University of Vicosa).

D. Sc Domingos Sarvio Magalhes Valente (Federal University of Vicosa). 

************************************************************
Citation:
Please cite our manuscript: 

Andre Luiz de Freitas Coelho, Daniel Marçal de Queiroz, Domingos Sárvio Magalhães Valente, Francisco de Assis de Carvalho Pinto,
An open-source spatial analysis system for embedded systems,
Computers and Electronics in Agriculture,
Volume 154,
2018,
Pages 289-295,
ISSN 0168-1699,
https://doi.org/10.1016/j.compag.2018.09.019.
(http://www.sciencedirect.com/science/article/pii/S0168169918302552)
Keywords: Precision agriculture; Ordinary kriging; Clustering analysis; Yield map; Data filtering

************************************************************
Version : 1.2
- The function for automatic calculate the contour was removed. The user need to enter with limite coordinates.

- Add option to define the delimitter of dataset file and if coordinate system is metric (UTM) or geodesic (Latitude/Longitude).

Dataset and limit file need to using the same delimitter and the same coordinate sytem (Latitude/Longitude or UTM). 

- Function to embedded matplotlib graphics in PyQt was inserted

- A bug in gauss function was adjusted.

- Pandas function was used to read and write files

- The clustering function was change to k-means, using scikit-learn python module


**********************************************************
Tips for use

-The result files it is save in a Res Directory.

-It is recommend create this directory, in the projetc directory

-Before click in select dataset and select limite file, the delimitter need to be selected

-It is necessary informe if dataset coordinate already it is UTM

-The limite file and dataset file need to have a header.

-In limite file, first column need to be longitude (or coordinate e-w) 

-In limite file,second column need to be latitude (or coordinate n-s)

**********************************************************
This software was tested in:

- Laptop with Windows 7 64 bits

- Laptop with Xubuntu 16.04 64 bits

- Beablebone Black Rev B with Debian 8.6

*********************************************************
This repository contains:

- Source Code : This folder contains the files of source code in Python 3.6

- limite_spad.csv and limite_soybean.csv : Limite files for two dataset, used in tests

- spad.csv and soybean.csv : Two dataset used in test,  using metric coordinate system (UTM) and geodesic coordinate system


**************************************************************
The results files, in csv format, for filtering, kriging and clustering is save in /res folder
***************************************************************
Instructions for use the source code:

-In Windows OS it is necessary to install a python compiler. You can use:

Python IDE :  https://www.python.org/downloads/

Anaconda : https://anaconda.org/anaconda/python

Or antoher python compiler

-In linux the python is native. But you can use a python IDE:

Oficial python IDE : apt-get install idle3

-----------------------------------------------
In Windows and Linux OS is necessary install specific python modules for use the software developed. For this, 
it is necessary to use the pip comand of python. 

In linux, for install pip in Python3: apt-get install python3-pip

For use pip : python3 -m pip install (module name)

In Anaconda, you can use the pip comand on terminal or use the conda comand.

-PyQt5

apt-get install pyqt5-dev-tools (in linux)
or 

apt-get install python3-qt5 (in linux)
or 

conda install -c dsdale24 pyqt5  (in Anaconda)
or

https://www.riverbankcomputing.com/software/pyqt/download5  (In Windows)

-Numpy

pip install numpy (in Windows)
or

python3 -m pip install numpy (Linux)
or

conda install -c anaconda numpy  (in Anaconda)

-Scipy

pip install scipy (in Windows)
or

python3 -m pip install scipy (Linux)
or

conda install -c anaconda scipy  (in Anaconda)

-MatplotLib

pip install matplotlib (in Windows)
or

python3 -m pip install matplotlib (Linux)
or

conda install -c conda-forge matplotlib  (in Anaconda)

-Pandas

pip install pandas (in Windows)
or

python3 -m pip install pandas (Linux)
or

conda install -c anaconda pandas  (in Anaconda)

-Utm

pip install utm (in Windows)
or

python3 -m pip install utm (Linux)
or

conda install -c conda-forge utm  (in Anaconda)


-Scikit-learn
pip install  scikit-learn

or
conda install scikit-learn

-----------------------------------------------
Description of files in Source Code Folder:

DAGApy.py 

It is the main script. It is responsible for generating a graphical user interface (GUI), 
using PyQt5. In GUI, the user can perform data acquisition, and perform
data analysis : execute filtering, kriging and clustering
In the tab Parameter, is possible to modify the default parameters
This script automatically import the other script file as a python module.

All operations it is performed using metric coordinate system

RunCluster Function executes k-means clustering.
By defaults, 2 classes are generated. 
The scikit-learn module it is used 

-----------------------------------------------
DAGApy.ipynb 

File for used Anaconda Jupyter Notebook

-----------------------------------------------
filtering.py

Software that executes filtering of data.
It is implemented three filtering processes:
- Pre filtering : the user input minimum and maximum acceptable values. The software deletes all the data that is outside of the set range.
- Outlier filtering : The software deletes data outside the range defined by mean 
and standard deviation of the data set. This is, the range is defined by mean - x * sd and mean + x * sd. The default x factor, in this software is 2.5, but can be modified by the user. 
Before the data removal, the software compares the candidate point with its neighbours. By default, all neighbours within a distance lower than 25 m is used. 
The algorithm calculates the median for the selected neighbhours, if the median value differs a lot of the analyzed point value, then the analyzed point is considered an error point and will be deleted from the data set.
The admissible range is (1-x) * median and (1 + x) * median. The default x factor is 0.25, and cannot be modified by user.
- Inlier filtering. The software compares each point with its neighbor.By default, all 
neighbors within a distance lower than 25 m is used. For the neighbor point, the median of data is calculate. If median of neighbor is much different of candidate point. This point it's really false and is deleted. The admissible range is (1-x) * median and (1 + x) * median. The default x factor is 0.25, and cannot be modified by user.
-----------------------------------------------
kriging.py

Software that executes ordinary kriging interpolation.
The default distance for the grid is 10 m, but can be modified by user.
The grid it is generate based on limite points informing by user.Three semivariogram models are fit : spherical, exponential and gauss. 
The range, sill and nuuget automatic calculates for the three models. The selection criteria is to minimize the sum of square errors.
Automatically, the software chooses the model that results in the lower sum of square errors.


--------------------------------------------
gui.py 

This file contains all the information regarding the GUI. This script is imported by main.py software.
It is generated by yield.ui, using the command pyuic5 -x guiv11.ui -o gui.py
The yield.ui is created using Qt Designer 4 software. 
For more information, read about PyQt5.

---------------------------------------------------
icon.png

It's the icon software

-------------------------------------------------

mplwidget.py

It is a python file that allow embedded a matplotlib graphics in a Widget Element of PyQt.
In Qt Designer it is necessary insert a Widget element and promote to mplwidget.
