grab-tag-graph (GTG) Introduction
=================================
Grab 'em, Tag 'em, Graph 'em (GTG) is a feature detection, evolution and feature characterization algorithm. It was built for the weather/ climate application of identifying mesoscale convective complexes in highly resolved temporal and geospatial remote-sensed datasets and characterizing the features using various data sources. 

The algorithm is implemented in Python. The user provides a set of inputs that define the feature to be identified as well as the smallest feature that could evolve in time to be that larger feature of interest. Data is read from netCDF files into arrays with the dimensions time,latitude,longitude,value. The data is searched for areas of interest that evolved in time to make the feature, and these areas are stored in [Networkx](https://networkx.github.io/) graph objects. The graphs are traversed via graph algorithms such as Dijkstra's shortest path to determine the core feature and the feature. 

Using GTG
=========
1. Download the source code tarball. The tarball contains some data to get you started (baselinDataDir.zip)
2. Once you unzip the tarball, unzip and change directory to **code/**
3. `python mccSearchUI.py`
mccSearchUI.py is the use of GTG for the weather application of finding mesoscale convective complexes (MCCs) in
When you run the command line, mccSearchUI.py, there are a number of inputs that will be required.
There are the main inputs you will have to supply (as well as their variable names in the code - especially useful for editing mainProg.py or tracing the code:
* mainDirStr : this is the working directory where you wish all the output files –images, textfiles, clipped netCDF files- to be store
* There is an option to preprocess data. The program 
* CEoriDirName : this is the directory where the original MERG data in netCDF format is stored
* TRMMdirName : this is the directory where the original TRMM data in netCDF format is stored
* Start date and time in the format yyyymmddhr
* End date and time in the format yyyymmddhr

The following assumptions are made:
* input data are in one folder. For MERG data this is CEoriDirName and for the TRMM data this is TRMMdirName in mainProg.py.  These directories cannot be the same.
* THERE IS NO FILE CHECKING. So please ensure ALL your files are there in netCDF format. 

4. Wait a bit!
5. Once everything went well, the directory you indicated where outputs should be stored when prompted in the CL will be generated, and four folders should appear in it. 
* image/  stores all images generated from plots.
* textFiles/ stores all the text files generated during the run, e.g. cloudElementsUserFile.txt that contains information about each cloud element identified
* MERGnetcdfCEs/ contains the infrared data in masked netCDF files that have been generated for each cloud element identified
* TRMMnetcdfCEs/ contains the precipitation data in clipped netCDF files that have been generated for each cloud element identified.


Project website
===============
The project wiki can be found at [https://github.com/kwhitehall/grab-tag-graph/wiki](https://github.com/kwhitehall/grab-tag-graph/wiki). Please visit here to learn more about the codebase structure and the program. Details of installation and usage are also available here. 


Community Mailing List
======================
Please visit the [Google Group](https://groups.google.com/d/forum/gtg-users) to see what's going on this project. It is a forum for persons users (not necessarily developers to discuss problems compiling & running the code, and suggestions for upgrades.

Subscribe at [gtg-users@googlegroups.com](gtg-users@googlegroups.com)
