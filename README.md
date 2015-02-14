grab-tag-graph (GTG) Introduction
=================================
Grab 'em, Tag 'em, Graph 'em (GTG) is a feature dectection, evolution and feature characterization algorithm. It was built for the weather/ climate aplication of identifying mesoscale convective complexes in highly resolved temporal and geospatial remote-sensed datasets and characterizing the features using various data sources. 

The algorithm is implemented in Python. The user provides a set of inputs that define the feature to be identified as well as the smallest feature that could evolve in time to be that larger feature of interest. Data is read from netCDF files into arrays with the dimensions time,latitude,longitude,value. The data is searched for areas of interest that evolved in time to make the feature, and these areas are stored in Networkx graph objects. The graphs are traversed via graph algorithms such as Dijkstra's shortest path to determine the core feature and the feature. 


Project website
===============
The project wiki can be found at [https://github.com/kwhitehall/grab-tag-graph/wiki](https://github.com/kwhitehall/grab-tag-graph/wiki). Please visit here to learn more about the codebase structure and the program. Details of installation and usage are also available here. 


Community Mailing List
======================
Please visit the [Google Group](https://groups.google.com/d/forum/gtg-users) to see what's going on this project. It is a forum for persons users (not necessarily developers to discuss problems compiling & running the code, and suggestions for upgrades.

Subscribe at [gtg-users@googlegroups.com](gtg-users@googlegroups.com)
