# HFPx3D
3D BEM solver for pressurized fractures

Volume Control version

developed by Geo-Energy Laboratory, École Polytechnique Fédérale de Lausanne
(EPFL), Switzerland.

![GEL-EPFL logo](http://gel.epfl.ch/files/content/sites/gel/files/Pictures/LOGOGEL-final-right-01.png?v=2&s=200)

Please note that the project is a work in progress.
Updates will be followed up with a short description here
or in more details on the **wiki** page
(https://github.com/GeoEnergyLab-EPFL/HFPx3D/wiki).

## Important notes
**MeshFiles/** folder contains surface mesh (triangulation) data.
There are two **numpy** binary files for each surface:
- **Elems_... .npy** for connectivity matrix (node numbers of each element)
- **Nodes_... .npy** for coordinates of the nodes

32- and 64-bit floating point data are marked as **_32.npy** and **_64.npy**.

**Save_m_as_npy** sctipt can be used to save data from MATLAB .mat files in
**numpy** format.
##
This code uses the **InsideLoop library** (https://github.com/insideloop/InsideLoop)

![InsideLoop icon](http://www.insideloop.io/wp-content/uploads/2014/09/inside-loop-logo-front.png)

Although the source codes of the library are provided, it is strongly suggested
to use the latest version of the library. Please clone/download the source
files from https://github.com/insideloop/InsideLoop
and copy the **il/** folder to the project directory (not to **src/**!).

##
**CMakeLists.txt** may require modifications depending on the math library used.
If OpenBlas is used instead of MKL (typically, in Windows), comment and
un-comment the corresponding blocks.