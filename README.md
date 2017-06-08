# HFPx3D
3D BEM
Volume Control solver for pressurized fractures

Developed by Geo-Energy Laboratory, École Polytechnique Fédérale de Lausanne
(EPFL),
Switzerland.

![GEL-EPFL logo](http://gel.epfl.ch/files/content/sites/gel/files/Pictures/LOGOGEL-final-right-01.png)

## Important notes

**MeshFiles/** contains surface mesh (triangulation) data.
There are two binary files for each surface:
- **Elems_... .npy** for connectivity matrix (node numbers of each element)
- **Nodes_... .npy** for coordinates of the nodes

32- and 64-bit floating point data are marked with **_32** and **_64**.

**Save_m_as_npy** sctipt can be used to save data from MATLAB .mat files in
.npy format.
##
**Important:** This code uses the **InsideLoop library** (https://github.com/insideloop/InsideLoop)

![InsideLoop icon](http://www.insideloop.io/wp-content/uploads/2014/09/inside-loop-logo-front.png)

It is strongly suggested to use the latest version of the library. Download
and copy the **il/** folder to the project directory (not to **src/**!).

**Notice (temporary):** To resolve an issue with `il::I`, in **il/math.h** replace "I" with "ii"