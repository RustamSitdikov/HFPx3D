# HFPx3D
3D BEM
Volume Control solver for pressurized fractures

![GEL-EPFL logo](http://gel.epfl.ch/files/content/sites/gel/files/Pictures/LOGOGEL-final-right-01.png)

**MeshFiles/** contains surface mesh (triangulation) data.
There are two files for each surface:
- **Elems_.npy** for connectivity matrix (node numbers of each element)
- **Nodes_.npy** for coordinates of the nodes

**Save_m_as_npy** sctipt can be used to save data from MATLAB .mat files in
.npy format.

**Important:** This code uses the **InsideLoop library** (https://github.com/insideloop/InsideLoop)

![InsideLoop icon](http://www.insideloop.io/wp-content/uploads/2014/09/inside-loop-logo-front.png)

Copy the **il/** folder to the project directory (not to **src/**).

**Notice (temporary):** To resolve an issue with `il::I`, in **il/math.h** replace "I" with "ii"