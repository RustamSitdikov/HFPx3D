# HFPx3D
3D BEM solver for pressurized fractures

Volume Control version

developed by Geo-Energy Laboratory, École Polytechnique Fédérale de Lausanne
(EPFL), Switzerland.

![GEL-EPFL logo](http://gel.epfl.ch/files/content/sites/gel/files/Pictures/LOGOGEL-final-right-01.png?v=2&s=50)

Please note that the project is a work in progress.
Updates will be followed up with a short description here
or in more details on the **wiki** page
(https://github.com/GeoEnergyLab-EPFL/HFPx3D/wiki).

## Important notes
The code requires the configuration file **config.toml** containing the
following data:
- Input/output directories/files (including surface mesh data)
- Material properties
- Load parameters (in-situ stress, injection)
- Numerical model parameters
##
Mesh (triangulation) input options:
1. **numpy** binary files, there are two for the same surface mesh:
    - **Elems_... .npy** for connectivity matrix (node numbers of each element)
    - **Nodes_... .npy** for coordinates of the nodes

  Both are stored and parsed as 3-row arrays (or >3 if additional data,
e.g. on crack tip nodes, are provided)

  32- and 64-bit floating point data are marked as **_32.npy** and **_64.npy**.

**Save_m_as_npy** sctipt can be used to save data from MATLAB .mat files in
**numpy** format.

2. **.csv** text files, also two for the same surface mesh:
    - **Elems_... .csv** for connectivity matrix (node numbers of each element)
    - **Nodes_... .csv** for coordinates of the nodes

  Note: In **.csv** files, both are stored as 3-column arrays,
but parsed as 3-row arrays.
##
Output options:
- nodal points (fill set for 2nd order approximation)
& solution (surface displacements/displacement discontinuities)
at these points
    - mark ```do_save_solution = true``` in **config.toml**
- stresses (with positive tension convention) at given observation points
(provided in an input file, either **.npy** or **.csv**,
in a way similar to nodes of the mesh)
    - mark ```do_postprocess = true``` and specify
  the observation points file (+ directory and file format) in **config.toml**
- the elastisity "influence" matrix
    - mark ```do_save_matrix = true``` in **config.toml**
##
**Mesh_Files/** folder contains some surface mesh (triangulation) examples.

**Test_Output/** folder contains solution and postprocessing examples +
observation points files.

##
This code uses the **InsideLoop library** (https://github.com/insideloop/InsideLoop)

![InsideLoop icon](http://www.insideloop.io/wp-content/uploads/2014/09/inside-loop-logo-front.png)

Although the source codes of the library are provided, it is strongly suggested
to use the latest version of the library. Please clone/download the source
files from https://github.com/insideloop/InsideLoop
and copy the **il/** folder to the project directory (not to **src/**!).

##
**CMakeLists.txt** may require modifications depending on the math library used.
If OpenBlas is used instead of MKL (typically, in Windows), set the
corresponding variables as shown:
```
set(IL_OPENBLAS 1)
set(IL_MKL 0)
```
or comment and un-comment the corresponding blocks.