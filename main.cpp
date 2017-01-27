#include <cstdio>
#include <il/Array2D.h>
#include <il/StaticArray.h>
#include <il/StaticArray2D.h>
//#include <il/io/numpy.h>
#include "mesh_file_io.h"
#include "matrix_asm.h"
#include "ele_base.h"

// main.cpp will be used for testing the code parts under development

int main() {

/*
    il::Array2D<il::int_t> mesh_conn(3, 1);
    il::Array2D<double> nodes_crd(3, 3);
    for (il::int_t k = 0; k < 3; ++k) {
        mesh_conn(k, 0) = k;
    }
    nodes_crd(0, 0) = 0.0;
    nodes_crd(1, 0) = 0.1;
    nodes_crd(2, 0) = 0.0;

    nodes_crd(0, 1) = 1.8;
    nodes_crd(1, 1) = 0.0;
    nodes_crd(2, 1) = 0.0;

    nodes_crd(0, 2) = 1.2;
    nodes_crd(1, 2) = 1.8;
    nodes_crd(2, 2) = 0.0;
*/

    std::string src_directory{"C:/Users/nikolski/.spyder-py3/3DBEM/"};
    std::string mesh_conn_fname{"Elems_pennymesh24el.npy"};
    std::string nodes_crd_fname{"Nodes_pennymesh24el.npy"};

    il::Array2D<il::int_t> mesh_conn;
    il::Array2D<double> nodes_crd;
    il::StaticArray<il::int_t, 2> nums = hfp3d::load_mesh_from_numpy
            (src_directory, mesh_conn_fname, nodes_crd_fname, true,
             il::io, mesh_conn, nodes_crd);

    il::int_t num_elems = mesh_conn.size(1), num_dof = 18 * num_elems;

    double mu = 1.0, nu = 0.35;

    il::Array2D<double> bem_matrix(num_dof, num_dof);
    bem_matrix = hfp3d::make_3dbem_matrix_s(mu, nu, 0.25, mesh_conn, nodes_crd);

    std::string work_directory{"C:/Users/nikolski/ClionProjects/3D-bem"
                                       "/Test_Output/"};
    std::string of_name{"test_assembly_1_ele.csv"};

    hfp3d::save_matrix_to_csv(bem_matrix, work_directory, of_name);

    //il::save(IM_1, "C:/Users/nikolski/.spyder-py3/3DBEM/matrix_24_el.npy",
    // il::io, status);
    //status.abort_on_error();

    return 0;
}