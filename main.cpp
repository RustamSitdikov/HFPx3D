#include <cstdio>
#include <il/Array2D.h>
#include <il/StaticArray.h>
#include <il/StaticArray2D.h>
//#include <il/io/numpy.h>
#include "mesh_file_io.h"
#include "matrix_asm.h"
#include "ele_base.h"
//#include "tensor_oper.h"
#include <complex>
#include <il/StaticArray3D.h>
#include "h_potential.h"

// main.cpp will be used for testing the code parts under development

int main() {

    double mu = 1.0, nu = 0.35;

    std::string work_directory{"C:/Users/nikolski/ClionProjects/3D-bem"
                                       "/Test_Output/"};

    // Matrix assembly for a penny-shaped crack (24 elements)

    std::string src_directory{"C:/Users/nikolski/.spyder-py3/3DBEM/"};
    std::string mesh_conn_fname{"Elems_pennymesh24el.npy"};
    std::string nodes_crd_fname{"Nodes_pennymesh24el.npy"};

    std::string of_name{"test_assembly_24_ele.csv"};

    il::Array2D<il::int_t> mesh_conn;
    il::Array2D<double> nodes_crd;
    //il::StaticArray<il::int_t, 2> nums =
            hfp3d::load_mesh_from_numpy
            (src_directory, mesh_conn_fname, nodes_crd_fname, true,
             il::io, mesh_conn, nodes_crd);

    il::int_t num_elems = mesh_conn.size(1), num_dof = 18 * num_elems;

    il::Array2D<double> bem_matrix(num_dof, num_dof);
    bem_matrix = hfp3d::make_3dbem_matrix_s(mu, nu, 0.25, mesh_conn, nodes_crd);

/*
    // Testing kernel integration for one element

    il::StaticArray2D<il::int_t, 3, 1> mesh_conn;
    il::StaticArray2D<double, 3, 3> nodes_crd;
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
    il::StaticArray<double, 3> m_p_crd;
    m_p_crd[0] = 1.0;
    m_p_crd[1] = 0.4;
    m_p_crd[2] = 1.21;
    il::StaticArray2D<double, 3, 3> r_tensor;
    il::StaticArray2D<std::complex<double>, 6, 6> sfm =
            hfp3d::make_el_sfm_uniform(nodes_crd, il::io, r_tensor);
    il::StaticArray<std::complex<double>, 3> tau =
            hfp3d::make_el_tau_crd(nodes_crd, r_tensor);
    hfp3d::HZ hz = hfp3d::make_el_pt_hz(nodes_crd, m_p_crd, r_tensor);
    il::StaticArray2D<double, 6, 18> bem_matrix =
            hfp3d::make_local_3dbem_submatrix(1, mu, nu, hz.h, hz.z, tau, sfm);
    std::string of_name{"test_integration_1_ele.csv"};
    std::string sf_name{"test_SFM_1_ele.csv"};

    hfp3d::save_data_to_csv(sfm, work_directory, sf_name);
*/

/*
    // Testing matrix assembly for one element

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
    std::string of_name{"test_assembly_1_ele.csv"};
*/

    hfp3d::save_data_to_csv(bem_matrix, work_directory, of_name);

/*
    il::save(IM_1, "C:/Users/nikolski/.spyder-py3/3DBEM/matrix_24_el.npy",
             il::io, status);
    status.abort_on_error();
*/

    return 0;
}