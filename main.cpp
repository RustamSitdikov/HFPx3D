// main.cpp will be used for testing the code parts under development

//#include <cstdio>
#include <iostream>
#include <il/Array.h>
#include <il/Array2D.h>
#include <il/linear_algebra/dense/factorization/LU.h>
#include <il/linear_algebra.h>
//#include <il/linear_algebra/dense/factorization/linear_solve.h>
#include "mesh_file_io.h"
#include "system_assembly.h"
#include "element_utilities.h"
//#include <complex>
//#include <il/StaticArray.h>
//#include <il/StaticArray2D.h>
//#include "ele_base.h"

int main() {

    double mu = 1.0, nu = 0.35;

    std::string src_directory{""}; // add path
    std::string mesh_conn_fname{"Elems_pennymesh121el_32.npy"};
    std::string nodes_crd_fname{"Nodes_pennymesh121el_32.npy"};

    std::string work_directory{""}; // add path
    std::string mf_name{"test_assembly_121_ele.csv"};
    std::string of_name{"test_solution_121_ele.csv"};

    il::Array2D<il::int_t> mesh_conn;
    il::Array2D<double> nodes_crd;
    hfp3d::load_mesh_from_numpy_32
            (src_directory, mesh_conn_fname, nodes_crd_fname, true,
             il::io, mesh_conn, nodes_crd);

    hfp3d::Dof_Handle dof_hndl;
    il::Array2D<double> bem_matrix;
    bem_matrix = hfp3d::make_3dbem_matrix_s
            (mu, nu, 0.25, mesh_conn, nodes_crd, 1, il::io, dof_hndl);
    il::int_t num_elems = mesh_conn.size(1);
    il::int_t num_dof = dof_hndl.n_dof;
    std::cout << "Full No of DOF = " << 18 * num_elems << std::endl;
    std::cout << "No of used DOF = " << num_dof << std::endl;

    hfp3d::save_data_to_csv(bem_matrix, work_directory, mf_name);

    il::Array<double> rhs(num_dof);
    for (il::int_t j = 0; j < num_elems; ++j) {
        for (int k = 0; k < 6; ++k) {
            il::int_t n = j * 6 + k;
            for (int l = 0; l < 3; ++l) {
                //il::int_t dof = n * 3 + l;
                il::int_t dof = dof_hndl.dof_h(n, l);
                if (dof != -1) {
                    rhs[dof] = (l != 1 ? 1.0 : 0.0);
                }
            }
        }
    }
    il::Status status{};
    il::Array<double> dd_v;
    il::LU<il::Array2D<double>> lu_decomposition(bem_matrix, il::io, status);
    // if (!status.ok()) {
    //     // The matrix is singular to the machine precision. You should deal with
    //     // the error.
    // }
    status.abort_on_error();
    // double cnd = lu_decomposition.condition_number(il::Norm::L2, );
    // std::cout << cnd << std::endl;
    dd_v = lu_decomposition.solve(rhs);
    // dd_v = il::linear_solve(bem_matrix, rhs, il::io, status);

    il::Array2D<double> dd(6 * num_elems, 6);
    for (il::int_t j = 0; j < num_elems; ++j) {
        il::StaticArray2D<double, 3, 3> el_vert;
        for (il::int_t k = 0; k < 3; ++k) {
            il::int_t n = mesh_conn(k, j);
            for (il::int_t l = 0; l < 3; ++l) {
                el_vert(l, k) = nodes_crd(l, n);
            }
        }
        il::StaticArray<il::StaticArray<double, 3>, 6> el_np;
        el_np = hfp3d::el_cp_uniform(el_vert, 0.0);
        // el_np = hfp3d::el_cp_nonuniform(el_vert, v_wts, 0.0);
        for (il::int_t k = 0; k < 6; ++k) {
            il::int_t n = j * 6 + k;
            for (il::int_t l = 0; l < 3; ++l) {
                //il::int_t dof = n * 3 + l;
                il::int_t dof = dof_hndl.dof_h(n, l);
                dd(n, l) = el_np[k][l];
                if (dof != -1) {
                    //dd(n, l) = dd_v[dof];
                    dd(n, l + 3) = dd_v[dof];
                } else {
                    //dd(n, l) = 0.0;
                    dd(n, l + 3) = 0.0;
                }
            }
        }
    }
    hfp3d::save_data_to_csv(dd, work_directory, of_name);

/*
    std::string npy_of_name{"/test_solution_24_ele.npy"};
    // il::Status status{};
    il::save(dd, out_directory + npy_of_name, il::io, status);
    status.abort_on_error();
*/

    return 0;
}