// main.cpp will be used for testing the code parts under development

//#include <cstdio>
#include <iostream>
//#include <cassert>
//#include <sys/stat.h>
//#include <complex>
//#include <ittnotify.h>

#include <il/Timer.h>
#include <il/Toml.h>
#include <il/String.h>
#include <il/Array.h>
#include <il/Array2D.h>
//#include <il/StaticArray.h>
//#include <il/StaticArray2D.h>
#include <il/linear_algebra/dense/factorization/LU.h>
#include <il/linear_algebra.h>
//#include <il/linear_algebra/dense/factorization/linear_solve.h>

#include "src/mesh_file_io.h"
#include "src/system_assembly.h"
#include "src/surface_mesh_utilities.h"
#include "src/element_utilities.h"
#include "src/tensor_utilities.h"

int main() {

    // source files directory (containing main.cpp)
    std::string src_f = __FILE__;
    while (src_f.find("\\")!=std::string::npos) {
        src_f.replace(src_f.find("\\"),1,"/");
    }
    std::string src_dir = src_f.substr(0, src_f.rfind("/"));
    // std::string src_dir{"C:/Users/nikolski/ClionProjects/HFPx3D_VC"};
    // std::string src_dir{"/home/nikolski/Documents/HFPx3D"};
    // std::string src_dir{"/home/lecampio/Documents/HFPx3D"};

    std::string default_f_name = src_dir + "/config.toml";
    il::String config_f_name(default_f_name.c_str());
    // il::String config_f_name =
    // "C:/Users/nikolski/ClionProjects/HFPx3D_VC/config.toml";

    std::string default_input_dir{src_dir + "/Mesh_Files/"};
    // std::string mesh_conn_fname{"Elems_pennymesh24el_32.npy"};
    // std::string nodes_crd_fname{"Nodes_pennymesh24el_32.npy"};

    std::string default_output_dir{src_dir + "/Test_Output/"};

    il::String mf_name = "test_assembly_24_ele.csv";
    il::String of_name = "test_solution_24_ele.csv";
    // std::string mf_name{"test_assembly_24_ele.csv"};
    // std::string of_name{"test_solution_24_ele.csv"};

    //il::String src_f_name(src_dir.c_str());
    il::String in_dir_name((src_dir + "/").c_str());
    il::String out_dir_name((src_dir + "/").c_str());
    il::String m_c_f_name;
    il::String m_n_f_name;
    il::String i_f_format;

    il::Status status{};

    // reading configuration & parameters
    auto config =
            il::load<il::MapArray<il::String, il::Dynamic>>
                    (config_f_name, il::io, status);

    status.abortOnError();

    // reading input (triangulation) files
    il::int_t pos = config.search("input_directory");
    if (config.found(pos) && config.value(pos).isString()) {
        in_dir_name.append(config.value(pos).asString());
    } else {
        in_dir_name = il::String(default_input_dir.c_str());
    }
    in_dir_name.append("/");
    pos = config.search("mesh_conn_fname");
    if (config.found(pos) && config.value(pos).isString()) {
        m_c_f_name.append(config.value(pos).asString());
    } else {
        std::cout << "Can't find the input file" << std::endl;
        abort();
    }
    pos = config.search("nodes_crd_fname");
    if (config.found(pos) && config.value(pos).isString()) {
        m_n_f_name.append(config.value(pos).asString());
    } else {
        std::cout << "Can't find the input file" << std::endl;
        abort();
    }
    pos = config.search("input_format");
    if (config.found(pos) && config.value(pos).isString()) {
        i_f_format = config.value(pos).asString();
    } else {
        i_f_format = "npy32";
    }
    pos = config.search("output_directory");
    if (config.found(pos) && config.value(pos).isString()) {
        out_dir_name.append(config.value(pos).asString());
    } else {
        out_dir_name = il::String(default_output_dir.c_str());
    }


    // material properties (default)
    //hfp3d::Properties_T properties_list;
    //properties_list.mu = il::Array<double>{1, 1.0};
    //properties_list.nu = il::Array<double>{1, 0.35};
    double mu = 1.0, nu = 0.35;

    // numerical simulation parameters
    hfp3d::Num_Param_T n_par; // default: beta = 0.125; tip_type = 1; DD in global

    hfp3d::Load_T load;
    // stress at infinity (in-situ stress)
    load.s_inf = il::StaticArray<double, 6> {0.0};
    // load.s_inf[2] = 1.0; load.s_inf[4] = 1.0;
    pos = config.search("S_xx");
    if (config.found(pos) && config.value(pos).isDouble()) {
        load.s_inf[0] = (config.value(pos).toDouble());
    }
    pos = config.search("S_yy");
    if (config.found(pos) && config.value(pos).isDouble()) {
        load.s_inf[1] = (config.value(pos).toDouble());
    }
    pos = config.search("S_zz");
    if (config.found(pos) && config.value(pos).isDouble()) {
        load.s_inf[2] = (config.value(pos).toDouble());
    }
    pos = config.search("S_xy");
    if (config.found(pos) && config.value(pos).isDouble()) {
        load.s_inf[3] = (config.value(pos).toDouble());
    }
    pos = config.search("S_xz");
    if (config.found(pos) && config.value(pos).isDouble()) {
        load.s_inf[4] = (config.value(pos).toDouble());
    }
    pos = config.search("S_yz");
    if (config.found(pos) && config.value(pos).isDouble()) {
        load.s_inf[5] = (config.value(pos).toDouble());
    }

    // mesh
    hfp3d::Mesh_Geom_T mesh;

    // loading the mesh from files
    if (i_f_format == "csv") {
        hfp3d::load_mesh_from_csv
                (in_dir_name, m_c_f_name, m_n_f_name, true, il::io, mesh);
    } else if (i_f_format == "npy64") {
        hfp3d::load_mesh_from_numpy_64
                (in_dir_name, m_c_f_name, m_n_f_name, true, il::io, mesh);
    } else { // treat as 32-bit numpy by default
        hfp3d::load_mesh_from_numpy_32
                (in_dir_name, m_c_f_name, m_n_f_name, true, il::io, mesh);
    }

    //hfp3d::load_mesh_from_numpy_32
    // (input_dir, mesh_conn_fname, nodes_crd_fname,
    // true, il::io, mesh);

    // resetting the timer
    il::Timer timer{};
    timer.start();

//    __itt_resume();

    //hfp3d::DoF_Handle_T dof_handle;
    //dof_handle = hfp3d::make_dof_h_crack(mesh ,2, 1);
    //il::Array2D<double> bem_matrix;
    //bem_matrix = hfp3d::make_3dbem_matrix_s(mu, nu, mesh, n_par, il::io, dof_handle);
    //il::Array<double> rhs{num_dof, 0.0};

    il::int_t num_elems = mesh.conn.size(1);
    hfp3d::Mesh_Data_T mesh_data;
    mesh_data.mesh = mesh;
    //mesh_data.mat_id = hfp3d::make_mat_id_triv(mesh, 2);
    mesh_data.dof_h_dd = hfp3d::make_dof_h_crack(mesh, 2, 1);
    il::int_t num_dof = mesh_data.dof_h_dd.n_dof;
    mesh_data.dd = il::Array2D<double> {num_elems * 6, 3, 0.0};

    il::Array2D<il::int_t> ip(0, 7, 0);
    hfp3d::Mesh_Data_T mdf = hfp3d::init_mesh_data_p_fault(mesh, 2, ip);

    hfp3d::SAE_T sae;
    sae.matrix = hfp3d::make_3dbem_matrix_s
            (mu, nu, mesh, n_par, il::io, mesh_data.dof_h_dd);

    sae.rhs_v = il::Array<double> {num_dof, 0.0};
    for (il::int_t el = 0; el < num_elems; ++el) {
        // element vertices' coordinates
        il::StaticArray2D<double, 3, 3> el_vert;
        for (il::int_t j = 0; j < 3; ++j) {
            il::int_t n = mesh.conn(j, el);
            for (il::int_t k = 0; k < 3; ++k) {
                el_vert(k, j) = mesh.nods(k, n);
            }
        }

        // Rotation tensor for the element
        il::StaticArray2D<double, 3, 3> r_tensor =
                hfp3d::make_el_r_tensor(el_vert);

        // Normal vector at collocation point (x)
        il::StaticArray<double, 3> norm_v;
        for (int j = 0; j < 3; ++j) {
            norm_v[j] = -r_tensor(2, j);
        }
        // induced traction
        il::StaticArray<double, 3> trac_inf =
                hfp3d::nv_dot_sim(norm_v, load.s_inf);

        // setting the RHS of the system (loads)
        for (int ln = 0; ln < 6; ++ln) {
            //il::int_t n = el * 6 + ln;
            for (int l = 0; l < 3; ++l) {
                il::int_t ldof = ln * 3 + l;
                //il::int_t dof = dof_handle.dof_h(el, ldof);
                il::int_t dof = mesh_data.dof_h_dd.dof_h(el, ldof);
                if (dof != -1) {
                    //rhs[dof] = (l != 1 ? 1.0 : 0.0);
                    sae.rhs_v[dof] = -trac_inf[l];
                }
            }
        }
    }
    il::Array<double> dd_v;

//    __itt_pause();
    timer.stop();

    std::cout << "Assembly: " << timer.elapsed() << "s" << std::endl;
    std::cout << 18 * num_elems << " DoF Total" << std::endl;
    std::cout << 18 * num_elems - num_dof << " Fixed DoF" << std::endl;

    timer.reset();
    timer.start();

    // dd_v = il::linear_solve(bem_matrix, rhs, il::io, status);
    //il::LU<il::Array2D<double>> lu_decomposition(bem_matrix, il::io, status);
    il::LU<il::Array2D<double>> lu_decomposition(sae.matrix, il::io, status);
    // if (!status.ok()) {
    //     // The matrix is singular to the machine precision. You should deal with
    //     // the error.
    // }
    status.abortOnError();
    // double cnd = lu_decomposition.condition_number(il::Norm::L2, );
    // std::cout << cnd << std::endl;
    //dd_v = lu_decomposition.solve(rhs);
    dd_v = lu_decomposition.solve(sae.rhs_v);

    timer.stop();
    std::cout << "Solution: " << timer.elapsed() << "s" << std::endl;

    hfp3d::write_dd_vector_to_md
            (dd_v, mesh_data.dof_h_dd,
             false, mesh_data.dof_h_pp,
             il::io, mesh_data);

    // saving matrix to a .CSV file
    hfp3d::save_data_to_csv(sae.matrix, out_dir_name, mf_name);

    // the 2D array for nodal points' coordinates and DD - initialization
    il::Array2D<double> out_dd(6 * num_elems, 6);

    // saving the solution (nodes + DD) to a .CSV file
    for (il::int_t j = 0; j < num_elems; ++j) {
        il::StaticArray2D<double, 3, 3> el_vert;
        for (il::int_t k = 0; k < 3; ++k) {
            il::int_t n = mesh.conn(k, j);
            for (il::int_t l = 0; l < 3; ++l) {
                el_vert(l, k) = mesh.nods(l, n);
            }
        }
        il::StaticArray<il::StaticArray<double, 3>, 6> el_np;
        el_np = hfp3d::el_cp_uniform(el_vert, 0.0);
        // el_np = hfp3d::el_cp_nonuniform(el_vert, v_wts, 0.0);
        for (il::int_t k = 0; k < 6; ++k) {
            il::int_t n = j * 6 + k;
            for (il::int_t l = 0; l < 3; ++l) {
                //il::int_t dof = n * 3 + l;
                il::int_t l_dof = k * 3 + l;
                il::int_t g_dof = mesh_data.dof_h_dd.dof_h(j, l_dof);
                out_dd(n, l) = el_np[k][l];
                if (g_dof != -1) {
                    //out_dd(n, l) = dd_v[dof];
                    out_dd(n, l + 3) = dd_v[g_dof];
                } else {
                    //out_dd(n, l) = 0.0;
                    out_dd(n, l + 3) = 0.0;
                }
            }
        }
    }

    hfp3d::save_data_to_csv(out_dd, out_dir_name, of_name);

/*
    std::string npy_of_name{"/test_solution_24_ele.npy"};
    // il::Status status{};
    il::save(out_dd, output_dir + npy_of_name, il::io, status);
    status.abortOnError();
*/

    return 0;
}