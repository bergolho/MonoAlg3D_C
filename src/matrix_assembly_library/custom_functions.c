//
// Created by bergolho and jennyhelyanwe on 03/06/24.
//

#include <stdbool.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "../3dparty/sds/sds.h"
#include "../3dparty/stb_ds.h"
#include "../alg/grid/grid.h"
#include "../config/assembly_matrix_config.h"
#include "../libraries_common/common_data_structures.h"
#include "../utils/file_utils.h"
#include "../utils/utils.h"
#include "../domains_library/custom_mesh_info_data.h"
#include "purkinje_coupling_matrix_assembly.c"

#ifdef COMPILE_CUDA
#include "../gpu_utils/gpu_utils.h"
#endif

// Matrix assembly function to configure the conductivities for the NORMAL, DENSE and SPARSE regions.
// Works with all DTI meshes
// Function used for the T-wave personalisation paper.
ASSEMBLY_MATRIX(anisotropic_sigma_assembly_matrix_with_dense_and_sparse_region_twave) {
    
    if(the_grid->adaptive) {
        log_error_and_exit("anisotropic_sigma_assembly_matrix_with_dense_and_sparse_region_twave function does not support mesh adaptivity yet!. Aborting!\n");
    }

    uint32_t num_active_cells = the_grid->num_active_cells;
    struct cell_node **ac = the_grid->active_cells;

    initialize_diagonal_elements(the_solver, the_grid);

    //      D tensor    //
    // | sx    sxy   sxz |
    // | sxy   sy    syz |
    // | sxz   syz   sz  |
    real_cpu D[3][3];
    real_cpu D_dense[3][3];
    real_cpu D_sparse[3][3];
    int i;

    real_cpu sigma_l = 0.0;
    real_cpu sigma_t = 0.0;
    real_cpu sigma_n = 0.0;
    real_cpu sigma_l_dense = 0.0;
    real_cpu sigma_t_dense = 0.0;
    real_cpu sigma_n_dense = 0.0;
    real_cpu sigma_l_sparse = 0.0;
    real_cpu sigma_t_sparse = 0.0;
    real_cpu sigma_n_sparse = 0.0;

    char *fiber_file = NULL;
    GET_PARAMETER_STRING_VALUE_OR_USE_DEFAULT(fiber_file, config, "fibers_file");

    bool fibers_in_mesh = false;
    GET_PARAMETER_BOOLEAN_VALUE_OR_USE_DEFAULT(fibers_in_mesh, config, "fibers_in_mesh");

    struct fiber_coords *fibers = NULL;

    // Normal myocardium conductivities
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real_cpu, sigma_l, config, "sigma_l");
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real_cpu, sigma_t, config, "sigma_t");
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real_cpu, sigma_n, config, "sigma_n");

    // Dense region myocardium conductivities
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real_cpu, sigma_l_dense, config, "sigma_l_dense");
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real_cpu, sigma_t_dense, config, "sigma_t_dense");
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real_cpu, sigma_n_dense, config, "sigma_n_dense");

    // Sparse myocardium conductivities
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real_cpu, sigma_l_sparse, config, "sigma_l_sparse");
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real_cpu, sigma_t_sparse, config, "sigma_t_sparse");
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real_cpu, sigma_n_sparse, config, "sigma_n_sparse");

    real_cpu *f = NULL;
    real_cpu *s = NULL;
    real_cpu *n = NULL;

    if(fiber_file) {
        log_info("Loading mesh fibers\n");
        fibers = read_fibers(fiber_file, false);
    }
    else if(!fibers_in_mesh) {
        GET_PARAMETER_VECTOR_VALUE_OR_USE_DEFAULT(f, config, "f", 3);
        GET_PARAMETER_VECTOR_VALUE_OR_USE_DEFAULT(s, config, "s", 3);
        GET_PARAMETER_VECTOR_VALUE_OR_USE_DEFAULT(n, config, "n", 3);

        if(!f) {
            f = malloc(sizeof(real_cpu)*3);
            f[0] = 1.0;
            f[1] = 0.0;
            f[2] = 0.0;
        }

        if(!s) {
            s = malloc(sizeof(real_cpu)*3);
            s[0] = 0.0;
            s[1] = 1.0;
            s[2] = 0.0;
        }

        if(!n) {
            n = malloc(sizeof(real_cpu)*3);
            n[0] = 0.0;
            n[1] = 0.0;
            n[2] = 1.0;
        }
    }

    OMP(parallel for private(D))
    for(i = 0; i < num_active_cells; i++) {

        if(fibers) {
            int fiber_index = ac[i]->original_position_in_file;

            if(fiber_index == -1) {
                log_error_and_exit("fiber_index should not be -1, but it is for cell in index %d - %lf, %lf, %lf\n", i, ac[i]->center.x, ac[i]->center.y, ac[i]->center.z);
            }

            if(sigma_t == sigma_n) {
                calc_tensor2(D, fibers[fiber_index].f, sigma_l, sigma_t);
                calc_tensor2(D_dense, fibers[fiber_index].f, sigma_l_dense, sigma_t_dense);
                calc_tensor2(D_sparse, fibers[fiber_index].f, sigma_l_sparse, sigma_t_sparse);
            }
            else {
                calc_tensor(D, fibers[fiber_index].f, fibers[fiber_index].s, fibers[fiber_index].n, sigma_l, sigma_t, sigma_n);
                calc_tensor(D_dense, fibers[fiber_index].f, fibers[fiber_index].s, fibers[fiber_index].n, sigma_l_dense, sigma_t_dense, sigma_n_dense);
                calc_tensor(D_sparse, fibers[fiber_index].f, fibers[fiber_index].s, fibers[fiber_index].n, sigma_l_sparse, sigma_t_sparse, sigma_n_sparse);
            }
            ac[i]->sigma.fibers = fibers[fiber_index];
        }
        else if(fibers_in_mesh) {
            if(sigma_t == sigma_n) {
                calc_tensor2(D, ac[i]->sigma.fibers.f, sigma_l, sigma_t);
                calc_tensor2(D_dense, ac[i]->sigma.fibers.f, sigma_l_dense, sigma_t_dense);
                calc_tensor2(D_sparse, ac[i]->sigma.fibers.f, sigma_l_sparse, sigma_t_sparse);
            }
            else {
                calc_tensor(D, ac[i]->sigma.fibers.f, ac[i]->sigma.fibers.s, ac[i]->sigma.fibers.n, sigma_l, sigma_t, sigma_n);
                calc_tensor(D_dense, ac[i]->sigma.fibers.f, ac[i]->sigma.fibers.s, ac[i]->sigma.fibers.n, sigma_l_dense, sigma_t_dense, sigma_n_dense);
                calc_tensor(D_sparse, ac[i]->sigma.fibers.f, ac[i]->sigma.fibers.s, ac[i]->sigma.fibers.n, sigma_l_sparse, sigma_t_sparse, sigma_n_sparse);
            }

        }
        else {
            if(sigma_t == sigma_n) {
                calc_tensor2(D, f, sigma_l, sigma_t);
                calc_tensor2(D_dense, f, sigma_l_dense, sigma_t_dense);
                calc_tensor2(D_sparse, f, sigma_l_sparse, sigma_t_sparse);
            }
            else {
                calc_tensor(D, f, s, n, sigma_l, sigma_t, sigma_n);
                calc_tensor(D_dense, f, s, n, sigma_l_dense, sigma_t_dense, sigma_n_dense);
                calc_tensor(D_sparse, f, s, n, sigma_l_sparse, sigma_t_sparse, sigma_n_sparse);
            }
        }

        // Check if the cell is tagged with an endocardium region tag
        // 2 = DENSE region || 3 = SPARSE region
        int tag = DTI_MESH_FAST_ENDO(ac[i]);

        // Cell is inside the DENSE endocardium region
        if (tag == 2) {
            ac[i]->sigma.x = D_dense[0][0];
            ac[i]->sigma.y = D_dense[1][1];
            ac[i]->sigma.z = D_dense[2][2];

            ac[i]->sigma.xy = D_dense[0][1];
            ac[i]->sigma.xz = D_dense[0][2];
            ac[i]->sigma.yz = D_dense[1][2];
        }
        // Cell is inside the SPARSE endocardium region
        else if (tag == 3) {
            ac[i]->sigma.x = D_sparse[0][0];
            ac[i]->sigma.y = D_sparse[1][1];
            ac[i]->sigma.z = D_sparse[2][2];

            ac[i]->sigma.xy = D_sparse[0][1];
            ac[i]->sigma.xz = D_sparse[0][2];
            ac[i]->sigma.yz = D_sparse[1][2];
        }
        // Otherwise, NORMAL type of cell
        else {
            ac[i]->sigma.x = D[0][0];
            ac[i]->sigma.y = D[1][1];
            ac[i]->sigma.z = D[2][2];

            ac[i]->sigma.xy = D[0][1];
            ac[i]->sigma.xz = D[0][2];
            ac[i]->sigma.yz = D[1][2];
        }
    }

    OMP(parallel for)
    for(i = 0; i < num_active_cells; i++) {
        fill_discretization_matrix_elements_aniso(ac[i]);
    }

    free(f);
    free(s);
    free(n);
}

// Matrix assembly function to configure the conductivities for the NORMAL, DENSE and SPARSE regions.
// + Works with Purkinje coupled simulations.
// Works with all DTI meshes
// Function used for the T-wave personalisation paper.
ASSEMBLY_MATRIX(purkinje_coupling_with_anisotropic_sigma_assembly_matrix_with_dense_and_sparse_region_twave) {

    if(the_grid->adaptive) {
        log_error_and_exit("purkinje_coupling_with_anisotropic_sigma_assembly_matrix_with_dense_and_sparse_region_twave function does not support mesh adaptivity yet!. Aborting!\n");
    }

// [MYOCARDIUM SETUP]
    uint32_t num_active_cells = the_grid->num_active_cells;
    struct cell_node **ac = the_grid->active_cells;

    initialize_diagonal_elements(the_solver, the_grid);

    //      D tensor    //
    // | sx    sxy   sxz |
    // | sxy   sy    syz |
    // | sxz   syz   sz  |
    real_cpu D[3][3];
    real_cpu D_dense[3][3];
    real_cpu D_sparse[3][3];
    int i;

    real_cpu sigma_l = 0.0;
    real_cpu sigma_t = 0.0;
    real_cpu sigma_n = 0.0;
    real_cpu sigma_l_dense = 0.0;
    real_cpu sigma_t_dense = 0.0;
    real_cpu sigma_n_dense = 0.0;
    real_cpu sigma_l_sparse = 0.0;
    real_cpu sigma_t_sparse = 0.0;
    real_cpu sigma_n_sparse = 0.0;
    real_cpu sigma_purkinje = 0.0;

    char *fiber_file = NULL;
    GET_PARAMETER_STRING_VALUE_OR_USE_DEFAULT(fiber_file, config, "fibers_file");

    bool fibers_in_mesh = false;
    GET_PARAMETER_BOOLEAN_VALUE_OR_USE_DEFAULT(fibers_in_mesh, config, "fibers_in_mesh");

    struct fiber_coords *fibers = NULL;

    // Normal myocardium conductivities
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real_cpu, sigma_l, config, "sigma_l");
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real_cpu, sigma_t, config, "sigma_t");
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real_cpu, sigma_n, config, "sigma_n");

    // Dense region myocardium conductivities
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real_cpu, sigma_l_dense, config, "sigma_l_dense");
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real_cpu, sigma_t_dense, config, "sigma_t_dense");
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real_cpu, sigma_n_dense, config, "sigma_n_dense");

    // Sparse myocardium conductivities
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real_cpu, sigma_l_sparse, config, "sigma_l_sparse");
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real_cpu, sigma_t_sparse, config, "sigma_t_sparse");
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real_cpu, sigma_n_sparse, config, "sigma_n_sparse");

    // Purkinje conductivity
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real, sigma_purkinje, config, "sigma_purkinje");

    real_cpu *f = NULL;
    real_cpu *s = NULL;
    real_cpu *n = NULL;

    if(fiber_file) {
        log_info("Loading mesh fibers\n");
        fibers = read_fibers(fiber_file, false);
    }
    else if(!fibers_in_mesh) {
        GET_PARAMETER_VECTOR_VALUE_OR_USE_DEFAULT(f, config, "f", 3);
        GET_PARAMETER_VECTOR_VALUE_OR_USE_DEFAULT(s, config, "s", 3);
        GET_PARAMETER_VECTOR_VALUE_OR_USE_DEFAULT(n, config, "n", 3);

        if(!f) {
            f = malloc(sizeof(real_cpu)*3);
            f[0] = 1.0;
            f[1] = 0.0;
            f[2] = 0.0;
        }

        if(!s) {
            s = malloc(sizeof(real_cpu)*3);
            s[0] = 0.0;
            s[1] = 1.0;
            s[2] = 0.0;
        }

        if(!n) {
            n = malloc(sizeof(real_cpu)*3);
            n[0] = 0.0;
            n[1] = 0.0;
            n[2] = 1.0;
        }
    }

    OMP(parallel for private(D))
    for(i = 0; i < num_active_cells; i++) {

        if(fibers) {
            int fiber_index = ac[i]->original_position_in_file;

            if(fiber_index == -1) {
                log_error_and_exit("fiber_index should not be -1, but it is for cell in index %d - %lf, %lf, %lf\n", i, ac[i]->center.x, ac[i]->center.y, ac[i]->center.z);
            }

            if(sigma_t == sigma_n) {
                calc_tensor2(D, fibers[fiber_index].f, sigma_l, sigma_t);
                calc_tensor2(D_dense, fibers[fiber_index].f, sigma_l_dense, sigma_t_dense);
                calc_tensor2(D_sparse, fibers[fiber_index].f, sigma_l_sparse, sigma_t_sparse);
            }
            else {
                calc_tensor(D, fibers[fiber_index].f, fibers[fiber_index].s, fibers[fiber_index].n, sigma_l, sigma_t, sigma_n);
                calc_tensor(D_dense, fibers[fiber_index].f, fibers[fiber_index].s, fibers[fiber_index].n, sigma_l_dense, sigma_t_dense, sigma_n_dense);
                calc_tensor(D_sparse, fibers[fiber_index].f, fibers[fiber_index].s, fibers[fiber_index].n, sigma_l_sparse, sigma_t_sparse, sigma_n_sparse);
            }
            ac[i]->sigma.fibers = fibers[fiber_index];
        }
        else if(fibers_in_mesh) {
            if(sigma_t == sigma_n) {
                calc_tensor2(D, ac[i]->sigma.fibers.f, sigma_l, sigma_t);
                calc_tensor2(D_dense, ac[i]->sigma.fibers.f, sigma_l_dense, sigma_t_dense);
                calc_tensor2(D_sparse, ac[i]->sigma.fibers.f, sigma_l_sparse, sigma_t_sparse);
            }
            else {
                calc_tensor(D, ac[i]->sigma.fibers.f, ac[i]->sigma.fibers.s, ac[i]->sigma.fibers.n, sigma_l, sigma_t, sigma_n);
                calc_tensor(D_dense, ac[i]->sigma.fibers.f, ac[i]->sigma.fibers.s, ac[i]->sigma.fibers.n, sigma_l_dense, sigma_t_dense, sigma_n_dense);
                calc_tensor(D_sparse, ac[i]->sigma.fibers.f, ac[i]->sigma.fibers.s, ac[i]->sigma.fibers.n, sigma_l_sparse, sigma_t_sparse, sigma_n_sparse);
            }

        }
        else {
            if(sigma_t == sigma_n) {
                calc_tensor2(D, f, sigma_l, sigma_t);
                calc_tensor2(D_dense, f, sigma_l_dense, sigma_t_dense);
                calc_tensor2(D_sparse, f, sigma_l_sparse, sigma_t_sparse);
            }
            else {
                calc_tensor(D, f, s, n, sigma_l, sigma_t, sigma_n);
                calc_tensor(D_dense, f, s, n, sigma_l_dense, sigma_t_dense, sigma_n_dense);
                calc_tensor(D_sparse, f, s, n, sigma_l_sparse, sigma_t_sparse, sigma_n_sparse);
            }
        }

        // Check if the cell is tagged with an endocardium region tag
        // 2 = DENSE region || 3 = SPARSE region
        int tag = DTI_MESH_FAST_ENDO(ac[i]);

        // Cell is inside the DENSE endocardium region
        if (tag == 2) {
            ac[i]->sigma.x = D_dense[0][0];
            ac[i]->sigma.y = D_dense[1][1];
            ac[i]->sigma.z = D_dense[2][2];

            ac[i]->sigma.xy = D_dense[0][1];
            ac[i]->sigma.xz = D_dense[0][2];
            ac[i]->sigma.yz = D_dense[1][2];
        }
        // Cell is inside the SPARSE endocardium region
        else if (tag == 3) {
            ac[i]->sigma.x = D_sparse[0][0];
            ac[i]->sigma.y = D_sparse[1][1];
            ac[i]->sigma.z = D_sparse[2][2];

            ac[i]->sigma.xy = D_sparse[0][1];
            ac[i]->sigma.xz = D_sparse[0][2];
            ac[i]->sigma.yz = D_sparse[1][2];
        }
        // Otherwise, NORMAL type of cell
        else {
            ac[i]->sigma.x = D[0][0];
            ac[i]->sigma.y = D[1][1];
            ac[i]->sigma.z = D[2][2];

            ac[i]->sigma.xy = D[0][1];
            ac[i]->sigma.xz = D[0][2];
            ac[i]->sigma.yz = D[1][2];
        }
    }

    OMP(parallel for)
    for(i = 0; i < num_active_cells; i++) {
        fill_discretization_matrix_elements_aniso(ac[i]);
    }

    free(f);
    free(s);
    free(n);

// [PURKINJE SETUP]
    static bool sigma_purkinje_initialized = false;

    uint32_t num_purkinje_active_cells = the_grid->purkinje->num_active_purkinje_cells;
    struct cell_node **ac_purkinje = the_grid->purkinje->purkinje_cells;
    struct node *pk_node = the_grid->purkinje->network->list_nodes;
    bool has_point_data = the_grid->purkinje->network->has_point_data;

    initialize_diagonal_elements_purkinje(the_solver, the_grid);

    if(!sigma_purkinje_initialized) {
        // Check if the Purkinje network file has the POINT_DATA section
        if (has_point_data) {
            struct node *tmp = the_grid->purkinje->network->list_nodes;
            uint32_t i = 0;
            while (tmp != NULL)
            {
                // Copy the prescribed conductivity from the Purkinje network file into the ALG Purkinje cell structure
                ac_purkinje[i]->sigma.x = tmp->sigma;

                tmp = tmp->next; i++;
            }
        }
        // Otherwise, initilize the conductivity of all cells homogenously with the value from the configuration file
        else {
            OMP(parallel for)
            for (uint32_t i = 0; i < num_active_cells; i++) {
                ac[i]->sigma.x = sigma_purkinje;
            }
        }
        sigma_purkinje_initialized = true;
    }
    fill_discretization_matrix_elements_purkinje(has_point_data,sigma_purkinje,ac_purkinje,num_purkinje_active_cells,pk_node);
}
