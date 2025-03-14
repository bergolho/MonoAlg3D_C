#include "../../gpu_utils/sycl_utils.h"
#include <cstdio>
#include <sycl/sycl.hpp>
#include <dpct/dpct.hpp>

#include "../../common_types/common_types.h"

#include "ten_tusscher_3_RS.h"

extern "C" SET_ODE_INITIAL_CONDITIONS_SYCL(set_model_initial_conditions_sycl)
{
    char *cell_type;
    #ifdef ENDO
    cell_type = strdup("ENDO");
    #endif

    #ifdef EPI
    cell_type = strdup("EPI");
    #endif

    #ifdef MCELL
    cell_type = strdup("MCELL");
    #endif

    printf("Using TenTusscher3 With Extra Parameters 2004 SYCL model\n");
    printf("TT3 model configured with %s cell\n", cell_type);
    dpct::device_ext &dev_ct1 = dpct::get_current_device();
    sycl::queue &q_ct1 = dev_ct1.default_queue();
    printf("Running on '%s'\n", q_ct1.get_device().get_info<sycl::info::device::name>().c_str());

    uint32_t num_cells = solver->original_num_cells;
    
    solver->sv = sycl::malloc_device<real>(num_cells*NEQ, q_ct1);
    
    if (solver->sv) {
        const int BLOCK_SIZE = 32;
        const int GRID = (num_cells + BLOCK_SIZE - 1) / BLOCK_SIZE;
        const int TOTAL_THREADS = BLOCK_SIZE*GRID;

        try {
            q_ct1.submit([&](sycl::handler& h) {
                real *sv = solver->sv;
                h.parallel_for(
                    sycl::nd_range<1>(sycl::range<1>(TOTAL_THREADS), sycl::range<1>(BLOCK_SIZE)),
                    [=](sycl::nd_item<1> item) {
                        int i = item.get_global_id(0);
                        if (i < num_cells) {
                            sv[0*num_cells+i] = -86.2f;   // V;       millivolt
                            sv[1*num_cells+i] = 0.0f; //M
                            sv[2*num_cells+i] = 0.75; //H
                            sv[3*num_cells+i] = 0.75; //J
                            sv[4*num_cells+i] = 0.0f; //Xr1
                            sv[5*num_cells+i] = 0.0f; //Xs
                            sv[6*num_cells+i] = 1.0f; //S
                            sv[7*num_cells+i] = 1.0f; //F
                            sv[8*num_cells+i] = 1.0f; //F2
                            sv[9*num_cells+i] = 0.0; //D_INF
                            sv[10*num_cells+i] = 0.0; //R_INF
                            sv[11*num_cells+i] = 0.0; //Xr2_INF
                        }
                    }
                );
            }).wait();
        } catch (sycl::exception &e) {
            printf("SYCL exception: %s\n", e.what());
        }
    }
}

extern "C" SOLVE_MODEL_ODES(solve_model_odes_sycl) {

    size_t num_cells_to_solve = ode_solver->num_cells_to_solve;
    uint32_t * cells_to_solve = ode_solver->cells_to_solve;
    real *sv = ode_solver->sv;
    real dt = ode_solver->min_dt;
    uint32_t num_steps = ode_solver->num_steps;

    dpct::device_ext &dev_ct1 = dpct::get_current_device();
    sycl::queue &q_ct1 = dev_ct1.default_queue();

    // Using Unifed Shared Memory (USM) instead of access buffer to improve performance
    real *d_stim = sycl::malloc_device<real>(num_cells_to_solve, q_ct1);
    uint32_t *d_cells_to_solve = NULL;
    if (cells_to_solve) {
        d_cells_to_solve = sycl::malloc_device<uint32_t>(num_cells_to_solve, q_ct1);
        q_ct1.memcpy(d_cells_to_solve, cells_to_solve, num_cells_to_solve*sizeof(uint32_t)).wait();
    }
    
    // Copy initial data to device
    q_ct1.memcpy(d_stim, stim_currents, num_cells_to_solve*sizeof(real)).wait();

    // Copy the extra data from host memory
    real *fibs = NULL;
    int num_extra_parameters = 8;
    real extra_par[num_extra_parameters];
    size_t extra_parameters_size = num_extra_parameters * sizeof(real);

    real *d_fibs;
    real *d_extra_par;
    real fibs_size = num_cells_to_solve*sizeof(real);

    struct extra_data_for_fibrosis* extra_data_from_cpu = (struct extra_data_for_fibrosis*)ode_solver->ode_extra_data;

    bool deallocate = false;
    if(ode_solver->ode_extra_data) {
        fibs = extra_data_from_cpu->fibrosis;
        extra_par[0] = extra_data_from_cpu->atpi;
        extra_par[1] = extra_data_from_cpu->Ko;
        extra_par[2] = extra_data_from_cpu->Ki;
        extra_par[3] = extra_data_from_cpu->Vm_modifier;
        extra_par[4] = extra_data_from_cpu->GNa_multiplicator;
        extra_par[5] = extra_data_from_cpu->GCaL_multiplicator;
        extra_par[6] = extra_data_from_cpu->INaCa_multiplicator;
        extra_par[7] = extra_data_from_cpu->Ikatp_multiplicator;
    }
    else {
        extra_par[0] = 6.8f;
        extra_par[1] = 5.4f;
        extra_par[2] = 138.3f;
        extra_par[3] = 0.0;
        extra_par[4] = 1.0f;
        extra_par[5] = 1.0f;
        extra_par[6] = 1.0f;
        extra_par[7] = 1.0f;

        fibs = (real*)malloc(fibs_size);

		for(uint64_t i = 0; i < num_cells_to_solve; i++) {
			fibs[i] = 1.0;
		}
        deallocate = true;
    }

    d_extra_par = sycl::malloc_device<real>(num_extra_parameters, q_ct1);
    q_ct1.memcpy(d_extra_par, extra_par, num_extra_parameters*sizeof(real)).wait();

    d_fibs = sycl::malloc_device<real>(num_cells_to_solve, q_ct1);
    q_ct1.memcpy(d_fibs, fibs, num_cells_to_solve*sizeof(real)).wait();

    // Define block and grid sizes
    const int BLOCK_SIZE = 32; 
    const int GRID = (num_cells_to_solve + BLOCK_SIZE - 1) / BLOCK_SIZE;

    try {
        q_ct1.submit([&](sycl::handler& hand) {
                // num_cells_to_solve -> 'i'
                hand.parallel_for(
                sycl::nd_range<1>(sycl::range<1>(GRID * BLOCK_SIZE), sycl::range<1>(BLOCK_SIZE)),
                [=](sycl::nd_item<1> item) [[sycl::reqd_work_group_size(BLOCK_SIZE)]] {
                    
                    size_t i = item.get_global_id(0);
                    if (i >= num_cells_to_solve) return;

                    int sv_id = (d_cells_to_solve) ? d_cells_to_solve[i] : i;

                    // Private memory for fast access
                    real rDY[NEQ];

                    for (int n = 0; n < num_steps; ++n) {
                        RHS_sycl(sv, d_stim[i], rDY, sv_id, dt, num_cells_to_solve, d_fibs[i], d_extra_par);

                        // Update state variables
                        sv[i] += dt * rDY[0];
                        for (int j = 1; j < NEQ; j++) {
                            sv[j * num_cells_to_solve + i] = rDY[j];
                        }
                    }
                });
        }).wait();

        // Free device memory
        sycl::free(d_stim, q_ct1);
        sycl::free(d_extra_par, q_ct1);
        sycl::free(d_fibs, q_ct1);
        if (d_cells_to_solve) sycl::free(d_cells_to_solve, q_ct1);
        if (deallocate) free(fibs);

    } catch (sycl::exception &e) {
        printf("SYCL exception: %s\n", e.what());
    }
}

inline void RHS_sycl(real *sv, real stim_current, real *rDY_, int sv_id, real dt, int num_cells, real fibrosis, real *extra_parameters) {
    //fibrosis = 0 means that the cell is fibrotic, 1 is not fibrotic. Anything between 0 and 1 means border zone
    real svolt      = sv[0 * num_cells + sv_id];
    real sm         = sv[1 * num_cells + sv_id];
    real sh         = sv[2 * num_cells + sv_id];
    real sj         = sv[3 * num_cells + sv_id];
    real sxr1       = sv[4 * num_cells + sv_id];
    real sxs        = sv[5 * num_cells + sv_id];
    real ss         = sv[6 * num_cells + sv_id];
    real sf         = sv[7 * num_cells + sv_id];
    real sf2        = sv[8 * num_cells + sv_id];
    real D_INF      = sv[9 * num_cells + sv_id];
    real R_INF      = sv[10 * num_cells + sv_id];
    real Xr2_INF    = sv[11 * num_cells + sv_id];

    #include "ten_tusscher_3_RS_common.inc"
}
