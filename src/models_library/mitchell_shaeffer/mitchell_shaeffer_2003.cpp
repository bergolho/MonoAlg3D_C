#include "../../gpu_utils/sycl_utils.h"
#include <cstdio>
#include <sycl/sycl.hpp>
#include <dpct/dpct.hpp>

#include "../../common_types/common_types.h"

#include "mitchell_shaeffer_2003.h"


extern "C" SET_ODE_INITIAL_CONDITIONS_SYCL(set_model_initial_conditions_sycl)
{
    printf("Using Mitchell-Shaeffer 2003 SYCL model\n");
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
                            sv[0*num_cells+i] = 0.0;  // V
                            sv[1*num_cells+i] = 0.88; // h
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

    real *d_stim = sycl::malloc_device<real>(num_cells_to_solve, q_ct1);
    uint32_t *d_cells_to_solve = NULL;
    if (cells_to_solve) {
        d_cells_to_solve = sycl::malloc_device<uint32_t>(num_cells_to_solve, q_ct1);
        q_ct1.memcpy(d_cells_to_solve, cells_to_solve, num_cells_to_solve*sizeof(uint32_t));
    }
    
    // Copy initial data to device
    q_ct1.memcpy(d_stim, stim_currents, num_cells_to_solve*sizeof(real)).wait();

    // Define block and grid sizes
    const int BLOCK_SIZE = 32;
    const int GRID = (num_cells_to_solve + BLOCK_SIZE - 1) / BLOCK_SIZE;

    try {
        q_ct1.submit([&](sycl::handler& hand) {
                // num_cells_to_solve -> 'i'
                hand.parallel_for(
                sycl::nd_range<1>(sycl::range<1>(GRID * BLOCK_SIZE), sycl::range<1>(BLOCK_SIZE)),
                [=](sycl::nd_item<1> item) [[sycl::reqd_work_group_size(32)]] {
                    
                    size_t i = item.get_global_id(0);
                    if (i >= num_cells_to_solve) return;

                    int sv_id = (d_cells_to_solve) ? d_cells_to_solve[i] : i;

                    // Private memory for fast access
                    real rDY[NEQ];

                    for (int n = 0; n < num_steps; ++n) {
                        RHS_sycl(sv, d_stim[i], rDY, sv_id, num_cells_to_solve);

                        // Update state variables
                        #pragma unroll
                        for (int j = 0; j < NEQ; j++) {
                            sv[j * num_cells_to_solve + i] += dt * rDY[j];
                        }
                    }
                });
        }).wait();

        // Free device memory
        sycl::free(d_stim, q_ct1);
        if (d_cells_to_solve) sycl::free(d_cells_to_solve, q_ct1);

    } catch (sycl::exception &e) {
        printf("SYCL exception: %s\n", e.what());
    }
}

inline void RHS_sycl(real *Y, real stim_current, real *dY, int sv_id, int num_cells) {

    // Load state variables
    real V = Y[0 * num_cells + sv_id];
    real h = Y[1 * num_cells + sv_id];
    //real stim_current = d_stim[i];

    // Constants (computed outside the loop for efficiency)
    constexpr real tau_in = 0.3;
    constexpr real tau_out = 6.0;
    constexpr real V_gate = 0.13;
    constexpr real tau_open = 120.0;
    constexpr real tau_close = 150.0;

    // Algebraics
    real J_stim = stim_current;
    real J_in = (h * (pow(V, 2.00000) * (1.00000 - V))) / tau_in;
    real J_out = - (V / tau_out);

    // Compute rates
    dY[0] = J_out + J_in + J_stim;
    dY[1] = (V < V_gate) ? (1.00000 - h) / tau_open : -h / tau_close;

}
