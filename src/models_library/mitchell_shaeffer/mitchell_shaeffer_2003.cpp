#include "../../gpu_utils/sycl_utils.h"
#include <cstdio>
#include <sycl/sycl.hpp>
#include <dpct/dpct.hpp>

#include "../../common_types/common_types.h"

#include "mitchell_shaeffer_2003.h"

// KERNELS
// Set the initial condition for the ODE system
//void SetInitialConditionsKernel (real *sv, const sycl::nd_item<3> &item_ct1, int num_volumes) {
//    int threadID = item_ct1.get_global_id(2);
//    if (threadID < num_volumes) {
//        // Directly access as a flat array (no pitched memory)
//        sv[threadID * NEQ] = 0.00000820413566106744f;  // V
//        sv[threadID * NEQ + 1] = 0.8789655121804799f;  // h
//    }
//}

extern "C" SET_ODE_INITIAL_CONDITIONS_SYCL(set_model_initial_conditions_sycl)
{
    printf("Using Mitchell-Shaeffer 2003 SYCL model\n");
    dpct::device_ext &dev_ct1 = dpct::get_current_device();
    sycl::queue &q_ct1 = dev_ct1.default_queue();
    printf("Running on '%s'\n", q_ct1.get_device().get_info<sycl::info::device::name>().c_str());

    uint32_t num_cells = solver->original_num_cells;
    
    // TODO: Try to allocate the 'sv' array using 'malloc_device()'
    solver->sv = new real[num_cells*NEQ]();
    sycl::buffer<real, 2> sv_buf(solver->sv, sycl::range<2>(num_cells,NEQ));

    if (solver->sv) {
        const int BLOCK_SIZE = 32;
        const int GRID = (num_cells + BLOCK_SIZE - 1) / BLOCK_SIZE;

        try {
            q_ct1.submit([&](sycl::handler& h) {
                auto sv = sv_buf.get_access<sycl::access::mode::write>(h);
                h.parallel_for(sycl::range<1>(num_cells), [=](sycl::id<1> i) {
                    sv[i][0] = 0.00000820413566106744;      // V
                    sv[i][1] = 0.8789655121804799;          // h
                });
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

    sycl::buffer<real, 2> sv_buf(ode_solver->sv, sycl::range<2>(num_cells_to_solve, NEQ));
    sycl::buffer<real, 1> stim_buf(stim_currents, sycl::range<1>(num_cells_to_solve));

    //const int BLOCK_SIZE = 32;
    //const int GRID = (num_cells_to_solve + BLOCK_SIZE - 1) / BLOCK_SIZE;

    try {
        q_ct1.submit([&](sycl::handler& hand) {
                auto sv_sycl = sv_buf.get_access<sycl::access::mode::read_write>(hand);
                auto stim = stim_buf.get_access<sycl::access::mode::read>(hand);

                // num_cells_to_solve -> 'i'
                hand.parallel_for(sycl::range<1>(num_cells_to_solve), [=](sycl::id<1> i) {
                    
                    int sv_id;
                    if(cells_to_solve)
                        sv_id = cells_to_solve[i];
                    else
                        sv_id = i;

                    real rDY[NEQ];

                    for (int n = 0; n < num_steps; ++n) {

                        //State variables
                        const real V = sv_sycl[i][0];
                        const real h = sv_sycl[i][1];
                        const real stim_current = stim[i];

                        // Constants
                        const real tau_in = 0.3;
                        const real tau_out = 6.0;
                        const real V_gate = 0.13;
                        const real tau_open = 120.0;
                        const real tau_close = 150.0;

                        // Algebraics
                        real J_stim = stim_current;
                        real J_in = ( h*( sycl::pow(V, 2.00000)*(1.00000 - V)))/tau_in;
                        real J_out = - (V/tau_out);

                        // Rates
                        rDY[0] = J_out + J_in + J_stim;
                        rDY[1] = (V < V_gate ? (1.00000 - h)/tau_open : - h/tau_close);
                        
                        // neq -> 'j'
                        for(int j = 0; j < NEQ; j++) {
                            sv_sycl[i][j] = dt * rDY[j] + sv_sycl[i][j]; 
                        }            
                    }
                });
        }).wait();
    } catch (sycl::exception &e) {
        printf("SYCL exception: %s\n", e.what());
    }
}

// Remember to use: sycl::pow, sycl::exp, ...
void RHS_sycl(real *sv, real *rDY_, real stim_current) {

    //State variables
    const real V = sv[0];
    const real h = sv[1];

    // Constants
    const real tau_in = 0.3;
    const real tau_out = 6.0;
    const real V_gate = 0.13;
    const real tau_open = 120.0;
    const real tau_close = 150.0;

    // Algebraics
    real J_stim = stim_current;
    real J_in = ( h*( sycl::pow(V, 2.00000)*(1.00000 - V)))/tau_in;
    real J_out = - (V/tau_out);

    // Rates
    rDY_[0] = J_out + J_in + J_stim;
    rDY_[1] = (V < V_gate ? (1.00000 - h)/tau_open : - h/tau_close);

}
