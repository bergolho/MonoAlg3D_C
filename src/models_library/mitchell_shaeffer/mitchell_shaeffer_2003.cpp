#include <cstdio>
#include <sycl/sycl.hpp>
#include <dpct/dpct.hpp>

#include "../../common_types/common_types.h"

#include "mitchell_shaeffer_2003.h"

#define BLOCK_SIZE 32
static dpct::constant_memory<size_t, 0> pitch;
size_t pitch_h;

// TODO: Move this struct to the "gpu_utils/sycl_utils.h" header
struct sycl_manager {
    sycl::queue queue;
    dpct::device_ext *dev_ct1;

    // TODO: Find a better way to set the device.
    // The user may want an specific GPU, maybe use the 'gpu_id' parameter as well.
    void set_queue (const bool use_gpu) {
        int device_id;
        // GPU
        if (use_gpu) {
            device_id = 0;
        }
        // CPU
        else {
            device_id = 1;
        }   
        dev_ct1 = &dpct::get_device(device_id);
        queue = dev_ct1->in_order_queue();
    }
};

// TODO: I think we will need one "struct sycl_manager" per cellular model.
//      For the linear system and ECG as well...
struct sycl_manager *sycl_config;

// KERNELS
// Set the initial condition for the ODE system
void SetInitialConditionsKernel (real *sv, const sycl::nd_item<3> &item_ct1, int num_volumes) {
    int threadID = item_ct1.get_global_id(2);
    if (threadID < num_volumes) {
        // Directly access as a flat array (no pitched memory)
        sv[threadID * NEQ] = 0.00000820413566106744f;  // V
        sv[threadID * NEQ + 1] = 0.8789655121804799f;  // h
    }
}

extern "C" SET_ODE_INITIAL_CONDITIONS_SYCL(set_model_initial_conditions_sycl)
{
    
    printf("Using Mitchell-Shaeffer 2003 SYCL model\n");
    sycl_config = (struct sycl_manager*)malloc(sizeof(struct sycl_manager));    // TODO: Remember to free me after the end of the main loop ...
    sycl_config->set_queue(solver->gpu);
    printf("Running on '%s'\n", sycl_config->queue.get_device().get_info<sycl::info::device::name>().c_str());

    uint32_t num_cells = solver->original_num_cells;
    solver->sv = malloc_device<real>(num_cells * NEQ, sycl_config->queue);
    real *sv_tmp = malloc_device<real>(num_cells * NEQ, sycl_config->queue);

    const int GRID = (num_cells + BLOCK_SIZE - 1) / BLOCK_SIZE;

    try {
        sycl_config->queue.submit([&](sycl::handler &cgh) {
        cgh.parallel_for(
            sycl::nd_range<3>(sycl::range<3>(1, 1, GRID * BLOCK_SIZE), sycl::range<3>(1, 1, BLOCK_SIZE)),
            [=](sycl::nd_item<3> item_ct1) {
                SetInitialConditionsKernel(sv_tmp, item_ct1, num_cells);
            });
        }).wait();
    }
    catch (sycl::exception &e) {
        printf("SYCL exception: %s\n", e.what());
    }

    sycl_config->queue.memcpy(solver->sv, sv_tmp, num_cells * NEQ * sizeof(real)).wait();

    sycl::free(sv_tmp, sycl_config->queue);
}

extern "C" SOLVE_MODEL_ODES(solve_model_odes_sycl) {

    uint32_t sv_id;

    size_t num_cells_to_solve = ode_solver->num_cells_to_solve;
    uint32_t * cells_to_solve = ode_solver->cells_to_solve;
    real *sv = ode_solver->sv;
    real dt = ode_solver->min_dt;
    uint32_t num_steps = ode_solver->num_steps;
    
    if (ode_solver->sv) {
        printf("sv is not NULL!\n");
    }
    else {
        printf("sv is NULL!\n");
    }
    exit(1);

    OMP(parallel for private(sv_id))
    for (uint32_t i = 0; i < num_cells_to_solve; i++)
    {
        if(cells_to_solve)
            sv_id = cells_to_solve[i];
        else
            sv_id = (uint32_t )i;

        for (int j = 0; j < num_steps; ++j)
        {
            solve_model_ode_sycl(dt, sv + (sv_id * NEQ), stim_currents[i]);
        }
    }

}

void solve_model_ode_sycl(real dt, real *sv, real stim_current)
{

    real rY[NEQ], rDY[NEQ];

    for(int i = 0; i < NEQ; i++)
        rY[i] = sv[i];

    RHS_sycl(rY, rDY, stim_current);

    for(int i = 0; i < NEQ; i++)
        sv[i] = dt*rDY[i] + rY[i];
}

// Remember to use: sycl::pow, sycl::exp, ...
void RHS_sycl(const real *sv, real *rDY_, real stim_current)
{

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
