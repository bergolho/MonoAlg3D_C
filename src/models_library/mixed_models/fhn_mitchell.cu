#include <stddef.h>
#include "../../monodomain/constants.h"
#include "../model_gpu_utils.h"

#include "fhn_mitchell.h"

extern "C" SET_ODE_INITIAL_CONDITIONS_GPU(set_model_initial_conditions_gpu) {

    log_to_stdout_and_file("Using mixed version of modified FHN 1961 + Mitchell-Shaeffer 2003 GPU model\n");

    uint32_t num_volumes = solver->original_num_cells;

    // execution configuration
    const int GRID  = (num_volumes + BLOCK_SIZE - 1)/BLOCK_SIZE;

    size_t size = num_volumes*sizeof(real);
    size_t extra_data_bytes_size = num_volumes*sizeof(uint32_t);

    // TODO: Think what to do when the number of equations are different between the cellular models ...
    check_cuda_error(cudaMallocPitch((void **) &(solver->sv), &pitch_h, size, (size_t )NEQ_1));
    check_cuda_error(cudaMemcpyToSymbol(pitch, &pitch_h, sizeof(size_t)));

    // Get the mapping array
    uint32_t *mapping = NULL;
    uint32_t *mapping_device = NULL;
    if(solver->ode_extra_data) {
        mapping = (uint32_t*)(solver->ode_extra_data);
        check_cuda_error(cudaMalloc((void **)&mapping_device, extra_data_bytes_size));
        check_cuda_error(cudaMemcpy(mapping_device, mapping, extra_data_bytes_size, cudaMemcpyHostToDevice));
    }
    else {
        log_to_stderr_and_file_and_exit("You need to specify a mask function when using a mixed model!\n");
    }

    kernel_set_model_inital_conditions <<<GRID, BLOCK_SIZE>>>(solver->sv, NULL, mapping_device, num_volumes);

    check_cuda_error( cudaPeekAtLastError() );
    cudaDeviceSynchronize();
    return pitch_h;

}

extern "C" SOLVE_MODEL_ODES(solve_model_odes_gpu) {

    size_t num_cells_to_solve = ode_solver->num_cells_to_solve;
    uint32_t * cells_to_solve = ode_solver->cells_to_solve;
    real *sv = ode_solver->sv;
    real dt = ode_solver->min_dt;
    uint32_t num_steps = ode_solver->num_steps;

    // execution configuration
    const int GRID  = ((int)num_cells_to_solve + BLOCK_SIZE - 1)/BLOCK_SIZE;

    size_t stim_currents_size = sizeof(real)*num_cells_to_solve;
    size_t cells_to_solve_size = sizeof(uint32_t)*num_cells_to_solve;

    real *stims_currents_device;
    check_cuda_error(cudaMalloc((void **) &stims_currents_device, stim_currents_size));
    check_cuda_error(cudaMemcpy(stims_currents_device, stim_currents, stim_currents_size, cudaMemcpyHostToDevice));

    uint32_t *cells_to_solve_device = NULL;
    if(cells_to_solve != NULL) {
        check_cuda_error(cudaMalloc((void **) &cells_to_solve_device, cells_to_solve_size));
        check_cuda_error(cudaMemcpy(cells_to_solve_device, cells_to_solve, cells_to_solve_size, cudaMemcpyHostToDevice));
    }

    // Get the mapping array
    uint32_t *mapping = NULL, *mapping_device = NULL;
    if (ode_solver->ode_extra_data) {
        mapping = (uint32_t*)(ode_solver->ode_extra_data);
        check_cuda_error(cudaMalloc((void **)&mapping_device, cells_to_solve_size));
        check_cuda_error(cudaMemcpy(mapping_device, mapping, cells_to_solve_size, cudaMemcpyHostToDevice));
    }
    else {
        log_to_stderr_and_file_and_exit("You need to specify a mask function when using a mixed model!\n");
    }

    solve_gpu<<<GRID, BLOCK_SIZE>>>(dt, sv, stims_currents_device, cells_to_solve_device, num_cells_to_solve, num_steps, mapping_device);

    check_cuda_error( cudaPeekAtLastError() );

    check_cuda_error(cudaFree(stims_currents_device));
    if(cells_to_solve_device) check_cuda_error(cudaFree(cells_to_solve_device));
    if(mapping_device) check_cuda_error(cudaFree(mapping_device));

}

__global__ void kernel_set_model_inital_conditions(real *sv, real *IC, uint32_t *mapping, int num_volumes) {

    int threadID = blockDim.x * blockIdx.x + threadIdx.x;

    if (threadID < num_volumes) {
        
        // FHN
        if (mapping[threadID] == 0) {
            *((real * )((char *) sv + pitch * 0) + threadID) = 0.0f;   // Vm millivolt
            *((real * )((char *) sv + pitch * 1) + threadID) = 0.0f;   // h dimensionless
        }
        // Mitchell
        else if (mapping[threadID] == 1) {
            *((real * )((char *) sv + pitch * 0) + threadID) = 0.00000820413566106744f; // Vm millivolt
            *((real * )((char *) sv + pitch * 1) + threadID) = 0.8789655121804799f;     // h dimensionless
        }
    }
}

// Solving the model for each cell in the tissue matrix ni x nj
__global__ void solve_gpu(real dt, real *sv, real* stim_currents,
                          uint32_t *cells_to_solve, uint32_t num_cells_to_solve,
                          int num_steps,
                          uint32_t *mapping) {

    int threadID = blockDim.x * blockIdx.x + threadIdx.x;
    int sv_id;

    // Each thread solves one cell 
    if(threadID < num_cells_to_solve) {

        if(cells_to_solve)
            sv_id = cells_to_solve[threadID];
        else
            sv_id = threadID;

        real rDY[NEQ_1];

        for (int n = 0; n < num_steps; ++n) {

            // FHN
            if (mapping[sv_id] == 0) {

                RHS_gpu_fhn(sv, rDY, stim_currents[threadID], sv_id);

                for(int i = 0; i < NEQ_1; i++) {
                    *((real *) ((char *) sv + pitch * i) + sv_id) = dt * rDY[i] + *((real *) ((char *) sv + pitch * i) + sv_id);
                }
            }
            // Mitchell
            else if (mapping[sv_id] == 1) {

                RHS_gpu_mitchell(sv, rDY, stim_currents[threadID], sv_id);

                for (int i = 0; i < NEQ_2; i++) {
                    *((real *) ((char *) sv + pitch * i) + sv_id) = dt * rDY[i] + *((real *) ((char *) sv + pitch * i) + sv_id);
                }
            }

        }
    }
}

inline __device__ void RHS_gpu_fhn (real *sv_, real *rDY_, real stim_current, int threadID_) 
{

    //State variables
    const real u = *((real*)((char*)sv_ + pitch * 0) + threadID_);
    const real v = *((real*)((char*)sv_ + pitch * 1) + threadID_);

    const real a = 0.2f;
    const real b = 0.5f;
    const real k = 36.0;
    const real epsilon  =  0.000150;

    rDY_[0] = k*(u*(1.0f - u)*(u - a) - u*v) + stim_current;
    rDY_[1] = k*epsilon*(b*u - v);

}

inline __device__ void RHS_gpu_mitchell(real *sv_, real *rDY_, real stim_current, int threadID_)
{
    //State variables
    const real V = *((real*)((char*)sv_ + pitch * 0) + threadID_);
    const real h = *((real*)((char*)sv_ + pitch * 1) + threadID_);

    // Constants
    const real tau_in = 0.3;
    const real tau_out = 6.0;
    const real V_gate = 0.13;
    const real tau_open = 120.0;
    const real tau_close = 150.0;

    // Algebraics
    real J_stim = stim_current;
    real J_in = ( h*( powf(V, 2.00000)*(1.00000 - V)))/tau_in;
    real J_out = - (V/tau_out);

    // Rates
    rDY_[0] = J_out + J_in + J_stim;
    rDY_[1] = (V < V_gate ? (1.00000 - h)/tau_open : - h/tau_close);
}
