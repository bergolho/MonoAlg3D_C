#ifndef MONOALG3D_MODEL_TORORD_FKATP_MIXED_ENDO_MID_EPI_H
#define MONOALG3D_MODEL_TORORD_FKATP_MIXED_ENDO_MID_EPI_H

// TOMEK, Jakub et al. Development, calibration, and validation of a novel human ventricular myocyte model in health, disease, and drug block. 
//  Elife, v. 8, p. e48890, 2019.

#include "../model_common.h"
#include "../../extra_data_library/helper_functions.h"

#define NEQ 43
#define INITIAL_V (-88.6369922306458)

#define ENDO 0.0
#define MID  1.0
#define EPI  2.0

// CPU macros
#define SOLVE_EQUATION_EULER_CPU(id) sv[id] = dt * rDY[id] + rY[id]

#define SOLVE_EQUATION_RUSH_LARSEN_CPU(id) sv[id] = (abs(a[id]) < TOLERANCE) ? rY[id] + dt * (rY[id] * a[id] + b[id]) : \
                                        exp(a[id] * dt)*(rY[id] + (b[id] / a[id])) - (b[id] / a[id] )

// GPU macros
#define SOLVE_EQUATION_EULER_GPU(id) *((real *)((char *)sv + pitch * id) + sv_id) = *((real *)((char *)sv + pitch * id) + sv_id) + dt * rDY[id]

#define SOLVE_EQUATION_RUSH_LARSEN_GPU(id) *((real *)((char *)sv + pitch * id) + sv_id) = (abs(a[id]) < TOLERANCE) ? *((real *)((char *)sv + pitch * id) + sv_id) + dt * ( *((real *)((char *)sv + pitch * id) + sv_id) * a[id] + b[id]) : \
                                                                                    exp(a[id] * dt)*(*((real *)((char *)sv + pitch * id) + sv_id) + (b[id] / a[id])) - (b[id] / a[id]) 

// SYCL macros
#define SOLVE_EQUATION_EULER_SYCL(id) sv[id * num_cells_to_solve + sv_id] += dt * rDY[id];

#define SOLVE_EQUATION_RUSH_LARSEN_SYCL(id) sv[id * num_cells_to_solve + sv_id] = (abs(a[id]) < TOLERANCE) ? sv[id * num_cells_to_solve + sv_id] + dt * ( sv[id * num_cells_to_solve + sv_id] * a[id] + b[id]) : \
                                                                                    exp(a[id] * dt)*(sv[id * num_cells_to_solve + sv_id] + (b[id] / a[id])) - (b[id] / a[id])

#ifdef __CUDACC__

#include "../../gpu_utils/gpu_utils.h"

__global__ void kernel_set_model_initial_conditions(real *sv, int num_volumes, size_t pitch, bool use_adpt_dt, real min_dt);
__global__ void kernel_set_model_initial_conditions_endo_mid_epi(real *sv, int num_volumes, size_t pitch, bool use_adpt_dt, real min_dt,\
                                                real *initial_endo, real *initial_epi, real *initial_mid, real *transmurality, real *sf_Iks);

__global__ void solve_endo_mid_epi_gpu(real cur_time, real dt, real *sv, real *stim_currents, uint32_t *cells_to_solve, real *transmurality, real *sf_Iks, real *extra_params,\
                          uint32_t num_cells_to_solve, int num_steps, size_t pitch, bool use_adpt, real abstol, real reltol, real max_dt, bool has_extra_params);

inline __device__ void RHS_gpu(real *sv, real *rDY_, real stim_current, real transmurality, real sf_Iks, real *extra_params, int threadID_, real dt, size_t pitch, bool use_adpt_dt);
inline __device__ void RHS_RL_gpu(real *a_, real *b_, real *sv, real *rDY_, real stim_current, real transmurality, real sf_Iks, real *extra_params, int threadID_, real dt, size_t pitch, bool use_adpt_dt);
inline __device__ void solve_forward_euler_gpu_adpt(real *sv, real stim_curr, real transmurality, real sf_Iks, real *extra_params, real final_time, int thread_id, size_t pitch, real abstol, real reltol, real min_dt, real max_dt);

#endif

#ifdef COMPILE_SYCL

void kernel_set_model_initial_conditions(real *sv, int num_volumes, bool adaptive, real min_dt);
void kernel_set_model_initial_conditions_endo_mid_epi(real *sv, int num_volumes, bool adaptive, real min_dt,\
                                                            real *initial_endo, real *initial_epi, real *initial_mid,\
                                                            real *transmurality, real *sf_Iks);
inline void RHS_sycl (real *sv, real *rDY_, real stim_current, real transmurality, real sf_Iks, real *extra_params, real dt, int sv_id, int num_volumes, bool use_adpt_dt);
inline void RHS_RL_sycl (real *a_, real *b_, real *sv, real *rDY_, real stim_current, real transmurality, real sf_Iks, real *extra_params, real dt, int sv_id, int num_volumes, bool use_adpt_dt);
inline void solve_forward_euler_sycl_adpt(real *sv, real stim_current, real transmurality, real sf_Iks, real *extra_params, real final_time, int sv_id, int num_volumes, real abstol, real reltol, real min_dt, real max_dt);

#endif

void RHS_cpu(const real *sv, real *rDY_, real stim_current, real dt, real transmurality, real sf_Iks, real const *extra_params);
void RHS_RL_cpu (real *a_, real *b_, const real *sv, real *rDY_, real stim_current, real dt, real transmurality, real sf_Iks, real const *extra_params);
void solve_forward_euler_cpu_adpt(real *sv, real stim_curr, real transmurality, real sf_Iks, real final_time, int sv_id, struct ode_solver *solver, real const *extra_params);
void solve_model_ode_cpu(real dt, real *sv, real stim_current, real transmurality, real sf_Iks, real const *extra_params);

#endif //MONOALG3D_MODEL_TORORD_FKATP_MIXED_ENDO_MID_EPI_H

