#include "../../gpu_utils/sycl_utils.h"
#include <cstdio>
#include <sycl/sycl.hpp>
#include <dpct/dpct.hpp>

#include "../../common_types/common_types.h"

#include "trovato_2020.h"

extern "C" SET_ODE_INITIAL_CONDITIONS_SYCL(set_model_initial_conditions_sycl)
{
    log_info("Using Trovato_2020 SYCL model\n");
    dpct::device_ext &dev_ct1 = dpct::get_current_device();
    sycl::queue &q_ct1 = dev_ct1.default_queue();
    log_info("Running on '%s'\n", q_ct1.get_device().get_info<sycl::info::device::name>().c_str());

    uint32_t num_volumes = solver->original_num_cells;

    // Allocate memory for state vector
    bool adaptive = solver->adaptive;
    if (adaptive) {
        log_error_and_exit("[-] Adaptive timestep is not implemented for this model yet!\n");
        //log_info("Using Adaptive timestep to solve the ODEs\n");
        //solver->sv = sycl::malloc_device<real>(num_volumes*(NEQ+3), q_ct1);
    }
    else {
        log_info("Using Fixed timestep to solve the ODEs\n");
        solver->sv = sycl::malloc_device<real>(num_volumes*NEQ, q_ct1);
    }

    kernel_set_model_initial_conditions(solver->sv, num_volumes, adaptive, solver->min_dt);
    
}

extern "C" SOLVE_MODEL_ODES(solve_model_odes_sycl) {

    const real TOLERANCE = 1e-8;
    size_t num_cells_to_solve = ode_solver->num_cells_to_solve;
    uint32_t * cells_to_solve = ode_solver->cells_to_solve;
    real *sv = ode_solver->sv;
    real dt = ode_solver->min_dt;
    real cur_time = current_t;
    real max_dt = ode_solver->max_dt;
    real abstol = ode_solver->abs_tol;
    real reltol = ode_solver->rel_tol;
    uint32_t num_steps = ode_solver->num_steps;
    bool has_extra_params = (ode_solver->ode_extra_data) ? true : false;
    bool adaptive = ode_solver->adaptive;

    dpct::device_ext &dev_ct1 = dpct::get_current_device();
    sycl::queue &q_ct1 = dev_ct1.default_queue();

    real *d_stim = sycl::malloc_device<real>(num_cells_to_solve, q_ct1);
    uint32_t *d_cells_to_solve = NULL;
    if (cells_to_solve) {
        d_cells_to_solve = sycl::malloc_device<uint32_t>(num_cells_to_solve, q_ct1);
        q_ct1.memcpy(d_cells_to_solve, cells_to_solve, num_cells_to_solve*sizeof(uint32_t)).wait();
    }
    q_ct1.memcpy(d_stim, stim_currents, num_cells_to_solve*sizeof(real)).wait();

    // Explicitly copy host memory to device
    int num_extra_parameters = 29;
    real extra_par[num_extra_parameters];
    real *d_extra_par = NULL;

    // Get the extra data array if exists
    if(has_extra_params) {
        struct extra_data_for_trovato *extra_data = (struct extra_data_for_trovato*)ode_solver->ode_extra_data;
        extra_par[0]     = extra_data->GNa_Multiplier;
        extra_par[1]     = extra_data->GNaL_Multiplier;
        extra_par[2]     = extra_data->GCaT_Multiplier;
        extra_par[3]     = extra_data->Gto_Multiplier;
        extra_par[4]     = extra_data->Gsus_Multiplier;
        extra_par[5]     = extra_data->Gkr_Multiplier;
        extra_par[6]     = extra_data->Gks_Multiplier;
        extra_par[7]     = extra_data->GfNa_Multiplier;
        extra_par[8]     = extra_data->GfK_Multiplier;
        extra_par[9]     = extra_data->GK1_Multiplier;
        extra_par[10]    = extra_data->GNCX_Multiplier;
        extra_par[11]    = extra_data->GNaK_Multiplier;
        extra_par[12]    = extra_data->INa_Multiplier;
        extra_par[13]    = extra_data->ICaL_Multiplier;
        extra_par[14]    = extra_data->ICaNa_Multiplier;
        extra_par[15]    = extra_data->ICaK_Multiplier;
        extra_par[16]    = extra_data->Ito_Multiplier;
        extra_par[17]    = extra_data->INaL_Multiplier;
        extra_par[18]    = extra_data->IKr_Multiplier;
        extra_par[19]    = extra_data->IKs_Multiplier;
        extra_par[20]    = extra_data->IK1_Multiplier;
        extra_par[21]    = extra_data->INaCa_Multiplier;
        extra_par[22]    = extra_data->INaK_Multiplier;
        extra_par[23]    = extra_data->INab_Multiplier;
        extra_par[24]    = extra_data->ICab_Multiplier;
        extra_par[25]    = extra_data->ICaT_Multiplier;
        extra_par[26]    = extra_data->Isus_Multiplier;
        extra_par[27]    = extra_data->If_Multiplier;
        extra_par[28]    = extra_data->IpCa_Multiplier;

        d_extra_par = sycl::malloc_device<real>(num_extra_parameters, q_ct1);
        q_ct1.memcpy(d_extra_par, extra_par, num_extra_parameters*sizeof(real)).wait();
    }
    // No [extra_data] section, we consider all cells ENDO!
    else {

        // Default: initialize all current modifiers
        for (int i = 0; i < num_extra_parameters; i++) {
            extra_par[i] = 1.0;
        }

        d_extra_par = sycl::malloc_device<real>(num_extra_parameters, q_ct1);
        q_ct1.memcpy(d_extra_par, extra_par, num_extra_parameters*sizeof(real)).wait();
    }

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

                    // Fixed timestep
                    if (!adaptive) {
                        // Private memory for fast access
                        real rDY[NEQ];
                        real a[NEQ], b[NEQ];

                        for(int n = 0; n < num_steps; ++n) {

                            RHS_RL_sycl(a, b, sv, rDY, d_stim[i], d_extra_par, dt, sv_id, num_cells_to_solve, false);

                            // Solve variables based on its type:
                            //  Non-linear = Euler
                            //  Hodkin-Huxley = Rush-Larsen || Euler (if 'a' coefficient is too small)
                            SOLVE_EQUATION_EULER_SYCL(0);        // v        
                            SOLVE_EQUATION_EULER_SYCL(1);        // CaMKt    
                            SOLVE_EQUATION_EULER_SYCL(2);        // cass 
                            SOLVE_EQUATION_EULER_SYCL(3);        // nai  
                            SOLVE_EQUATION_EULER_SYCL(4);        // nasl 
                            SOLVE_EQUATION_EULER_SYCL(5);        // nass 
                            SOLVE_EQUATION_EULER_SYCL(6);        // ki 
                            SOLVE_EQUATION_EULER_SYCL(7);        // kss
                            SOLVE_EQUATION_EULER_SYCL(8);        // ksl
                            SOLVE_EQUATION_EULER_SYCL(9);        // cai
                            SOLVE_EQUATION_EULER_SYCL(10);       // casl
                            SOLVE_EQUATION_EULER_SYCL(11);       // cansr
                            SOLVE_EQUATION_EULER_SYCL(12);       // cajsr
                            SOLVE_EQUATION_EULER_SYCL(13);       // cacsr
                            SOLVE_EQUATION_EULER_SYCL(14);       // Jrel1
                            SOLVE_EQUATION_EULER_SYCL(15);       // Jrel2
                            SOLVE_EQUATION_RUSH_LARSEN_SYCL(16); // m
                            SOLVE_EQUATION_RUSH_LARSEN_SYCL(17); // hf
                            SOLVE_EQUATION_RUSH_LARSEN_SYCL(18); // hs
                            SOLVE_EQUATION_RUSH_LARSEN_SYCL(19); // j
                            SOLVE_EQUATION_RUSH_LARSEN_SYCL(20); // hsp
                            SOLVE_EQUATION_RUSH_LARSEN_SYCL(21); // jp
                            SOLVE_EQUATION_RUSH_LARSEN_SYCL(22); // mL
                            SOLVE_EQUATION_RUSH_LARSEN_SYCL(23); // hL
                            SOLVE_EQUATION_RUSH_LARSEN_SYCL(24); // hLp
                            SOLVE_EQUATION_RUSH_LARSEN_SYCL(25); // a
                            SOLVE_EQUATION_RUSH_LARSEN_SYCL(26); // i1
                            SOLVE_EQUATION_RUSH_LARSEN_SYCL(27); // i2
                            SOLVE_EQUATION_RUSH_LARSEN_SYCL(28); // d
                            SOLVE_EQUATION_RUSH_LARSEN_SYCL(29); // ff
                            SOLVE_EQUATION_RUSH_LARSEN_SYCL(30); // fs
                            SOLVE_EQUATION_RUSH_LARSEN_SYCL(31); // fcaf
                            SOLVE_EQUATION_RUSH_LARSEN_SYCL(32); // fcas
                            SOLVE_EQUATION_RUSH_LARSEN_SYCL(33); // jca
                            SOLVE_EQUATION_RUSH_LARSEN_SYCL(34); // ffp
                            SOLVE_EQUATION_RUSH_LARSEN_SYCL(35); // fcafp
                            SOLVE_EQUATION_EULER_SYCL(36);       // nca
                            SOLVE_EQUATION_RUSH_LARSEN_SYCL(37); // b
                            SOLVE_EQUATION_RUSH_LARSEN_SYCL(38); // g
                            SOLVE_EQUATION_RUSH_LARSEN_SYCL(39); // xrf
                            SOLVE_EQUATION_RUSH_LARSEN_SYCL(40); // xrs
                            SOLVE_EQUATION_RUSH_LARSEN_SYCL(41); // xs1
                            SOLVE_EQUATION_RUSH_LARSEN_SYCL(42); // xs2
                            SOLVE_EQUATION_RUSH_LARSEN_SYCL(43); // y
                            SOLVE_EQUATION_RUSH_LARSEN_SYCL(44); // xk1
                            SOLVE_EQUATION_EULER_SYCL(45);       // u
                        }
                    }
                    // Adaptive timestep
                    else {
                        solve_forward_euler_sycl_adpt(sv, d_stim[i], d_extra_par, cur_time + max_dt, sv_id, num_cells_to_solve, abstol, reltol, dt, max_dt);
                    }
                });
        }).wait();

        // Free device memory
        sycl::free(d_stim, q_ct1);
        if (d_cells_to_solve) sycl::free(d_cells_to_solve, q_ct1);
        if (d_extra_par) sycl::free(d_extra_par, q_ct1);

    } catch (sycl::exception &e) {
        printf("SYCL exception: %s\n", e.what());
    }
}

void kernel_set_model_initial_conditions(real *_sv, int num_volumes, bool adaptive, real min_dt) {
    
    dpct::device_ext &dev_ct1 = dpct::get_current_device();
    sycl::queue &q_ct1 = dev_ct1.default_queue();
    if (_sv) {
        const int BLOCK_SIZE = 32;
        const int GRID = (num_volumes + BLOCK_SIZE - 1) / BLOCK_SIZE;
        const int TOTAL_THREADS = BLOCK_SIZE*GRID;

        try {
            q_ct1.submit([&](sycl::handler& h) {
                real *sv = _sv;
                h.parallel_for(
                    sycl::nd_range<1>(sycl::range<1>(TOTAL_THREADS), sycl::range<1>(BLOCK_SIZE)),
                    [=](sycl::nd_item<1> item) {
                        int i = item.get_global_id(0);
                        if (i < num_volumes) {
                            sv[0*num_volumes+i] = -8.668819e+01;
                            sv[1*num_volumes+i] = 4.962000e-03;
                            sv[2*num_volumes+i] = 1.020000e-04;
                            sv[3*num_volumes+i] = 8.240764e+00;
                            sv[4*num_volumes+i] = 8.240456e+00;
                            sv[5*num_volumes+i] = 8.240464e+00;
                            sv[6*num_volumes+i] = 1.437586e+02;
                            sv[7*num_volumes+i] = 1.437590e+02;
                            sv[8*num_volumes+i] = 1.437590e+02;
                            sv[9*num_volumes+i] = 4.400000e-05;
                            sv[10*num_volumes+i] = 1.020000e-04;
                            sv[11*num_volumes+i] = 1.260273e+00;
                            sv[12*num_volumes+i] = 1.245054e+00;
                            sv[13*num_volumes+i] = 1.261940e+00;
                            sv[14*num_volumes+i] = 1.530000e-04;
                            sv[15*num_volumes+i] = 0.000000e+00;
                            sv[16*num_volumes+i] = 6.321000e-03;
                            sv[17*num_volumes+i] = 7.887830e-01;
                            sv[18*num_volumes+i] = 7.887200e-01;
                            sv[19*num_volumes+i] = 7.906340e-01;
                            sv[20*num_volumes+i] = 5.799350e-01;
                            sv[21*num_volumes+i] = 7.911100e-01;
                            sv[22*num_volumes+i] = 2.420000e-04;
                            sv[23*num_volumes+i] = 4.638110e-01;
                            sv[24*num_volumes+i] = 2.404730e-01;
                            sv[25*num_volumes+i] = 2.730000e-04;
                            sv[26*num_volumes+i] = 6.497620e-01;
                            sv[27*num_volumes+i] = 9.899710e-01;
                            sv[28*num_volumes+i] = 0.000000e+00;
                            sv[29*num_volumes+i] = 1.000000e+00;
                            sv[30*num_volumes+i] = 9.268520e-01;
                            sv[31*num_volumes+i] = 1.000000e+00;
                            sv[32*num_volumes+i] = 1.000000e+00;
                            sv[33*num_volumes+i] = 9.999790e-01;
                            sv[34*num_volumes+i] = 1.000000e+00;
                            sv[35*num_volumes+i] = 1.000000e+00;
                            sv[36*num_volumes+i] = 5.431000e-03;
                            sv[37*num_volumes+i] = 3.040000e-04;
                            sv[38*num_volumes+i] = 9.942210e-01;
                            sv[39*num_volumes+i] = 3.290000e-04;
                            sv[40*num_volumes+i] = 5.707310e-01;
                            sv[41*num_volumes+i] = 1.909580e-01;
                            sv[42*num_volumes+i] = 2.230000e-04;
                            sv[43*num_volumes+i] = 2.334660e-01;
                            sv[44*num_volumes+i] = 9.970830e-01;
                            sv[45*num_volumes+i] = 4.655710e-01;
                        }
                        if (adaptive) {
                            sv[46*num_volumes+i] = min_dt;
                            sv[47*num_volumes+i] = 0.0;
                            sv[48*num_volumes+i] = 0.0;
                        }
                    }
                );
            }).wait();
        } catch (sycl::exception &e) {
            printf("SYCL exception: %s\n", e.what());
        }
    }
}

// TODO: Segmentation fault.
inline void solve_forward_euler_sycl_adpt(real *sv, real stim_current, real *extra_params, real final_time, int sv_id, int num_volumes, real abstol, real reltol, real min_dt, real max_dt) {
    #define DT sv[(NEQ)*num_volumes+sv_id]
    #define TIME_NEW sv[(NEQ+1)*num_volumes+sv_id]
    #define PREVIOUS_DT sv[(NEQ+2)*num_volumes+sv_id]

    real rDY[NEQ];

    real _tolerances_[NEQ];
    real _aux_tol = 0.0;
    real dt = DT;
    real time_new = TIME_NEW;
    real previous_dt = PREVIOUS_DT;

    real edos_old_aux_[NEQ];
    real edos_new_euler_[NEQ];
    real _k1__[NEQ];
    real _k2__[NEQ];
    real _k_aux__[NEQ];
    real sv_local[NEQ];

    const real _beta_safety_ = 0.8;

    const real __tiny_ = pow(abstol, 2.0);

    if(time_new + dt > final_time) {
        dt = final_time - time_new;
    }

    for(int i = 0; i < NEQ; i++) {
        sv_local[i] = sv[i*num_volumes+sv_id];
    }

    RHS_sycl(sv_local, rDY, stim_current, extra_params, dt, sv_id, num_volumes, true);
    time_new += dt;

    for(int i = 0; i < NEQ; i++) {
        _k1__[i] = rDY[i];
    }

	while(1) {

		for(int i = 0; i < NEQ; i++) {
			// stores the old variables in a vector
			edos_old_aux_[i] = sv_local[i];
			// computes euler method
			edos_new_euler_[i] = _k1__[i] * dt + edos_old_aux_[i];
			// steps ahead to compute the rk2 method
			sv_local[i] = edos_new_euler_[i];
		}

		time_new += dt;

		RHS_sycl(sv_local, rDY, stim_current, extra_params, dt, sv_id, num_volumes, true);
		time_new -= dt; // step back

		real greatestError = 0.0, auxError = 0.0;
		
		for(int i = 0; i < NEQ; i++) {

			// stores the new evaluation
			_k2__[i] = rDY[i];
			_aux_tol = fabs(edos_new_euler_[i]) * reltol;
			_tolerances_[i] = (abstol > _aux_tol) ? abstol : _aux_tol;

			// finds the greatest error between  the steps
			auxError = fabs(((dt / 2.0) * (_k1__[i] - _k2__[i])) / _tolerances_[i]);

			greatestError = (auxError > greatestError) ? auxError : greatestError;
		}

		/// adapt the time step
		greatestError += __tiny_;
		previous_dt = dt;

		/// adapt the time step
		dt = _beta_safety_ * dt * sqrt(1.0f / greatestError);

		if(dt < min_dt) {
			dt = min_dt;
		}
		else if(dt > max_dt) {
			dt = max_dt;
		}

		if(time_new + dt > final_time) {
			dt = final_time - time_new;
		}

		// it doesn't accept the solution or accept and risk a NaN
		if(greatestError >= 1.0f && dt > min_dt) {
			// restore the old values to do it again
			for(int i = 0; i < NEQ; i++) {
				sv_local[i] = edos_old_aux_[i];
			}
		
		} else {
			for(int i = 0; i < NEQ; i++) {
				_k_aux__[i] = _k2__[i];
				_k2__[i] = _k1__[i];
				_k1__[i] = _k_aux__[i];
			}

			for(int i = 0; i < NEQ; i++) {
				sv_local[i] = edos_new_euler_[i];
			}

			if(time_new + previous_dt >= final_time) {
				if(final_time == time_new) {
					break;
				} else if(time_new < final_time) {
					dt = previous_dt = final_time - time_new;
					time_new += previous_dt;
					break;
				} 	
			} else {
				time_new += previous_dt;
			}
		}
	}

    for(int i = 0; i < NEQ; i++) {
        sv[i*num_volumes+sv_id] = sv_local[i];
    }

    DT = dt;
    TIME_NEW = time_new;
    PREVIOUS_DT = previous_dt;
}

inline void RHS_sycl (real *sv, real *rDY_, real stim_current, real *extra_params, real dt, int sv_id, int num_volumes, bool use_adpt_dt) {

    // Current modifiers
    real GNa_Multiplier     =  extra_params[0];
    real GNaL_Multiplier    =  extra_params[1];
    real GCaT_Multiplier    =  extra_params[2];
    real Gto_Multiplier     =  extra_params[3];
    real Gsus_Multiplier    =  extra_params[4];
    real Gkr_Multiplier     =  extra_params[5];
    real Gks_Multiplier     =  extra_params[6];
    real GfNa_Multiplier    =  extra_params[7];
    real GfK_Multiplier     =  extra_params[8];
    real GK1_Multiplier     =  extra_params[9];
    real GNCX_Multiplier    =  extra_params[10];
    real GNaK_Multiplier    =  extra_params[11];
    real INa_Multiplier     =  extra_params[12]; 
    real ICaL_Multiplier    =  extra_params[13];
    real ICaNa_Multiplier   =  extra_params[14];
    real ICaK_Multiplier    =  extra_params[15];
    real Ito_Multiplier     =  extra_params[16];
    real INaL_Multiplier    =  extra_params[17];
    real IKr_Multiplier     =  extra_params[18]; 
    real IKs_Multiplier     =  extra_params[19]; 
    real IK1_Multiplier     =  extra_params[20]; 
    real INaCa_Multiplier   =  extra_params[21];
    real INaK_Multiplier    =  extra_params[22];  
    real INab_Multiplier    =  extra_params[23];  
    real ICab_Multiplier    =  extra_params[24];
    real ICaT_Multiplier    =  extra_params[25];
    real Isus_Multiplier    =  extra_params[26];
    real If_Multiplier      =  extra_params[27];
    real IpCa_Multiplier    =  extra_params[28];

    // Get the stimulus current from the current cell
    real calc_I_stim = stim_current;

    //State variables
    real v              = sv[0*num_volumes+sv_id];
    real CaMKt          = sv[1*num_volumes+sv_id];
    real cass           = sv[2*num_volumes+sv_id];
    real nai            = sv[3*num_volumes+sv_id];
    real nasl           = sv[4*num_volumes+sv_id];
    real nass           = sv[5*num_volumes+sv_id];
    real ki             = sv[6*num_volumes+sv_id];
    real kss            = sv[7*num_volumes+sv_id];
    real ksl            = sv[8*num_volumes+sv_id];
    real cai            = sv[9*num_volumes+sv_id];
    real casl           = sv[10*num_volumes+sv_id];
    real cansr          = sv[11*num_volumes+sv_id];
    real cajsr          = sv[12*num_volumes+sv_id];
    real cacsr          = sv[13*num_volumes+sv_id];
    real Jrel1          = sv[14*num_volumes+sv_id];
    real Jrel2          = sv[15*num_volumes+sv_id];
    real m              = sv[16*num_volumes+sv_id];
    real hf             = sv[17*num_volumes+sv_id];
    real hs             = sv[18*num_volumes+sv_id];
    real j              = sv[19*num_volumes+sv_id];
    real hsp            = sv[20*num_volumes+sv_id];
    real jp             = sv[21*num_volumes+sv_id];
    real mL             = sv[22*num_volumes+sv_id];
    real hL             = sv[23*num_volumes+sv_id];
    real hLp            = sv[24*num_volumes+sv_id];
    real a              = sv[25*num_volumes+sv_id];
    real i1             = sv[26*num_volumes+sv_id];
    real i2             = sv[27*num_volumes+sv_id];
    real d              = sv[28*num_volumes+sv_id];
    real ff             = sv[29*num_volumes+sv_id];
    real fs             = sv[30*num_volumes+sv_id];
    real fcaf           = sv[31*num_volumes+sv_id];
    real fcas           = sv[32*num_volumes+sv_id];
    real jca            = sv[33*num_volumes+sv_id];
    real ffp            = sv[34*num_volumes+sv_id];
    real fcafp          = sv[35*num_volumes+sv_id];
    real nca            = sv[36*num_volumes+sv_id];
    real b              = sv[37*num_volumes+sv_id];
    real g              = sv[38*num_volumes+sv_id];
    real xrf            = sv[39*num_volumes+sv_id];
    real xrs            = sv[40*num_volumes+sv_id];
    real xs1            = sv[41*num_volumes+sv_id];
    real xs2            = sv[42*num_volumes+sv_id];
    real y              = sv[43*num_volumes+sv_id];
    real xk1            = sv[44*num_volumes+sv_id];
    real u              = sv[45*num_volumes+sv_id];

    #include "trovato_2020_common.inc.c"    
}

inline void RHS_RL_sycl (real *a_, real *b_, real *sv, real *rDY_, real stim_current, real *extra_params, real dt, int sv_id, int num_volumes, bool use_adpt_dt) {
    
    // Current modifiers
    real GNa_Multiplier     =  extra_params[0];
    real GNaL_Multiplier    =  extra_params[1];
    real GCaT_Multiplier    =  extra_params[2];
    real Gto_Multiplier     =  extra_params[3];
    real Gsus_Multiplier    =  extra_params[4];
    real Gkr_Multiplier     =  extra_params[5];
    real Gks_Multiplier     =  extra_params[6];
    real GfNa_Multiplier    =  extra_params[7];
    real GfK_Multiplier     =  extra_params[8];
    real GK1_Multiplier     =  extra_params[9];
    real GNCX_Multiplier    =  extra_params[10];
    real GNaK_Multiplier    =  extra_params[11];
    real INa_Multiplier     =  extra_params[12]; 
    real ICaL_Multiplier    =  extra_params[13];
    real ICaNa_Multiplier   =  extra_params[14];
    real ICaK_Multiplier    =  extra_params[15];
    real Ito_Multiplier     =  extra_params[16];
    real INaL_Multiplier    =  extra_params[17];
    real IKr_Multiplier     =  extra_params[18]; 
    real IKs_Multiplier     =  extra_params[19]; 
    real IK1_Multiplier     =  extra_params[20]; 
    real INaCa_Multiplier   =  extra_params[21];
    real INaK_Multiplier    =  extra_params[22];  
    real INab_Multiplier    =  extra_params[23];  
    real ICab_Multiplier    =  extra_params[24];
    real ICaT_Multiplier    =  extra_params[25];
    real Isus_Multiplier    =  extra_params[26];
    real If_Multiplier      =  extra_params[27];
    real IpCa_Multiplier    =  extra_params[28];

    // Get the stimulus current from the current cell
    real calc_I_stim = stim_current;

    //State variables
    real v              = sv[0*num_volumes+sv_id];
    real CaMKt          = sv[1*num_volumes+sv_id];
    real cass           = sv[2*num_volumes+sv_id];
    real nai            = sv[3*num_volumes+sv_id];
    real nasl           = sv[4*num_volumes+sv_id];
    real nass           = sv[5*num_volumes+sv_id];
    real ki             = sv[6*num_volumes+sv_id];
    real kss            = sv[7*num_volumes+sv_id];
    real ksl            = sv[8*num_volumes+sv_id];
    real cai            = sv[9*num_volumes+sv_id];
    real casl           = sv[10*num_volumes+sv_id];
    real cansr          = sv[11*num_volumes+sv_id];
    real cajsr          = sv[12*num_volumes+sv_id];
    real cacsr          = sv[13*num_volumes+sv_id];
    real Jrel1          = sv[14*num_volumes+sv_id];
    real Jrel2          = sv[15*num_volumes+sv_id];
    real m              = sv[16*num_volumes+sv_id];
    real hf             = sv[17*num_volumes+sv_id];
    real hs             = sv[18*num_volumes+sv_id];
    real j              = sv[19*num_volumes+sv_id];
    real hsp            = sv[20*num_volumes+sv_id];
    real jp             = sv[21*num_volumes+sv_id];
    real mL             = sv[22*num_volumes+sv_id];
    real hL             = sv[23*num_volumes+sv_id];
    real hLp            = sv[24*num_volumes+sv_id];
    real a              = sv[25*num_volumes+sv_id];
    real i1             = sv[26*num_volumes+sv_id];
    real i2             = sv[27*num_volumes+sv_id];
    real d              = sv[28*num_volumes+sv_id];
    real ff             = sv[29*num_volumes+sv_id];
    real fs             = sv[30*num_volumes+sv_id];
    real fcaf           = sv[31*num_volumes+sv_id];
    real fcas           = sv[32*num_volumes+sv_id];
    real jca            = sv[33*num_volumes+sv_id];
    real ffp            = sv[34*num_volumes+sv_id];
    real fcafp          = sv[35*num_volumes+sv_id];
    real nca            = sv[36*num_volumes+sv_id];
    real b              = sv[37*num_volumes+sv_id];
    real g              = sv[38*num_volumes+sv_id];
    real xrf            = sv[39*num_volumes+sv_id];
    real xrs            = sv[40*num_volumes+sv_id];
    real xs1            = sv[41*num_volumes+sv_id];
    real xs2            = sv[42*num_volumes+sv_id];
    real y              = sv[43*num_volumes+sv_id];
    real xk1            = sv[44*num_volumes+sv_id];
    real u              = sv[45*num_volumes+sv_id];

    #include "trovato_2020_RL_common.inc.c"
}
