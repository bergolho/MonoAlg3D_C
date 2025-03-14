#include "../../gpu_utils/sycl_utils.h"
#include <cstdio>
#include <sycl/sycl.hpp>
#include <dpct/dpct.hpp>

#include "../../common_types/common_types.h"

#include "ToRORd_fkatp_mixed_endo_mid_epi.h"

extern "C" SET_ODE_INITIAL_CONDITIONS_SYCL(set_model_initial_conditions_sycl)
{
    log_info("Using ToRORd_fkatp_2019 SYCL model\n");
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

    // Get initial condition from extra_data
    real *initial_conditions_endo = NULL;
    real *initial_conditions_epi = NULL;
    real *initial_conditions_mid = NULL;
    real *transmurality = NULL;
	real *sf_Iks = NULL;
    real *d_initial_conditions_endo = NULL;
    real *d_initial_conditions_epi = NULL;
    real *d_initial_conditions_mid = NULL;
    real *d_transmurality = NULL;
	real *d_sf_Iks = NULL;
    if(solver->ode_extra_data) {
        log_info("[INFO] Applying a mask function using the defined extra data function!\n");
        struct extra_data_for_torord_gksgkrtjca_twave *extra_data = (struct extra_data_for_torord_gksgkrtjca_twave*)solver->ode_extra_data;
        initial_conditions_endo = extra_data->initial_ss_endo;
        initial_conditions_epi = extra_data->initial_ss_epi;
        initial_conditions_mid = extra_data->initial_ss_mid;
        transmurality = extra_data->transmurality;
		sf_Iks = extra_data->sf_IKs;

        d_initial_conditions_endo = sycl::malloc_device<real>(NEQ, q_ct1);
        q_ct1.memcpy(d_initial_conditions_endo, initial_conditions_endo, NEQ*sizeof(real));
        d_initial_conditions_epi = sycl::malloc_device<real>(NEQ, q_ct1);
        q_ct1.memcpy(d_initial_conditions_epi, initial_conditions_epi, NEQ*sizeof(real));
        d_initial_conditions_mid = sycl::malloc_device<real>(NEQ, q_ct1);
        q_ct1.memcpy(d_initial_conditions_mid, initial_conditions_mid, NEQ*sizeof(real));
        d_transmurality = sycl::malloc_device<real>(num_volumes, q_ct1);
        q_ct1.memcpy(d_transmurality, transmurality, num_volumes*sizeof(real));
        d_sf_Iks = sycl::malloc_device<real>(num_volumes, q_ct1);
        q_ct1.memcpy(d_sf_Iks, sf_Iks, num_volumes*sizeof(real));
    }
    else {
        log_info("[INFO] You should supply a mask function to tag the cells when using this mixed model!\n");
        log_info("[INFO] Considering all cells ENDO!\n");
    }

    if (solver->ode_extra_data) {
        kernel_set_model_initial_conditions_endo_mid_epi(solver->sv, num_volumes, adaptive, solver->min_dt,\
                                                            d_initial_conditions_endo, d_initial_conditions_epi, d_initial_conditions_mid,\
                                                            d_transmurality, d_sf_Iks);
    }
    else {
        kernel_set_model_initial_conditions(solver->sv, num_volumes, adaptive, solver->min_dt);
    }

    if (d_initial_conditions_endo) sycl::free(d_initial_conditions_endo, q_ct1);
    if (d_initial_conditions_epi) sycl::free(d_initial_conditions_epi, q_ct1);
    if (d_initial_conditions_mid) sycl::free(d_initial_conditions_mid, q_ct1);
    if (d_transmurality) sycl::free(d_transmurality, q_ct1);
    if (d_sf_Iks) sycl::free(d_sf_Iks, q_ct1);
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
        q_ct1.memcpy(d_cells_to_solve, cells_to_solve, num_cells_to_solve*sizeof(uint32_t));
    }
    q_ct1.memcpy(d_stim, stim_currents, num_cells_to_solve*sizeof(real)).wait();

    // Explicitly Copy host memory to device
    int num_extra_parameters = 20;
    real *transmurality = NULL;
	real *sf_Iks = NULL;
	real *d_transmurality = NULL;
    real *d_sf_Iks = NULL;
    real extra_par[num_extra_parameters];
    real *d_extra_par = NULL;

    // Get the extra data array if exists
    // 'Transmurality' and 'sf_IKs' mapping defined on 'extra_data' function
    if(has_extra_params) {
        struct extra_data_for_torord_gksgkrtjca_twave *extra_data = (struct extra_data_for_torord_gksgkrtjca_twave*)ode_solver->ode_extra_data;
        extra_par[0]  = extra_data->INa_Multiplier; 
        extra_par[1]  = extra_data->INaL_Multiplier;
        extra_par[2]  = extra_data->INaCa_Multiplier;
        extra_par[3]  = extra_data->INaK_Multiplier;
        extra_par[4]  = extra_data->INab_Multiplier; 
        extra_par[5]  = extra_data->Ito_Multiplier;
        extra_par[6]  = extra_data->IKr_Multiplier; 
        extra_par[7]  = extra_data->IKs_Multiplier; 
        extra_par[8]  = extra_data->IK1_Multiplier;
        extra_par[9]  = extra_data->IKb_Multiplier;
        extra_par[10]  = extra_data->IKCa_Multiplier;
        extra_par[11] = extra_data->ICaL_Multiplier;  
        extra_par[12] = extra_data->ICab_Multiplier;  
        extra_par[13] = extra_data->IpCa_Multiplier;
        extra_par[14] = extra_data->ICaCl_Multiplier; 
        extra_par[15] = extra_data->IClb_Multiplier;
        extra_par[16] = extra_data->Jrel_Multiplier;
        extra_par[17] = extra_data->Jup_Multiplier;
        extra_par[18] = extra_data->aCaMK_Multiplier;
        extra_par[19] = extra_data->taurelp_Multiplier;
        sf_Iks = extra_data->sf_IKs;
        transmurality = extra_data->transmurality;

        d_transmurality = sycl::malloc_device<real>(num_cells_to_solve, q_ct1);
        q_ct1.memcpy(d_transmurality, transmurality, num_cells_to_solve*sizeof(real));
        d_sf_Iks = sycl::malloc_device<real>(num_cells_to_solve, q_ct1);
        q_ct1.memcpy(d_sf_Iks, sf_Iks, num_cells_to_solve*sizeof(real));
        d_extra_par = sycl::malloc_device<real>(num_extra_parameters, q_ct1);
        q_ct1.memcpy(d_extra_par, extra_par, num_extra_parameters*sizeof(real));
    }
    // No [extra_data] section, we consider all cells ENDO!
    else {

        // Default: initialize all current modifiers
        for (uint32_t i = 0; i < num_extra_parameters; i++) {
            if (i == 10)
                extra_par[i] = 0.0;
            else 
                extra_par[i] = 1.0;
        }
        d_extra_par = sycl::malloc_device<real>(num_extra_parameters, q_ct1);
        q_ct1.memcpy(d_extra_par, extra_par, num_extra_parameters*sizeof(real));
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

                            // Transmurality and sf_IKs are given by the columns of the ALG mesh file 
                            if (has_extra_params) {
                                RHS_RL_sycl(a, b, sv, rDY, d_stim[i], d_transmurality[i], d_sf_Iks[i], d_extra_par, dt, sv_id, num_cells_to_solve, false);
                            }
                            // Default: all cells are ENDO and sf_IKs equals to 1
                            else {
                                RHS_RL_sycl(a, b, sv, rDY, d_stim[i], 0.0, 1.0, d_extra_par, dt, sv_id, num_cells_to_solve, false);
                            }

                            // Solve variables based on its type:
                            //  Non-linear = Euler
                            //  Hodkin-Huxley = Rush-Larsen || Euler (if 'a' coefficient is too small)
                            SOLVE_EQUATION_EULER_SYCL(0);        // v        
                            SOLVE_EQUATION_EULER_SYCL(1);        // nai
                            SOLVE_EQUATION_EULER_SYCL(2);        // nass
                            SOLVE_EQUATION_EULER_SYCL(3);        // ki   
                            SOLVE_EQUATION_EULER_SYCL(4);        // kss
                            SOLVE_EQUATION_EULER_SYCL(5);        // cai
                            SOLVE_EQUATION_EULER_SYCL(6);        // cass 
                            SOLVE_EQUATION_EULER_SYCL(7);        // cansr
                            SOLVE_EQUATION_EULER_SYCL(8);        // cajsr
                            SOLVE_EQUATION_RUSH_LARSEN_SYCL(9);  // m
                            SOLVE_EQUATION_RUSH_LARSEN_SYCL(10); // hp
                            SOLVE_EQUATION_RUSH_LARSEN_SYCL(11); // h
                            SOLVE_EQUATION_RUSH_LARSEN_SYCL(12); // j
                            SOLVE_EQUATION_RUSH_LARSEN_SYCL(13); // jp
                            SOLVE_EQUATION_RUSH_LARSEN_SYCL(14); // mL
                            SOLVE_EQUATION_RUSH_LARSEN_SYCL(15); // hL
                            SOLVE_EQUATION_RUSH_LARSEN_SYCL(16); // hLp
                            SOLVE_EQUATION_RUSH_LARSEN_SYCL(17); // a
                            SOLVE_EQUATION_RUSH_LARSEN_SYCL(18); // iF
                            SOLVE_EQUATION_RUSH_LARSEN_SYCL(19); // iS
                            SOLVE_EQUATION_RUSH_LARSEN_SYCL(20); // ap
                            SOLVE_EQUATION_RUSH_LARSEN_SYCL(21); // iFp
                            SOLVE_EQUATION_RUSH_LARSEN_SYCL(22); // iSp
                            SOLVE_EQUATION_RUSH_LARSEN_SYCL(23); // d
                            SOLVE_EQUATION_RUSH_LARSEN_SYCL(24); // ff
                            SOLVE_EQUATION_RUSH_LARSEN_SYCL(25); // fs
                            SOLVE_EQUATION_RUSH_LARSEN_SYCL(26); // fcaf
                            SOLVE_EQUATION_RUSH_LARSEN_SYCL(27); // fcas
                            SOLVE_EQUATION_RUSH_LARSEN_SYCL(28); // jca
                            SOLVE_EQUATION_EULER_SYCL(29);       // nca
                            SOLVE_EQUATION_EULER_SYCL(30);       // nca_i
                            SOLVE_EQUATION_RUSH_LARSEN_SYCL(31); // ffp
                            SOLVE_EQUATION_RUSH_LARSEN_SYCL(32); // fcafp
                            SOLVE_EQUATION_RUSH_LARSEN_SYCL(33); // xs1
                            SOLVE_EQUATION_RUSH_LARSEN_SYCL(34); // xs2
                            SOLVE_EQUATION_RUSH_LARSEN_SYCL(35); // Jrel_np
                            SOLVE_EQUATION_EULER_SYCL(36);        // CaMKt 
                            SOLVE_EQUATION_EULER_SYCL(37);       // ikr_c0
                            SOLVE_EQUATION_EULER_SYCL(38);       // ikr_c1
                            SOLVE_EQUATION_EULER_SYCL(39);       // ikr_c2
                            SOLVE_EQUATION_EULER_SYCL(40);       // ikr_o
                            SOLVE_EQUATION_EULER_SYCL(41);       // ikr_i
                            SOLVE_EQUATION_RUSH_LARSEN_SYCL(42); // Jrel_p
                        }
                    }
                    // Adaptive timestep
                    else {
                        solve_forward_euler_sycl_adpt(sv, d_stim[i], d_transmurality[i], d_sf_Iks[i], d_extra_par, cur_time + max_dt, sv_id, num_cells_to_solve, abstol, reltol, dt, max_dt);
                    }
                });
        }).wait();

        // Free device memory
        sycl::free(d_stim, q_ct1);
        if (d_cells_to_solve) sycl::free(d_cells_to_solve, q_ct1);
        if (d_transmurality) sycl::free(d_transmurality, q_ct1);
        if (d_sf_Iks) sycl::free(d_sf_Iks, q_ct1);
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
                            sv[0*num_volumes+i]  = -88.6369922306458;
                            sv[1*num_volumes+i]  = 11.8973412949238;
                            sv[2*num_volumes+i]  = 11.897661047085;
                            sv[3*num_volumes+i]  = 141.234464714982;
                            sv[4*num_volumes+i]  = 141.234423402713;
                            sv[5*num_volumes+i]  = 7.26747296460659e-05;
                            sv[6*num_volumes+i]  = 6.33786975780735e-05;
                            sv[7*num_volumes+i]  = 1.5326530637197;
                            sv[8*num_volumes+i]  = 1.53394579180493;
                            sv[9*num_volumes+i]  = 0.000828007761976018;
                            sv[10*num_volumes+i] = 0.666527193684116;
                            sv[11*num_volumes+i] = 0.826020806005678;
                            sv[12*num_volumes+i] = 0.826055985895856;
                            sv[13*num_volumes+i] = 0.825850881115628;
                            sv[14*num_volumes+i] = 0.000166868626513013;
                            sv[15*num_volumes+i] = 0.522830604669169;
                            sv[16*num_volumes+i] = 0.285969584294187;
                            sv[17*num_volumes+i] = 0.000959137028030184;
                            sv[18*num_volumes+i] = 0.999601150012565;
                            sv[19*num_volumes+i] = 0.5934016398361;
                            sv[20*num_volumes+i] = 0.000488696137242056;
                            sv[21*num_volumes+i] = 0.999601147267179;
                            sv[22*num_volumes+i] = 0.654668660159696;
                            sv[23*num_volumes+i] = 9.50007519781516e-32;
                            sv[24*num_volumes+i] = 0.999999992317577;
                            sv[25*num_volumes+i] = 0.939258048397962;
                            sv[26*num_volumes+i] = 0.999999992317557;
                            sv[27*num_volumes+i] = 0.999898379647465;
                            sv[28*num_volumes+i] = 0.99997825156004;
                            sv[29*num_volumes+i] = 0.000444816183420527;
                            sv[30*num_volumes+i] = 0.000755072490632667;
                            sv[31*num_volumes+i] = 0.999999992318446;
                            sv[32*num_volumes+i] = 0.999999992318445;
                            sv[33*num_volumes+i] = 0.24240468344952;
                            sv[34*num_volumes+i] = 0.000179537726989804;
                            sv[35*num_volumes+i] = -6.88308558109975e-25;
                            sv[36*num_volumes+i] = 0.0111749845355653;
                            sv[37*num_volumes+i] = 0.998036620213316;
                            sv[38*num_volumes+i] = 0.000858801779013532;
                            sv[39*num_volumes+i] = 0.000709744678350176;
                            sv[40*num_volumes+i] = 0.000381261722195702;
                            sv[41*num_volumes+i] = 1.35711566929992e-05;
                            sv[42*num_volumes+i] = 2.30252452954649e-23;
                        }
                        if (adaptive) {
                            sv[43*num_volumes+i] = min_dt;
                            sv[44*num_volumes+i] = 0.0;
                            sv[45*num_volumes+i] = 0.0;
                        }
                    }
                );
            }).wait();
        } catch (sycl::exception &e) {
            printf("SYCL exception: %s\n", e.what());
        }
    }
}

void kernel_set_model_initial_conditions_endo_mid_epi(real *_sv, int num_volumes, bool adaptive, real min_dt,\
                                                            real *initial_endo, real *initial_epi, real *initial_mid,\
                                                            real *transmurality, real *sf_Iks) {
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
                            
                            for (int j = 0; j < NEQ; j++) {
                                if (transmurality[i] == ENDO) {
                                    sv[j*num_volumes+i] = initial_endo[j];
                                }
                                else if (transmurality[i] == EPI) {
                                    sv[j*num_volumes+i] = initial_epi[j];
                                }
                                else if (transmurality[i] == MID) {
                                    sv[j*num_volumes+i] = initial_mid[j];
                                }
                            }
                            if(adaptive) {
                                sv[43*num_volumes+i] = min_dt; // dt
                                sv[44*num_volumes+i] = 0.0;    // time_new
                                sv[45*num_volumes+i] = 0.0;    // previous dt
                            }
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
inline void solve_forward_euler_sycl_adpt(real *sv, real stim_current, real transmurality, real sf_Iks, real *extra_params, real final_time, int sv_id, int num_volumes, real abstol, real reltol, real min_dt, real max_dt) {
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

    RHS_sycl(sv_local, rDY, stim_current, transmurality, sf_Iks, extra_params, dt, sv_id, num_volumes, true);
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

		RHS_sycl(sv_local, rDY, stim_current, transmurality, sf_Iks, extra_params, dt, sv_id, num_volumes, true);
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

inline void RHS_sycl (real *sv, real *rDY_, real stim_current, real transmurality, real sf_Iks, real *extra_params, real dt, int sv_id, int num_volumes, bool use_adpt_dt) {

    // Current modifiers
    real INa_Multiplier = extra_params[0];   
    real INaL_Multiplier = extra_params[1];  
    real INaCa_Multiplier = extra_params[2];  
    real INaK_Multiplier = extra_params[3];  
    real INab_Multiplier = extra_params[4];   
    real Ito_Multiplier = extra_params[5];  
    real IKr_Multiplier = extra_params[6];   
    real IKs_Multiplier = extra_params[7];   
    real IK1_Multiplier = extra_params[8];  
    real IKb_Multiplier = extra_params[9];  
    real IKCa_Multiplier = extra_params[10];  
    real ICaL_Multiplier = extra_params[11];   
    real ICab_Multiplier = extra_params[12];   
    real IpCa_Multiplier = extra_params[13]; 
    real ICaCl_Multiplier = extra_params[14];  
    real IClb_Multiplier = extra_params[15]; 
    real Jrel_Multiplier = extra_params[16]; 
    real Jup_Multiplier = extra_params[17]; 
    real aCaMK_Multiplier = extra_params[18]; 
    real taurelp_Multiplier = extra_params[19];

    // Get the celltype for the current cell
    real celltype = transmurality;

    // Get the stimulus current from the current cell
    real calc_I_stim = stim_current;

    // State variables
    real v       = sv[0 * num_volumes + sv_id];
    real nai     = sv[1 * num_volumes + sv_id];
    real nass    = sv[2 * num_volumes + sv_id];
    real ki      = sv[3 * num_volumes + sv_id];
    real kss     = sv[4 * num_volumes + sv_id];
    real cai     = sv[5 * num_volumes + sv_id];
    real cass    = sv[6 * num_volumes + sv_id];
    real cansr   = sv[7 * num_volumes + sv_id];
    real cajsr   = sv[8 * num_volumes + sv_id];
    real m       = sv[9 * num_volumes + sv_id];
    real hp      = sv[10 * num_volumes + sv_id];
    real h       = sv[11 * num_volumes + sv_id];
    real j       = sv[12 * num_volumes + sv_id];
    real jp      = sv[13 * num_volumes + sv_id];
    real mL      = sv[14 * num_volumes + sv_id];
    real hL      = sv[15 * num_volumes + sv_id];
    real hLp     = sv[16 * num_volumes + sv_id];
    real a       = sv[17 * num_volumes + sv_id];
    real iF      = sv[18 * num_volumes + sv_id];
    real iS      = sv[19 * num_volumes + sv_id];
    real ap      = sv[20 * num_volumes + sv_id];
    real iFp     = sv[21 * num_volumes + sv_id];
    real iSp     = sv[22 * num_volumes + sv_id];
    real d       = sv[23 * num_volumes + sv_id];
    real ff      = sv[24 * num_volumes + sv_id];
    real fs      = sv[25 * num_volumes + sv_id];
    real fcaf    = sv[26 * num_volumes + sv_id];
    real fcas    = sv[27 * num_volumes + sv_id];
    real jca     = sv[28 * num_volumes + sv_id];
    real nca     = sv[29 * num_volumes + sv_id];
    real nca_i   = sv[30 * num_volumes + sv_id];
    real ffp     = sv[31 * num_volumes + sv_id];
    real fcafp   = sv[32 * num_volumes + sv_id];
    real xs1     = sv[33 * num_volumes + sv_id];
    real xs2     = sv[34 * num_volumes + sv_id];
    real Jrel_np = sv[35 * num_volumes + sv_id];
    real CaMKt   = sv[36 * num_volumes + sv_id];
    real ikr_c0  = sv[37 * num_volumes + sv_id];
    real ikr_c1  = sv[38 * num_volumes + sv_id];
    real ikr_c2  = sv[39 * num_volumes + sv_id];
    real ikr_o   = sv[40 * num_volumes + sv_id];
    real ikr_i   = sv[41 * num_volumes + sv_id];
    real Jrel_p  = sv[42 * num_volumes + sv_id];

    #include "ToRORd_fkatp_mixed_endo_mid_epi.common.c"    
}

inline void RHS_RL_sycl (real *a_, real *b_, real *sv, real *rDY_, real stim_current, real transmurality, real sf_Iks, real *extra_params, real dt, int sv_id, int num_volumes, bool use_adpt_dt) {
    
        // Current modifiers
    real INa_Multiplier = extra_params[0];   
    real INaL_Multiplier = extra_params[1];  
    real INaCa_Multiplier = extra_params[2];  
    real INaK_Multiplier = extra_params[3];  
    real INab_Multiplier = extra_params[4];   
    real Ito_Multiplier = extra_params[5];  
    real IKr_Multiplier = extra_params[6];   
    real IKs_Multiplier = extra_params[7];   
    real IK1_Multiplier = extra_params[8];  
    real IKb_Multiplier = extra_params[9];  
    real IKCa_Multiplier = extra_params[10];  
    real ICaL_Multiplier = extra_params[11];   
    real ICab_Multiplier = extra_params[12];   
    real IpCa_Multiplier = extra_params[13]; 
    real ICaCl_Multiplier = extra_params[14];  
    real IClb_Multiplier = extra_params[15]; 
    real Jrel_Multiplier = extra_params[16]; 
    real Jup_Multiplier = extra_params[17]; 
    real aCaMK_Multiplier = extra_params[18]; 
    real taurelp_Multiplier = extra_params[19];

    // Get the celltype for the current cell
    real celltype = transmurality;

    // Get the stimulus current from the current cell
    real calc_I_stim = stim_current;

    // State variables
    real v       = sv[0 * num_volumes + sv_id];
    real nai     = sv[1 * num_volumes + sv_id];
    real nass    = sv[2 * num_volumes + sv_id];
    real ki      = sv[3 * num_volumes + sv_id];
    real kss     = sv[4 * num_volumes + sv_id];
    real cai     = sv[5 * num_volumes + sv_id];
    real cass    = sv[6 * num_volumes + sv_id];
    real cansr   = sv[7 * num_volumes + sv_id];
    real cajsr   = sv[8 * num_volumes + sv_id];
    real m       = sv[9 * num_volumes + sv_id];
    real hp      = sv[10 * num_volumes + sv_id];
    real h       = sv[11 * num_volumes + sv_id];
    real j       = sv[12 * num_volumes + sv_id];
    real jp      = sv[13 * num_volumes + sv_id];
    real mL      = sv[14 * num_volumes + sv_id];
    real hL      = sv[15 * num_volumes + sv_id];
    real hLp     = sv[16 * num_volumes + sv_id];
    real a       = sv[17 * num_volumes + sv_id];
    real iF      = sv[18 * num_volumes + sv_id];
    real iS      = sv[19 * num_volumes + sv_id];
    real ap      = sv[20 * num_volumes + sv_id];
    real iFp     = sv[21 * num_volumes + sv_id];
    real iSp     = sv[22 * num_volumes + sv_id];
    real d       = sv[23 * num_volumes + sv_id];
    real ff      = sv[24 * num_volumes + sv_id];
    real fs      = sv[25 * num_volumes + sv_id];
    real fcaf    = sv[26 * num_volumes + sv_id];
    real fcas    = sv[27 * num_volumes + sv_id];
    real jca     = sv[28 * num_volumes + sv_id];
    real nca     = sv[29 * num_volumes + sv_id];
    real nca_i   = sv[30 * num_volumes + sv_id];
    real ffp     = sv[31 * num_volumes + sv_id];
    real fcafp   = sv[32 * num_volumes + sv_id];
    real xs1     = sv[33 * num_volumes + sv_id];
    real xs2     = sv[34 * num_volumes + sv_id];
    real Jrel_np = sv[35 * num_volumes + sv_id];
    real CaMKt   = sv[36 * num_volumes + sv_id];
    real ikr_c0  = sv[37 * num_volumes + sv_id];
    real ikr_c1  = sv[38 * num_volumes + sv_id];
    real ikr_c2  = sv[39 * num_volumes + sv_id];
    real ikr_o   = sv[40 * num_volumes + sv_id];
    real ikr_i   = sv[41 * num_volumes + sv_id];
    real Jrel_p  = sv[42 * num_volumes + sv_id];

    #include "ToRORd_fkatp_mixed_endo_mid_epi_RL.common.c"
}
