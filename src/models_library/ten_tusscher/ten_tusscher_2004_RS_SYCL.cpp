#include "../../gpu_utils/sycl_utils.h"
#include <cstdio>
#include <sycl/sycl.hpp>
#include <dpct/dpct.hpp>

#include "../../common_types/common_types.h"

#include "ten_tusscher_2004.h"

extern "C" SET_ODE_INITIAL_CONDITIONS_SYCL(set_model_initial_conditions_sycl)
{
    printf("Using TenTusscher3 2004 SYCL model\n");
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
                            sv[0*num_cells+i] = INITIAL_V; // V millivolt
                            sv[1*num_cells+i] = 0.f;       // M
                            sv[2*num_cells+i] = 0.75;      // H
                            sv[3*num_cells+i] = 0.75f;     // J
                            sv[4*num_cells+i] = 0.f;       // Xr1
                            sv[5*num_cells+i] = 1.f;       // Xr2
                            sv[6*num_cells+i] = 0.f;       // Xs
                            sv[7*num_cells+i] = 1.f;       // S
                            sv[8*num_cells+i] = 0.f;       // R
                            sv[9*num_cells+i] = 0.f;       // D
                            sv[10*num_cells+i] = 1.f;      // F
                            sv[11*num_cells+i] = 1.f;      // FCa
                            sv[12*num_cells+i] = 1.f;      // G
                            sv[13*num_cells+i] = 0.0002;   // Cai
                            sv[14*num_cells+i] = 0.2f;     // CaSR
                            sv[15*num_cells+i] = 11.6f;    // Nai
                            sv[16*num_cells+i] = 138.3f;   // Ki
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
        q_ct1.memcpy(d_cells_to_solve, cells_to_solve, num_cells_to_solve*sizeof(uint32_t));
    }
    
    // Copy initial data to device
    q_ct1.memcpy(d_stim, stim_currents, num_cells_to_solve*sizeof(real));

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

                    #pragma unroll
                    for (int n = 0; n < num_steps; ++n) {
                        RHS_sycl(sv, d_stim[i], rDY, sv_id, dt, num_cells_to_solve);

                        // Update state variables
                        #pragma unroll
                        for (int j = 0; j < NEQ; j++) {
                            sv[j * num_cells_to_solve + i] = rDY[j];
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

inline void RHS_sycl(real *Y, real stim_current, real *dY, int sv_id, real dt, int num_cells) {
    // State variables
    real svolt = Y[0 * num_cells + sv_id];
    real sm    = Y[1 * num_cells + sv_id];
    real sh    = Y[2 * num_cells + sv_id];
    real sj    = Y[3 * num_cells + sv_id];
    real sxr1  = Y[4 * num_cells + sv_id];
    real sxr2  = Y[5 * num_cells + sv_id];
    real sxs   = Y[6 * num_cells + sv_id];
    real ss    = Y[7 * num_cells + sv_id];
    real sr    = Y[8 * num_cells + sv_id];
    real sd    = Y[9 * num_cells + sv_id];
    real sf    = Y[10 * num_cells + sv_id];
    real sfca  = Y[11 * num_cells + sv_id];
    real sg    = Y[12 * num_cells + sv_id];
    real Cai   = Y[13 * num_cells + sv_id];
    real CaSR  = Y[14 * num_cells + sv_id];
    real Nai   = Y[15 * num_cells + sv_id];
    real Ki    = Y[16 * num_cells + sv_id];

    //External concentrations
    real Ko=5.4;
    real Cao=2.0;
    real Nao=140.0;

    //Intracellular volumes
    real Vc=0.016404;
    real Vsr=0.001094;

    //Calcium dynamics
    real Bufc=0.15f;
    real Kbufc=0.001f;
    real Bufsr=10.f;
    real Kbufsr=0.3f;
    real taufca=2.f;
    real taug=2.f;
    real Vmaxup=0.000425f;
    real Kup=0.00025f;

//Constants
    const real R = 8314.472f;
    const real F = 96485.3415f;
    const real T =310.0f;
    real RTONF   =(R*T)/F;

//Cellular capacitance
    real CAPACITANCE=0.185;

//Parameters for currents
//Parameters for IKr
    real Gkr=0.096;
//Parameters for Iks
    real pKNa=0.03;
#ifdef EPI
    real Gks=0.245;
#endif
#ifdef ENDO
    real Gks=0.245;
#endif
#ifdef MCELL
    real Gks=0.062;
#endif
//Parameters for Ik1
    real GK1=5.405;
//Parameters for Ito
#ifdef EPI
    real Gto=0.294;
#endif
#ifdef ENDO
    real Gto=0.073;
#endif
#ifdef MCELL
    real Gto=0.294;
#endif
//Parameters for INa
    real GNa=14.838;
//Parameters for IbNa
    real GbNa=0.00029;
//Parameters for INaK
    real KmK=1.0;
    real KmNa=40.0;
    real knak=1.362;
//Parameters for ICaL
    real GCaL=0.000175;
//Parameters for IbCa
    real GbCa=0.000592;
//Parameters for INaCa
    real knaca=1000;
    real KmNai=87.5;
    real KmCa=1.38;
    real ksat=0.1;
    real n=0.35;
//Parameters for IpCa
    real GpCa=0.825;
    real KpCa=0.0005;
//Parameters for IpK;
    real GpK=0.0146;


    real IKr;
    real IKs;
    real IK1;
    real Ito;
    real INa;
    real IbNa;
    real ICaL;
    real IbCa;
    real INaCa;
    real IpCa;
    real IpK;
    real INaK;
    real Irel;
    real Ileak;


    real dNai;
    real dKi;
    real dCai;
    real dCaSR;

    real A;
//    real BufferFactorc;
//    real BufferFactorsr;
    real SERCA;
    real Caisquare;
    real CaSRsquare;
    real CaCurrent;
    real CaSRCurrent;


    real fcaold;
    real gold;
    real Ek;
    real Ena;
    real Eks;
    real Eca;
    real CaCSQN;
    real bjsr;
    real cjsr;
    real CaBuf;
    real bc;
    real cc;
    real Ak1;
    real Bk1;
    real rec_iK1;
    real rec_ipK;
    real rec_iNaK;
    real AM;
    real BM;
    real AH_1;
    real BH_1;
    real AH_2;
    real BH_2;
    real AJ_1;
    real BJ_1;
    real AJ_2;
    real BJ_2;
    real M_INF;
    real H_INF;
    real J_INF;
    real TAU_M;
    real TAU_H;
    real TAU_J;
    real axr1;
    real bxr1;
    real axr2;
    real bxr2;
    real Xr1_INF;
    real Xr2_INF;
    real TAU_Xr1;
    real TAU_Xr2;
    real Axs;
    real Bxs;
    real Xs_INF;
    real TAU_Xs;
    real R_INF;
    real TAU_R;
    real S_INF;
    real TAU_S;
    real Ad;
    real Bd;
    real Cd;
    real TAU_D;
    real D_INF;
    real TAU_F;
    real F_INF;
    real FCa_INF;
    real G_INF;

    real inverseVcF2=1/(2*Vc*F);
    real inverseVcF=1./(Vc*F);
    real Kupsquare=Kup*Kup;
//    real BufcKbufc=Bufc*Kbufc;
//    real Kbufcsquare=Kbufc*Kbufc;
//    real Kbufc2=2*Kbufc;
//    real BufsrKbufsr=Bufsr*Kbufsr;
//    const real Kbufsrsquare=Kbufsr*Kbufsr;
//    const real Kbufsr2=2*Kbufsr;
    const real exptaufca=sycl::exp(-dt/taufca);
    const real exptaug=sycl::exp(-dt/taug);

    real sItot;

    //Needed to compute currents
    Ek=RTONF*(log((Ko/Ki)));
    Ena=RTONF*(log((Nao/Nai)));
    Eks=RTONF*(log((Ko+pKNa*Nao)/(Ki+pKNa*Nai)));
    Eca=0.5*RTONF*(log((Cao/Cai)));
    Ak1=0.1/(1.+sycl::exp(0.06*(svolt-Ek-200)));
    Bk1=(3.*sycl::exp(0.0002*(svolt-Ek+100))+
         sycl::exp(0.1*(svolt-Ek-10)))/(1.+sycl::exp(-0.5*(svolt-Ek)));
    rec_iK1=Ak1/(Ak1+Bk1);
    rec_iNaK=(1./(1.+0.1245*sycl::exp(-0.1*svolt*F/(R*T))+0.0353*sycl::exp(-svolt*F/(R*T))));
    rec_ipK=1./(1.+sycl::exp((25-svolt)/5.98));


    //Compute currents
    INa=GNa*sm*sm*sm*sh*sj*(svolt-Ena);
    ICaL=GCaL*sd*sf*sfca*4*svolt*(F*F/(R*T))*
         (sycl::exp(2*svolt*F/(R*T))*Cai-0.341*Cao)/(sycl::exp(2*svolt*F/(R*T))-1.);
    Ito=Gto*sr*ss*(svolt-Ek);
    IKr=Gkr*sycl::sqrt(Ko/5.4)*sxr1*sxr2*(svolt-Ek);
    IKs=Gks*sxs*sxs*(svolt-Eks);
    IK1=GK1*rec_iK1*(svolt-Ek);
    INaCa=knaca*(1./(KmNai*KmNai*KmNai+Nao*Nao*Nao))*(1./(KmCa+Cao))*
          (1./(1+ksat*sycl::exp((n-1)*svolt*F/(R*T))))*
          (sycl::exp(n*svolt*F/(R*T))*Nai*Nai*Nai*Cao-
           sycl::exp((n-1)*svolt*F/(R*T))*Nao*Nao*Nao*Cai*2.5);
    INaK=knak*(Ko/(Ko+KmK))*(Nai/(Nai+KmNa))*rec_iNaK;
    IpCa=GpCa*Cai/(KpCa+Cai);
    IpK=GpK*rec_ipK*(svolt-Ek);
    IbNa=GbNa*(svolt-Ena);
    IbCa=GbCa*(svolt-Eca);


    //Determine total current
    (sItot) = IKr    +
              IKs   +
              IK1   +
              Ito   +
              INa   +
              IbNa  +
              ICaL  +
              IbCa  +
              INaK  +
              INaCa +
              IpCa  +
              IpK   +
              stim_current;


    //update concentrations
    Caisquare=Cai*Cai;
    CaSRsquare=CaSR*CaSR;
    CaCurrent=-(ICaL+IbCa+IpCa-2.0f*INaCa)*inverseVcF2*CAPACITANCE;
    A=0.016464f*CaSRsquare/(0.0625f+CaSRsquare)+0.008232f;
    Irel=A*sd*sg;
    Ileak=0.00008f*(CaSR-Cai);
    SERCA=Vmaxup/(1.f+(Kupsquare/Caisquare));
    CaSRCurrent=SERCA-Irel-Ileak;
    CaCSQN=Bufsr*CaSR/(CaSR+Kbufsr);
    dCaSR=dt*(Vc/Vsr)*CaSRCurrent;
    bjsr=Bufsr-CaCSQN-dCaSR-CaSR+Kbufsr;
    cjsr=Kbufsr*(CaCSQN+dCaSR+CaSR);
    CaSR=(sycl::sqrt(bjsr*bjsr+4.*cjsr)-bjsr)/2.;
    CaBuf=Bufc*Cai/(Cai+Kbufc);
    dCai=dt*(CaCurrent-CaSRCurrent);
    bc=Bufc-CaBuf-dCai-Cai+Kbufc;
    cc=Kbufc*(CaBuf+dCai+Cai);
    Cai=(sycl::sqrt(bc*bc+4*cc)-bc)/2;



    dNai=-(INa+IbNa+3*INaK+3*INaCa)*inverseVcF*CAPACITANCE;
    Nai+=dt*dNai;

    dKi=-(stim_current+IK1+Ito+IKr+IKs-2*INaK+IpK)*inverseVcF*CAPACITANCE;
    Ki+=dt*dKi;

    //compute steady state values and time constants
    AM=1./(1.+sycl::exp((-60.-svolt)/5.));
    BM=0.1/(1.+sycl::exp((svolt+35.)/5.))+0.10/(1.+sycl::exp((svolt-50.)/200.));
    TAU_M=AM*BM;
    M_INF=1./((1.+sycl::exp((-56.86-svolt)/9.03))*(1.+sycl::exp((-56.86-svolt)/9.03)));
    if (svolt>=-40.)
    {
        AH_1=0.;
        BH_1=(0.77/(0.13*(1.+sycl::exp(-(svolt+10.66)/11.1))));
        TAU_H= 1.0/(AH_1+BH_1);
    }
    else
    {
        AH_2=(0.057*sycl::exp(-(svolt+80.)/6.8));
        BH_2=(2.7*sycl::exp(0.079*svolt)+(3.1e5)*sycl::exp(0.3485*svolt));
        TAU_H=1.0/(AH_2+BH_2);
    }
    H_INF=1./((1.+sycl::exp((svolt+71.55)/7.43))*(1.+sycl::exp((svolt+71.55)/7.43)));
    if(svolt>=-40.)
    {
        AJ_1=0.;
        BJ_1=(0.6*sycl::exp((0.057)*svolt)/(1.+sycl::exp(-0.1*(svolt+32.))));
        TAU_J= 1.0/(AJ_1+BJ_1);
    }
    else
    {
        AJ_2=(((-2.5428e4)*sycl::exp(0.2444*svolt)-(6.948e-6)*
                                             sycl::exp(-0.04391*svolt))*(svolt+37.78)/
              (1.+sycl::exp(0.311*(svolt+79.23))));
        BJ_2=(0.02424*sycl::exp(-0.01052*svolt)/(1.+sycl::exp(-0.1378*(svolt+40.14))));
        TAU_J= 1.0/(AJ_2+BJ_2);
    }
    J_INF=H_INF;

    Xr1_INF=1./(1.+sycl::exp((-26.-svolt)/7.));
    axr1=450./(1.+sycl::exp((-45.-svolt)/10.));
    bxr1=6./(1.+sycl::exp((svolt-(-30.))/11.5));
    TAU_Xr1=axr1*bxr1;
    Xr2_INF=1./(1.+sycl::exp((svolt-(-88.))/24.));
    axr2=3./(1.+sycl::exp((-60.-svolt)/20.));
    bxr2=1.12/(1.+sycl::exp((svolt-60.)/20.));
    TAU_Xr2=axr2*bxr2;

    Xs_INF=1./(1.+sycl::exp((-5.-svolt)/14.));
    Axs=1100./(sycl::sqrt(1.+sycl::exp((-10.-svolt)/6)));
    Bxs=1./(1.+sycl::exp((svolt-60.)/20.));
    TAU_Xs=Axs*Bxs;

#ifdef EPI
    R_INF=1./(1.+sycl::exp((20-svolt)/6.));
    S_INF=1./(1.+sycl::exp((svolt+20)/5.));
    TAU_R=9.5*sycl::exp(-(svolt+40.)*(svolt+40.)/1800.)+0.8;
    TAU_S=85.*sycl::exp(-(svolt+45.)*(svolt+45.)/320.)+5./(1.+sycl::exp((svolt-20.)/5.))+3.;
#endif
#ifdef ENDO
    R_INF=1./(1.+sycl::exp((20-svolt)/6.));
    S_INF=1./(1.+sycl::exp((svolt+28)/5.));
    TAU_R=9.5*sycl::exp(-(svolt+40.)*(svolt+40.)/1800.)+0.8;
    TAU_S=1000.*sycl::exp(-(svolt+67)*(svolt+67)/1000.)+8.;
#endif
#ifdef MCELL
    R_INF=1./(1.+sycl::exp((20-svolt)/6.));
    S_INF=1./(1.+sycl::exp((svolt+20)/5.));
    TAU_R=9.5*sycl::exp(-(svolt+40.)*(svolt+40.)/1800.)+0.8;
    TAU_S=85.*sycl::exp(-(svolt+45.)*(svolt+45.)/320.)+5./(1.+sycl::exp((svolt-20.)/5.))+3.;
#endif


    D_INF=1./(1.+sycl::exp((-5-svolt)/7.5));
    Ad=1.4/(1.+sycl::exp((-35-svolt)/13))+0.25;
    Bd=1.4/(1.+sycl::exp((svolt+5)/5));
    Cd=1./(1.+sycl::exp((50-svolt)/20));
    TAU_D=Ad*Bd+Cd;
    F_INF=1./(1.+sycl::exp((svolt+20)/7));
    TAU_F=1125*sycl::exp(-(svolt+27)*(svolt+27)/300)+80+165/(1.+sycl::exp((25-svolt)/10));


    FCa_INF=(1./(1.+sycl::pow((Cai/0.000325),8))+
             0.1/(1.+sycl::exp((Cai-0.0005)/0.0001))+
             0.20/(1.+sycl::exp((Cai-0.00075)/0.0008))+
             0.23 )/1.46;
    if(Cai<0.00035)
        G_INF=1./(1.+sycl::pow((Cai/0.00035),6));
    else
        G_INF=1./(1.+sycl::pow((Cai/0.00035),16));

    //Update gates
    dY[1]  = M_INF-(M_INF-sm)*sycl::exp(-dt/TAU_M);
    dY[2]  = H_INF-(H_INF-sh)*sycl::exp(-dt/TAU_H);
    dY[3]  = J_INF-(J_INF-sj)*sycl::exp(-dt/TAU_J);
    dY[4]  = Xr1_INF-(Xr1_INF-sxr1)*sycl::exp(-dt/TAU_Xr1);
    dY[5]  = Xr2_INF-(Xr2_INF-sxr2)*sycl::exp(-dt/TAU_Xr2);
    dY[6]  = Xs_INF-(Xs_INF-sxs)*sycl::exp(-dt/TAU_Xs);
    dY[7]  = S_INF-(S_INF-ss)*sycl::exp(-dt/TAU_S);
    dY[8]  = R_INF-(R_INF-sr)*sycl::exp(-dt/TAU_R);
    dY[9]  = D_INF-(D_INF-sd)*sycl::exp(-dt/TAU_D);
    dY[10] = F_INF-(F_INF-sf)*sycl::exp(-dt/TAU_F);
    fcaold= sfca;
    sfca = FCa_INF-(FCa_INF-sfca)*exptaufca;
    if(sfca>fcaold && (svolt)>-37.0)
        sfca = fcaold;
    gold = sg;
    sg = G_INF-(G_INF-sg)*exptaug;

    if(sg>gold && (svolt)>-37.0)
        sg=gold;

    //update voltage
    dY[0] = svolt + dt*(-sItot);
    dY[11] = sfca;
    dY[12] = sg;
    dY[13] = Cai;
    dY[14] = CaSR;
    dY[15] = Nai;
    dY[16] = Ki;
}
