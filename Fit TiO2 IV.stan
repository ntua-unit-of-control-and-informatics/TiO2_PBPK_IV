functions{
    real [] pbpk(real t,
                real[] M,
                real[] theta,
                real[] rdata,
                int[] idata) {


    real dMdt[33];
    real Q_total; real V_blood; real Vven; real Vart; real Wm_ven; real Wm_art;

    real W_st; real W_lu; real W_spl; real W_li; real W_ki;
    real W_ht; real W_br; real W_ut; real W_skel;

    real Vtis_st; real Vcap_st; real Wm_st; real Q_st;
    real Vtis_ht; real Vcap_ht; real Wm_ht; real Q_ht;
    real Vtis_ki; real Vcap_ki; real Wm_ki; real Q_ki;
    real Vtis_br; real Vcap_br; real Wm_br; real Q_br;
    real Vtis_spl; real Vcap_spl; real Wm_spl; real Q_spl;
    real Vtis_lu; real Vcap_lu; real Wm_lu; real Q_lu;
    real Vtis_li; real Vcap_li; real Wm_li; real Q_li;
    real Vtis_ut; real Vcap_ut; real Wm_ut; real Q_ut;
    real Vtis_skel; real Vcap_skel; real Wm_skel; real Q_skel;

    real Ccap_lu; real Ctis_lu; real Cm_lu;
    real Ccap_spl; real Ctis_spl; real Cm_spl;
    real Ccap_li; real Ctis_li; real Cm_li;
    real Ccap_ki; real Ctis_ki; real Cm_ki;
    real Ccap_ht; real Ctis_ht; real Cm_ht;
    real Ccap_br; real Ctis_br; real Cm_br;
    real Ccap_ut; real Ctis_ut; real Cm_ut;
    real Ccap_skel; real Ctis_skel; real Cm_skel;
    real Ccap_st; real Ctis_st; real Cm_st;

    real C_art; real C_ven;

    real xfast; real xslow; real xbr; real Pup_max; real uptake_st;
    real k_de; real P_li; real P_st; real P_rest; real CLE_f; real CLE_u;
    real uptake; real uptake_spl; real uptake_li;

    real PA_lu; real PA_spl; real PA_li; real PA_ki; real PA_ht;
    real PA_br; real PA_ut; real PA_skel; real PA_st;

    real P_lu; real P_spl; real P_ki; real P_ht;
    real P_br; real P_ut; real P_skel;

    real Pup_lu; real Pup_spl; real Pup_li; real Pup_ki; real Pup_ht;
    real Pup_br; real Pup_ut; real Pup_skel; real Pup_st;
    real Pup_ven; real Pup_art;

    real x_lu; real x_spl; real x_li; real x_ki; real x_ht; real x_br;
    real x_ut; real x_skel; real x_st;

    real uptake_ki; real uptake_ht; real uptake_br; real uptake_ut;
    real uptake_skel; real uptake_blood; real uptake_lu;

    Q_total = rdata[1];
    V_blood = rdata[2];
    Vven = rdata[3];
    Vart = rdata[4];
    Wm_ven = rdata[5];
    Wm_art = rdata[6];

    W_st = rdata[7];
    W_ht = rdata[8];
    W_ki = rdata[9];
    W_br = rdata[10];
    W_spl = rdata[11];
    W_lu = rdata[12];
    W_li = rdata[13];
    W_ut = rdata[14];
    W_skel = rdata[15];

    Vtis_st = rdata[16];
    Vtis_ht = rdata[17];
    Vtis_ki = rdata[18];
    Vtis_br = rdata[19];
    Vtis_spl = rdata[20];
    Vtis_lu = rdata[21];
    Vtis_li = rdata[22];
    Vtis_ut = rdata[23];
    Vtis_skel = rdata[24];

    Vcap_st = rdata[25];
    Vcap_ht = rdata[26];
    Vcap_ki = rdata[27];
    Vcap_br = rdata[28];
    Vcap_spl = rdata[29];
    Vcap_lu = rdata[30];
    Vcap_li = rdata[31];
    Vcap_ut = rdata[32];
    Vcap_skel = rdata[33];

    Wm_st = rdata[34];
    Wm_ht = rdata[35];
    Wm_ki = rdata[36];
    Wm_br = rdata[37];
    Wm_spl = rdata[38];
    Wm_lu = rdata[39];
    Wm_li = rdata[40];
    Wm_ut = rdata[41];
    Wm_skel = rdata[42];

    Q_st = rdata[43];
    Q_ht = rdata[44];
    Q_ki = rdata[45];
    Q_br = rdata[46];
    Q_spl = rdata[47];
    Q_lu = rdata[48];
    Q_li = rdata[49];
    Q_ut = rdata[50];
    Q_skel = rdata[51];

    xfast = rdata[52];
    xslow = rdata[53];
    xbr = rdata[54];
    Pup_max = rdata[55];
    uptake_st = rdata[56];
    k_de = rdata[57];

    P_st = rdata[58];
    CLE_u = rdata[59];
    P_li = rdata[60];

    ///////////params///////////
    // set up the stohastic parameters to be estimated
    uptake = theta[1];
    uptake_spl = theta[2];
    uptake_li = theta[3];
    P_rest = theta[4];
    CLE_f = theta[5];



    //Capillary concentrations
   Ccap_lu = M[1]/Vcap_lu;
   Ccap_spl = M[4]/Vcap_spl;
   Ccap_li = M[7]/Vcap_li;
   Ccap_ki = M[10]/Vcap_ki;
   Ccap_ht = M[13]/Vcap_ht;
   Ccap_br = M[16]/Vcap_br;
   Ccap_ut = M[19]/Vcap_ut;
   Ccap_skel = M[22]/Vcap_skel;
   Ccap_st = M[25]/Vcap_st;

   //Tissue concentration
   Ctis_lu = M[2]/Vtis_lu;
   Ctis_spl = M[5]/Vtis_spl;
   Ctis_li = M[8]/Vtis_li;
   Ctis_ki = M[11]/Vtis_ki;
   Ctis_ht = M[14]/Vtis_ht;
   Ctis_br = M[17]/Vtis_br;
   Ctis_ut = M[20]/Vtis_ut;
   Ctis_skel = M[23]/Vtis_skel;
   Ctis_st = M[26]/Vtis_st;

   C_ven = M[28]/Vven;
   C_art = M[30]/Vart;

   PA_lu = xslow*Q_total;
   PA_spl = xfast*Q_spl;
   PA_li = xfast*Q_li;
   PA_ki = xslow*Q_ki;
   PA_ht = xslow*Q_ht;
   PA_br = xbr*Q_br;
   PA_ut = xslow*Q_ut;
   PA_skel = xfast*Q_skel;
   PA_st = xslow*Q_st;

   P_lu = P_rest;
   P_spl = P_rest;
   //P_li = P_li;
   P_ki = P_rest;
   P_ht = P_rest;
   P_br = P_rest;
   P_ut = P_rest;
   P_skel = P_rest;
   //P_st = P_st;

   Pup_lu = Pup_max*(1-(M[3]/(Wm_lu*uptake)));
   Pup_spl = Pup_max*(1-(M[6]/(Wm_spl*uptake_spl)));
   Pup_li = Pup_max*(1-(M[9]/(Wm_li*uptake_li)));
   Pup_ki = Pup_max*(1-(M[12]/(Wm_ki*uptake)));
   Pup_ht = Pup_max*(1-(M[15]/(Wm_ht*uptake)));
   Pup_br = Pup_max*(1-(M[18]/(Wm_br*uptake)));
   Pup_ut = Pup_max*(1-(M[21]/(Wm_ut*uptake)));
   Pup_skel = Pup_max*(1-(M[24]/(Wm_skel*uptake)));
   Pup_st = Pup_max*(1-(M[27]/(Wm_st*uptake_st)));
   Pup_ven = Pup_max*(1-(M[29]/(Wm_ven*uptake)));
   Pup_art = Pup_max*(1-(M[31]/(Wm_art*uptake)));

   //Lungs
   dMdt[1] = Q_total*C_ven - Q_total*Ccap_lu - PA_lu*Ccap_lu + PA_lu*Ctis_lu/P_lu;
   dMdt[2] = PA_lu*Ccap_lu - PA_lu*Ctis_lu/P_lu - Pup_lu*Vtis_lu*Ctis_lu + k_de*M[3];
   dMdt[3]   = Pup_lu*Vtis_lu*Ctis_lu - k_de*M[3];

   //Spleen
   dMdt[4] = Q_spl*C_art - Q_spl*Ccap_spl - PA_spl*Ccap_spl + PA_spl*Ctis_spl/P_spl;
   dMdt[5] = PA_spl*Ccap_spl - PA_spl*Ctis_spl/P_spl - Pup_spl*Vtis_spl*Ctis_spl + k_de*M[6];
   dMdt[6]   = Pup_spl*Vtis_spl*Ctis_spl - k_de*M[6];

   //Liver
   dMdt[7] = Q_li*C_art + Q_spl*Ccap_spl - (Q_li+Q_spl)*Ccap_li - PA_li*Ccap_li + PA_li*Ctis_li/P_li;
   dMdt[8] = PA_li*Ccap_li - PA_li*Ctis_li/P_li - Pup_li*Vtis_li*Ctis_li + k_de*M[9] - CLE_f*M[8];
   dMdt[9] = Pup_li*Vtis_li*Ctis_li - k_de*M[9];
   dMdt[32] = CLE_f*M[8];

   //Kidneys
   dMdt[10] = Q_ki*C_art - Q_ki*Ccap_ki - PA_ki*Ccap_ki + PA_ki*Ctis_ki/P_ki - CLE_u*M[10];
   dMdt[11] = PA_ki*Ccap_ki - PA_ki*Ctis_ki/P_ki - Pup_ki*Vtis_ki*Ctis_ki + k_de*M[12];
   dMdt[12] = Pup_ki*Vtis_ki*Ctis_ki - k_de*M[12];
   dMdt[33] = CLE_u*M[10];

   //Heart
   dMdt[13] = Q_ht*C_art - Q_ht*Ccap_ht - PA_ht*Ccap_ht + PA_ht*Ctis_ht/P_ht;
   dMdt[14] = PA_ht*Ccap_ht - PA_ht*Ctis_ht/P_ht - Pup_ht*Vtis_ht*Ctis_ht + k_de*M[15];
   dMdt[15] = Pup_ht*Vtis_ht*Ctis_ht - k_de*M[15];

   //Brain
   dMdt[16] = Q_br*C_art - Q_br*Ccap_br - PA_br*Ccap_br + PA_br*Ctis_br/P_br;
   dMdt[17] = PA_br*Ccap_br - PA_br*Ctis_br/P_br - Pup_br*Vtis_br*Ctis_br + k_de*M[18];
   dMdt[18] = Pup_br*Vtis_br*Ctis_br - k_de*M[18];

   //Uterus
   dMdt[19] = Q_ut*C_art - Q_ut*Ccap_ut - PA_ut*Ccap_ut + PA_ut*Ctis_ut/P_ut;
   dMdt[20] = PA_ut*Ccap_ut - PA_ut*Ctis_ut/P_ut - Pup_ut*Vtis_ut*Ctis_ut + k_de*M[21];
   dMdt[21] = Pup_ut*Vtis_ut*Ctis_ut - k_de*M[21];

   //Skeleton
   dMdt[22] = Q_skel*C_art - Q_skel*Ccap_skel - PA_skel*Ccap_skel + PA_skel*Ctis_skel/P_skel;
   dMdt[23] = PA_skel*Ccap_skel - PA_skel*Ctis_skel/P_skel - Pup_skel*Vtis_skel*Ctis_skel + k_de*M[24];
   dMdt[24] = Pup_skel*Vtis_skel*Ctis_skel - k_de*M[24];

   //Soft tissue
   dMdt[25] = Q_st*C_art - Q_st*Ccap_st - PA_st*Ccap_st + PA_st*Ctis_st/P_st;
   dMdt[26] = PA_st*Ccap_st - PA_st*Ctis_st/P_st - Pup_st*Vtis_st*Ctis_st + k_de*M[27];
   dMdt[27] = Pup_st*Vtis_st*Ctis_st - k_de*M[27];

   //Veins
   dMdt[28] = - Q_total*C_ven + (Q_li+Q_spl)*Ccap_li + Q_ki*Ccap_ki + Q_ht*Ccap_ht +
              Q_br*Ccap_br + Q_ut*Ccap_ut + Q_skel*Ccap_skel +
              Q_st*Ccap_st - Pup_ven*C_ven + k_de*M[29];
   dMdt[29] = Pup_ven*C_ven - k_de*M[29];

   //Arteries
   dMdt[30] = Q_total*Ccap_lu - Q_spl*C_art - Q_li*C_art - Q_ki*C_art -
              Q_ht*C_art - Q_br*C_art - Q_ut*C_art - Q_skel*C_art - Q_st*C_art -
              Pup_art*C_art + k_de*M[31];
   dMdt[31] = Pup_art*C_art - k_de*M[31];

   return dMdt;
   }
}
//////////////////////////////////////////////////////////////////////////
// Input data from R code
data{

        int<lower=0> N_param;                // Number of parameters to be estimated
        int<lower=0> N_compart;              //number of observed compartments
        int<lower=0> N_diff;                 // number of differential equations
        real  time[5];
        real  urine_time[11];
        real  feces_time[4];
        real  mass[5,N_compart];
        real  feces[4];
        real  urine[11];
        real  m0[N_diff];           // Initial concentration in compartments
        real  t_init;                  // Initial time
        real<lower=0>  eta_tr[N_param];
        real<lower=0>  eta_tr_std[N_param];
        real  params[60];      // Matrix containing the individual parameters
        real  rel_tol;
        real  abs_tol;
        real  max_num_steps;

}
//////////////////////////////////////////////////////////////////
transformed data {
      real rdata[0];
      int idata[0];
      //vector[N_param]  eta_tr_std ;
      vector[N_param]  eta_std ;
      vector[N_param]  eta ;
      vector [N_param] H;                //covariance matrix


// Tranform from normal to lognormal distribution
      for (i in 1:N_param){
              //eta_tr_std[i] = eta_tr[i]*0.5;
              eta_std[i] = sqrt(log(((eta_tr_std[i]^2)/(eta_tr[i])^2)+1));
              eta[i]=log(((eta_tr[i])^2)/sqrt((eta_tr_std[i]^2)+(eta_tr[i])^2));
              H[i] = eta_std[i];
      }
// Remove the comments below in order to test if the ODEs give right results
// considering the mean values of the parameters

    //  print(integrate_ode_bdf(pbpk, m0, t_init, time,
    //             to_array_1d(eta_tr[:]), params[:], idata,
    //             rel_tol, abs_tol, max_num_steps))
}
//////////////////////////////////////////////////////////////////
parameters{
// Declare the number of statistical errors and the stohastic parametric vector
// and their lower (and upper, if needed) bounds.

        real<lower=0>  sigma[3];
        vector<lower=0> [N_param] theta;
}

////////////////////////////////////////////////////////////////////
model{
real m_hat[5,N_diff];
real total_m_hat[5,N_compart];
real feces_excr[4,N_diff];
real urine_excr[11,N_diff];
real feces_hat[4];
real urine_hat[11];

//priors distributions
sigma[1] ~ normal(0,0.1);
sigma[2] ~ normal(0,0.01);
sigma[3] ~ normal(0,0.001);
theta[:] ~ lognormal(eta[:], H[:]);

//likelihood

// Solve the ODEs system for each sampling
// We have three different solutions because the time points of the experimental
// data for the feces and urine are different.

m_hat[:,:] = integrate_ode_bdf(pbpk, m0, t_init, time,
            to_array_1d(theta[:]), params[:], idata,
            rel_tol, abs_tol, max_num_steps);


feces_excr[:,:] = integrate_ode_bdf(pbpk, m0, t_init, feces_time,
            to_array_1d(theta[:]), params[:], idata,
            rel_tol, abs_tol, max_num_steps);

urine_excr[:,:] = integrate_ode_bdf(pbpk, m0, t_init, urine_time,
            to_array_1d(theta[:]), params[:], idata,
            rel_tol, abs_tol, max_num_steps);

for (i in 1:5){
            //Total amount of NPs in each organ

            //Amount in Liver
            total_m_hat[i,1] = m_hat[i,8] + m_hat[i,9];

            //Amount in Spleen
            total_m_hat[i,2] = m_hat[i,5] + m_hat[i,6];

            //Amount in Kidneys
            total_m_hat[i,3] =  m_hat[i,11] + m_hat[i,12];

            //Ammount in Lungs
            total_m_hat[i,4] =  m_hat[i,2] + m_hat[i,3];

            //Amount in Heart
            total_m_hat[i,5] =  m_hat[i,14] + m_hat[i,15];

            //Amount in Brain
            total_m_hat[i,6] =  m_hat[i,17] + m_hat[i,18];

            //Amount in Uterus
            total_m_hat[i,7] =  m_hat[i,20] + m_hat[i,21];

            //Amount in Blood
            total_m_hat[i,8] = m_hat[i,28] + m_hat[i,29] +
                               m_hat[i,30] + m_hat[i,31];

            //Amount in Skeleton
            total_m_hat[i,9] =  m_hat[i,23] + m_hat[i,24];

            //Amount in Soft Tissue
            total_m_hat[i,10] = m_hat[i,26] + m_hat[i,27];
}

            feces_hat[:] = feces_excr[:,32];
            urine_hat[:] = urine_excr[:,33];

//Liver
to_vector(mass[:,1]) ~ normal(to_vector(total_m_hat[:,1]),sigma[1]);
//Spleen
to_vector(mass[:,2]) ~ normal(to_vector(total_m_hat[:,2]),sigma[1]);
//Kidneys
to_vector(mass[:,3]) ~ normal(to_vector(total_m_hat[:,3]),sigma[2]);
//Lungs
to_vector(mass[:,4]) ~ normal(to_vector(total_m_hat[:,4]),sigma[2]);
//Heart
to_vector(mass[:,5]) ~ normal(to_vector(total_m_hat[:,5]),sigma[3]);
//Brain
to_vector(mass[:,6]) ~ normal(to_vector(total_m_hat[:,6]),sigma[3]);
//Uterus
to_vector(mass[:,7]) ~ normal(to_vector(total_m_hat[:,7]),sigma[3]);
//Blood
to_vector(mass[:,8]) ~ normal(to_vector(total_m_hat[:,8]),sigma[1]);
//Skeleton
to_vector(mass[:,9]) ~ normal(to_vector(total_m_hat[:,9]),sigma[1]);
//Soft Tissue
to_vector(mass[:,10]) ~ normal(to_vector(total_m_hat[:,10]),sigma[2]);

to_vector(feces[:]) ~ normal(to_vector(feces_hat[:]),sigma[2]);
to_vector(urine[:]) ~ normal(to_vector(urine_hat[:]),sigma[2]);
}
////////////////////////////////////////////////////////////////////

generated quantities{

}
