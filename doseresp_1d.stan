// AP. Woodward, University of Georgia, 2025.
// Implementation of a dose-response model for a single agent (4-parameter log-logistic).
// This uses a multilevel structure to incorporate between-subject and between-plate variation.
//      For the PD parameters, there is partial pooling at the plate level, but not the subject level.
// A hierarchical observation model handles information from plate-specific controls.
//      Conditional variability may also vary at the subject level (subject-varying observation model).
// The model specification is adapted from (https://doi.org/10.1371/journal.pone.0146021).

functions {
  // 'DR_model()' contains the structural model, and returns predicted response for a single set of parameters.
  vector DR_model(vector dose, real Emin, real Emax, real EC50, real H){
    // Initialize internal variables.
    int nd;
    vector [size(dose)] result;
    nd = size(dose);
    for (i in 1:nd){
      if (dose[i] == 0){
        result[i] = Emin;
      }else{
        result[i] = (Emax + ((Emin-Emax)/((1+exp(H*(log(dose[i])-log(EC50)))))));
      }
    }
    return result;
  }
  // 'DR_caller()' is a high-level function that is called in the transformed parameters block.
  // This calls the total structural model on each plate (uniform set of parameters). 
  // It accepts from the data a count of observations P, and vectors of plate and subject labels.
  // The return is an (P) vector containing the model predictions to match the data.
  vector DR_caller(vector dose, real [] Emin, real [] Emax, real [] EC50, real [] H, int P, int R, int [] plate_ind, int [] subject_ind, int [] first_ind, int [] last_ind){
      vector [P] raw_predictions;
       // for each plate-subject;
      for (L in 1:R){
          raw_predictions[first_ind[L]:last_ind[L]] = DR_model(dose[first_ind[L]:last_ind[L]], Emin[L], Emax[L], EC50[L], H[L]);
      }
    return raw_predictions;
  }
}

data{
  int<lower=1> Nsub;               // the number of subjects.
  int<lower=1> Nplate;             // the number of plates.
  int<lower=1> P;                  // the number of experimental observations (non-controls).
  real response_exp[P];            // response variable for the experimental observations.
  int<lower=1> plate_exp[P];       // plate index of the experimental observations.
  int<lower=1> subject_exp[P];     // subject index of the experimental observations.
  vector [P] dose;                 // dose/concentration for drug 1.
  int<lower=1> plate_ind[Nplate];  // The plate labels in the plate index.
  int<lower=1> subj_ind[Nplate];   // The subject labels in the plate index.
  int<lower=1> first_ind[Nplate];  // First row in the experimental dataset for each plate.
  int<lower=1> last_ind[Nplate];   // Last row in the experimental dataset for each plate.
  int<lower=1> D;                  // The number of dead control observations.
  real response_dead[D];           // response variable for the dead-control observations.
  int<lower=1> subject_dead[D];    // subject index of the dead-control observations.
  int<lower=1> plate_dead[D];      // plate index of the dead-control observations.
}

parameters{
  // Observation model parameters.
  real dead_resp_mu[Nsub];
  real <lower = 0> dead_resp_sd[Nsub];
  // Locations of the PD parameters (for each of 2 subjects).
  real Emin_mu[Nsub];     
  real Emax_mu[Nsub];
  real EC50_mu[Nsub];  
  real H_mu[Nsub];  
  //
  real <lower = 0> dead_resp_mu_sd;
  real <lower = 0> Emin_sd;
  real <lower = 0> Emax_sd;
  real <lower = 0> EC50_sd;
  real <lower = 0> H_sd;
  //
  matrix [5,Nplate] plate_ZS_mat;
  cholesky_factor_corr[5] L_plate;
}

transformed parameters{
  // Initialize the internal variables.
  real <lower = 0> dead_VEC [Nplate];
  real <lower = 0> Emin_VEC [Nplate];
  real <lower = 0> Emax_VEC [Nplate];
  real <lower = 0> EC50_VEC [Nplate];
  real <lower = 0> H_VEC [Nplate];
  vector <lower = 0> [P] resp_pred;
  // Initialize the matrix of 'random' effects at the plate level.
  matrix[5, Nplate] plate_ZS_ranef;
  plate_ZS_ranef   = diag_pre_multiply([1,1,1,1,1], L_plate)*plate_ZS_mat;
  // Assign the plate-level pharmacodynamic model parameters on the centered scale.
  //    This parameterization constrains the dose-response relationship to be inhibitory.
  //    The response is constrained to be positive and is defined as a multiple of the response for the killed control.
  for (E in 1:Nplate){
    dead_VEC[E]   = exp(dead_resp_mu[subj_ind[E]] + (plate_ZS_ranef[1,E] * dead_resp_mu_sd));
    Emin_VEC[E]   = dead_VEC[E] * (1+exp(Emin_mu[subj_ind[E]]    + (plate_ZS_ranef[2,E] * Emin_sd)));
    Emax_VEC[E]   = dead_VEC[E] + ((Emin_VEC[E]-dead_VEC[E]) * inv_logit(Emax_mu[subj_ind[E]] + (plate_ZS_ranef[3,E] * Emax_sd)));
    EC50_VEC[E]   = exp(EC50_mu[subj_ind[E]] + (plate_ZS_ranef[4,E] * EC50_sd));
    H_VEC[E]      = exp(H_mu[subj_ind[E]]    + (plate_ZS_ranef[5,E] * H_sd));
  }
  // Generate the model-predicted concentrations at the plate level, by calling to the structural model.
  resp_pred = DR_caller(dose, Emin_VEC, Emax_VEC, EC50_VEC, H_VEC, P, Nplate, plate_ind, subj_ind, first_ind, last_ind);
}

model{
  // Priors for the observation model.
  dead_resp_mu    ~ normal(-1,1);
  dead_resp_sd    ~ normal(0,0.5);
  // Priors for the pharmacodynamic model.
  Emin_mu         ~ normal(-1,1);     
  Emax_mu         ~ normal(0,2);
  EC50_mu         ~ normal(-5,3);    
  H_mu            ~ normal(0,1);       
  //
  dead_resp_mu_sd ~ normal(0,0.5);
  Emin_sd         ~ normal(0,0.2);
  Emax_sd         ~ normal(0,0.2);
  EC50_sd         ~ normal(0,0.2);
  H_sd            ~ normal(0,0.2);
  //
  to_vector(plate_ZS_mat) ~ normal(0,1);
  L_plate ~ lkj_corr_cholesky(2);
  // The observation model for the dead controls.
  response_dead ~ lognormal(dead_resp_mu[subject_dead], dead_resp_sd[subject_dead]);
  // The observation model for the responses.
  response_exp ~ lognormal(log(resp_pred), dead_resp_sd[subj_ind[plate_exp]]);
}

generated quantities{
  real<lower = 0> Emin_pop[Nsub];
  real<lower = 0, upper = 1> Emax_logit[Nsub];
  real<lower = 0> Emax_pop[Nsub];
  real<lower = 0> EC50_pop[Nsub];
  real<lower = 0> H_pop[Nsub];
  for (G in 1:Nsub){
    Emin_pop[G]   = exp(dead_resp_mu[G]) * (1+exp(Emin_mu[G]));
    Emax_logit[G] = inv_logit(Emax_mu[G]);
    Emax_pop[G]   = exp(dead_resp_mu[G]) + ((Emin_pop[G] - exp(dead_resp_mu[G])) * Emax_logit[G]);
    EC50_pop[G]   = exp(EC50_mu[G]);
    H_pop[G]      = exp(H_mu[G]);
  }
}
