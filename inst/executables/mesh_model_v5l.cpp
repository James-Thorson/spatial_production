#include <TMB.hpp>

// Posfun
template<class Type>
Type posfun(Type x, Type lowerlimit, Type &pen){
  pen += CppAD::CondExpLt(x,lowerlimit,Type(0.01)*pow(x-lowerlimit,2),Type(0));
  return CppAD::CondExpGe(x,lowerlimit,x,lowerlimit/(Type(2)-x/lowerlimit));
}

// Function for detecting NAs
template<class Type>
bool isNA(Type x){
  return R_IsNA(asDouble(x));
}

// function for logistic transform
template<class Type>
Type plogis(Type x){
  return 1.0 / (1.0 + exp(-x));
}

// dlognorm
template<class Type>
Type dlognorm(Type x, Type meanlog, Type sdlog, int give_log=false){
  Type Return;
  if(give_log==false) Return = dnorm( log(x), meanlog, sdlog, false) / x;
  if(give_log==true) Return = dnorm( log(x), meanlog, sdlog, true) - log(x);
  return Return;
}

// dzinflognorm
template<class Type>
Type dzinflognorm(Type x, Type meanlog, Type encounter_prob, Type sdlog, int give_log=false){
  Type Return;
  if(x==0){
    if(give_log==false) Return = 1.0 - encounter_prob;
    if(give_log==true) Return = log(1.0 - encounter_prob);
  }else{
    if(give_log==false) Return = encounter_prob * dlognorm( x, meanlog, sdlog, false );
    if(give_log==true) Return = log(encounter_prob) + dlognorm( x, meanlog, sdlog, true );
  } 
  return Return;
}

// dzinfgamma, shape = 1/CV^2, scale = mean*CV^2
template<class Type>
Type dzinfgamma(Type x, Type posmean, Type encounter_prob, Type cv, int give_log=false){
  Type Return;
  if(x==0){
    if(give_log==false) Return = 1.0 - encounter_prob;
    if(give_log==true) Return = log(1.0 - encounter_prob);
  }else{
    if(give_log==false) Return = encounter_prob * dgamma( x, pow(cv,-2), posmean*pow(cv,2), false );
    if(give_log==true) Return = log(encounter_prob) + dgamma( x, pow(cv,-2), posmean*pow(cv,2), true );
  } 
  return Return;
}

// dzinfnorm
template<class Type>
Type dzinfnorm(Type x, Type posmean, Type encounter_prob, Type cv, int give_log=false){
  Type Return;
  if(x==0){
    if(give_log==false) Return = 1.0 - encounter_prob;
    if(give_log==true) Return = log(1.0 - encounter_prob);
  }else{
    if(give_log==false) Return = encounter_prob * dnorm( x, posmean, posmean*cv, false );
    if(give_log==true) Return = log(encounter_prob) + dnorm( x, posmean, posmean*cv, true );
  } 
  return Return;
}

// Main function
template<class Type>
Type objective_function<Type>::operator() ()
{
  // Options
  DATA_FACTOR( Options_vec );
  // Slot 0: Include movement
  // Slot 1: Which dynamical model (0=Gompertz, 1=Ricker)
  // Slot 2: what data type (0=Poisson, 1=delta-lognormal)
  // Slot 3: what transform for occupancy prob (0=logistic, 1=asympotic-negative-exponential)
  // Slot 4: Parameterization for east-north-west-south movement log_mvec: 0=Cardinal directions(0=east,1=north,2=west,3=south); 1=Relative(0=east,1=north-relative-to-east,2=west-relative-to-est,3=south-relative-to-east)
  
  // Dimensions
  DATA_INTEGER(n_i);         // Number of observations (stacked across all years)
  DATA_INTEGER(n_r);        // Number of triangles
  DATA_INTEGER(n_g);        // Number of vertices in GMRF approximation to triangle values
  DATA_INTEGER(n_t);         // Number of real "strata" (i.e., verticies containing data)
  DATA_INTEGER(n_tdiv);

  // Data vectors
  DATA_VECTOR( c_i );
  DATA_FACTOR( r_i );
  DATA_FACTOR( t_i );
  DATA_VECTOR( km2_i );
  DATA_VECTOR( km2_r );

  // effort and catch data
  DATA_VECTOR( catch_t );
  DATA_MATRIX( effortdens_rt );
  
  // Movement matrices
  DATA_SPARSE_MATRIX( M1 );
  DATA_SPARSE_MATRIX( M2 );
  DATA_SPARSE_MATRIX( M3 );
  DATA_SPARSE_MATRIX( M4 );

  // SPDE objects
  DATA_SPARSE_MATRIX(G0);
  DATA_SPARSE_MATRIX(G1);
  DATA_SPARSE_MATRIX(G2);

  // life history parameters
  DATA_MATRIX( bio_priors );
  
  // Fixed effects
  PARAMETER( log_beta );
  PARAMETER( alpha );
  PARAMETER_VECTOR( log_mvec );   // Movement parameters (east, north, west, south)
  PARAMETER( logkappa );
  PARAMETER( logSigmaU );
  PARAMETER( logSigmaO );
  PARAMETER( logmeanu0 );
  PARAMETER_VECTOR( densvar_z );

  // Fishing mortality
  PARAMETER( ln_sigmacatch );
  PARAMETER_VECTOR( ln_F_t );

  // -- Gaussian random fields
  PARAMETER_ARRAY( ln_u_gt );            // State-matrix
  PARAMETER_VECTOR( Omegainput_g );      // Spatial variation in growth rates

  // Biological parameters
  PARAMETER_VECTOR( ln_biovec_z );

  // Global parameters
  using namespace density;
  Type pi = 3.141592;
  Type beta = exp( log_beta );
  Type rho = 1.0 - beta; 
  Type pos_penalty = 0;
  
  // Objective function
  vector<Type> NLL_c(6);    // 0:ln_u; 1:Epsilon; 2:Likelihood
  NLL_c.setZero();
  vector<Type> NLL_u_t(n_t);
  NLL_u_t.setZero();
  Type NLL = 0;                // Objective function
  vector<Type> NLL_i(n_i);
  
  // Derived parameters
  Type SigmaU = exp(logSigmaU);
  Type SigmaO = exp(logSigmaO);
  Type kappa_pow2 = exp(2.0*logkappa);
  Type kappa_pow4 = kappa_pow2*kappa_pow2;
  Type tauU = 1 / sqrt(4*pi*exp(2*logSigmaU)*exp(2*logkappa));
  Type tauO = 1 / sqrt(4*pi*exp(2*logSigmaO)*exp(2*logkappa));
  Type Range_raw = sqrt(8) / exp( logkappa );
  vector<Type> F_t(n_t);
  F_t = exp( ln_F_t );
  vector<Type> mvec(4);
  mvec.setZero();
 
  // Movement matrix
  Eigen::SparseMatrix<Type> M_sparse(n_r,n_r);
  Eigen::SparseMatrix<Type> Mdiv_sparse(n_r,n_r);
  Eigen::SparseMatrix<Type> Identity_sparse(n_r,n_r);
  if( Options_vec(0)==0 ){
    M_sparse.setIdentity();
    Identity_sparse.setIdentity();
    Mdiv_sparse.setIdentity();
  }
  if( Options_vec(0)!=0 ){
    if(Options_vec(4)==0) mvec = exp(log_mvec);
    if(Options_vec(4)==1){
      mvec(0) = exp(log_mvec(0));
      mvec(1) = exp(log_mvec(0)) * exp(log_mvec(1));
      mvec(2) = exp(log_mvec(0)) * exp(log_mvec(2));
      mvec(3) = exp(log_mvec(0)) * exp(log_mvec(3));
    }
    M_sparse = mvec(0)*M1 + mvec(1)*M2 + mvec(2)*M3 + mvec(3)*M4;
    Identity_sparse.setIdentity();
    Mdiv_sparse = M_sparse / Type(n_tdiv) + Identity_sparse; // Euler approximation
  }
  
  // Random field probability
  Eigen::SparseMatrix<Type> Q;
  Q = kappa_pow4*G0 + Type(2.0)*kappa_pow2*G1 + G2;
  GMRF_t<Type> gmrf_nll = GMRF(Q);

  // Derived fields
  vector<Type> Omega_g(n_g);
  Omega_g = Omegainput_g / tauO;

  // Track estimates and predictions of effort (dimensions n_r x n_t, i.e., density for each triangle in the domain)
  vector<Type> catchpred_t(n_t);
  catchpred_t.setZero();
  array<Type> catchpred_rt(n_r,n_t);
  catchpred_rt.setZero();
  array<Type> u_rt(n_r,n_t);
  array<Type> OFL_rt(n_r,n_t);
  array<Type> upred_rt(n_r,n_t); 
  // Track estimates and predictions of effort (dimensions n_gmrf x n_t, i.e., padded with zeros)
  array<Type> ln_uhat_gt(n_g,n_t);

  // Estimated density and catches
  for(int t=0; t<n_t; t++){  
    for(int r=0; r<n_r; r++){
      // Book-keeping
      u_rt(r,t) = exp( ln_u_gt(r,t) );
      catchpred_rt(r,t) = u_rt(r,t) * (1-exp(-1 * effortdens_rt(r,t) * F_t(t)));
    }
    // Calculate predicted total catch
    catchpred_t(t) = catchpred_rt.col(t).sum();
  }
  
  // Predict dynamics of state-vector
    // Survey is prior to effort in a given year
    // i.e., effort and catch in year t control change from u(t) to u(t+1)
    // Sequence: density(t) -> Survey(t) -> catch(t) -> movement -> production -> Process error(t) -> density(t+1)
  // First year -- Approximate equilibrium upred(r,0)
    // Assumes effortdens_rt(r,0) is equal to equilibrium effort
    // This implies that upred(r,0) ~= upred(r,1) for all r
  for(int r=0; r<n_r; r++){
    // Approximation to expected density
    if(Options_vec(1)==0) upred_rt(r,0) = km2_r(r) * exp( logmeanu0 + Omega_g(r)/beta );
    if(Options_vec(1)==1) upred_rt(r,0) = km2_r(r) * (exp(logmeanu0) + Omega_g(r))/beta;
    // Survival
    upred_rt(r,0) = upred_rt(r,0) * exp(-1 * effortdens_rt(r,0) * F_t(0));
  }  
  // Movement during initialization (Euler approximation)
  for(int tdev=0; tdev<n_tdiv; tdev++){
    upred_rt.col(0) = Mdiv_sparse * upred_rt.col(0).matrix();
  }
  // Check that it density is positive for all triangles (necessary for Ricker dynamics)
  for(int r=0; r<n_r; r++){
    ln_uhat_gt(r,0) = log( posfun(upred_rt(r,0), Type(1e-10), pos_penalty) );
  }
  // Fill out empty slots
  for(int g=n_r; g<n_g; g++) ln_uhat_gt(g,0) = 0.0;
  
  // Later years -- Calculate upred(r,t) given u_rt(r,t-1)
  for(int t=1; t<n_t; t++){
    // Survival
    for(int r=0; r<n_r; r++){
      upred_rt(r,t) = u_rt(r,t-1) * exp(-1 * effortdens_rt(r,t-1) * F_t(t-1));
    }
    // Movement (Euler approximation)
    for(int tdev=0; tdev<n_tdiv; tdev++){
      upred_rt.col(t) = Mdiv_sparse * upred_rt.col(t).matrix();
    }
    // Density-dependent function
    for(int r=0; r<n_r; r++){
      if(Options_vec(1)==0) ln_uhat_gt(r,t) = log( upred_rt(r,t) * exp(alpha + Omega_g(r) - beta*log( upred_rt(r,t)/km2_r(r) )) );
      if(Options_vec(1)==1) ln_uhat_gt(r,t) = log( upred_rt(r,t) * exp(alpha + Omega_g(r) - beta*upred_rt(r,t)/km2_r(r)) );
    }
    // Fill out empty slots
    for(int g=n_r; g<n_g; g++) ln_uhat_gt(g,t) = 0.0;
  }
  
  // Probability
  for(int t=0; t<n_t; t++){
    NLL_u_t(t) += SCALE(gmrf_nll, 1.0/tauU)( ln_u_gt.col(t)-ln_uhat_gt.col(t) );
  }
  NLL_c(0) = NLL_u_t.sum();

  // Epsilon probability
  NLL_c(1) = gmrf_nll(Omegainput_g);

  // Likelihood
  vector<Type> chat_i(n_i);
  vector<Type> encounterprob_i(n_i);
  for(int i=0; i<n_i; i++){
    // chat_i is the expectation for compound process including both encounter probability and positive catch rate components
    chat_i(i) = exp( ln_u_gt(r_i(i),t_i(i)) ) * km2_i(i) / km2_r(r_i(i));
    // Calculate deviance
    if( Options_vec(2)==0 ) NLL_i(i) = -1 * dpois(c_i(i), chat_i(i), true);
    if( Options_vec(2)!=0 ){
      // Calculate encounter probability (only used for delta-lognormal model)
      if( Options_vec(3)==0 ) encounterprob_i(i) = plogis(densvar_z(1) + densvar_z(2)*log(chat_i(i)) );
      if( Options_vec(3)==1 ) encounterprob_i(i) = plogis(densvar_z(1)) * ( 1.0 - exp(-1 * chat_i(i) * exp(densvar_z(2))) );
      if( Options_vec(3)==2 ) encounterprob_i(i) = plogis(densvar_z(3)) * plogis( densvar_z(1) + densvar_z(2) * log(chat_i(i)) );
      // probability of data
      if( Options_vec(2)==1 ) NLL_i(i) = -1 * dzinflognorm(c_i(i), log(chat_i(i)/encounterprob_i(i)), encounterprob_i(i), exp(densvar_z(0)), true);
      if( Options_vec(2)==2 ) NLL_i(i) = -1 * dzinfgamma(c_i(i), chat_i(i)/encounterprob_i(i), encounterprob_i(i), exp(densvar_z(0)), true);
      if( Options_vec(2)==3 ) NLL_i(i) = -1 * dzinfnorm(c_i(i), chat_i(i)/encounterprob_i(i), encounterprob_i(i), chat_i(i)/encounterprob_i(i)*exp(densvar_z(0)), true);
    }
  }
  NLL_c(2) = NLL_i.sum();
  
  // Catch penalty
  for(int t=0; t<n_t; t++){
    if( !isNA(catch_t(t)) ) NLL_c(3) -= dnorm( log(catch_t(t)), log(catchpred_t(t)), exp(ln_sigmacatch), true );
  }

  // Priors on biological parameters
  vector<Type> biovec_z(2);
  biovec_z = exp( ln_biovec_z );
  if( !isNA(bio_priors(0,2)) ) NLL_c(4) -= dlognorm( biovec_z(0), bio_priors(0,2), bio_priors(0,3), true ) + ln_biovec_z(0); // M
  if( !isNA(bio_priors(1,2)) ) NLL_c(4) -= dlognorm( biovec_z(1), bio_priors(1,2), bio_priors(1,3), true ) + ln_biovec_z(1); // M
  
  // Calculate OFL
  for(int r=0; r<n_r; r++){
  for(int t=0; t<n_t; t++){
    OFL_rt(r,t) = u_rt(r,t) * exp( ln_biovec_z(0) ) * exp( ln_biovec_z(1) );
  }}
  
  // Calculate derived quantities
  vector<Type> OFL_t(n_t);
  vector<Type> u_t(n_t);
  vector<Type> ln_OFL_t(n_t);
  vector<Type> ln_u_t(n_t);
  for(int t=0; t<n_t; t++){
    OFL_t(t) = OFL_rt.col(t).sum();
    u_t(t) = u_rt.col(t).sum();
  }
  ln_OFL_t = log( OFL_t );
  ln_u_t = log( u_t );
  
  // Total likelihood
  NLL_c(5) = pos_penalty;
  NLL = NLL_c.sum();

  // Diagnostic output
  REPORT( chat_i );
  REPORT( encounterprob_i );
  //REPORT( encounterprob_rt );
  REPORT( F_t );
  REPORT( catchpred_t );
  REPORT( catchpred_rt );
  REPORT( NLL_u_t );
  REPORT( u_rt );
  REPORT( ln_u_gt );
  REPORT( upred_rt );
  REPORT( ln_uhat_gt );
  REPORT( NLL_i );
  REPORT( Omega_g );
  REPORT( Range_raw );
  REPORT( SigmaU );
  REPORT( SigmaO );
  REPORT( tauU );
  REPORT( tauO );
  REPORT( NLL_c );
  REPORT( NLL );
  REPORT( rho );
  REPORT( log_mvec );
  REPORT( mvec );
  REPORT( densvar_z );
  REPORT( bio_priors );
  REPORT( OFL_rt );
  REPORT( alpha );
  REPORT( beta );
  REPORT( pos_penalty );  
  
  // Movement stuff
  REPORT( M_sparse );
  REPORT( Mdiv_sparse );

  // Unbiased predictors
  ADREPORT( u_rt );
  ADREPORT( OFL_rt );
  ADREPORT( biovec_z );

  // standard errors 
  ADREPORT( OFL_t );
  ADREPORT( u_t );  
  ADREPORT( rho );
  ADREPORT( alpha );
  ADREPORT( beta );
  ADREPORT( mvec );
  ADREPORT( ln_OFL_t );
  ADREPORT( ln_u_t );  

  return NLL;
}
