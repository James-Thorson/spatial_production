
#  Matches mesh_model_v5j
SpatProdSimulator_OriginalTiming_Fn = function( MoveMat, SD_omega=1, SD_epsilon=1, SD_effort=1, CV_obs=1, effort_par=c(0.2,0.5), sizepar=c(1,0.5), Scale, Dynamical_Model, n_s, n_t, r_s, n_r, loc_r, logmeanu0, alpha, beta, km2_r ){
  # Load library
  require( RandomFields )

  # Simulate
  RF_omega = RMgauss(var=SD_omega^2, scale=Scale)
  RF_epsilon = RMgauss(var=SD_epsilon^2, scale=Scale)
  RF_effort = RMgauss(var=SD_effort^2, scale=Scale)

  # Simulate effort
  effortdens_t = rlnorm( n_t, meanlog=log(effort_par[1])-effort_par[2]^2/2, sdlog=effort_par[2] )
  effortdens_r = exp( RFsimulate(model=RF_effort, x=loc_r[,1], y=loc_r[,2])@data[,1] - SD_effort^2/2 )
  effortdens_rt = outer(effortdens_r, rep(1,n_t)) * outer(rep(1,n_r), effortdens_t)

  # Simulate density for each triangle
  # Gompertz: u(t+1) = u(t) * exp( alpha - beta*log(u(t)) )
  # Moran-Ricker: u(t+1) = u(t) * exp( alpha - beta*u(t) )
  catch_rt = upred_rt = u_rt = Epsilon_rt = matrix(NA, ncol=n_t+n_years_burnin, nrow=n_r)
  Omega_r = RFsimulate(model=RF_omega, x=loc_r[,1], y=loc_r[,2])@data[,1] - SD_omega^2/2
  for(t in 1:n_t){
    Epsilon_rt[,t] = RFsimulate(model=RF_epsilon, x=loc_r[,1], y=loc_r[,2])@data[,1]
    if(t==1){
      # Approximate stationary density
      if( Dynamical_Model=="Gompertz" ) u_rt[,t] = km2_r * exp( logmeanu0 + Omega_r )
      if( Dynamical_Model=="Ricker" ) u_rt[,t] = km2_r * exp( logmeanu0 ) + Omega_r
      # Fishing mortality and catch
      catch_rt[,t] = (1 - exp(-effortdens_rt[,t]) ) * u_rt[,t]
      u_rt[,t] = exp(-effortdens_rt[,t]) * u_rt[,t]
      # Movement
      u_rt[,t] = as.vector( MoveMat %*% u_rt[,t] )
      # Process error
      u_rt[,t] = u_rt[,t] * exp( Epsilon_rt[,t] )
    }
    if(t>=2){
      # Fishing effort
      catch_rt[,t] = (1 - exp(-effortdens_rt[,t]) ) * u_rt[,t-1]
      upred_rt[,t] = exp(-effortdens_rt[,t]) * u_rt[,t-1]
      # Movement
      upred_rt[,t] = as.vector( MoveMat %*% upred_rt[,t] )
      # Production
      if( Dynamical_Model=="Gompertz" ) u_rt[,t] = upred_rt[,t] * exp(alpha + Omega_r - beta*log(upred_rt[,t]/km2_r) + Epsilon_rt[,t])
      if( Dynamical_Model=="Ricker" ) u_rt[,t] = upred_rt[,t] * exp(alpha + Omega_r - beta*(upred_rt[,t]/km2_r) + Epsilon_rt[,t])
    }
  }

  # Simulate samples for each site and year
  DF = expand.grid("s_i"=1:n_s, "t_i"=1:n_t)
  DF = cbind( DF, "r_i"=r_s[DF[,'s_i']], "km2_i"=1 )
  DF = cbind( DF, "cexp_i"=u_rt[ as.matrix(DF[,c('r_i','t_i')]) ] * DF[,'km2_i']/km2_r[DF[,'r_i']] )
  DF = cbind( DF, "c_i"=rpois(n_s*n_t, lambda=DF[,'cexp_i']) )
  DF = cbind( DF, "encounterprob_i"=1-exp(-DF[,'cexp_i']) )
  logSD_obs = sqrt( log(CV_obs^2+1) )
  DF = cbind( DF, "zinflognorm_i"=ifelse(DF[,'c_i']>0,1,0) * rlnorm(nrow(DF), meanlog=log(DF[,'cexp_i']/DF[,'encounterprob_i']), sdlog=logSD_obs) )
  DF = cbind( DF, "zinfgamma_i"=ifelse(DF[,'c_i']>0,1,0) * rgamma(nrow(DF), shape=CV_obs^(-2), scale=DF[,'cexp_i']/DF[,'encounterprob_i']*CV_obs^2) )

  # Total catches
  catch_t = colSums( catch_rt )

  # Return stuff
  SettingsList = list( MoveMat=MoveMat, SD_omega=SD_omega, SD_epsilon=SD_epsilon, SD_effort=SD_effort, CV_obs=CV_obs, effort_par=effort_par, sizepar=sizepar, Scale=Scale, Dynamical_Model=Dynamical_Model, n_s=n_s, n_t=n_t, r_s=r_s, n_r=n_r, loc_r=loc_r, logmeanu0=logmeanu0, alpha=alpha, beta=beta, km2_r=km2_r )
  Return = list("SettingsList"=SettingsList, "DF"=DF, "catch_t"=catch_t, "catch_rt"=catch_rt, "effortdens_rt"=effortdens_rt, "MoveMat"=MoveMat, "upred_rt"=upred_rt, "u_rt"=u_rt, "Epsilon_rt"=Epsilon_rt, "Omega_r"=Omega_r)
  return( Return )
}

