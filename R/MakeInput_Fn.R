
MakeInput_Fn = function( Version, Options_vec, DF, n_r, n_g, n_t, n_tdiv, MoveList, Catch_DF, km2_r, EffortDens_rt, spde_gmrf, Bio_Priors, n_years_burnin ){

  # Data inputs
  if(Version=="mesh_model_v1") Data = list("n_i"=nrow(DF), "n_r"=n_r, "n_g"=n_g, "n_t"=n_t, "n_tdiv"=n_tdiv, "c_i"=DF[,'c_i'], "r_i"=DF[,'r_i']-1, "t_i"=DF[,'t_i']-1, "M1"=inla.as.dgTMatrix(MoveList[["M1"]]), "M2"=inla.as.dgTMatrix(MoveList[["M2"]]), "M3"=inla.as.dgTMatrix(MoveList[["M3"]]), "M4"=inla.as.dgTMatrix(MoveList[["M4"]]), "G0"=spde_gmrf$param.inla$M0, "G1"=spde_gmrf$param.inla$M1, "G2"=spde_gmrf$param.inla$M2)
  if(Version%in%c("mesh_model_v3b","mesh_model_v3","mesh_model_v2")) Data = list("Options_vec"=Options_vec, "n_i"=nrow(DF), "n_r"=n_r, "n_g"=n_g, "n_t"=length(unique(DF[,'t_i'])), "n_tdiv"=n_tdiv, "c_i"=DF[,'c_i'], "r_i"=DF[,'r_i']-1, "t_i"=DF[,'t_i']-min(DF[,'t_i']), "M1"=inla.as.dgTMatrix(MoveList[["M1"]]), "M2"=inla.as.dgTMatrix(MoveList[["M2"]]), "M3"=inla.as.dgTMatrix(MoveList[["M3"]]), "M4"=inla.as.dgTMatrix(MoveList[["M4"]]), "G0"=spde_gmrf$param.inla$M0, "G1"=spde_gmrf$param.inla$M1, "G2"=spde_gmrf$param.inla$M2)
  if(Version%in%c("mesh_model_v4")) Data = list("Options_vec"=Options_vec, "n_i"=nrow(DF), "n_r"=n_r, "n_g"=n_g, "n_t"=length(unique(DF[,'t_i'])), "n_tdiv"=n_tdiv, "c_i"=DF[,'c_i'], "r_i"=DF[,'r_i']-1, "t_i"=DF[,'t_i']-min(DF[,'t_i']), "catch_t"=Catch_DF[,'TWL_catch'], "effort_rt"=Effort_rt, "M1"=inla.as.dgTMatrix(MoveList[["M1"]]), "M2"=inla.as.dgTMatrix(MoveList[["M2"]]), "M3"=inla.as.dgTMatrix(MoveList[["M3"]]), "M4"=inla.as.dgTMatrix(MoveList[["M4"]]), "G0"=spde_gmrf$param.inla$M0, "G1"=spde_gmrf$param.inla$M1, "G2"=spde_gmrf$param.inla$M2)
  if(Version%in%c("mesh_model_v4b")) Data = list("Options_vec"=Options_vec, "n_i"=nrow(DF), "n_r"=n_r, "n_g"=n_g, "n_t"=length(unique(DF[,'t_i'])), "n_tdiv"=n_tdiv, "c_i"=DF[,'c_i'], "r_i"=DF[,'r_i']-1, "t_i"=DF[,'t_i']-min(DF[,'t_i']), "km2_i"=DF[,'km2_i'], "km2_r"=km2_r, "catch_t"=Catch_DF[,'TWL_catch'], "effortdens_rt"=EffortDens_rt, "M1"=inla.as.dgTMatrix(MoveList[["M1"]]), "M2"=inla.as.dgTMatrix(MoveList[["M2"]]), "M3"=inla.as.dgTMatrix(MoveList[["M3"]]), "M4"=inla.as.dgTMatrix(MoveList[["M4"]]), "G0"=spde_gmrf$param.inla$M0, "G1"=spde_gmrf$param.inla$M1, "G2"=spde_gmrf$param.inla$M2)
  if(Version%in%c("mesh_model_v5l","mesh_model_v5k","mesh_model_v5j","mesh_model_v5i","mesh_model_v5h","mesh_model_v5g","mesh_model_v5f","mesh_model_v5e","mesh_model_v5d","mesh_model_v5c","mesh_model_v5b","mesh_model_v5")) Data = list("Options_vec"=Options_vec, "n_i"=nrow(DF), "n_r"=n_r, "n_g"=n_g, "n_t"=length(unique(DF[,'t_i'])), "n_tdiv"=n_tdiv, "c_i"=DF[,'c_i'], "r_i"=DF[,'r_i']-1, "t_i"=DF[,'t_i']-min(DF[,'t_i']), "km2_i"=DF[,'km2_i'], "km2_r"=km2_r, "catch_t"=Catch_DF[,'TWL_catch'], "effortdens_rt"=EffortDens_rt, "M1"=inla.as.dgTMatrix(MoveList[["M1"]]), "M2"=inla.as.dgTMatrix(MoveList[["M2"]]), "M3"=inla.as.dgTMatrix(MoveList[["M3"]]), "M4"=inla.as.dgTMatrix(MoveList[["M4"]]), "G0"=spde_gmrf$param.inla$M0, "G1"=spde_gmrf$param.inla$M1, "G2"=spde_gmrf$param.inla$M2, "bio_priors"=Bio_Priors)
  # Add burn-in-years if necessary
  if( n_years_burnin>0 ){
    Data[["n_t"]] = Data[["n_t"]] + n_years_burnin
    Data[["t_i"]] = Data[["t_i"]] + n_years_burnin
    Data[["catch_t"]] = c("burnin"=rep(Data[["catch_t"]][1],n_years_burnin), "observed"=Data[["catch_t"]])
    Data[["effortdens_rt"]] = cbind("burnin"=outer(Data[["effortdens_rt"]][,1],rep(1,n_years_burnin)), "observed"=Data[["effortdens_rt"]])
  }

  # Parameter inputs
    # Start beta very low so Ricker model starts at permissible location
  if(Version=="mesh_model_v1") Params = list("rho"=rho, "alpha"=alpha, "log_mvec"=log(mvec), "logkappa"=log(1), "logSigmaU"=sqrt(Var), "logSigmaO"=sqrt(Var), "logmeanu0"=logmeanu0, "ln_u_gt"=matrix(logmeanu0, nrow=n_g, ncol=n_t), "Omegainput_g"=rnorm(n_g))
  if(Version=="mesh_model_v2") Params = list("log_beta"=log(0.1), "alpha"=log(1), "log_mvec"=log(rep(0.2,4)), "logkappa"=log(1), "logSigmaU"=log(1), "logSigmaO"=log(1), "logmeanu0"=log(1), "ln_u_gt"=matrix( log(mean(Data$c_i)), nrow=n_g, ncol=n_t), "Omegainput_g"=rnorm(Data$n_g))
  if(Version%in%c("mesh_model_v3")) Params = list("log_beta"=log(0.1), "alpha"=log(1), "log_mvec"=log(rep(0.2,4)), "logkappa"=log(1), "logSigmaU"=log(1), "logSigmaO"=log(1), "logmeanu0"=log(1), "densvar_z"=c(log(1),0,0), "ln_u_gt"=matrix( log(mean(Data$c_i))+rnorm(Data$n_g*Data$n_t,mean=0,sd=0.1), nrow=Data$n_g, ncol=Data$n_t), "Omegainput_g"=rnorm(Data$n_g,mean=0,sd=0.1))
  if(Version%in%c("mesh_model_v3b")) Params = list("log_beta"=log(0.1), "alpha"=log(1), "log_mvec"=log(rep(0.2,4)), "logkappa"=log(1), "logSigmaU"=log(1), "logSigmaO"=log(1), "logmeanu0"=log(1), "densvar_z"=c(log(1),0,0,0), "ln_u_gt"=matrix( log(mean(Data$c_i))+rnorm(Data$n_g*Data$n_t,mean=0,sd=0.1), nrow=Data$n_g, ncol=Data$n_t), "Omegainput_g"=rnorm(Data$n_g,mean=0,sd=0.1))
  if(Version%in%c("mesh_model_v4")) Params = list("log_beta"=log(0.1), "alpha"=log(1), "log_mvec"=log(rep(0.2,4)), "logkappa"=log(1), "logSigmaU"=log(1), "logSigmaO"=log(1), "logmeanu0"=log(1), "densvar_z"=c(log(1),0,0,0), "ln_sigmacatch"=log(0.01), "ln_F_t"=log(rep(0.1,Data$n_t)), "ln_u_gt"=matrix( log(mean(Data$c_i))+rnorm(Data$n_g*Data$n_t,mean=0,sd=0.1), nrow=Data$n_g, ncol=Data$n_t), "Omegainput_g"=rnorm(Data$n_g,mean=0,sd=0.1))
  if(Version%in%c("mesh_model_v4b")) Params = list("log_beta"=log(0.01), "alpha"=log(1), "log_mvec"=log(rep(0.2,4)), "logkappa"=log(1), "logSigmaU"=log(1), "logSigmaO"=log(1), "logmeanu0"=log(1), "densvar_z"=c(log(1),0,0,0), "ln_sigmacatch"=log(0.01), "ln_F_t"=log(rep(0.1,Data$n_t)), "ln_u_gt"=matrix( log(mean(Data$c_i)*mean(Data$km2_r)/mean(Data$km2_i))+rnorm(Data$n_g*Data$n_t,mean=0,sd=0.1), nrow=Data$n_g, ncol=Data$n_t), "Omegainput_g"=rnorm(Data$n_g,mean=0,sd=0.1))
  if(Version%in%c("mesh_model_v5l","mesh_model_v5k","mesh_model_v5j","mesh_model_v5i","mesh_model_v5h","mesh_model_v5g","mesh_model_v5f","mesh_model_v5e","mesh_model_v5d","mesh_model_v5c","mesh_model_v5b","mesh_model_v5")) Params = list("log_beta"=log(0.01), "alpha"=log(1), "log_mvec"=log(rep(0.2,4)), "logkappa"=log(1), "logSigmaU"=log(1), "logSigmaO"=log(1), "logmeanu0"=log(1), "densvar_z"=c(log(100),0,0,0), "ln_sigmacatch"=log(0.01), "ln_F_t"=log(rep(0.1,Data$n_t)), "ln_u_gt"=matrix( log(mean(Data$c_i)*mean(Data$km2_r)/mean(Data$km2_i))+0.1*rnorm(Data$n_g*Data$n_t), nrow=Data$n_g, ncol=Data$n_t), "Omegainput_g"=0.1*rnorm(Data$n_g), "ln_biovec_z"=Bio_Priors[c('M','FoverM'),'meanlog'])
  if(Dynamical_Model=="Ricker") Params$logmeanu0 = 2

  # Define random effects
  Random = c("ln_u_gt", "Omegainput_g") # NULL #
  if( Version%in%c("mesh_model_v5l","mesh_model_v5k","mesh_model_v5j","mesh_model_v5i","mesh_model_v5h","mesh_model_v5g","mesh_model_v5f","mesh_model_v5e","mesh_model_v5d","mesh_model_v5c","mesh_model_v5b","mesh_model_v5","mesh_model_v4b","mesh_model_v4") ) Random = c( Random, "ln_F_t" )
  if( Version%in%c("mesh_model_v5l","mesh_model_v5k","mesh_model_v5j","mesh_model_v5i","mesh_model_v5h","mesh_model_v5g","mesh_model_v5f","mesh_model_v5e","mesh_model_v5d","mesh_model_v5c","mesh_model_v5b","mesh_model_v5") ) Random = c( Random, "ln_biovec_z" )

  # Make "Map" -- Turn off params
  Map = list()
  if( Version%in%c("mesh_model_v5l","mesh_model_v5k","mesh_model_v5j","mesh_model_v5i","mesh_model_v5h","mesh_model_v5g","mesh_model_v5f","mesh_model_v5e","mesh_model_v5d","mesh_model_v5c","mesh_model_v5b","mesh_model_v5","mesh_model_v4b","mesh_model_v4") ){
    Map[["ln_sigmacatch"]] = factor( NA )
  }
  # Turn off asymptote for encounter probability
  if( Options_vec["Encounter_fn"]!=2 & Version%in%c("mesh_model_v5l","mesh_model_v5j","mesh_model_v5i","mesh_model_v5h","mesh_model_v5g","mesh_model_v5f","mesh_model_v5e","mesh_model_v5d","mesh_model_v5c","mesh_model_v5b","mesh_model_v5","mesh_model_v4b","mesh_model_v4","mesh_model_v3b") ){
    Map[["densvar_z"]] = factor( c(1,2,3,NA) )
    Params[["densvar_z"]] = c( log(1),0,0,0 )
  }
  # Deal with movement options
  if( Options_vec["Include_movement"]==0 ){
    Map[["log_mvec"]] = factor( rep(NA,length(Params[["log_mvec"]])) )
    Params[["log_mvec"]] = rep(-1e20,length(Params[["log_mvec"]]))
  }
  if( Options_vec["Include_movement"]==1 ){
    Map[["log_mvec"]] = factor( c(1,1,1,1) )
    Params[["log_mvec"]] = rep(log(0.1),length(Params[["log_mvec"]]))
  }
  if( Options_vec["Include_movement"]==3 ){
    Map[["log_mvec"]] = factor( c(1,2,1,2) )
    Params[["log_mvec"]] = rep(log(0.1),length(Params[["log_mvec"]]))
  }
  if( Options_vec["Include_movement"]==4 ){
    Map[["log_mvec"]] = factor( rep(NA,4) )
    residprob_fn = function( log_mvec, residprob_params, log_mvec_parameterization, outputtype="objective" ){
      if(log_mvec_parameterization==0) Params[["log_mvec"]] = rep( log_mvec,length(Params[["log_mvec"]]))
      if(log_mvec_parameterization==1) Params[["log_mvec"]] = c( log_mvec, log(residprob_params[[2]]), 0, log(residprob_params[[2]]))
      setwd( TmbFile )
      dyn.load( dynlib(Version) )
      Obj = MakeADFun( data=Data, parameters=Params, map=Map, random=Random, inner.control=list(maxit=2500)) #, inner.method="Nelder-Mead"
      Report = Obj$report()
      M = matpow( Report$Mdiv_sparse, Data$n_tdiv )
      residprob = mean( diag(as.matrix(M)) )
      print( paste0("mvec=c(",paste(round(Report[["mvec"]],3),collapse=","),") residprob=",round(residprob,3)) )
      if( outputtype=="objective" ) Return = abs(residprob-residprob_params[[1]])
      if( outputtype=="log_mvec" ) Return = Params[["log_mvec"]]
      return( Return )
    }
    #residprob_fn( log_mvec=log(0.2), residprob_params=Target_Residency_Prob, log_mvec_parameterization=Options_vec[["Movement_parameterization"]])
    Opt = optimize( f=residprob_fn, interval=log(c(0.01,1)), residprob_params=Target_Residency_Prob, log_mvec_parameterization=Options_vec[["Movement_parameterization"]] )
    Params[["log_mvec"]] = residprob_fn( log_mvec=Opt$minimum, outputtype="log_mvec", residprob_params=Target_Residency_Prob, log_mvec_parameterization=Options_vec[["Movement_parameterization"]])
  }

  # Return stuff
  Return = list( "Data"=Data, "Params"=Params, "Random"=Random, "Map"=Map )
  return( Return )
}
