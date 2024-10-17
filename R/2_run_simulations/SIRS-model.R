#' function to simulate forced SIRS model 
#' (vectorized to simulate multiple parameter sets simultaneously)
#' @param N number of weeks to simulate
#' @param parms matrix of parameter values (e.g. from posterior distribution)
#'        dimensions should be simulations x parameters
#' @param pop population size
#' @param npi vector of NPI forcing; 1 if no NPI forcing
#'            dimensions should be simulations x parameters
#' @param rsvforce matrix of RSV forcing; all 1 if no RSV forcing
#' @param testing_scale vector of testing values (between 0 and 1)
SIRS_force_discrete <- function(N, parms, pop, npi = 1, rsvforce = 1, testing_scalar){
  if(sum(npi) == 1){npi = rep(1, N)}
  if(sum(rsvforce) == 1){rsvforce = matrix(0, nrow = N, ncol = nrow(parms)); parms <- cbind(parms, data.frame("c" = 1))}
  with(as.list(c(parms)),{
    S <- matrix(NA, nrow = N, ncol = nrow(parms))
    I <- matrix(NA, nrow = N, ncol = nrow(parms))
    R <- matrix(NA, nrow = N, ncol = nrow(parms))
    Ipred <- matrix(NA, nrow = N, ncol = nrow(parms))
    beta <- matrix(NA, nrow = N, ncol = nrow(parms))
    S[1,] = parms[,"S0"] * pop
    I[1,] = parms[,"I0"] * pop
    R[1,] = (1 - parms[,"S0"] - parms[,"I0"]) * pop
    Ipred[1,] = parms[,"I0"] * pop * parms[,"rho"]
    for (i in 2:N) {
      beta[i,] = npi[i]*parms[,"b"]*parms[,"r"]*(1 + parms[,"c"]*rsvforce[i,])*(parms[,"a"]*cos(2*pi*((i-40)/52 - parms[,"p"])) + 1)
      foi = beta[i,]*I[i-1,]/pop 
      # S transitions
      Soutall = (1 - exp(-foi - parms[,"mu"]))*S[i-1,] 
      Sdeath = Soutall*(parms[,"mu"]/(parms[,"mu"] + foi))
      Sout = Soutall - Sdeath 
      birth = (1-exp(-parms[,"mu"]))*pop
      # I
      Ioutall = (1 - exp(-parms[,"gamma"] - parms[,"mu"]))*I[i-1, ] 
      Ideath = Ioutall*(parms[,"mu"]/(parms[,"mu"] + parms[,"gamma"]))
      Iout = Ioutall - Ideath 
      # R
      Routall = (1 - exp(-parms[,"omega"] - parms[,"mu"]))*R[i-1, ] 
      Rdeath = Routall*(parms[,"mu"]/(parms[,"mu"] + parms[,"omega"])) 
      Rout = Routall - Rdeath
      # implement transitions
      S[i, ] = S[i-1, ] + birth - Sout + Rout - Sdeath
      I[i, ] = I[i-1, ] + Sout - Iout - Ideath
      R[i, ] = R[i-1, ] + Iout - Rout - Rdeath
      Ipred[i, ] = Sout * parms[,"rho"] * testing_scalar[i] 
    }
    return(bind_rows(melt(S) %>% mutate(variable = "S"), 
                     melt(I) %>% mutate(variable = "I"), 
                     melt(Ipred) %>% mutate(variable = "Ipred"), 
                     melt(R) %>% mutate(variable = "R"), 
                     melt(beta) %>% mutate(variable = "beta")) %>%
             rename(t = Var1, draw_id = Var2))
  })
}

#' function to simulate RSV and HPMV using RSV output as forcing for HMPV
#' @param N number of weeks to simulate
#' @param parms_rsv matrix of RSV parameters (e.g. from posterior distribution)
#'        dimensions should be simulations x parameters
#' @param parms_mpv matrix of HMPV parameters
#' @param pop double population size
#' @param testing_scalar_rsv vector of testing scalar for RSV (between 0 and 1)
#' @param testing_scalar_mpv vector of testing scalars for HMPV (between 0 and 1)
#' @param npi vector of NPI forcing; 1 if no NPI forcing
#' @param rsvscaling_cutoff integer representing t before which RSV incidence should
#'        be scaled as forcing for HMPV tranmsission; if NA, use all RSV incidence
SIRS_force_discrete_both = function(N, parms_rsv, parms_mpv, pop, 
                                     testing_scalar_rsv, testing_scalar_mpv,
                                     npi = 1, rsvscaling_cuttoff = NA){
  rsv_sim = SIRS_force_discrete(N = N, parms = parms_rsv, pop = pop, 
                                 npi = npi, rsvforce = 1, testing_scalar = testing_scalar_rsv)
  rsv_I = as.matrix(
    rsv_sim %>% 
      filter(variable == "I") %>%
      dcast(t ~ draw_id, value.var = "value") %>%
      select(-t)
  )
  if(is.na(rsvscaling_cuttoff)){
    if(ncol(rsv_I) == 1){rsv_max = max(rsv_I)}
    else{rsv_max = apply(rsv_I, 2, max)}
  }
  else{
    if(ncol(rsv_I) == 1){rsv_max = max(rsv_I)}
    else{rsv_max = apply(rsv_I[1:rsvscaling_cuttoff, ], 2, max)}
  }
  rsv_force_sim = t(t(rsv_I)/rsv_max)
  mpv_sim = SIRS_force_discrete(N = N, parms = parms_mpv, pop = pop, 
                                 npi = npi, rsvforce = rsv_force_sim, testing_scalar = testing_scalar_mpv)
  return(
    bind_rows(
      rsv_sim %>% mutate(pathogen = "rsv"), 
      mpv_sim %>% mutate(pathogen = "mpv"), 
    )
  )
}

