#' @param rsv_force list of functions for age-specific rsv_force values
seirs_age_structured <- function(t, x, parms, compartments, age_classes, mort, 
                                 fert, waifw, rsv_force = NA, 
                                 adjust_beta_flag = FALSE, print_warnings_flag = FALSE){
  with(as.list(parms),{
    # if(any(x < 0)){print("X NEG!")}
    nage = length(age_classes)
    ncomp = length(compartments)
    aging <- 1/diff(c(0,age_classes))#/52 # for weekly time steps
    x_mat = matrix(x, ncol = ncomp, nrow = nage, dimnames = list(age_classes, compartments), byrow = TRUE)
    # calculate age-specific transmission rates
    beta <- b*r*(a*cos(2*pi*((t-40)/52-p))+1)
    lambda <- waifw%*%((beta/sum(x)*(x_mat[, "I"])))
    if(adjust_beta_flag){
      beta_hat = sum(lambda * x_mat[, "S"])/ (sum(x_mat[, "I"]) * sum(x_mat[, "S"])/sum(x))
      lambda = lambda * (beta)/beta_hat
    }
    # adjust lambda for interactions
    if(any(is.na(rsv_force))){
      f = rep(0, length(age_classes))
    }
    else{
      f = c * unlist(lapply(rsv_force, function(i){i(t)})) / rowSums(x_mat)
      if(any(f < -1)){
        if(print_warnings_flag){print(paste0("some lambdas less than 0, t: ", t))}
        f[which(f < -1)] = -1+1e-6
      }
    }
    lambda = lambda * (1 + f)
    # fertility
    Fmat <- buildFMatrix(age.classes = age_classes, fert = fert, ncompartments = ncomp, time.step = 1)
    N_fert <- Fmat %*% x
    N_fert <- matrix(N_fert, ncol = ncomp, nrow = nage, dimnames = list(age_classes, compartments))
    N_age <- x_mat * aging
    if(length(age_classes) == 1){
      N_age_in = matrix(0, 1, ncomp, dimnames = list(age_classes, compartments))
      N_age_out = matrix(0, 1, ncomp, dimnames = list(age_classes, compartments))
    }
    else{
      N_age_in <- rbind(rep(0, ncomp), N_age[1:(nrow(N_age)-1), ])
      N_age_out <- rbind(N_age[1:(nrow(N_age)-1),], rep(0, ncomp))
    }
    # calculate age-specific derivatives
    dS <- - (lambda + mort) * x_mat[, "S"] + omega * x_mat[, "R"] +
      N_fert[, "S"] + N_age_in[, "S"] - N_age_out[, "S"]
    dE <- lambda*x_mat[, "S"] - (mort + sigma) * x_mat[, "E"] +
      N_fert[, "E"] + N_age_in[, "E"] - N_age_out[, "E"]
    dI <- sigma * x_mat[, "E"] - (mort + gamma) * x_mat[, "I"] +
      N_fert[, "I"] + N_age_in[, "I"] - N_age_out[, "I"]
    dR <- gamma * x_mat[, "I"] - (mort + omega) * x_mat[, "R"] +
      N_fert[, "R"] + N_age_in[, "R"] - N_age_out[, "R"]
    der <- c(dS, dE, dI, dR)
    names(der) <- sapply(compartments, function(i){paste0(i, "_", age_classes)})
    der <- matrix(c(dS, dE, dI, dR), 
                  nrow = nage, ncol = ncomp)
    der <- c(t(der))
    names(der) <- sapply(age_classes, function(i){paste0(compartments, "_", i)})
    return(list(der))
  })
} 

#' max_t in weeks
run_ode <- function(age_classes, mort, fert, start_pop, 
                    waifw = NA, IC_type, max_t, params_rsv, params_mpv, 
                    beep_flag = FALSE, adjust_beta_flag = FALSE, print_warnings_flag = FALSE,
                    plot_flag = FALSE, plot_title = NA){
  IC = setup_IC(start_pop, age_classes, compartments, mort, fert, IC_type)
  if(any(is.na(waifw))){
    waifw = matrix(1, length(age_classes), length(age_classes)) 
  }
  # run model
  times = seq(1, max_t, 1)
  rslts_rsv <- as.data.frame(
    ode(
      y = IC,
      times = times,
      func = seirs_age_structured,
      compartments = compartments,
      age_classes = age_classes,
      mort = mort, 
      fert = fert,
      waifw = waifw,
      parms = params_rsv, 
      rsv_force = NA, 
      adjust_beta_flag = adjust_beta_flag, 
      print_warnings_flag = print_warnings_flag
    ))
  rsv_force = lapply(seq(4, 4*length(age_classes), 4), 
                     function(i){approxfun(rslts_rsv[,1], rslts_rsv[,i], yright = 0)})
  rslts_mpv <- as.data.frame(
    ode(
      y = IC,
      times = times,
      func = seirs_age_structured,
      compartments = compartments,
      age_classes = age_classes,
      mort = mort, 
      fert = fert,
      waifw = waifw,
      parms = params_mpv, 
      rsv_force = rsv_force, 
      adjust_beta_flag = adjust_beta_flag, 
      print_warnings_flag = print_warnings_flag
    ))
  if(beep_flag){beep()}
  return(process_results(rslts_rsv, rslts_mpv, plot_flag, plot_title, max_t))
}


buildFMatrix <- function(age.classes=c(1:60, seq(72,120,by=12), seq(180,600,by=60)),  
                         fert =  c(rep(0,66), rep(0.1,7)),   	#one for every age class. start reproducing at age 20, and here assume ~ flat
                         ncompartments, time.step, maternal_immunity_flag = FALSE){
  nage <- length(age.classes)
  Fmat <- matrix(0,ncompartments*nage,ncompartments*nage)
  birth_compartment = 1#ifelse(maternal_immunity_flag, 1, 2)
  for (j in 1:nage) {
    # same fertility for all compartments, all babies have maternal immunity or born into susceptible (based on birth_compartment)
    Fmat[birth_compartment,((j-1)*ncompartments+1):(j*ncompartments)] <- rep(fert[j]*time.step, ncompartments)    
  }
  return(Fmat)
}

findStableStruct <- function(age.classes=c(1:60,seq(72,120,by=12),seq(180,600,by=60)), 
                             mort=c(rep(1e-9,72),1), 
                             fert =  c(rep(0,66),rep(0.1,7)), time.step = 1){
  nage <- length(age.classes)
  # aging.rate <- time.step/diff(c(0,age.classes))
  aging.rate <- 1/diff(c(0,age_classes))/365 # for daily time steps
  Fmat <- Tmat <- matrix(0,nage,nage)
  for (j in 1:(nage-1)) { 
    Tmat[j,j] <- (1-mort[j]*time.step)*(1-aging.rate[j])
    Tmat[j+1,j] <- (1-mort[j]*time.step)*aging.rate[j]
    Tmat[j,j] <- (1-mort[j])*(1-aging.rate[j])
    Tmat[j+1,j] <- (1-mort[j])*aging.rate[j]
  }
  j <- nage	
  Tmat[j,j] <- (1-mort[j]*time.step)
  Fmat[1,] <- fert*time.step
  # calculate equilibrium values
  stable.age <- Re(eigen(Tmat+Fmat)$vector[,1])
  stable.age <- stable.age/sum(stable.age)
  lambda <- Re(eigen(Tmat+Fmat)$value[1])
  reprod.value <- Re(eigen(Tmat+Fmat)$vector[1,])
  return(list(stable.age = stable.age, lambda = lambda,
              reprod.value = reprod.value, age.classes = age.classes))
}

# std give 2000 infected individuals across all PA and PB
setup_IC <- function(start_pop, age_classes, compartments, mort, fert, type = "std"){
  # indexing - the rows for maternal, susceptible, etc
  indx_comp = rep(compartments, length(age_classes))
  if(type == "std"){
    # setup initial conditions
    IC <- rep(0,length(age_classes)*length(compartments))
    names(IC) = paste0(rep(compartments, length(age_classes)), "_", 
                       sort(rep(age_classes, length(compartments))))
    IC[which(indx_comp == "S")] = start_pop/length(age_classes) 			 # susceptibles
    IC[which(indx_comp == "I")] = 500/length(age_classes)
    IC[which(indx_comp == "S")] = IC[which(indx_comp == "S")] - 500/length(age_classes)
    return(IC)
  }
  if(type == "stable-age"){
    IC <- rep(0,length(age_classes)*length(compartments))
    names(IC) = paste0(rep(compartments, length(age_classes)), "_", 
                       sort(rep(age_classes, length(compartments))))
    expected_stable <- findStableStruct(age.classes = age_classes, mort = mort, fert = fert, time.step = 1)
    IC[which(indx_comp == "S")] = start_pop*expected_stable$stable.age 			 # susceptibles
    IC[which(indx_comp == "I")] = IC[which(indx_comp == "S")]*0.02
    IC[which(indx_comp == "S")] = IC[which(indx_comp == "S")]*0.98
    if(any(IC < 0)){browser()}
    if(abs(sum(IC) - start_pop) > 1e-6){browser()}
    return(IC)
  }
}

process_results <- function(rslts_rsv, rslts_mpv, plot_flag = FALSE, 
                            plot_title = NA, max_t){
  rslts_long <- bind_rows(
    melt(rslts_rsv, c("time")) %>%
      tidytable::separate(variable, into = c("variable", "age")) %>% 
      mutate(pathogen = "rsv"), 
    melt(rslts_mpv, c("time")) %>%
      tidytable::separate(variable, into = c("variable", "age")) %>% 
      mutate(pathogen = "mpv")
  )
  if(plot_flag){
    rslts_tot <- rslts_long %>% 
      filter(variable == "I") %>%
      summarize(value = sum(value), .by = c("variable", "time", "pathogen"))
    p <- ggplot(data = rslts_tot %>% filter(time > max_t-20*52), 
                aes(x = time, y = value, color = pathogen)) + 
      geom_line() +
      facet_wrap(vars(variable), scales = "free") + 
      theme_bw()
    if(!is.na(plot_title)){
      p <- p + ggtitle(plot_title)
    }
    print(p)
  }
  return(rslts_long)
}

create_polymod_matrix = function(age_classes, plot_flag = FALSE, 
                                 age_classes_to_label = c(seq(12, 60, 12), seq(120, 840, 120))){
  # create polymod using code from Bjornstad book
  data(polymod)
  x = y = polymod$contactor[1:30]
  z = matrix(polymod$contact.rate, ncol = 30, nrow = 30)
  n = length(x)
  # symmetrize
  z2 = (z + t(z))/2
  z3 = as.vector(z2)
  xy = data.frame(x = rep(x[1:n], n), y = rep(y[1:n], each = n))
  polysmooth = Tps(xy, z3, df = 100)
  # surface(polysmooth, xlab = "", ylab = "", col = gray((12:32)/32))
  # annualize & symmetrize
  ps = predict(polysmooth, x = expand.grid(age_classes/12, age_classes/12))
  ps2 = matrix(ps, ncol = length(age_classes))
  ps2 = ps2 + t(ps2)
  W = ps2/mean(ps2)
  if(plot_flag){
    # plot W matrix
    p <- ggplot(data = melt(W), aes(x = Var1, y = Var2, fill = value)) + 
      geom_tile() + 
      scale_fill_viridis_c() + 
      scale_x_continuous(expand = c(0,0),
                         breaks = which(age_classes %in% age_classes_to_label),
                         labels = age_classes_to_label/12,
                         name = "age (years)") +
      scale_y_continuous(expand = c(0,0),
                         breaks = which(age_classes %in% age_classes_to_label),
                         labels = age_classes_to_label/12,
                         name = "age (years)")
    print(p)
  }
  return(W)
}

