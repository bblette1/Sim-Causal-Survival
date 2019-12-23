##############################################################################
### simCausalSurvival()
### Description: Function to simulate data (as desribed in Young et al. 2010)
###              for which the assumptions of three causal models hold: a
###              Structural Nested Accelerated Failure Time Model, a Marginal
###              Structural Cox Model, and a Structural Nested Cumulative
###              Failure Time Model
###   Arguments:
###          n : number of subjects in simulated data set
###    tpoints : number of follow-up time points assuming no events
###        psi : treatment effect parameter value
###     lambda : hazard rate for potential time-to-event given never exposed
###   gamma0-3 : parameter values for confounder model
###   alpha0-3 : parameter values for treatment model
###
###     Returns: the function returns a data frame including columns for
###         ID : subject identifier variable
###      Delta : event indicator
###          A : time-varying treatment indicator
###          L : time-varying confounder indicator
### Delta_prev : indicator of event by previous time point
###     A_prev : treatment at previous time point
###     L_prev : confounder at previous time point
###       Time : event time
##############################################################################
simCausalSurvival <- function(n = 1000, tpoints = 10, psi = 0.3,
                              lambda = 0.05, gamma0 = 0, gamma1 = 0.5,
                              gamma2 = 0.2, gamma3 = 0.3, alpha0 = 0.1,
                              alpha1 = 0.3, alpha2 = 0.5, alpha3 = 0.2) {
  #######################################################################
  
  # Simulate potential time under no exposure
  T_0 <- rexp(n, lambda)
  
  # Initialize confounder, treatment, and outcome arrays
  L <- array(NA, dim = c(n, tpoints))
  A <- array(NA, dim = c(n, tpoints))
  Y <- array(NA, dim = c(n, tpoints))
  L_prev <- array(NA, dim = c(n, tpoints))
  A_prev <- array(NA, dim = c(n, tpoints))
  Y_prev <- array(NA, dim = c(n, tpoints))
  check_times <- array(NA, dim = c(n, tpoints))
  Times <- array(NA, dim = c(n, tpoints))

  # Do first measurements separately
  L_prev[, 1] <- rep(0, n)
  A_prev[, 1] <- rep(0, n)
  Y_prev[, 1] <- rep(0, n)
  
  L_logit <- gamma0 + gamma1*(T_0 < 0.5)
  L_probs <- exp(L_logit) / (1 + exp(L_logit))
  L[, 1] <- rbinom(n, 1, L_probs)
  
  A_logit <- alpha0 + alpha1*L[, 1]
  A_probs <- exp(A_logit) / (1 + exp(A_logit))
  A[, 1] <- rbinom(n, 1, A_probs)

  check_times[, 1] <- exp(psi*A[, 1])
  Y[, 1] <- as.numeric(T_0 <= check_times[, 1])
  Times[T_0 <= check_times[, 1], 1] <- 
    T_0[T_0 <= check_times[, 1]]*exp(-psi*A[T_0 <= check_times[, 1], 1])
  
  
  # Now generate second through m-th measurements in a loop
  for (j in 2:tpoints) {
    
    Y_prev[, j] <- Y[, (j - 1)]
    
    L_logit <- gamma0 + gamma1*(T_0 < 0.5) + gamma2*A[, (j - 1)] +
               gamma3*L[, (j - 1)]
    L_probs <- exp(L_logit) / (1 + exp(L_logit))
    L[, j] <- rbinom(n, 1, L_probs)
    
    A_logit <- alpha0 + alpha1*L[, j] + alpha2*L[, (j - 1)] +
               alpha3*A[, (j - 1)]
    A_probs <- exp(A_logit) / (1 + exp(A_logit))
    A[, j] <- rbinom(n, 1, A_probs)
    
    L[, j][Y_prev[, j] == 1] <- 0
    A[, j][Y_prev[, j] == 1] <- 0
    
    check_times[, j] <- check_times[, (j - 1)] + exp(psi*A[, j])
    
    Y[, j] <- as.numeric(T_0 <= check_times[, j])
    
    group1 <- (Y[, j] == 1 & Y_prev[, j] == 1)
    group2 <- (Y[, j] == 1 & Y_prev[, j] == 0)
    
    Times[, j][group1] <-
      Times[, (j - 1)][group1]
    
    Times[, j][group2] <-
      (j - 1 + ((T_0[group2] - check_times[, (j-1)][group2])*
                 exp(-psi*A[, j][group2])))
    
    #Times[, j] <- 
      #Times[, (j - 1)]*Y[, j]*Y_prev[, j] +
      #(j - 1 + ((T_0 - check_times[, (j-1)])*exp(-psi*A[,(j-1)])))*
      #(1 - Y[, (j - 1)])
    
    L_prev[, j] <- L[, (j - 1)] 
    A_prev[, j] <- A[, (j - 1)] 
    
  }
  
  # Make data set in counting process format
  ID <- rep(1:n, each = tpoints)
  #Time[is.na(Time)] <- tpoints
  #Start <- rep(1:tpoints-1, n)
  #Start[Start > c(t(Time)) & c(t(Time)) != 0] <- 
    #c(t(Time))[Start > c(t(Time)) & c(t(Time)) != 0]
  #Stop <- rep(1:tpoints, n)
  #Stop[c(t(Y)) == 1 & c(t(Time)) != 0] <- 
    #c(t(Time))[c(t(Y)) == 1 & c(t(Time)) != 0]
  
  data <- data.frame(ID, c(t(Y)), c(t(A)), c(t(L)), c(t(Y_prev)),
                     c(t(A_prev)), c(t(L_prev)), c(t(Times)))
  colnames(data) <- c("ID", "Delta", "A", "L", "Delta_prev", "A_prev",
                      "L_prev", "Time")
  
  return(data)
  #return(data[data$Start != data$Stop, ])
  
}
