#--- Simulations for g-estimation of structural mean models ---
#--------------- An application to MR analysis ---------------- 

# Loading packages
if (!require('tidyverse')) install.packages('tidyverse'); library('tidyverse')


# ----------------------- Appendix 3.1  -----------------------
# Identifying the point effect using Mendelian randomization
#  - Causal estimand of interest: point effect
#  - Number of exposure measurements considered in the model: one
#  - Instrument-exposure relationship changes over time? Yes
#  - Effect of exposure modified by previous exposure? No
#  - Presence of time-varying outcome-exposure confounding? No
# -------------------------------------------------------------

set.seed(2543789)
n <- 10000

  # Specifying parameters
  alpha_ZA <- expand.grid(Z1=c(0.25, 0), Z2=0.25, Z3=c(0.25, 0))
  alpha_U <- 0.3
  alpha_AA <- expand.grid(A1_A2=c(0.3), A1_A3=c(0.3, 0), A2_A3=c(0.3))
  beta_AY <- expand.grid(B1=c(-0.5, 0), B2=0.1, B3=c(1, 0))

  # Simulation function
  simulation1 <- function(b){
    
    # Looping over different values of beta_AY
    f.beta_AY <- function(ay){
      
      #Looping over different values of alpha_AA
      f.alpha_AA <- function(aa){
        
        # Looping over different values of alpha_ZA
        f.alpha_ZA <- function(za){
          
          # Generating variables
          U <- rnorm(n, 0, 1)
          Z <- rbinom(n, size=2, prob=0.3)
          
          A1 <- Z*alpha_ZA[za,"Z1"] + U*alpha_U + rnorm(n, 0, 1)
          A2 <- Z*alpha_ZA[za,"Z2"] + U*alpha_U + A1*alpha_AA[aa,"A1_A2"] + rnorm(n, 0, 1)
          A3 <- Z*alpha_ZA[za,"Z3"] + U*alpha_U + A1*alpha_AA[aa,"A1_A3"] + A2*alpha_AA[aa,"A2_A3"] + rnorm(n, 0, 1)
          A <- A2
          Y <-  A1*beta_AY[ay, "B1"] + A2*beta_AY[ay, "B2"] + A3*beta_AY[ay, "B3"] + U*alpha_U + rnorm(n, 0, 1)
          
          # Analyzing data using SMMs
          predZ <- predict(lm(Z~1), type="response")
          smm.results <- t(Y)%*%(Z - predZ)%*%t(Z-predZ)%*%A%*%solve(t(A)%*%(Z-predZ)%*%t(Z-predZ)%*%A)
          names(smm.results) <- c("psi.hat")
          
          # Results
          return(unlist(c(alpha_ZA[za,], beta_AY[ay,], alpha_AA[aa,], smm.results)))
        }
        return(lapply(1:nrow(alpha_ZA), f.alpha_ZA) %>% do.call(rbind,.))
      }
      return(lapply(1:nrow(alpha_AA), f.alpha_AA) %>% do.call(rbind,.))
    }
    return(lapply(1:nrow(beta_AY), f.beta_AY) %>% do.call(rbind,.))
  }

  # Obtaining simulation results
  sim1.results <- lapply(1:1000, simulation1) %>% do.call(rbind,.)
  
  # Analyzing results
  sim1 <- sim1.results %>% 
    data.frame() %>%
    group_by(Z1, Z2, Z3, B1, B2, B3, A1_A2, A1_A3, A2_A3) %>%
    summarise_all(mean) %>%
    mutate(psi.true=B2+A2_A3*B3) %>%
    mutate_at(vars(c("psi.hat")), ~sprintf("%.3f", .)) %>%
    select(Z1, Z2, Z3, B1, B2, B3, A1_A2, A1_A3, A2_A3, psi.true, psi.hat)

  
# ----------------------- Appendix 3.2 -----------------------
# Identifying period effects using MR analysis: all exposure time
# time points are measured, and exposure-outcome confounders are 
# time-fixed
#  - Causal estimand of interest: period effect
#  - Number of exposure measurements considered in the model: all (3)
#  - Instrument-exposure relationship changes over time? Yes
#  - Effect of exposure modified by previous exposure? No
#  - Presence of time-varying outcome-exposure confounding? No
# -------------------------------------------------------------
  
set.seed(123543)
n <- 10000

  # Specifying parameters
  alpha_U <- 0.5
  alpha_AA <- expand.grid(A1_A2=c(0.3, 0), A1_A3=c(0.3, 0), A2_A3=c(0.3, 0))
  alpha_ZA <- cbind(c(0.5, 0.3, 0.1),
                    c(0.1, 0.4, 0.2),
                    c(0.2, 0.5, 0.3),
                    c(0.3, 0.1, 0.4),
                    c(0.4, 0.2, 0.5))
  beta_AY <- c(-0.5, 0.1, 1)

  # Data Generation
  simulation2 <- function(i){
    
    # Looping over different values of alpha_AA
    f.alpha_AA <- function(aa){
      
      # Generating variables
      Z <- lapply(1:5,function(x) {
        return(rbinom(n = n, size = 2, prob = 0.3))
      }) %>% do.call(cbind,.)
      
      U <- rnorm(n, 0, 1)
      A1 <- Z %*% alpha_ZA[1,] + alpha_U*U + rnorm(n, 0, 1)
      A2 <- Z %*% alpha_ZA[2,] + alpha_U*U + rnorm(n, 0, 1) + alpha_AA[aa,"A1_A2"]*A1
      A3 <- Z %*% alpha_ZA[3,] + alpha_U*U + rnorm(n, 0, 1) + alpha_AA[aa,"A1_A3"]*A1 + alpha_AA[aa,"A2_A3"]*A2
      A <- cbind(A1, A2, A3)
      
      Y <- A1*beta_AY[1] + A2*beta_AY[2] + A3*beta_AY[3] + U*alpha_U + rnorm(n, 0, 1)
      
      # Analyzing data using SMMs
      predZ <- apply(Z, 2, function(i){predict(lm(i~1), type="response")})
      smm.results <- t(Y)%*%(Z - predZ)%*%t(Z-predZ)%*%A%*%solve(t(A)%*%(Z-predZ)%*%t(Z-predZ)%*%A)
      names(smm.results) <- c("psi1", "psi2", "psi3")
      
      # Analyzing data using multivariable MR
      predA <- predict(lm(A~Z), type="response")
      multimr.results <- lm(Y~predA)$coef[-1]
      names(multimr.results) <- c("delta1", "delta2", "delta3")
      
      # Returning results
      return(unlist(c(alpha_AA[aa,], smm.results, multimr.results)))
    }
    return(lapply(1:nrow(alpha_AA), f.alpha_AA) %>% do.call(rbind,.))
  }

  # Obtaining simulation results
  sim2.results <- lapply(1:1000, simulation2) %>% do.call(rbind,.)
  
  # Analyzing results
  sim2 <- sim2.results %>% 
    data.frame() %>%
    mutate(psi_sum=psi1+psi2+psi3) %>%
    mutate(delta_sum=delta1+delta2+delta3) %>%
    group_by(A1_A2, A1_A3, A2_A3) %>%
    summarise_all(mean) %>%
    mutate_at(vars(c(contains("psi"), contains("delta"))), ~sprintf("%.3f", .)) %>%
    relocate(psi_sum, .after=psi3)

  
# ----------------------- Appendix 3.3 -----------------------
# Identifying period effects using MR analysis: all exposure time
# time points are measured, and exposure-outcome confounders are 
# time-varying
#  - Causal estimand of interest: period effect
#  - Number of exposure measurements considered in the model: all (3)
#  - Instrument-exposure relationship changes over time? Yes
#  - Effect of exposure modified by previous exposure? No
#  - Presence of time-varying outcome-exposure confounding? Yes
# -------------------------------------------------------------

set.seed(123543)
n <- 10000

  # Specifying parameters
  alpha_U <- 0.5
  alpha_AA <- expand.grid(A1_A2=c(0.3, 0), A1_A3=c(0.3, 0), A2_A3=c(0.3, 0))
  alpha_ZA <- cbind(c(0.5, 0.3, 0.1),
                    c(0.1, 0.4, 0.2),
                    c(0.2, 0.5, 0.3),
                    c(0.3, 0.1, 0.4),
                    c(0.4, 0.2, 0.5))
  beta_AY <- c(-0.5, 0.1, 1)

  # Simulation function
  simulation3 <- function(i){
    
    # Looping over different values of alpha_AA
    f.alpha_AA <- function(aa){
      
      # Generating variables
      Z <- lapply(1:5,function(x) {
        return(rbinom(n = n, size = 2, prob = 0.3))
      }) %>% do.call(cbind,.)
      
      U1 <- rnorm(n, 0, 1)
      A1 <- Z %*% alpha_ZA[1,] + alpha_U*U1 + rnorm(n, 0, 1)
      U2 <- alpha_U*A1 + rnorm(n, 0, 1)
      A2 <- Z %*% alpha_ZA[2,] + alpha_U*U2 + rnorm(n, 0, 1) + alpha_AA[aa,"A1_A2"]*A1
      U3 <- alpha_U*A2 + rnorm(n, 0, 1)
      A3 <- Z %*% alpha_ZA[3,] + alpha_U*U3 + rnorm(n, 0, 1) + alpha_AA[aa,"A1_A3"]*A1 + alpha_AA[aa,"A2_A3"]*A2
      A <- cbind(A1, A2, A3)
      U <- cbind(U1, U2, U3)
      
      Y <- U %*% rep(alpha_U, 3) + beta_AY[1]*A1 + beta_AY[2]*A2 + beta_AY[3]*A3 + rnorm(n, 0, 1)
      
      # Analyzing data using SMM
      predZ <- apply(Z, 2, function(i){predict(lm(i~1), type="response")})
      smm.results <- t(Y)%*%(Z - predZ)%*%t(Z-predZ)%*%A%*%solve(t(A)%*%(Z-predZ)%*%t(Z-predZ)%*%A)
      names(smm.results) <- c("psi1", "psi2", "psi3")
      
      # Analyzing data using multivariable MR
      predA <- predict(lm(A~Z), type="response")
      multimr.results <- lm(Y~predA)$coef[-1]
      names(multimr.results) <- c("delta1", "delta2", "delta3")
      
      # Returning results
      return(unlist(c(alpha_AA[aa,], smm.results, multimr.results)))
    }
    return(lapply(1:nrow(alpha_AA), f.alpha_AA) %>% do.call(rbind,.))
  }

  # Obtaining simulation results
  sim3.results <- lapply(1:1000, simulation3) %>% do.call(rbind,.)
  
  # Analyzing results
  sim3 <- sim3.results %>% 
    data.frame() %>%
    mutate(psi_sum=psi1+psi2+psi3) %>%
    mutate(delta_sum=delta1+delta2+delta3) %>%
    group_by(A1_A2, A1_A3, A2_A3) %>%
    summarise_all(mean) %>%
    mutate_at(vars(c(contains("psi"), contains("delta"))), ~sprintf("%.3f", .)) %>%
    relocate(psi_sum, .after=psi3)

  
# ----------------------- Appendix 3.4 -----------------------
# Identifying period effects using MR analysis: interaction
# between exposures
#  - Causal estimand of interest: period effect
#  - Number of exposure measurements considered in the model: all (3)
#  - Instrument-exposure relationship changes over time? Yes
#  - Effect of exposure modified by previous exposure? Yes
#  - Presence of time-varying outcome-exposure confounding? No
# -------------------------------------------------------------

set.seed(653146)
n <- 10000

  # Specifying parameters
  alpha_U <- 0.5
  alpha_AA <- expand.grid(A1_A2=c(0.3, 0), A1_A3=c(0.3, 0), A2_A3=c(0.3, 0))
  alpha_ZA <- cbind(c(0.5, 0.3, 0.1),
                    c(0.1, 0.4, 0.2),
                    c(0.2, 0.5, 0.3),
                    c(0.3, 0.1, 0.4),
                    c(0.4, 0.2, 0.5))
  beta_AY <- c(-0.5, 0.1, 1, 0.3)
  
  # Simulation function
  simulation4 <- function(i){
    
    # Looping over different values of alpha_AA
    f.alpha_AA <- function(aa){
      
      # Generating variables
      Z <- lapply(1:5,function(x) {
        return(rbinom(n = n, size = 2, prob = 0.3))
      }) %>% do.call(cbind,.)
      
      U <- rnorm(n, 0, 1)
      A1 <- Z %*% alpha_ZA[1,] + alpha_U*U + rnorm(n, 0, 1)
      A2 <- Z %*% alpha_ZA[2,] + alpha_U*U + rnorm(n, 0, 1) + alpha_AA[aa,"A1_A2"]*A1
      A3 <- Z %*% alpha_ZA[3,] + alpha_U*U + rnorm(n, 0, 1) + alpha_AA[aa,"A1_A3"]*A1 + alpha_AA[aa,"A2_A3"]*A2
      A <- cbind(A1, A2, A3, A1*A2)
      
      Y <- A1*beta_AY[1] + A2*beta_AY[2] + A3*beta_AY[3] + A1*A2*beta_AY[4] + U*alpha_U + rnorm(n, 0, 1)
      
      # Analyzing data using SMMs
      predZ <- apply(Z, 2, function(i){predict(lm(i~1), type="response")})
      smm.results <- t(Y)%*%(Z - predZ)%*%t(Z-predZ)%*%A%*%solve(t(A)%*%(Z-predZ)%*%t(Z-predZ)%*%A)
      names(smm.results) <- c("psi1", "psi2", "psi3", "psi4")
      
      # Analyzing data using multivariable MR
      predA <- predict(lm(A~Z), type="response")
      multimr.results <- lm(Y~predA)$coef[-1]
      names(multimr.results) <- c("delta1", "delta2", "delta3", "delta4")
      
      # Returning results
      return(unlist(c(alpha_AA[aa,], smm.results, multimr.results)))
    }
    return(lapply(1:nrow(alpha_AA), f.alpha_AA) %>% do.call(rbind,.))
  }
  
  # Obtaining simulation results
  sim4.results <- lapply(1:1000, simulation4) %>% do.call(rbind,.)
  
  # Analyzing results
  sim4 <- sim4.results %>% 
    data.frame() %>%
    mutate(psi_sum=psi1+psi2+psi3+psi4) %>%
    mutate(delta_sum=delta1+delta2+delta3+delta4) %>%
    group_by(A1_A2, A1_A3, A2_A3) %>%
    summarise_all(mean) %>%
    mutate_at(vars(c(contains("psi"), contains("delta"))), ~sprintf("%.3f", .)) %>%
    relocate(psi_sum, .after=psi4)

  
# ----------------------- Appendix 3.5 -----------------------
# Identifying period effects using MR analysis of a single 
# measurement of the exposure
#  - Causal estimand of interest: period effect
#  - Number of exposure measurements considered in the model: one
#  - Instrument-exposure relationship changes over time? Yes
#  - Effect of exposure modified by previous exposure? No
#  - Presence of time-varying outcome-exposure confounding? No
# -------------------------------------------------------------

set.seed(123543)
n <- 10000
  
  # Specifying parameters
  alpha_U <- 0.5
  alpha_AA <- expand.grid(A1_A2=c(0.3, 0), A1_A3=c(0.3, 0), A2_A3=c(0.3, 0))
  alpha_ZA <- cbind(c(0.5, 0.3, 0.1),
                    c(0.1, 0.4, 0.2),
                    c(0.2, 0.5, 0.3),
                    c(0.3, 0.1, 0.4),
                    c(0.4, 0.2, 0.5))
  beta_AY <- c(-0.5, 0.1, 1)
  
  # Simulation function
  simulation5 <- function(i){
    
    # Looping over different alpha_AA
    f.alpha_AA <- function(aa){
      
      # Looping over exposure time point being analyzed 
      f.time <- function(t){
        
        # Generating variables
        Z <- lapply(1:5,function(x) {
          return(rbinom(n = n, size = 2, prob = 0.3))
        }) %>% do.call(cbind,.)
        
        U <- rnorm(n, 0, 1)
        A1 <- Z %*% alpha_ZA[1,] + alpha_U*U + rnorm(n, 0, 1)
        A2 <- Z %*% alpha_ZA[2,] + alpha_U*U + rnorm(n, 0, 1) + alpha_AA[aa,"A1_A2"]*A1
        A3 <- Z %*% alpha_ZA[3,] + alpha_U*U + rnorm(n, 0, 1) + alpha_AA[aa,"A1_A3"]*A1 + alpha_AA[aa,"A2_A3"]*A2
        
        Y <- A1*beta_AY[1] + A2*beta_AY[2] + A3*beta_AY[3] + U*alpha_U + rnorm(n, 0, 1)
        
        if (t==1){A <- A1}
        if (t==2){A <- A2}
        if (t==3){A <- A3}
        if (t==4){
          randA <- cut(runif(n), c(0, 1/3, 2/3, 1), labels=c(1,2,3))
          A <- ifelse(randA==1, A1, ifelse(randA==2, A2, A3))
        }
        
        # Analyzing data using SMMs
        predZ <- apply(Z, 2, function(i){predict(lm(i~1), type="response")})
        smm.results <- t(Y)%*%(Z - predZ)%*%t(Z-predZ)%*%A%*%solve(t(A)%*%(Z-predZ)%*%t(Z-predZ)%*%A)
        names(smm.results) <- c("psi")
        
        # Returning results
        return(unlist(c(alpha_AA[aa,], exposure=t, smm.results)))
      }
      return(lapply(1:4, f.time) %>% do.call(rbind,.))
    }
    return(lapply(1:nrow(alpha_AA), f.alpha_AA) %>% do.call(rbind,.))
  }
  
  # Obtaining simulation results
  sim5.results <- lapply(1:1000, simulation5) %>% do.call(rbind,.)
  
  # Analyzing results
  sim5 <- sim5.results %>% 
    data.frame() %>%
    group_by(A1_A2, A1_A3, A2_A3, exposure) %>%
    summarise_all(mean) %>%
    mutate_at(vars(c(contains("psi"), contains("delta"))), ~sprintf("%.3f", .)) %>%
    pivot_wider(id_cols=c(A1_A2, A1_A3, A2_A3), names_from=exposure, values_from=psi)


# ----------------------- Appendix 3.6 -----------------------
# Identifying period effects using MR analysis of a single 
# measurement of the exposure and assuming the instrument-exposure
# relationship stays constant
#  - Causal estimand of interest: period effect
#  - Number of exposure measurements considered in the model: one
#  - Instrument-exposure relationship changes over time? No
#  - Effect of exposure modified by previous exposure? No
#  - Presence of time-varying outcome-exposure confounding? No
# -------------------------------------------------------------  
  
set.seed(123543)
n <- 10000
  
  # Specifying parameters
  alpha_U <- 0.5
  alpha_AA <- expand.grid(A1_A2=c(0.3, 0),
                          A1_A3=c(0.3, 0),
                          A2_A3=c(0.3, 0))
  beta_AY <- c(-0.5, 0.1, 1)
  
  # Data Generation
  simulation6 <- function(i){
    
    # Looping over different values of alpha_AA
    f.alpha_AA <- function(aa){
      
      # Looping over different exposure time points considered in the analysis
      f.time <- function(t){
        
        alpha_ZA1 <- seq(0.1, 0.5, 0.1)
        alpha_ZA2 <- alpha_ZA1-alpha_ZA1*alpha_AA[aa, "A1_A2"]  
        alpha_ZA3 <- alpha_ZA1-alpha_ZA1*alpha_AA[aa, "A1_A2"]*alpha_AA[aa, "A2_A3"]-
          alpha_ZA1*alpha_AA[aa, "A1_A3"]-alpha_ZA2*alpha_AA[aa, "A2_A3"]
        
        alpha_ZA <- rbind(alpha_ZA1, alpha_ZA2, alpha_ZA3)  
        
        # Generating variables
        Z <- lapply(1:5,function(x) {
          return(rbinom(n = n, size = 2, prob = 0.3))
        }) %>% do.call(cbind,.)
        
        U <- rnorm(n, 0, 1)
        A1 <- Z %*% alpha_ZA[1,] + alpha_U*U + rnorm(n, 0, 1)
        A2 <- Z %*% alpha_ZA[2,] + alpha_U*U + rnorm(n, 0, 1) + alpha_AA[aa,"A1_A2"]*A1
        A3 <- Z %*% alpha_ZA[3,] + alpha_U*U + rnorm(n, 0, 1) + alpha_AA[aa,"A1_A3"]*A1 + alpha_AA[aa,"A2_A3"]*A2
        Y <- A1*beta_AY[1] + A2*beta_AY[2] + A3*beta_AY[3] + U*alpha_U + rnorm(n, 0, 1)
        
        if (t==1){A <- A1}
        if (t==2){A <- A2}
        if (t==3){A <- A3}
        
        # Analyzing data using SMMs
        predZ <- apply(Z, 2, function(i){predict(lm(i~1), type="response")})
        smm.results <- t(Y)%*%(Z - predZ)%*%t(Z-predZ)%*%A%*%solve(t(A)%*%(Z-predZ)%*%t(Z-predZ)%*%A)
        names(smm.results) <- c("psi")
        
        # Returning results
        return(unlist(c(alpha_AA[aa,], exposure=t, smm.results)))
      }
      return(lapply(1:3, f.time) %>% do.call(rbind,.))
    }
    return(lapply(1:nrow(alpha_AA), f.alpha_AA) %>% do.call(rbind,.))
  }
  
  # Obtaining simulation results
  sim6.results <- lapply(1:1000, simulation6) %>% do.call(rbind,.)
  
  # Analyzing results
  sim6 <- sim6.results %>% 
    data.frame() %>%
    group_by(A1_A2, A1_A3, A2_A3, exposure) %>%
    summarise_all(mean) %>%
    mutate_at(vars(c(contains("psi"), contains("delta"))), ~sprintf("%.3f", .)) %>%
    pivot_wider(id_cols=c(A1_A2, A1_A3, A2_A3), names_from=exposure, values_from=psi)

  
# ----------------------- Appendix 3.7 -----------------------
# Identifying period effects using MR analysis of multiple
# exposure measurements
#  - Causal estimand of interest: period effect
#  - Number of exposure measurements considered in the model: subset (two)
#  - Instrument-exposure relationship changes over time? Yes
#  - Effect of exposure modified by previous exposure? No
#  - Presence of time-varying outcome-exposure confounding? No
# -------------------------------------------------------------  

set.seed(123543)
n <- 10000
  
  # Simulation Parameters
  alpha_U <- 0.5
  alpha_AA <- expand.grid(A1_A2=c(0.3), A1_A3=c(0.3), A2_A3=c(0.3))
  alpha_ZA <- cbind(c(0.5, 0.3, 0.1),
                    c(0.1, 0.4, 0.2),
                    c(0.2, 0.5, 0.3),
                    c(0.3, 0.1, 0.4),
                    c(0.4, 0.2, 0.5))
  beta_AY <- c(-0.5, 0.1, 1)
  
  # Data Generation
  simulation7 <- function(i){
    
    # Looping over which exposure measurements are considered in the analysis
    f.A_included <- function(a_in){
      
      # Generating variables
      Z <- lapply(1:5,function(x) {
        return(rbinom(n = n, size = 2, prob = 0.3))
      }) %>% do.call(cbind,.)
      
      U <- rnorm(n, 0, 1)
      A1 <- Z %*% alpha_ZA[1,] + alpha_U*U + rnorm(n, 0, 1)
      A2 <- Z %*% alpha_ZA[2,] + alpha_U*U + rnorm(n, 0, 1) + alpha_AA[1,"A1_A2"]*A1
      A3 <- Z %*% alpha_ZA[3,] + alpha_U*U + rnorm(n, 0, 1) + alpha_AA[1,"A1_A3"]*A1 + alpha_AA[1,"A2_A3"]*A2
      
      if (a_in==12){A <- cbind(A1, A2)}
      if (a_in==13){A <- cbind(A1, A3)}
      if (a_in==23){A <- cbind(A2, A3)}  
      
      Y <- A1*beta_AY[1] + A2*beta_AY[2] + A3*beta_AY[3] + U*alpha_U + rnorm(n, 0, 1)
      
      # Analyzing data using SMMs
      predZ <- apply(Z, 2, function(i){predict(lm(i~1), type="response")})
      smm.results <- t(Y)%*%(Z - predZ)%*%t(Z-predZ)%*%A%*%solve(t(A)%*%(Z-predZ)%*%t(Z-predZ)%*%A)
      names(smm.results) <- c("psi1", "psi2")
      
      # Returning results
      return(unlist(c(a_in, smm.results)))
    }
    return(lapply(c(12, 13, 23), f.A_included) %>% do.call(rbind,.))
  }
  
  # Obtaining simulation results
  sim7.results <- lapply(1:1000, simulation7) %>% do.call(rbind,.)
  
  # Analyzing results
  sim7 <- sim7.results %>% 
    data.frame() %>%
    mutate(smm.psi1=ifelse(V1==12|V1==13, psi1, NA)) %>%
    mutate(smm.psi2=ifelse(V1==12, psi2, ifelse(V1==23, psi1, NA))) %>%
    mutate(smm.psi3=ifelse(V1==13|V1==23, psi2, NA)) %>%
    mutate(smm.sum=psi1+psi2) %>%
    group_by(V1) %>%
    summarise_all(mean) %>%
    mutate_at(vars(c(contains("smm"))), ~ifelse(!is.na(.), sprintf("%.3f", .), "N/A")) %>%
    select(-psi1, -psi2)

  
# ----------------------- Appendix 3.8 -----------------------
# Identifying period effects using MR analysis of multiple
# exposure measurements and assuming the instrument-exposure relationship
# stays constant
#  - Causal estimand of interest: period effect
#  - Number of exposure measurements considered in the model: subset (two)
#  - Instrument-exposure relationship changes over time? No
#  - Effect of exposure modified by previous exposure? No
#  - Presence of time-varying outcome-exposure confounding? No
# -------------------------------------------------------------  
  
set.seed(123543)
n <- 10000
  
  # Specifying parameters
  alpha_U <- 0.5
  alpha_AA <- expand.grid(A1_A2=c(0.3), A1_A3=c(0.3), A2_A3=c(0.3))
  beta_AY <- c(-0.5, 0.1, 1)
  
  alpha_ZA1 <- seq(0.1, 0.5, 0.1)
  alpha_ZA2 <- alpha_ZA1-alpha_ZA1*alpha_AA[aa, "A1_A2"]  
  alpha_ZA3 <- alpha_ZA1-alpha_ZA1*alpha_AA[aa, "A1_A2"]*alpha_AA[aa, "A2_A3"]-
    alpha_ZA1*alpha_AA[aa, "A1_A3"]-alpha_ZA2*alpha_AA[aa, "A2_A3"]
  
  alpha_ZA <- rbind(alpha_ZA1, alpha_ZA2, alpha_ZA3)  
  
  # Simulation function
  simulation8 <- function(i){
    
    # Looping over different combination of exposure measurements to be included
    # in the analysis
    f.A_included <- function(a_in){
      
      # Generating variables
      Z <- lapply(1:5,function(x) {
        return(rbinom(n = n, size = 2, prob = 0.3))
      }) %>% do.call(cbind,.)
      
      U <- rnorm(n, 0, 1)
      A1 <- Z %*% alpha_ZA[1,] + alpha_U*U + rnorm(n, 0, 1)
      A2 <- Z %*% alpha_ZA[2,] + alpha_U*U + rnorm(n, 0, 1) + alpha_AA[1,"A1_A2"]*A1
      A3 <- Z %*% alpha_ZA[3,] + alpha_U*U + rnorm(n, 0, 1) + alpha_AA[1,"A1_A3"]*A1 + alpha_AA[1,"A2_A3"]*A2
      
      if (a_in==12){A <- cbind(A1, A2)}
      if (a_in==13){A <- cbind(A1, A3)}
      if (a_in==23){A <- cbind(A2, A3)}  
      
      Y <- A1*beta_AY[1] + A2*beta_AY[2] + A3*beta_AY[3] + U*alpha_U + rnorm(n, 0, 1)
      
      # Analyzing data using SMMs
      predZ <- apply(Z, 2, function(i){predict(lm(i~1), type="response")})
      smm.results <- t(Y)%*%(Z - predZ)%*%t(Z-predZ)%*%A%*%solve(t(A)%*%(Z-predZ)%*%t(Z-predZ)%*%A)
      names(smm.results) <- c("psi1", "psi2")
      
      # Returning results
      return(unlist(c(a_in, smm.results)))
    }
    return(lapply(c(12, 13, 23), f.A_included) %>% do.call(rbind,.))
  }
  
  # Obtaining simulation results
  sim8.results <- lapply(1:1000, simulation8) %>% do.call(rbind,.)
  
  # Analyzing results
  sim8 <- sim8.results %>% 
    data.frame() %>%
    mutate(smm.psi1=ifelse(V1==12|V1==13, psi1, NA)) %>%
    mutate(smm.psi2=ifelse(V1==12, psi2, ifelse(V1==23, psi1, NA))) %>%
    mutate(smm.psi3=ifelse(V1==13|V1==23, psi2, NA)) %>%
    mutate(smm.sum=psi1+psi2) %>%
    group_by(V1) %>%
    summarise_all(mean) %>%
    mutate_at(vars(c(contains("smm"))), ~ifelse(!is.na(.), sprintf("%.3f", .), "N/A")) %>%
    select(-psi1, -psi2)


# ----------------------- Appendix 3.9 -----------------------
# Identifying period effects using MR analysis of multiple
# exposure measurements and assuming the instrument-exposure relationship
# stays constant during certain periods of time
#  - Causal estimand of interest: period effect
#  - Number of exposure measurements considered in the model: subset (two)
#  - Instrument-exposure relationship changes over time? Yes (but only
#    over certain, not all, intervals)
#  - Effect of exposure modified by previous exposure? No
#  - Presence of time-varying outcome-exposure confounding? No
# -------------------------------------------------------------  

set.seed(123543)
n <- 10000

  # Specifying parameters
  alpha_U <- 0.5
  alpha_AA <- expand.grid(A1_A2=c(0.3), A1_A3=c(0.3), A2_A3=c(0.3))
  beta_AY <- c(-0.5, 0.1, 1)
  
  # Simulation function
  simulation9 <- function(i){
    
    # Looping over different combination of exposure measurements to be included
    # in the analysis    
    f.A_included <- function(a_in){
      
      # Looping over different combination of time intervals where z-a 
      # relationship is constant
      f.z_same <- function(z_same){
        
        # Generating variables
        Z <- lapply(1:5,function(x) {
          return(rbinom(n = n, size = 2, prob = 0.3))
        }) %>% do.call(cbind,.)
        
        if (z_same==12){
          alpha_ZA1 <- seq(0.1, 0.5, 0.1)
          alpha_ZA2 <- seq(0.1, 0.5, 0.1)-alpha_ZA1*alpha_AA[aa, "A1_A2"]  
          alpha_ZA3 <- c(0.35, 0.25, 0.35, 0.25, 0.35)-alpha_ZA1*alpha_AA[aa, "A1_A2"]*alpha_AA[aa, "A2_A3"]-
            alpha_ZA1*alpha_AA[aa, "A1_A3"]-alpha_ZA2*alpha_AA[aa, "A2_A3"]
          
          alpha_ZA <- rbind(alpha_ZA1, alpha_ZA2, alpha_ZA3)  
        }
        if (z_same==13){
          alpha_ZA1 <- seq(0.1, 0.5, 0.1)
          alpha_ZA2 <- c(0.35, 0.25, 0.35, 0.25, 0.35)-alpha_ZA1*alpha_AA[aa, "A1_A2"]  
          alpha_ZA3 <- seq(0.1, 0.5, 0.1)-alpha_ZA1*alpha_AA[aa, "A1_A2"]*alpha_AA[aa, "A2_A3"]-
            alpha_ZA1*alpha_AA[aa, "A1_A3"]-alpha_ZA2*alpha_AA[aa, "A2_A3"]
          
          alpha_ZA <- rbind(alpha_ZA1, alpha_ZA2, alpha_ZA3)  
        }
        if (z_same==23){
          alpha_ZA1 <- c(0.35, 0.25, 0.35, 0.25, 0.35)
          alpha_ZA2 <- seq(0.1, 0.5, 0.1)-alpha_ZA1*alpha_AA[aa, "A1_A2"]  
          alpha_ZA3 <- seq(0.1, 0.5, 0.1)-alpha_ZA1*alpha_AA[aa, "A1_A2"]*alpha_AA[aa, "A2_A3"]-
            alpha_ZA1*alpha_AA[aa, "A1_A3"]-alpha_ZA2*alpha_AA[aa, "A2_A3"]
          
          alpha_ZA <- rbind(alpha_ZA1, alpha_ZA2, alpha_ZA3)  
        }
        
        U <- rnorm(n, 0, 1)
        A1 <- Z %*% alpha_ZA[1,] + alpha_U*U + rnorm(n, 0, 1)
        A2 <- Z %*% alpha_ZA[2,] + alpha_U*U + rnorm(n, 0, 1) + alpha_AA[1,"A1_A2"]*A1
        A3 <- Z %*% alpha_ZA[3,] + alpha_U*U + rnorm(n, 0, 1) + alpha_AA[1,"A1_A3"]*A1 + alpha_AA[1,"A2_A3"]*A2
        
        if (a_in==12){A <- cbind(A1, A2)}
        if (a_in==13){A <- cbind(A1, A3)}
        if (a_in==23){A <- cbind(A2, A3)}  
        
        Y <- A1*beta_AY[1] + A2*beta_AY[2] + A3*beta_AY[3] + U*alpha_U + rnorm(n, 0, 1)
        
        # Analyzing data using SMMs
        predZ <- apply(Z, 2, function(i){predict(lm(i~1), type="response")})
        smm.results <- t(Y)%*%(Z - predZ)%*%t(Z-predZ)%*%A%*%solve(t(A)%*%(Z-predZ)%*%t(Z-predZ)%*%A)
        names(smm.results) <- c("psi1", "psi2")
        
        # Returning results
        return(unlist(c(z_same=z_same, a_in=a_in, smm.results)))
      }
      return(lapply(c(12,13,23), f.z_same) %>% do.call(rbind,.))
    }
    return(lapply(c(12, 13, 23), f.A_included) %>% do.call(rbind,.))
  }
  
  # Obtaining simulation results
  sim9.results <- lapply(1:1000, simulation9) %>% do.call(rbind,.)
  
  # Analyzing results
  sim9 <- sim9.results %>% 
    data.frame() %>%
    mutate(smm.psi1=ifelse(a_in==12|a_in==13, psi1, NA)) %>%
    mutate(smm.psi2=ifelse(a_in==12, psi2, ifelse(a_in==23, psi1, NA))) %>%
    mutate(smm.psi3=ifelse(a_in==13|a_in==23, psi2, NA)) %>%
    mutate(smm.sum=psi1+psi2) %>%
    group_by(z_same, a_in) %>%
    summarise_all(mean) %>%
    mutate_at(vars(c(contains("smm"))), ~ifelse(!is.na(.), sprintf("%.3f", .), "N/A")) %>%
    select(-psi1, -psi2)