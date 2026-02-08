#' Calculate estimates of variability explained by model component
#'
#' @details Calculate variability explained over the entire domain using the
#' standard R2 technique, separately by model component
#'
#' @param Samples, aligned and formated posterior samples
#' @param B, orthonormal basis matrix
#' @param Data, data used to fit the model
#' @param Type, model type/structure
#'
#' @return named list of variability explained by component
#'
Var_Exp <- function(Estimates, Data, Type){
  X = Data$X
  var_tot = var(c(Data$Y))
  output = list()

  FE_mod = X %*% t(Estimates$FE)
  output$FE = 1 - var(c(FE_mod - Data$Y))/var_tot

  if(Type %in% c("Single-level", "Multilevel")){
    L1_dev = (Estimates$S1 %*% t(Estimates$EF1))[Data$ID, ]
    L1_mod = FE_mod + L1_dev
    output$L1 = 1 - var(c(L1_mod - Data$Y))/var_tot

    if(Type == "Multilevel"){
      L2_dev = Estimates$S2 %*% t(Estimates$EF2)
      L2_mod = L1_mod + L2_dev
      output$L2 = 1 - var(c(L2_mod - Data$Y))/var_tot
    }
  }

  return(output)
}

#' Calculate R-Hat statistics for fixed effects functions (median and maximum
#' over the domain)
#'
#' @details Calculates Gelman-Rubin R-Hat statistics at all points along
#' the functional domain, summarizing by the median and maximum observed values
#'
#' @param FE, fixed effect posterior samples organized by-chain
#' @param Data, data used to fit the model
#' @param VarNames, variable names for the fixed effects
#'
#' @return named vector of R-Hat statistic summaries
#'
RHat_FE <- function(FE, Data, VarNames){

  # Format samples
  mu_chains = map(FE, function(x){
    return(abind(x, along = 3))
  }) %>% abind(along = 4)

  # Calculate R-Hats
  mu_rhats = matrix(0, nrow = Data$M, ncol = Data$P)
  for(m in 1:Data$M){
    for(p in 1:Data$P){
      mu_rhats[m,p] = Rhat(mu_chains[m,p,,])
    }
  }

  # Label outputs for readability
  median_rhats = apply(mu_rhats, 2, median)
  names(median_rhats) = VarNames

  max_rhats = apply(mu_rhats, 2, max)
  names(max_rhats) = VarNames

  return(list(Median = median_rhats,Max = max_rhats))
}

#' Calculate R-Hat statistics for eigenfunctions (median and maximum
#' over the domain) and scores after alignment
#'
#' @details Calculates Gelman-Rubin R-Hat statistics for eigenfunctions at all
#' points along the functional domain, and for all scores. Summarizes by the
#' median and maximum observed values within each FPC
#'
#' @param EF, eigenfunction posterior samples organized by-chain
#' @param Scores, list of score samples
#' @param EF_Estimate, matrix FPC estimate
#' @param Data, data used to fit the model
#' @param Level, level of the hierarchy to consider
#' @param Routine, method of posterior FPC alignment, options are:
#' * **sign**: simple alignment based on cross-correlation signs
#' * **procrustes**: alignment based on Procrustes analysis
#'
#' @return vector of R-Hat statistic summaries
#'
RHat_FPC_Score <- function(EF, Scores, EF_Estimate, Data, Level, Routine){
  if(Level == 1){
    K = Data$K1
    Nv = Data$I
  }
  else{
    K = Data$K2
    Nv = Data$N
  }

  # Align
  Align_Objs = map(1:length(EF), function(x){
    if(Routine == "sign"){
      aligned = sign_al(EF[[x]], Scores[[x]], EF_Estimate)
    }
    else if(Routine == "procrustes"){
      aligned = procrustes_al(EF[[x]], Scores[[x]], EF_Estimate)
    }
    else{
      stop("Routine not supported")
    }
    return(aligned)
  })

  # Format
  FPC_chains = map(Align_Objs, function(x){
    return(abind(x$EF, along = 3))
  }) %>% abind(along = 4)

  Score_chains = map(Align_Objs, function(x){
    return(abind(x$Score, along = 3))
  }) %>% abind(along = 4)

  print(dim(FPC_chains))
  print(dim(Score_chains))

  # Calculate R-Hats
  FPC_rhats = matrix(0, nrow = Data$M, ncol = K)
  for(m in 1:Data$M){
    for(k in 1:K){
      FPC_rhats[m,k] = Rhat(FPC_chains[m,k,,])
    }
  }

  Score_rhats = matrix(0, nrow = Nv, ncol = K)
  for(n in 1:Nv){
    for(k in 1:K){
      Score_rhats[n,k] = Rhat(Score_chains[n,k,,])
    }
  }

  # Label outputs for readability
  all_names = paste0("FPC ", 1:K)

  median_FPC_rhats = apply(FPC_rhats, 2, median)
  names(median_FPC_rhats) = all_names
  max_FPC_rhats = apply(FPC_rhats, 2, max)
  names(max_FPC_rhats) = all_names

  median_Score_rhats = apply(Score_rhats, 2, median)
  names(median_Score_rhats) = all_names
  max_Score_rhats = apply(Score_rhats, 2, max)
  names(max_Score_rhats) = all_names

  return(list(EF = list(Median = median_FPC_rhats, Max = max_FPC_rhats),
              Score = list(Median = median_Score_rhats, Max = max_Score_rhats)))
}

#' Calculate R-Hat statistics for parameters within the MCMC model (fixed
#' effects and any present FPCs/scores)
#'
#' @details Calculates Gelman-Rubin R-Hat statistics for functional
#' components at all points along the functional domain, and every score.
#' Summarizes by the median and maximum observed values within function
#'
#' @param Samples, formatted posterior samples organized by chain
#' @param Estimates, object containing posterior estimates
#' @param Data, data used to fit the model
#' @param VarNames, variable names for the fixed effects
#' @param Type, model type/structure
#' @param Routine, method of posterior FPC alignment, options are:
#' * **sign**: simple alignment based on cross-correlation signs
#' * **procrustes**: alignment based on Procrustes analysis
#'
#' @return vector of R-Hat statistic summaries
#'
RHat_FullMod <- function(Samples, Estimates, Data, VarNames, Type,
                         Routine = "sign"){

  output = list()
  output$Beta = RHat_FE(Samples$FE, Data, VarNames)

  if(Type %in% c("Single-level", "Multilevel")){
    level1 = RHat_FPC_Score(Samples$EF1, Samples$S1, Estimates$EF1,
                            Data, 1, Routine)
    output$EF1 = level1$EF
    output$S1 = level1$Score

    if(Type == "Multilevel"){
      level2 = RHat_FPC_Score(Samples$EF2, Samples$S2, Estimates$EF2,
                              Data, 2, Routine)
      output$EF2 = level2$EF
      output$S2 = level2$Score
    }
  }

  return(output)
}
