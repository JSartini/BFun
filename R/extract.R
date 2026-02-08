#' Extract permuted samples from output of cmdstanr
#'
#' @details Places posterior samples from cmdstanr into standard rstan format
#'
#' @param fit_obj, result of calling $sample()
#'
#' @return Named list of formated posterior samples
#'
cmd_extract_permuted <- function(fit_obj){
  vars = fit_obj$metadata()$stan_variables
  draws = posterior::as_draws_rvars(fit_obj$draws())

  lapply(vars, \(var_name){
    posterior::draws_of(draws[[var_name]], with_chains = FALSE)
  }) |> setNames(vars)
}

#' Extract by-chain samples from output of cmdstanr
#'
#' @details Places posterior samples from cmdstanr into standard rstan format
#'
#' @param fit_obj, result of calling $sample()
#'
#' @return Named array of posterior samples by chain
#'
cmd_extract_chains <- function(fit_obj){
  return(fit_obj$draws(format = "draws_array"))
}

#' Place posterior samples in usable format, interleaving the chains
#'
#' @details Places posterior samples into named lists of appropriate
#' dimension. Samples from multiple chains are interleaved in
#' random order.
#'
#' @param Samples, samples from call to STAN
#' @param Data, data used to fit the STAN model
#' @param B_out, orthonormal basis matrix
#' @param Type, model type/structure
#'
#' @return list containing the formatted samples
#'
extract_permuted <- function(Samples, Data, B_out, Type){
  # Constants
  n_samp = length(Samples$sigma2)

  # Fixed effects
  wBeta = map(1:n_samp, function(x){
    wB_s = Samples$w_mu[x,,]
    if(is.null(dim(wB_s))){
      dim(wB_s) = c(length(wB_s), 1)
    }
    return(wB_s)
  })
  Beta = map(wBeta, function(x){
    FE = B_out %*% x
    return(FE)
  })

  # FPCA structure
  if(Type %in% c("Single-level", "Multilevel")){
    S1 = map(1:n_samp, function(x){
      return(Samples$Scores1[x,,])
    })
    W1 = map(1:n_samp, function(x){
      return(Samples$Psi_1[x,,])
    })
    EF1 = map(W1, function(psi){
      EF_sample = B_out %*% psi
      return(EF_sample)
    })

    if(Type == "Multilevel"){
      S2 = map(1:n_samp, function(x){
        return(Samples$Scores2[x,,])
      })
      W2 = map(1:n_samp, function(x){
        return(Samples$Psi_2[x,,])
      })
      EF2 = map(W2, function(psi){
        EF_sample = B_out %*% psi
        return(EF_sample)
      })
    }
    else{
      S2 = NULL; W2 = NULL; EF2 = NULL
    }
  }
  else{
    S1 = NULL; W1 = NULL; EF1 = NULL
    S2 = NULL; W2 = NULL; EF2 = NULL
  }

  return(list(wFE = wBeta, FE = Beta,
              S1 = S1, W1 = W1, EF1 = EF1,
              S2 = S2, W2 = W2, EF2 = EF2))
}

#' Place posterior samples in usable format, separating the chains
#'
#' @details Places posterior samples into named lists of appropriate
#' dimension. Samples from multiple chains are each given their own list entry.
#'
#' @param Samples, samples from call to STAN, organized by-chain
#' @param Data, data used to fit the STAN model
#' @param B_out, orthonormal basis matrix
#' @param Type, model type/structure
#'
#' @return List of posterior draws organized by chain
#'
extract_chain <- function(Samples, Data, B_out, Type){
  # Constants
  dimen_names = dimnames(Samples)$variable
  n_samp = dim(Samples)[1]
  n_chain = dim(Samples)[2]

  # Fixed effects
  FE_idx = grepl("w_mu", dimen_names)
  FE = list()

  for(j in 1:n_chain){
    FE[[j]] = map(1:n_samp, function(i){
      FE_mat = Samples[i,j,FE_idx]
      dim(FE_mat) = c(Data$Q, Data$P)
      FE_sample = B_out %*% FE_mat
      return(FE_sample)
    })
  }

  # FPCA structure
  if(Type %in% c("Single-level", "Multilevel")){
    S1_idx = grepl("Scores1", dimen_names)
    W1_idx = grepl("Psi_1", dimen_names)
    S1 = list()
    W1 = list()
    EF1 = list()

    for(j in 1:n_chain){
      S1[[j]] = map(1:n_samp, function(i){
        Xi = Samples[i,j,S1_idx]
        dim(Xi) = c(Data$I, Data$K1)
        return(Xi)
      })
      W1[[j]] = map(1:n_samp, function(i){
        Psi = Samples[i,j,W1_idx]
        dim(Psi) = c(Data$Q, Data$K1)
        return(Psi)
      })
      EF1[[j]] = map(W1[[j]], function(psi){
        EF_sample = B_out %*% psi
        return(EF_sample)
      })
    }

    if(Type == "Multilevel"){
      S2_idx = grepl("Scores2", dimen_names)
      W2_idx = grepl("Psi_2", dimen_names)
      S2 = list()
      W2 = list()
      EF2 = list()

      for(j in 1:n_chain){
        S2[[j]] = map(1:n_samp, function(i){
          Xi = Samples[i,j,S2_idx]
          dim(Xi) = c(Data$N, Data$K2)
          return(Xi)
        })
        W2[[j]] = map(1:n_samp, function(i){
          Psi = Samples[i,j,W2_idx]
          dim(Psi) = c(Data$Q, Data$K2)
          return(Psi)
        })
        EF2[[j]] = map(W2[[j]], function(psi){
          EF_sample = B_out %*% psi
          return(EF_sample)
        })
      }
    }
    else{
      S2 = NULL; W2 = NULL; EF2 = NULL
    }
  }
  else{
    S1 = NULL; W1 = NULL; EF1 = NULL
    S2 = NULL; W2 = NULL; EF2 = NULL
  }

  return(list(S1 = S1, W1 = W1, EF1 = EF1,
              S2 = S2, W2 = W2, EF2 = EF2,
              FE = FE))
}
