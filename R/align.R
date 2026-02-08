#' Align posterior FPC/score samples via sign flips
#'
#' @details Changes signs of FPCs/scores to align with a fixed reference point
#' (according to cross-correlation between FPCs)
#'
#' @param FPCs, list of eigenfunction samples
#' @param Scores, list of score samples
#' @param Reference, FPC matrix (# observations x # FPCs)
#' @param Weights, list of eigenfunction spline weight samples, not required
#'
#' @return Named list with aligned FPC/Score samples
#'
sign_al <- function(FPCs, Scores, Reference, Weights = NULL){
  # Constants
  K = FPCs[[1]] %>% ncol()

  # Outputs
  out_FPC = list()
  out_Score = list()
  out_Weight = list()

  # Order/signs
  align_objs = map(1:length(FPCs), function(x){
    evals = apply(Scores[[x]], 2, var)
    fpc_order = sort(evals, decreasing = T, index.return = T)$ix

    Phi = FPCs[[x]][, fpc_order]
    dot_prods = colSums(Phi * Reference)
    fpc_signs = sign(dot_prods)
    fpc_signs[fpc_signs == 0] = 1 # Ensure no zeros

    return(list(order = fpc_order, signs = fpc_signs))
  })

  for(i in 1:length(align_objs)){
    ord = align_objs[[i]]$order
    sig = align_objs[[i]]$signs

    out_FPC[[i]] = FPCs[[i]][, ord] %*% diag(sig)
    out_Score[[i]] = Scores[[i]][, ord] %*% diag(sig)
  }

  if(!is.null(Weights)){
    for(i in 1:length(align_objs)){
      ord = align_objs[[i]]$order
      sig = align_objs[[i]]$signs

      out_Weight[[i]] = Weights[[i]][, ord] %*% diag(sig)
    }
  }
  else{
    out_Weight = NA
  }

  output = list(EF = out_FPC, Weight = out_Weight, Score = out_Score)

  return(output)
}

#' Align posterior FPC/score samples via Procrustes analysis
#'
#' @details Changes signs of FPCs/scores to align with a fixed reference point
#' (using Procrustes analysis to minimize the difference to the reference)
#'
#' @param FPCs, list of eigenfunction samples
#' @param Scores, list of score samples
#' @param Reference, FPC matrix (# observations x # FPCs)
#' @param Weights, list of eigenfunction spline weight samples, not required
#'
#' @return Named list with aligned FPC/Score samples
#'
procrustes_al <- function(FPCs, Scores, Reference, Weights = NULL){
  # Outputs
  out_FPC = list()
  out_Score = list()
  out_Weight = list()

  transforms = map(FPCs, function(Phi){
    R_comps = svd(t(Phi) %*% Reference)
    R = R_comps$u %*% t(R_comps$v)
    return(R)
  })

  for(i in 1:length(transforms)){
    out_FPC[[i]] = FPCs[[i]] %*% transforms[[i]]
    out_Score[[i]] = Scores[[i]] %*% transforms[[i]]
  }

  if(!is.null(Weights)){
    for(i in 1:length(transforms)){
      out_Weight[[i]] = Weights[[i]] %*% transforms[[i]]
    }
  }
  else{
    out_Weight = NA
  }

  output = list(EF = out_FPC, Score = out_Score, Weight = out_Weight)

  return(output)
}

#' Align posterior FPC/score samples via Procrustes analysis
#'
#' @details Changes signs of FPCs/scores to align with a fixed reference point
#' (using Procrustes analysis to minimize the difference to the reference)
#'
#' @param Samples, formated poseterior samples
#' @param Data, data used to fit the model
#' @param B, orthonormal basis matrix
#' @param Reference, named list of FPC matrices against which to align. Later
#' entries correspond to lower levels of the hierarchy
#' @param Type, model type/structure
#' @param Routine, method of posterior FPC alignment, options are:
#' * **sign**: simple alignment based on cross-correlation signs
#' * **procrustes**: alignment based on Procrustes analysis
#'
#' @return BFM_Samples object with extracted samples after alignment
#'
Align_Posterior <- function(Samples, Data, B, Type, Routine = "sign"){
  if(Type == "Fixed Effects"){
    stop("No samples to align")
  }

  # Set reference to central FPC sample
  Reference = list()
  Reference[[1]] = FPC_Est(Samples$W1, Samples$S1, B, Data$K1)
  if(Type == "Multilevel"){
    Reference[[2]] = FPC_Est(Samples$W2, Samples$S2, B, Data$K2)
  }

  # Constants
  output = list(FE = Samples$FE,
                wFE = Samples$wFE)

  if(Routine == "sign"){
    level1_align = sign_al(Samples$EF1, Samples$S1, Reference[[1]], Samples$W1)
    output$EF1 = level1_align$EF
    output$W1 = level1_align$Weight
    output$S1 = level1_align$Score
    output$EV1 = map(output$S1, function(s){
      EV = apply(s, 2, var)
      dim(EV) = c(length(EV), 1)
      return(EV)
    })

    if(Type == "Multilevel"){
      level2_align = sign_al(Samples$EF2, Samples$S2, Reference[[2]], Samples$W2)
      output$EF2 = level2_align$EF
      output$W2 = level2_align$Weight
      output$S2 = level2_align$Score
      output$EV2 = map(output$S2, function(s){
        EV = apply(s, 2, var)
        dim(EV) = c(length(EV), 1)
        return(EV)
      })
    }
  }
  else if(Routine == "procrustes"){
    level1_align = procrustes_al(Samples$EF1, Samples$S1, Reference[[1]], Samples$W1)
    output$EF1 = level1_align$EF
    output$W1 = level1_align$Weight
    output$S1 = level1_align$Score
    output$EV1 = map(output$S1, function(s){
      return(apply(s, 2, var))
    })

    if(Type == "Multilevel"){
      level2_align = procrustes_al(Samples$EF2, Samples$S2, Reference[[2]], Samples$W2)
      output$EF2 = level2_align$EF
      output$W2 = level2_align$Weight
      output$S2 = level2_align$Score
      output$EV2 = map(output$S2, function(s){
        return(apply(s, 2, var))
      })
    }
  }
  else{
    stop("Routine not supported")
  }
  return(output)
}
