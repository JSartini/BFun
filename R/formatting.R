#' Formats input data and constants for Bayesian Functional (Mixed) Models
#'
#' @details Places all inputs in the necessary format and record model
#' structure
#'
#' @param form, a formula indicating the functional response (in matrix form)
#' and the predictor covariates (scalar)
#' @param data, dataframe containing the data to be fit. Functional data
#' should be in matrix format.
#' @param id, column name containing subject-level identifiers (if available)
#' @param K1, number of FPCs used at the first level, default chosen empirically
#' to account for 95% variability
#' @param visit, column name containing visit-level identifiers (if available)
#' @param K2, number of FPCs used at the second level, default chosen empirically
#' to account for 95% variability
#' @param spline_basis, type of spline basis used. Several are supported:
#' * **B**: Simple degree 3 B-splines, the default
#' * **Splinets**: A localized, orthonormal basis
#' * **OB**: Orthonormalized B-splines
#' * **Vec**: Vector orthonormal basis created through SVD
#' @param spline_dim, dimension of the spline basis, defaults to 25
#' @param args, argument points at which the function is evaluated, defaults to
#' equally-spaced points between 0 and 1
#' @param out_args, output argument points at which functions are evaluated, defaults
#' to args
#'
#' @return list of formatted input data and model type labels
#'
model_input <- function(form, data, id, visit, K1, K2, spline_basis,
                        spline_dim, args, out_args){

  # Extract data
  form_list = as.character(form)
  Y_mat = as.matrix(data[,form_list[2]])
  X_mat = model.matrix(formula(paste("~", form_list[3], collapse = " ")), data = data)
  x_names = colnames(X_mat)

  # Handle input/output arguments
  if(is.null(args)){
    scaled_args = seq(0, 1, length.out = ncol(Y_mat))
  }
  else{
    scaled_args = (args - min(args))/(max(args) - min(args))
  }
  if(is.null(out_args)){
    scaled_out = scaled_args
  }
  else{
    scaled_out = (out_args - min(out_args))/(max(out_args) - min(out_args))
  }

  # Set up inputs
  spline_objs = calc_splines(spline_basis, spline_dim, scaled_args)
  B_out = calc_splines(spline_basis, spline_dim, scaled_out)$Basis
  if(is.null(id)){
    model_file = "FAST_FE"
    mod_type = "Fixed Effects"
    input_data = list(N = nrow(Y_mat), M = ncol(wider_Y), P = ncol(X_mat),
                      Q = spline_dim, Y = Y_mat, X = X_mat,
                      B = spline_objs$Basis, P_alpha = spline_objs$Penalty)
  }
  else{
    if(spline_basis == "B"){
      stop("Basis not supported when modeling random effects")
    }

    id_vec = data %>%
      group_by(across(all_of(id))) %>%
      mutate(new_id = cur_group_id()) %>%
      ungroup() %>%
      pull(new_id)

    input_data = list(N = nrow(Y_mat), I = max(id_vec), ID = id_vec,
                      M = ncol(Y_mat), P = ncol(X_mat), Q = spline_dim, Y = Y_mat,
                      X = X_mat, B = spline_objs$Basis, P_alpha = spline_objs$Penalty)

    FE_mod = fosr2s(Y_mat, X_mat, nbasis = spline_dim)
    FE_Fitted = X_mat %*% t(FE_mod$est.func)

    if(is.null(visit)){
      model_file = "FAST_SL"
      mod_type = "Single-level"

      if(is.null(K1)){
        freq_fpca = fpca.face(Y = Y_mat - FE_Fitted, knots = spline_dim, pve = 0.95)
        input_data$K1 = freq_fpca$npc
      }
      else{
        input_data$K1 = K1
      }
    }
    else{
      model_file = "FAST_ML"
      mod_type = "Multilevel"

      if(is.null(K1)){
        freq_fpca = mfpca.face(Y = Y_mat - FE_Fitted, id = data[, id],
                               knots = spline_dim,
                               pve = 0.95, pve2 = 0.95)
        input_data$K1 = freq_fpca$npc$level1
        input_data$K2 = freq_fpca$npc$level2
      }
      else{
        input_data$K1 = K1
        input_data$K2 = K2
      }
    }
  }

  return(list(input_data = input_data, file = model_file, type = mod_type,
              new_B = B_out, x_names = x_names))
}


#' Formats model output of Bayesian Functional (Mixed) Model fits to be
#' provided to the user
#'
#' @details Places all outputs in readable, accessible formats
#'
#' @param Estimates, list of posterior mean estimates
#' @param Inferences, list of posterior credible intervals
#' @param Type, model type/structure
#' @param VarNames, ordered names of predictor covariates
#'
#' @return list of formatted input data and model type labels
#'
model_output <- function(Estimates, Inferences, Type, VarNames){
  output = list(Predictors = VarNames)

  output$BetaHat = Estimates$FE
  output$YHat = Estimates$YHat

  output$Beta_CI = Inferences$FE
  output$YHat_CI = Inferences$YHat

  if(Type == "Single-Level"){
    output$Scores = Estimates$S1
    output$FPC = Estimates$EF1
    output$EVal = Estimates$EV1

    output$Scores_CI = Inferences$S1
    output$FPC_CI = Inferences$EF1
    output$EVal_CI = Inferences$EV1
  }
  else if(Type == "Multilevel"){
    output$Scores = list(Level1 = Estimates$S1,
                         Level2 = Estimates$S2)
    output$FPC = list(Level1 = Estimates$EF1,
                      Level2 = Estimates$EF2)
    output$EVal = list(Level1 = Estimates$EV1,
                       Level2 = Estimates$EV2)

    output$Scores_CI = list(Level1 = Inferences$S1,
                            Level2 = Inferences$S2)
    output$FPC_CI = list(Level1 = Inferences$EF1,
                      Level2 = Inferences$EF2)
    output$EVal_CI = list(Level1 = Inferences$EV1,
                          Level2 = Inferences$EV2)
  }

  return(output)
}
