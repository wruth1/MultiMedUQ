
# The following packages are used here, but nowhere else: MASS, boot, purrr


#### Simulate a sample dataset ####


# make_validation_data(50, 20, b_Y, theta_Y, b_M, theta_M, output_list = F)


# Utilities ####

#' Store all regression parameters in a single list
#'
#' @param beta_Y Vector of length 4 containing fixed-effects for the outcome model. Order of variables is: Intercept, M, X, C1, C2.
#' @param Gamma_Y Covariance matrix of size 3x3 containing random effects for the outcome model. Order of variables is: Intercept, M, X.
#' @param beta_M Vector of length 3 containing fixed-effects for the outcome model. Order is: Intercept, X, C1, C2.
#' @param Gamma_M Covariance matrix of size 2x2 containing random effects for the outcome model. Order of variables is: Intercept, X.
#'
#' @return A named list containing fixed and random effects parameters for both the outcome and mediator regression models. Names in list match names of arguments to this function.
#' @export
#'
#' @examples
#' make_all_reg_pars()
make_all_reg_pars <- function(beta_Y = c(1, -1.5, 1, -1.5, 1), Gamma_Y = diag(3), beta_M = c(-1, 1, -1, 1), Gamma_M = diag(2)){
  names(beta_Y) = c("Int", "M", "X")
  names(beta_M) = c("Int", "X")
  list(beta_Y = beta_Y, Gamma_Y = Gamma_Y, beta_M = beta_M, Gamma_M = Gamma_M)
}


#' Construct a single data frame by `rbind`-ing elements of a list. Add/Create group labels.
#'
#' @param X_list A list of data frames. Must all have the same number of columns.
#' @param group_labels Optionally, a vector of labels for each element of `X_list`. Must be either the same length as `X_list` or equal to the total number of rows among elements in `X_list`. If `NULL`, labels are `G1`, `G2`,...
#'
#' @return A data frame containing all elements of `X_list` with an extra column for labels.
#' @export
#'
#' @examples
#' data = as.data.frame(matrix(c(1,0,0,1), nrow = 2))
#' more_data = as.data.frame(matrix(c(0,1,1,0), nrow = 2))
#' data_list = list(data, more_data)
#'
#' list_2_data(data_list)
#'
#' data_names = c("Alice", "Bob")
#' list_2_data(data_list, data_names)
list_2_data <- function(X_list, group_labels = NULL){
  # Stack elements of X_list ----
  X_data = purrr::list_rbind(X_list)

  # Create group labels ----
  all_labels = list_2_data_make_labels(X_list, group_labels)

  X_data$group = all_labels

  return(X_data)
}






#' Create group labels
#'
#' @param X_list A list of data frames.
#' @param group_labels Optionally, a vector of labels for each element of `X_list`. Must be either the same length as `X_list` or equal to the total number of rows among elements in `X_list`. If `NULL`, labels are `G1`, `G2`,...
#'
#' @return A vector of group labels with length equal to the total number of rows among elements of `X_list`.
#' @export
#' @keywords internal
#'
#' @examples
#' data = as.data.frame(matrix(c(1,0,0,1), nrow = 2))
#' more_data = as.data.frame(matrix(c(0,1,1,0), nrow = 2))
#' data_list = list(data, more_data)
#'
#' list_2_data_make_labels(data_list)
#'
#' data_names = c("Alice", "Bob")
#' list_2_data_make_labels(data_list, data_names)
list_2_data_make_labels <- function(X_list, group_labels = NULL){
  group_sizes = purrr::map_int(X_list, nrow)

  if(is.null(group_labels)){                              # Nothing supplied. Construct names and labels.
    K = length(X_list)
    group_names = paste0("G", 1:K)
    all_labels = rep(group_names, times = group_sizes)
  } else if(length(group_labels) == length(X_list)){      # Group names supplied. Construct labels.
    all_labels = rep(group_labels, times = group_sizes)
  } else if(length(group_labels) == sum(group_sizes)){        # Labels supplied.
    all_labels = group_labels
  } else{                                                 # Non-conformable arguments.
    stop("In list_2_data: Non-conformable sizes between X_list and group_labels.")
  }

  return(all_labels)
}






#' Generate a dataset for validating methodology
#'
#' @description
#' Generate a validation dataset with `n` observations in each of `K` groups. Variables include one outcome (Y), one mediator (M), one exposure (X), and two confounders (C1 and C2). All variables are binary. Default values are available for regression parameters (fixed effects and covariances of random effects).
#'
#' Order of variables is `Y`, `M`, `X`, `C1`, `C2`
#'
#' @param n Number of observations in each group.
#' @param K Number of groups
#' @param b_Y,b_M Coefficient vectors for \eqn{Y}-model and \eqn{M}-model, respectively.
#' @param theta_Y,theta_M Covariance parameters of random effects in \eqn{Y}-model and \eqn{M}-model, respectively. See details.
#' @param output_list Should output be formatted as a list with one component per group (of size n-by-5) or a single data.frame of size (Kn)-by-6 with a column labelled `group`?
#' @param return_REs Should random effects also be returned? If TRUE, output is a list containing two components: `data` and `REs`.
#' @param which_REs Which random effects to include in the calculation. Default is all. See the \href{../vignettes/which_REs.Rmd}{vignette} for more details.
#'
#' @return A simulated dataset, optionally inside a list which also contains random effects.
#' @export
#'
#' @examples
#' n = 20
#' K = 3
#' all_reg_pars = make_all_reg_pars()
#'
#' # Format output as a list
#' make_validation_data(n, K, all_reg_pars, output_list = TRUE)
#'
#' # Format output as a data.frame
#' make_validation_data(n, K, all_reg_pars, output_list = FALSE)
make_validation_data <-
  function(n, K, b_Y, theta_Y, b_M, theta_M,
           output_list = FALSE,
           return_REs = FALSE,
           which_REs = c("Y.Int", "Y.X", "Y.M", "M.Int", "M.X")) {
    # all_reg_pars = make_all_reg_pars()

    # Check that the number of variance parameters supplied matches which_REs
    if(num_Y_REs(which_REs) != theta2num_REs(theta_Y)){
      stop("In make_validation_data: The number of variance parameters supplied for Y does not match which_REs.")
    }
    if(num_M_REs(which_REs) != theta2num_REs(theta_M)){
      stop("In make_validation_data: The number of variance parameters supplied for M does not match which_REs.")
    }


    all_reg_pars = list(beta_M = b_M, Gamma_M = theta2Sigma(theta_M),
                        beta_Y = b_Y, Gamma_Y = theta2Sigma(theta_Y))


    # Generate data as a list
    data_list_output = make_validation_data_list(n, K, all_reg_pars, return_REs, which_REs)
    if (!return_REs) {
      data_list = data_list_output
    } else{
      data_list = data_list_output[["data"]]
      all_REs = data_list_output[["all_REs"]]
    }


    # Optionally, stack groups into a single data.frame
    if (!output_list) {
      data_output = list_2_data(data_list)
    } else{
      data_output = data_list
    }

    # Return data, optionally alongside random effects
    if (!return_REs) {
      return(data_output)
    } else{
      return(list(data = data_output, all_REs = all_REs))
    }
  }


make_validation_data_list <-
  function(n, K, all_reg_pars, return_REs = FALSE, which_REs = c("Y.Int", "Y.X", "Y.M", "M.Int", "M.X")) {
    output_list = purrr::map(1:K,
                             ~ make_one_group_validation(n, all_reg_pars, return_REs, which_REs))

    if (!return_REs) {
      return(output_list)
    } else{
      data_list = purrr::map(output_list, "data")
      RE_list = purrr::map(output_list, "REs")
      return(list(data = data_list, all_REs = RE_list))
    }
  }


make_one_group_validation <-
  function(n, all_reg_pars, return_REs = FALSE, which_REs = c("Y.Int", "Y.X", "Y.M", "M.Int", "M.X")) {
    X = make_X_validation(n)
    all_Cs = make_C_validation(n)

    M_info = make_M_validation(X, all_Cs, all_reg_pars, return_REs, which_REs)

    if (!return_REs) {
      M = M_info
    } else{
      M = M_info[["M"]]
      REs_M = M_info[["REs"]]
    }



    Y_info = make_Y_validation(M, X, all_Cs, all_reg_pars, return_REs, which_REs)

    if (!return_REs) {
      Y = Y_info
    } else{
      Y = Y_info[["Y"]]
      REs_Y = Y_info[["REs"]]
    }

    output_data = data.frame(Y = Y, X = X, M = M, all_Cs)

    if (!return_REs) {
      return(output_data)
    } else{
      return(list(
        data = output_data,
        REs = list(M = REs_M, Y = REs_Y, data = output_data)
      ))
    }
  }

make_X_validation <- function(n) {
  return(stats::rbinom(n, 1, 0.5))
}

make_C_validation <- function(n) {
  C1 = stats::rbinom(n, 1, 0.5)
  C2 = stats::rbinom(n, 1, 0.5)

  return(data.frame(C1 = C1, C2 = C2))
}



make_M_validation <-
  function(X, all_Cs, all_reg_pars, return_REs = FALSE, which_REs = c("Y.Int", "Y.X", "Y.M", "M.Int", "M.X")) {
    beta_M = all_reg_pars$beta_M    # Coefficient vector for fixed effects
    Gamma_M = all_reg_pars$Gamma_M  # Covariance matrix of random effects

    n = length(X)
    if (nrow(all_Cs) != n)
      stop("In make_M_validation: X and all_Cs must have same number of rows.")

    # Organize X and all_Cs into datasets
    ## Fixed effects data
    data_fix = data.frame(X = X, all_Cs)

    ## Random effects data
    RE_names = expand_REs(which_REs)
    data_ran_raw = c()
    if("M.Int" %in% RE_names){
      data_ran_raw = cbind(data_ran_raw, rep(1, times=n))
    }
    if("M.X" %in% RE_names){
      data_ran_raw = cbind(data_ran_raw, X)
    }
    data_ran = data.frame(data_ran_raw)


    # Compute linear predictors
    lin_pred_info = get_lin_preds(
      data_fix,
      data_ran,
      beta_M,
      Gamma_M,
      add_intercept_fix = TRUE,
      add_intercept_ran = FALSE,
      return_REs = return_REs
    )
    if (!return_REs) {
      lin_preds = lin_pred_info
    } else{
      lin_preds = lin_pred_info[["lin_preds"]]
      REs = lin_pred_info[["REs"]]
    }

    all_probs = boot::inv.logit(lin_preds)

    M = stats::rbinom(n, 1, all_probs)

    if (!return_REs) {
      return(M)
    } else{
      return(list(M = M, REs = REs))
    }
  }


# hist(lin_preds)
# plot(lin_preds)

# hist(all_probs)
# plot(all_probs)
# plot(X, all_probs)
# plot(all_Cs$C1, all_probs)
# plot(all_Cs$C2, all_probs)

make_Y_validation <-
  function(M, X, all_Cs, all_reg_pars, return_REs = FALSE, which_REs = c("Y.Int", "Y.X", "Y.M", "M.Int", "M.X")) {
    beta_Y = all_reg_pars$beta_Y    # Coefficient vector for fixed effects
    Gamma_Y = all_reg_pars$Gamma_Y  # Covariance matrix of random effects

    n = length(M)
    if ((length(X) != n) ||
        (nrow(all_Cs) != n))
      stop("In make_Y_validation: M, X and all_Cs must have same number of rows.")

    # Organize M, X and all_Cs into datasets
    ## Fixed effects
    data_fix = data.frame(X = X, M = M, all_Cs)

    ## Random effects
    RE_names = expand_REs(which_REs)
    data_ran_raw = c()
    if("Y.Int" %in% RE_names){
      data_ran_raw = cbind(data_ran_raw, rep(1, times=n))
    }
    if("Y.X" %in% RE_names){
      data_ran_raw = cbind(data_ran_raw, X)
    }
    if("Y.M" %in% RE_names){
      data_ran_raw = cbind(data_ran_raw, M)
    }
    data_ran = data.frame(data_ran_raw)


    ## Compute linear predictors
    lin_pred_info = get_lin_preds(
      data_fix,
      data_ran,
      beta_Y,
      Gamma_Y,
      add_intercept_fix = TRUE,
      add_intercept_ran = FALSE,
      return_REs = return_REs
    )
    if (!return_REs) {
      lin_preds = lin_pred_info
    } else{
      lin_preds = lin_pred_info[["lin_preds"]]
      REs = lin_pred_info[["REs"]]
    }

    all_probs = boot::inv.logit(lin_preds)

    Y = stats::rbinom(n, 1, all_probs)

    if (!return_REs) {
      return(Y)
    } else{
      return(list(Y = Y, REs = REs))
    }
  }




#### Compute linear predictors ####


#' Compute the contribution to the linear predictor due to the provided dataset and coefficient vector
#'
#' @param data A dataset. Can be a data frame or a matrix (former is converted to latter internally)
#' @param beta A vector of coefficients. Length must either match number of columns of `data` or be 1 less, if `add_intercept` is `FALSE` or `TRUE` respectively. In the latter case, intercept should be the first element of beta.
#' @param add_intercept Should `data` be augmented with a column of `1`s to represent an intercept? If so, this column is added to the left of `data`.
#'
#' @return A vector of linear predictors with length equal to the number of rows in `data`.
#' @export
#' @keywords internal
#'
#' @examples
#' data = data.frame(X1 = c(1,0), X2 = c(0,1))
#'
#' # No intercept
#' beta1 = c(1,2)
#' lin_pred_contrib(data, beta1, add_intercept = FALSE)
#'
#' # With intercept
#' beta2 = c(1,2,3)
#' lin_pred_contrib(data, beta2, add_intercept = TRUE)
lin_pred_contrib <- function(data, beta, add_intercept = TRUE){
  # Validate dimensions of input
  p_data = ncol(data)
  p_beta = length(beta)
  if((p_data == p_beta) && (!add_intercept)){
    # Do nothing
  } else if((p_data == p_beta - 1) && (add_intercept)){
    # Do nothing
  } else{
    stop("In lin_pred_contrib(): Incompatible dimensions of data and beta.")
  }

  # Setup data for multiplication by beta
  data_mat = as.matrix(data)
  if(add_intercept) data_mat = cbind(1, data_mat)

  # Compute linear predictor
  eta_vec = data_mat %*% beta
  eta = as.numeric(eta_vec)
  return(eta)
}



make_REs <- function(Gamma){
  MASS::mvrnorm(1, rep(0, times = nrow(Gamma)), Gamma)
}




#' Compute the linear predictor based on provided datasets and parameters
#'
#' Note : Please specify either add_intercept or both add_intercept_fix and add_intercept_ran. Do not provide all three!
#'
#' @param data_fix A dataset containing covariates with fixed effects. Intercept can optionally be added later.
#' @param data_ran A dataset containing coveriates with random effects. Intercept can optionally be added later.
#' @param beta Vector of fixed effects coefficients.
#' @param Gamma Covariance matrix of the random effects coefficients.
#' @param add_intercept Should an intercept column be added to the datasets for both fixed and random effects?
#' @param add_intercept_fix,add_intercept_ran Should an intercept column be added to the datasets for fixed and random effects respectively?
#' @param return_REs Should generated random effects be returned?
#'
#' @return A vector containing a linear predictor for each observation in the provided datasets.
#' @export
#'
#' @examples
#' data_fix = data.frame(X1 = c(1,0), X2 = c(0,1))
#' data_ran = data_fix
#'
#' # Fixed and random effects for intercept
#' beta1 = c(1,2,3)
#' Gamma1 = diag(3)
#'
#' get_lin_preds(data_fix, data_ran, beta1, Gamma1,
#'   add_intercept = TRUE)
#'
#'
#' # Fixed effect only for intercept
#' beta2 = c(1,2,3)
#' Gamma2 = diag(2)
#'
#' get_lin_preds(data_fix, data_ran, beta2, Gamma2,
#'   add_intercept_fix = TRUE, add_intercept_ran = FALSE)
#'
#'
#' # Return generated random effects in the second example
#' get_lin_preds(data_fix, data_ran, beta2, Gamma2,
#'  add_intercept_fix = TRUE, add_intercept_ran = FALSE, return_REs = TRUE)
get_lin_preds <- function(data_fix, data_ran, beta, Gamma, add_intercept = NULL, add_intercept_fix = NULL, add_intercept_ran = NULL, return_REs = FALSE){
  # Check that fixed and random effects datasets have the same number of observations
  if(nrow(data_fix) != nrow(data_ran)) stop("In get_lin_preds(): Different number of observations in fixed and random effects datasets.")

  # Where should intercepts be added?
  if(is.logical(add_intercept) && is.null(add_intercept_fix) && is.null(add_intercept_ran)){
    add_intercept_fix = add_intercept
    add_intercept_ran = add_intercept
  } else if(is.null(add_intercept) && is.logical(add_intercept_fix) && is.logical(add_intercept_ran)){
    # Do nothing
  } else{
    stop("In get_lin_preds(): Please specify either add_intercept or both add_intercept_fix and add_intercept_ran, but not all three.")
  }


  # Compute fixed and random effects components of the linear predictor
  contrib_fix = lin_pred_contrib(data_fix, beta, add_intercept_fix)

  ran_effs = make_REs(Gamma)
  contrib_ran = lin_pred_contrib(data_ran, ran_effs, add_intercept_ran)

  lin_preds = contrib_fix + contrib_ran

  # Return linear predictors, optionally in a list alongside the random effects
  if(return_REs){
    output = list(lin_preds = lin_preds, REs = ran_effs)
    return(output)
  } else{
    return(lin_preds)
  }
}

