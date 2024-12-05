



#' Generate a dataset using the parametric bootstrap
#'
#' @description
#' Generate a bootstrap dataset with `n` observations in each of `K` groups. Variables include one outcome (Y), one mediator (M), one exposure (X), and two confounders (C1 and C2). All variables are binary. Default values are available for regression parameters (fixed effects and covariances of random effects).
#'
#' Order of variables is `Y`, `M`, `X`, `C1`, `C2`
#'
#' @param n Number of observations in each group.
#' @param K Number of groups
#' @param b_Y,b_M Coefficient vectors for \eqn{Y}-model and \eqn{M}-model, respectively.
#' @param theta_Y,theta_M Covariance parameters of random effects in \eqn{Y}-model and \eqn{M}-model, respectively. See details.
#' @param X_list A list of vectors of X values. One vector for each group.
#' @param all_Cs_list A list of data frames of confounders. One data frame for each group. Each data frame contains two columns, `C1` and `C2`.
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
#' make_bootstrap_data(n, K, all_reg_pars, output_list = TRUE)
#'
#' # Format output as a data.frame
#' make_bootstrap_data(n, K, all_reg_pars, output_list = FALSE)
make_bootstrap_data <-
  function(n, K, b_Y, theta_Y, b_M, theta_M,
           X_list, all_Cs_list,
           output_list = FALSE,
           return_REs = FALSE,
           which_REs = c("Y.Int", "Y.X", "Y.M", "M.Int", "M.X")) {
    # all_reg_pars = make_all_reg_pars()

    # Check that the number of variance parameters supplied matches which_REs
    if(num_Y_REs(which_REs) != theta2num_REs(theta_Y)){
      stop("In make_bootstrap_data: The number of variance parameters supplied for Y does not match which_REs.")
    }
    if(num_M_REs(which_REs) != theta2num_REs(theta_M)){
      stop("In make_bootstrap_data: The number of variance parameters supplied for M does not match which_REs.")
    }
    if(K != length(X_list)){
      stop("In make_bootstrap_data: The number of groups supplied does not match the length of X_list.")
    }
    if(K != length(all_Cs_list)){
      stop("In make_bootstrap_data: The number of groups supplied does not match the length of all_Cs_list.")
    }


    all_reg_pars = list(beta_M = b_M, Gamma_M = theta2Sigma(theta_M),
                        beta_Y = b_Y, Gamma_Y = theta2Sigma(theta_Y))


    # Generate data as a list
    data_list_output = make_bootstrap_data_list(n, K, all_reg_pars, X_list, all_Cs_list, return_REs, which_REs)
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


make_bootstrap_data_list <-
  function(n, K, all_reg_pars, X_list, all_Cs_list, return_REs = FALSE, which_REs = c("Y.Int", "Y.X", "Y.M", "M.Int", "M.X")) {
    output_list = purrr::map(1:K,
                             ~ make_one_group_bootstrap(n, all_reg_pars, X_list[[.x]], all_Cs_list[[.x]], return_REs, which_REs))

    if (!return_REs) {
      return(output_list)
    } else{
      data_list = purrr::map(output_list, "data")
      RE_list = purrr::map(output_list, "REs")
      return(list(data = data_list, all_REs = RE_list))
    }
  }


make_one_group_bootstrap <-
  function(n, all_reg_pars, X, all_Cs, return_REs = FALSE, which_REs = c("Y.Int", "Y.X", "Y.M", "M.Int", "M.X")) {

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
