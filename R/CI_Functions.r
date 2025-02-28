
cov_mat_2_SEs <- function(cov_mat){
    sqrt(diag(cov_mat))
}


# One parameter at a time
build_CIs_one_par <- function(Theta_hats, SE){
    lcls = Theta_hats - qnorm(0.975) * SE
    ucls = Theta_hats + qnorm(0.975) * SE

    return(list(lcl = lcls, ucl = ucls))
}


all_CIs_one_dataset <- function(Theta_hats, cov_mat){
    SEs = cov_mat_2_SEs(cov_mat)
    some_CIs = lapply(seq_along(Theta_hats), function(i) build_CIs_one_par(Theta_hats[i], SEs[i]))

    return(some_CIs)
}

build_many_CIs <- function(Theta_hats_data, SEs){
    lapply(seq_along(SEs), function(i){
        this_Theta_hats = Theta_hats_data[,i]
        this_SE = SEs[i]

        return(build_CIs_one_par(this_Theta_hats, this_SE))
    })
}

# One dataset at a time
build_one_CI_many_pars <- function(Theta_hats, SEs){
    lapply(seq_along(Theta_hats), function(i){
        this_Theta_hat = Theta_hats[i]
        this_SE = SEs[i]

        return(build_CIs_one_par(this_Theta_hat, this_SE))
    })
}

## SEs can either be a vector of SEs or a list of such vectors with length == nrow(Theta_hats_data)
build_many_CIs_by_dataset <- function(Theta_hats_data, SEs){
    ## Check inputs
    if(is.list(SEs)){
        if(length(SEs) != nrow(Theta_hats_data)){
            stop("Number of estimates does not match number of SE vectors.")
        } else{
            if(length(SEs[[1]]) != ncol(Theta_hats_data)){
                stop("Number of parameters does not match number of SEs.")
            }
        }
        SE_list = TRUE
    } else if(is.numeric(SEs)){
        if(length(SEs) != ncol(Theta_hats_data)){
            stop("Number of parameters does not match number of SEs.")
        }
        SE_list = FALSE
    } else{
        stop("SEs must be either a vector or a list of vectors.")
    }

    ## Build CIs
    lapply(seq_len(nrow(Theta_hats_data)), function(i){
        this_Theta_hats = Theta_hats_data[i,]
        if(SE_list){
            this_SEs = SEs[[i]]
        } else{
            this_SEs = SEs
        }

        return(build_one_CI_many_pars(this_Theta_hats, this_SEs))
    })
}


# Check coverage
coverage_checks_one_par <- function(Theta, CIs){
    some_lcls = CIs$lcl
    some_ucls = CIs$ucl

    some_checks = some_lcls < Theta & Theta < some_ucls

    return(some_checks)
}

many_coverage_checks <- function(Theta, CIs){
    lapply(seq_along(CIs), function(i){
        this_CIs = CIs[[i]]
        this_Theta = Theta[i]

        return(coverage_checks_one_par(this_Theta, this_CIs))
    })
}

coverage_checks_2_rates <- function(coverage_checks_list){
    sapply(coverage_checks_list, function(x) mean(x))
}

#' Compute coverage rates for Wald CIs using a single covariance matrix for each estimate
#'
#' @param Theta_hats_data A data frame containing one estimated parameter per column and one dataset per row
#' @param cov_mat A covariance matrix to use with each estimate
#' @param true_Thetas True values of each parameter for checking coverage
#'
#' @returns A vector of coverage rates, one component for each parameter
#' @export
get_coverage_rates <- function(Theta_hats_data, cov_mat, true_Thetas){
    SEs = cov_mat_2_SEs(cov_mat)
    CIs = build_many_CIs(Theta_hats_data, SEs)
    coverage_checks = many_coverage_checks(true_Thetas, CIs)

    return(coverage_checks_2_rates(coverage_checks))
}

#' Compute coverage rates for Wald CIs using a separate covariance matrix for each estimate
#'
#' @param Theta_hats_data A data frame containing one estimated parameter per column and one dataset per row
#' @param cov_mat_list A list of covariance matrices, one for each dataset
#' @param true_Thetas True values of each parameter for checking coverage
#'
#' @returns A vector of coverage rates, one component for each parameter
#' @export
get_coverage_rates_many_cov_mats <- function(Theta_hats_data, cov_mat_list, true_Thetas){
    SEs_list = lapply(cov_mat_list, cov_mat_2_SEs)

    # Nesting structure is [dataset][parameter][lcl/ucl]
    # Needs to be [parameter][lcl/ucl][dataset]
    CIs_raw = build_many_CIs_by_dataset(Theta_hats_data, SEs_list)

    CIs = lapply(seq_len(ncol(Theta_hats_data)), function(i){
        this_lcls = sapply(seq_len(nrow(Theta_hats_data)), function(j) CIs_raw[[j]][[i]]$lcl) %>% unname()
        this_ucls = sapply(seq_len(nrow(Theta_hats_data)), function(j) CIs_raw[[j]][[i]]$ucl) %>% unname()

        return(list(lcl = this_lcls, ucl = this_ucls))
    })

    coverage_checks = many_coverage_checks(true_Thetas, CIs)
    coverage_rates = coverage_checks_2_rates(coverage_checks)

    return(coverage_rates)

}



# ---------------------------------------------------------------------------- #
#                                    Widths                                    #
# ---------------------------------------------------------------------------- #

#' Compute CI widths from endpoints
#'
#' @param CIs_one_par A list with two elements: `lcl` and `ucl`. Each component is a vector of CI endpoints.
#'
#' @returns A vector of CI widths
#' @export
get_widths_one_par <- function(CIs_one_par){
    some_lcls = CIs_one_par$lcl
    some_ucls = CIs_one_par$ucl

    some_widths = some_ucls - some_lcls

    return(some_widths)
}


#' Compute CI widths from endpoints for multiple parameters
#'
#' @param CIs_many_pars A list of lists. Each sub-list contains two elements: `lcl` and `ucl`. Each of `lcl` and `ucl` are vectors of CI endpoints.
#'
#' @returns A list of vectors of CI widths. Each vector contains the widths of all intervals for one parameter
#' @export
many_widths <- function(CIs_many_pars){
    lapply(seq_along(CIs_many_pars), function(i){
        this_CIs = CIs_many_pars[[i]]

        return(get_widths_one_par(this_CIs))
    })
}



#' Construct CIs and compute their widths
#'
#' @param Theta_hats_data A data frame containing estimates of each parameter. Each column is a parameter, and each row is a dataset.
#' @param cov_mat An estimated sampling covariance matrix for the estimates.
#' @param cov_mat_list A list of estimated sampling covariance matrices, one for each dataset.
#'
#' @name CI_Widths
#'
#' @returns A list of vectors containing widths of all intervals for each parameter. The list has one component per parameter, and the vector has one component per dataset.
#' @export
get_widths_one_cov_mat <- function(Theta_hats_data, cov_mat){
    SEs = cov_mat_2_SEs(cov_mat)
    CIs = build_many_CIs(Theta_hats_data, SEs)
    all_widths = many_widths(CIs)

    mean_widths = sapply(all_widths, mean)

    return(all_widths)
}



#' @rdname CI_Widths
#' @export
get_widths_many_cov_mats <- function(Theta_hats_data, cov_mat_list){
    SEs_list = lapply(cov_mat_list, cov_mat_2_SEs)

    # Nesting structure is [dataset][parameter][lcl/ucl]
    # Needs to be [parameter][lcl/ucl][dataset]
    CIs_raw = build_many_CIs_by_dataset(Theta_hats_data, SEs_list)

    CIs = lapply(seq_len(ncol(Theta_hats_data)), function(i){
        this_lcls = sapply(seq_len(nrow(Theta_hats_data)), function(j) CIs_raw[[j]][[i]]$lcl) %>% unname()
        this_ucls = sapply(seq_len(nrow(Theta_hats_data)), function(j) CIs_raw[[j]][[i]]$ucl) %>% unname()

        return(list(lcl = this_lcls, ucl = this_ucls))
    })

    all_widths = many_widths(CIs)

    return(all_widths)

}




#' Mean CI width for each parameter
#'
#' @param Theta_hats_data A data frame containing estimates of each parameter. Each column is a parameter, and each row is a dataset.
#' @param cov_mat An estimated sampling covariance matrix for the estimates.
#' @param cov_mat_list A list of estimated sampling covariance matrices, one for each dataset.
#'
#' @name Mean_CI_Widths
#'
#' @returns A vector of mean CI widths, one component for each parameter.
#' @export
mean_widths_one_cov_mat <- function(Theta_hats_data, cov_mat){
    all_widths = get_widths_one_cov_mat(Theta_hats_data, cov_mat)
    mean_widths = sapply(all_widths, mean)

    return(mean_widths)
}


#' @rdname Mean_CI_Widths
#' @export
mean_widths_many_cov_mats <- function(Theta_hats_data, cov_mat_list){
    all_widths = get_widths_many_cov_mats(Theta_hats_data, cov_mat_list)
    mean_widths = sapply(all_widths, mean)

    return(mean_widths)
}



#' SD of CI widths for each parameter
#'
#' @param Theta_hats_data A data frame containing estimates of each parameter. Each column is a parameter, and each row is a dataset.
#' @param cov_mat An estimated sampling covariance matrix for the estimates.
#' @param cov_mat_list A list of estimated sampling covariance matrices, one for each dataset.
#'
#' @name SD_CI_Widths
#'
#' @returns A vector of standard deviations of CI widths, one component for each parameter.
#' @export
SD_widths_one_cov_mat <- function(Theta_hats_data, cov_mat){
    all_widths = get_widths_one_cov_mat(Theta_hats_data, cov_mat)
    SD_widths = sapply(all_widths, sd)

    return(SD_widths)
}


#' @rdname SD_CI_Widths
#' @export
SD_widths_many_cov_mats <- function(Theta_hats_data, cov_mat_list){
    all_widths = get_widths_many_cov_mats(Theta_hats_data, cov_mat_list)
    SD_widths = sapply(all_widths, sd)

    return(SD_widths)
}
