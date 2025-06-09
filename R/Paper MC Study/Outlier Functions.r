
#* Use Grubbs' Test to detect outliers

# Like negative indexing, but can handle inds == c()
remove_indices <- function(X, inds){
    if(is.null(inds)){
        return(X)
    } else{
        return(X[-inds])
    }
}

remove_rows <- function(mat, inds){
    if(is.null(inds)){
        return(mat)
    } else{
        return(mat[-inds,])
    }
}

# Return indices of all outliers as determined by Grubbs' Test
find_all_outliers <- function(X, sig_level = 1e-5){
    outlier_inds = c()
    done = F

    while(!done){
        X_smaller = remove_indices(X, outlier_inds)

        grubbs_info = grubbs.test(X_smaller)
        p_val = grubbs_info$p.value

        if(p_val > sig_level){
            done = T
        } else{
            this_outlier_val = max(X_smaller)
            this_outlier_ind = which(X == this_outlier_val)
            outlier_inds = c(outlier_inds, this_outlier_ind)
        }

        # hist(X_smaller)

        # print(p_val)
    }

    return(outlier_inds)
}
