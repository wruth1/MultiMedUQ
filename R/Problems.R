# I require a very particular random effect structure. Specifically, the model for Y must contain REs for the intercept, X and M, and the model for M must contain REs for the intercept and X. This isn't a conceptual requirement, it's just how I've been doing my analysis so far. The problem arises when I want to analyze a dataset which doesn't play nicely with having this much structure. Specifically, the following problems can arise:

# In computing the ENC, I compute a quadratic form with the estimated RE covariance matrix and a vector of covariate values. Currently, those covariate values are hard-coded, e.g. c(1,x,1). This breaks when there are only two REs (i.e. the meat of the sandwich is a 2x2 matrix). This will eventually need to be made more flexible, but for now I mostly want to get something out the door.




# When testing the functions which construct asymptotic covariance matrices, the output is often not positive definite. This might be a mistake in the formulas, but I'm mostly just using other packages.  I think that the problem is with a nearly-singular information matrix; particularly because the problem goes away if I reduce the number of REs. For now, I'm just skipping the relevant tests.
    # I think that the problem is with insufficient sample size to invoke limit theory guaranteeing positive definiteness. More precisely, I expect that the number of groups is what needs to go to infinity, not the number of observations. See Gill, Jeff, and Gary King. "Numerical issues involved in inverting hessian matrices." Numerical issues in statistical computing for the social scientist (2004): 143-176 (https://gking.harvard.edu/files/gking/files/numhess.pdf) for a discussion
    # I'm just going to use simulated data for now, so I can control the number of groups. The machinery used to do this simulation is based on code from an old version of this project (ModelMediation), but does work fine for here. There are just differences between how I name certain things.


# There are some extra packages required for Generate_Sample_Data.R, which are not used elsewhere. Specifically, MASS, boot and purrr


# There may be some issues with for loops going from 2:(i-1) with i=2. I want this to be an empty loop, but R treats it as c(2,1). Maybe this won't be an issue, but I'm not sure. This is also hard to search for.


# My documentation for the labels of the random effects is copied across many, many functions. It would be better to move this to a vignette, even a short one, which I can then link in the @parameter description of which_REs.
