library(abind)
get_lli = function(to_skew = 1,
									 yob_term = 1,
									 age_term = 1,
									 smooth_age = 1,
									 smooth_yob = 1,
									 age_order = 1,
									 yob_order = 1,
									 weightv = "none",
									 interval_censor = 0,
									 data, 
                   sample_size, 
                   ref_age = NULL,
                   ref = list(scale = rep(0.06204, 45), shape = 10, skew  = 1.5),
                   K = 10, S = 1000, verbose = F, check = F)
{

	if (weightv == "none") 
		data$weight = 1

	# initial parameters
	parameters = list(intercept     = -3,
										log_shape     = log(10),
										log_skew      = log(1),
										beta_yob      = .01,
										yob_rw2       = rep(0, length(unique(data$yob))),
										log_yob_rw2_e = log(.1),
										beta_age      = .01,
										age_rw2       = rep(0, length(unique(data$age))),
										log_age_rw2_e = log(.1)
	)

	ll = NULL

	for (svysmp in 1:K) 
	{
		cat(svysmp,"...")
    # Generate a survey dataset
    fitdt = data %>% 
      group_by(svy) %>% 
      mutate(take = id %in% sample(id, sample_size, prob = weight)) %>% 
      filter(take)  %>% 
      select(svy, afs, biased_afs, age, event, yob, weight)

		# all data and meta data
    idata = list(afs             = fitdt$biased_afs,
                 afs_u           = fitdt$biased_afs + 1,
                 afs_l           = fitdt$biased_afs,
                 age             = fitdt$age-min(fitdt$age),
                 event           = fitdt$event,
                 yob             = fitdt$yob-min(fitdt$yob),
                 svw             = fitdt$weight,
                 to_skew         = to_skew,
                 age_term 			 = age_term,
                 smooth_age      = smooth_age,
                 yob_term        = yob_term,
                 smooth_yob      = smooth_yob,
                 interval_censor = 0,
                 to_fit          = rep(1, nrow(fitdt)),
                 sd_beta         = c(0, 1),
                 sd_yob          = c(1, 1e-3),
                 sd_age          = c(1, 1e-3),
                 age_order       = age_order,
                 yob_order       = yob_order,
                 shape_prior     = c(0, 1),
                 skew_prior      = c(0, 1),
                 R_age           = INLA:::inla.rw(length(unique(fitdt$age)), 1, F, T),
                 R_yob           = INLA:::inla.rw(length(unique(fitdt$yob)), 1, F, T))

		# change this depending on model, fixed parameters that are not estimated
		which_ones = NULL

		if (!to_skew)
			which_ones = c(which_ones, "log_skew")

		if (!yob_term)
			which_ones = c(which_ones, char(beta_yob, yob_rw2, log_yob_rw2_e))
		else {
			if (smooth_yob)
				which_ones = c(which_ones, "beta_yob")
			else
				which_ones = c(which_ones, char(yob_rw2, log_yob_rw2_e))
		}

		if (!age_term)
			which_ones = c(which_ones, char(beta_age, age_rw2, log_age_rw2_e))
		else {
			if (smooth_age)
				which_ones = c(which_ones, "beta_age")
			else
				which_ones = c(which_ones, char(age_rw2, log_age_rw2_e))
		}

		# set random effect
		re = NULL
		fe_only = TRUE
		re_only = FALSE
		if (smooth_age | smooth_yob) 
		{
      # there is case smoothing one of the two
			re = char(intercept, log_shape, log_skew, beta_yob, beta_age, yob_rw2, age_rw2)
			fe_only = FALSE
			re_only = TRUE
		}

		# model configs
		obj = MakeADFun(data 			 = idata, parameters = parameters,
										random 		 = re, map = tmb_fixit(parameters, which_ones),
										silent 		 = !verbose, DLL = 'model')
		fit = nlminb(obj$par, obj$fn, obj$gr)
		if (check) 
      browser()
		if (svysmp == 1) 
			cat(head(names(fit$par), 10), "\n")
    
    # sample posterior samples
		smp = sample_tmb(list(obj = obj, fit = fit), S, fe_only, re_only)

    yob_eff = grep('yob_rw2$', colnames(smp))
    age_eff = grep('age_rw2$', colnames(smp))

    eta = smp[, 'intercept'] + smp[, yob_eff]
    # how to treat age effect?
    if (!is.null(ref_age) && ref_age == 0) # consider age_eff = 0 is no bias
    {  
      age_eta = 0
      eta = eta + age_eta
    }
    else if (!is.null(ref_age) && ref_age != 0) # consider at ref_age there is no bias
    {
      age_eta = smp[, age_eff[which(15:49 == ref_age)]]
      eta = sweep(eta, 1, age_eta, '+')
    }

    scale = exp(eta)
    shape = exp(smp[, grep('log_shape', colnames(smp))])
    skew  = exp(smp[, grep('log_skew', colnames(smp))])

    # Kullback Liebler Distance
    ages  = seq(1, 100, 1)
    # the true parameters
    the_P  = sapply(ages, dskewlogis, scale = ref$scale, shape = ref$shape, skew = ref$skew)
    the_KLD = sapply(1:S, 
                function(s) {
                  the_Q = sapply(ages, dskewlogis, scale = scale[s,], shape = shape[s], skew = skew[s])
                  sapply(1:45, function(y) KLD(the_P[y,], the_Q[y,])) 
                }) 
    ll = abind(ll, the_KLD, along = 3) 
	} # end a sample
	cat("\n")
	ll
}
