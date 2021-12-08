get_posterior <- function(to_skew = 1,
                          yob_term = 1,
                          age_term = 1,
                          smooth_age = 1,
                          smooth_yob = 1,
                          weightv = "none",
                          interval_censor = 0,
                          ar_scale = .1,
                          data,
                          sample_size,
                          K = 1000, S = 100, verbose = F, check = F) {

    if (weightv == "none")
        data$sv_weight = 1

    loop_fn <- function(svysmp) {
        cat(svysmp, "...")
        # Generate a survey dataset

        fitdt <- data %>%
            group_by(svy) %>%
            mutate(take = id %in% sample(id, sample_size, prob = sampling_weight)) %>%
            filter(take) %>%
            select(svy, afs, biased_afs, age, event, yob, sv_weight)
        
        # all data and meta data
        idata <- list(
            afs = fitdt$biased_afs,
            age = fitdt$age - min(fitdt$age),
            event = fitdt$event,
            yob = fitdt$yob - min(fitdt$yob),
            sd_beta = c(0, 1),
            sd_yob = c(5e-1, 1e-1),
            sd_age = c(5e-1, 1e-1),
            shape_prior = c(0, 1),
            skew_prior = c(0, 1)
        )

        # initial parameters
        parameters <- list(
            intercept = -3,
            log_shape = log(10),
            log_skew = log(1),
            yob_rw2 = rep(0, length(unique(fitdt$yob))),
            yob_phi = rep(0, 2), # AR2
            age_rw2 = rep(0, length(unique(fitdt$age))),
            age_phi = rep(0, 2),
            log_ar_precision_yob = log(1e3),
            log_ar_precision_age = log(1e3)
        )

        # model configs
        invisible(suppressWarnings(openmp(1)))
        options(tape.parallel = FALSE, DLL = "model")
        
        obj <- MakeADFun(
            data = idata, parameters = parameters,
            random = char(yob_rw2, age_rw2),
            # random = char(intercept, log_shape, log_skew, yob_rw2, age_rw2, age_phi, yob_phi),
            silent = !verbose, DLL = "model"
        )

        fit <- nlminb(obj$par, obj$fn, obj$gr)

        # sample posterior samples
        smp <- sample_tmb(list(obj = obj, fit = fit), S, FALSE, TRUE)

        attributes(smp)$yob <- sort(unique(fitdt$yob))
        attributes(smp)$age <- sort(unique(fitdt$age))
        attributes(smp)$age_id <- grep("age_rw2$", colnames(smp))
        attributes(smp)$yob_id <- grep("yob_rw2$", colnames(smp))
        
        if (check) {
            browser()
        }
        
        smp
    } # end a sample

    posterior <- parallel::mclapply(1:K, loop_fn)
    posterior
}
