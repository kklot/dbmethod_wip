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
            afs_u = fitdt$biased_afs + 1,
            afs_l = fitdt$biased_afs,
            age = fitdt$age - min(fitdt$age),
            event = fitdt$event,
            yob = fitdt$yob - min(fitdt$yob),
            svw = fitdt$sv_weight,
            to_skew = to_skew,
            age_term = age_term,
            smooth_age = smooth_age,
            yob_term = yob_term,
            smooth_yob = smooth_yob,
            interval_censor = 0,
            to_fit = rep(1, nrow(fitdt)),
            sd_beta = c(0, 1),
            sd_yob = c(5e-1, 1e-1),
            sd_age = c(5e-1, 1e-1),
            shape_prior = c(0, 1),
            skew_prior = c(0, 1),
            ar_scale = ar_scale
        )

        # initial parameters
        parameters <- list(
            intercept = -3,
            log_shape = log(10),
            log_skew = log(1),
            beta_yob = .01,
            yob_rw2 = rep(0, length(unique(fitdt$yob))),
            yob_phi = rep(0, 2), # AR2
            beta_age = .01,
            age_rw2 = rep(0, length(unique(fitdt$age))),
            age_phi = 0
        )

        # change this depending on model, fixed parameters that are not estimated
        which_ones <- NULL

        if (!to_skew) {
            which_ones <- c(which_ones, "log_skew")
        }

        if (!yob_term) {
            which_ones <- c(which_ones, char(beta_yob, yob_rw2, yob_phi))
        } else {
            if (smooth_yob) {
                which_ones <- c(which_ones, "beta_yob")
            } else {
                which_ones <- c(which_ones, char(yob_rw2, yob_phi))
            }
        }

        if (!age_term) {
            which_ones <- c(which_ones, char(beta_age, age_rw2, age_phi))
        } else {
            if (smooth_age) {
                which_ones <- c(which_ones, "beta_age")
            } else {
                which_ones <- c(which_ones, char(age_rw2, age_phi))
            }
        }

        # set random effect
        re <- NULL
        fe_only <- TRUE
        re_only <- FALSE
        if (smooth_age | smooth_yob) {
            # there is case smoothing one of the two
            re <- char(intercept, log_shape, log_skew, beta_yob, beta_age, yob_rw2, age_rw2)
            fe_only <- FALSE
            re_only <- TRUE
        }

        # model configs
        openmp(1)
        options(tape.parallel = FALSE, DLL = "model")
        
        obj <- MakeADFun(
            data = idata, parameters = parameters,
            random = re, map = tmb_fixit(parameters, which_ones),
            silent = !verbose, DLL = "model"
        )
  
        fit <- nlminb(obj$par, obj$fn, obj$gr)
  
        if (svysmp == 1) {
            cat(head(names(fit$par), 10), "\n")
        }

        # sample posterior samples
        smp <- sample_tmb(list(obj = obj, fit = fit), S, fe_only, re_only)

        attributes(smp)$yob <- sort(unique(fitdt$yob))
        attributes(smp)$age <- sort(unique(fitdt$age))
        attributes(smp)$age_id <- grep("age_rw2$", colnames(smp))
        attributes(smp)$yob_id <- grep("yob_rw2$", colnames(smp))
        
        if (check) {
            browser()
            # Check if the age distribution is ok
            fitdt %>% ggplot(aes(age), alpha = .5) +
                geom_histogram() +
                facet_wrap(~svy)
            scale_i <- smp[1, "intercept"] +
                smp[1, attributes(smp)$yob_id] +
                smp[1, attributes(smp)$age_id[20]]            
            med_i <- qskewlogis(
                .5, scale_i %>% exp(),
                smp[1, "log_shape"] %>% exp(),
                smp[1, "log_skew"] %>% exp()
            )
            est_i <- tibble(yob = attributes(smp)$yob, med = med_i)
            fitdt %>%
                ggplot() +
                # geom_point(aes(yob, afs, color = factor(svy), alpha = factor(sv_weight)), se = F) +
                geom_smooth(aes(yob, afs, color = factor(svy)), se = F) +
                geom_line(aes(yob, median), lwd = 2, pdata) +
                geom_line(aes(yob, med), est_i, col = "red") +
                geom_vline(xintercept = c(1970, 2005))
        }

        smp
    } # end a sample

    posterior <- parallel::mclapply(1:K, loop_fn)
    posterior
}
