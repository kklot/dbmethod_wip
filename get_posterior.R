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
            age_phi = 0,
            log_ar_precision = log(1)
        )

        # model configs
        invisible(suppressWarnings(openmp(1)))
        options(tape.parallel = FALSE, DLL = "model")
        
        obj <- MakeADFun(
            data = idata, parameters = parameters,
            random = char(intercept, log_shape, log_skew, yob_rw2, age_rw2, age_phi, yob_phi),
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
            # Check if the age distribution is ok
            fitdt %>% ggplot(aes(age), alpha = .5) +
                geom_histogram() +
                facet_wrap(~svy)

            sml <- obj$simulate(, 1)
            sml$y2 <- sapply(1:length(sml$scale), \(x) rskewlogis(1, sml$scale[x], sml$shape, sml$skew))

            scale_i <- smp[1, "intercept"] +
                smp[1, attributes(smp)$yob_id] +
                smp[1, attributes(smp)$age_id[20]]            
            med_i <- qskewlogis(
                .5, scale_i %>% exp(),
                smp[1, "log_shape"] %>% exp(),
                smp[1, "log_skew"] %>% exp()
            )
            est_i <- tibble(yob = attributes(smp)$yob, med = med_i)
            est_i %>% filter(yob %in% range(wanted_cohort))            
            fitdt %>% bind_cols(y=sml$y2) %>%
                group_by(svy, yob) %>%
                summarise(afs = median(afs), y = median(y)) %>%
                ggplot() +
                geom_rect(aes(xmin = 1970, xmax = 2005, ymin = -Inf, ymax = Inf), fill = "#f0e5d9", alpha = .1) +
                geom_smooth(aes(yob, afs, color = factor(svy)), lwd = .7, se = F) +
                geom_line(aes(yob, y, color = factor(svy)), lty="dashed", lwd = .7) +
                geom_line(aes(yob, median, color="True"), lwd = 2, alpha=.7, pdata %>% filter(yob %in% wanted_cohort)) +
                geom_line(aes(yob, med, color="Fitted AR2"), est_i) +
                scale_color_manual(values = gen_colors(n = 9)) +
                labs(title = "AR(2)") -> g            
            ggsave("AR2xxx.pdf", g, width = 7, height = 7)            
            open_file("AR2xxx.pdf")

            tibble(yhat = sml$y2, y = sml$afs, event = sml$event) %>%
                ggplot() +
                geom_point(aes(y, yhat)) +
                facet_wrap(~event)

            put(2, 2)
            plotp(sml$yob_rw2)
            plotp(smp[1, attributes(smp)$yob_id])
            plotp(sml$age_rw2)
            plotp(smp[1, attributes(smp)$age_id])
        }

        smp
    } # end a sample

    posterior <- parallel::mclapply(1:K, loop_fn)
    posterior
}
