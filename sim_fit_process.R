task_id <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

library(abind)
library(TMB)
library(data.table)
library(tidyverse)
library(ktools)

set.seed(123)
# short hand
X = attributes

file <- list.files(pattern = "*.rds")[task_id]

post <- readRDS(file)

params <- c()
params$nsv <- file %>% gsub(pattern = ".*nsv([0-9]{1})_.*", replacement = "\\1") %>% as.numeric()
params$sample_size <- file %>% gsub(pattern = ".*sample_size([0-9]{4})_.*", replacement = "\\1") %>% as.numeric()
params$bias <- file %>% gsub(pattern = ".*bias(.*?)_.*", replacement = "\\1")
params$trend <- file %>% gsub(pattern = ".*trend(.*?).rds", replacement = "\\1")
params$theK <- file %>% gsub(pattern = ".*K([0-9]{1,3})_.*", replacement = "\\1") %>% as.numeric()

## Pool population

birth_cohorts <- 1940:2005
wanted_cohort <- 1970:2005

# AFS parameters and sampling
# Following skew log-logistic distribution, at the begining, year 1900

ref = list(scale = 0.06204, shape = 10, skew  = 1.5)
qskewlogis(.5, ref$scale, ref$shape, ref$skew) # 17

min_scale = 0.07; qskewlogis(.5, min_scale, ref$shape, ref$skew) 
max_scale = 0.0555; qskewlogis(.5, max_scale, ref$shape, ref$skew) 

if (params$trend == "none") {
    pdata = tibble(yob = birth_cohorts, scale = ref$scale, skew = ref$skew, shape = ref$shape)
} else if (params$trend == "increase") {
    scalev = seq(min_scale, max_scale, length.out = length(birth_cohorts))
    pdata = tibble(yob = birth_cohorts, scale = scalev, skew = ref$skew, shape = ref$shape)
} else if (params$trend == "decrease") {
    scalev = seq(max_scale, min_scale, length.out = length(birth_cohorts))
    pdata = tibble(yob = birth_cohorts, scale = scalev, skew = ref$skew, shape = ref$shape)
}

pdata %<>% mutate(median = qskewlogis(.5, scale, shape, skew))

# Get some stats

# apply to one scenario with 1000 posterior samples
stats = function(x, agref = 1) {
    # skew and shape 
    realskew = x[, "log_skew"] %>% exp()
    realshape = x[, "log_shape"] %>% exp()
    # scale depending on age reference
    ageref = tibble(age = X(x)$age, id = X(x)$age_id) %>% mutate(agr = findInterval(age, seq(15, 50, 5)))
    agecoef = x[, ageref %>% filter(agr == agref) %>% pull(id)] %>% rowMeans()
    realscale = exp(sweep(x[, X(x)$yob_id], 1, x[, "intercept"], "+") + agecoef)
    # median differences
    getmed = function(y, x) {
        tibble(
            med = qskewlogis(.5, realscale[y, ], realshape[y], realskew[y]),
            yob = X(x)$yob
        ) %>%
            left_join(pdata, "yob") %>%
            filter(yob %in% range(wanted_cohort)) %>%
            mutate(diff = med - median) %>%
            {
                tibble(
                    first_diff = .$diff[1],
                    last_diff = .$diff[2],
                    trend_diff = .$med[2] - .$med[1],
                    real_diff = .$median[2] - .$median[1]
                )
            }
    }    

    st = lapply(1:nrow(x), getmed, x = x) %>% bind_rows()

    st %>%
        mutate(
            up = quantile(trend_diff, prob = .975),
            lo = quantile(trend_diff, prob = .025),
            inside95 = real_diff <= up & real_diff >= lo,
            up = quantile(trend_diff, prob = .75),
            lo = quantile(trend_diff, prob = .25),
            insideIQR = real_diff <= up & real_diff >= lo
        ) %>% summarise(
            ave_first = mean(first_diff),
            ave_last = mean(last_diff),
            ave_trend = mean(trend_diff),
            cov_95 = mean(inside95),
            cov_iqr = mean(insideIQR)
        )
}

# parallel::detectCores()
# apply to the number of survey sample of a scenario
getstats = function(aref) {
    parallel::mclapply(post, stats, agref = aref, mc.cores = 40) %>%
        bind_rows(.id = "sample") %>%
        mutate(
            ageref = aref, nsv = params$nsv, size = params$sample_size,
            trend = params$trend, bias = params$bias, nsample = params$theK
        )
}
	
out1 = getstats(1)
out2 = getstats(2)
out3 = getstats(3)
out4 = getstats(4)

out = bind_rows(out1, out2, out3, out4)

dir.create("processed", FALSE)
saveRDS(out, paste0('processed', "/", file))