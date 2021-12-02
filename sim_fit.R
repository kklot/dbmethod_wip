library(abind)
library(TMB)
library(data.table)
library(tidyverse)
library(ktools)
library(loo)

set.seed(123)
options(mc.cores = parallel::detectCores()-2)

# cluster stuffs
task_id <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
base_dir <- "/scratch/fuchs/fias/knguyen/"
save_to  <- paste0(base_dir, "db_mt_fix_weight_smoother/")
dir.create(save_to, FALSE)

scenarios <- crossing(
    nsv = 2:7,
    sample_size = 1000,
    bias = char(none, logis),
    trend = char(increase, decrease, none),
    theK = 100
)
scenarios

params <- as.list(scenarios[task_id, ])

myname <- unlist(params) %>%
    paste0(names(.), .) %>%
    paste0(collapse = "_") %>%
    paste0(".rds")
mypath <- paste0(save_to, myname)

## Pool population
N <- 10^6
birth_cohorts <- 1940:2005
wanted_cohort <- 1970:2005
elig_age      <- 15:49 # age eligible for including in the year of surveys

# take age distribution of SZ
sz = data.frame(
    age = c(
        15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32,
        33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49
    ),
    n = c(
        525, 566, 514, 469, 443, 464, 351, 385, 321, 316, 262, 301, 273, 241,
        222, 245, 211, 188, 184, 185, 151, 184, 164, 235, 123, 151, 136, 144,
        109, 116, 108, 136, 124, 136, 81
    )
)
sz %>% {approxfun(.[,1], .[,2]/sum(.[,2]), yleft=0, yright=0)} -> age_weight

my_pop <- data.table(id  = 1:N, yob = sample(birth_cohorts, N, T), afs=0)

# choose survey year 5 year apart starting backward from 2020
# maximum 7
chosen_svy <- crossing(svy = 1990:2020, bch = wanted_cohort) %>%
    mutate(age = svy - bch) %>%
    filter(age >= 15) %>%
    group_by(svy) %>%
    summarise(min = min(age), max = max(age)) %>%
    arrange(desc(svy)) %>%
    mutate(id = 1:nrow(.) %% 5) %>%
    filter(svy == 2020 | id == 0) %>%
    slice_head(n = params$nsv) %>%
    pull(svy)

chosen_svy

# Load dll
# we do multiple fit on parallel so avoiding using it in TMB
openmp(1)
compile("model.cpp")
dyn.load(dynlib("model"))
invisible(config(tape.parallel=FALSE, DLL='model'))
source("tmb_sampling.R")

# AFS parameters and sampling
# Following skew log-logistic distribution, at the begining, year 1900
ref = list(scale = 0.06204, shape = 10, skew  = 1.5)
qskewlogis(.5, ref$scale, ref$shape, ref$skew) # 17

min_scale = 0.07
qskewlogis(.5, min_scale, ref$shape, ref$skew)

max_scale = 0.0555
qskewlogis(.5, max_scale, ref$shape, ref$skew)

if (params$trend == "none") {
	pdata = tibble(yob = birth_cohorts, scale = ref$scale, skew = ref$skew, shape = ref$shape)
} else if (params$trend == "increase") {
	scalev = seq(min_scale, max_scale, length.out = length(birth_cohorts))
	pdata  = tibble(yob = birth_cohorts, scale = scalev, skew = ref$skew, shape = ref$shape)
} else if (params$trend == "decrease") {
	scalev = seq(max_scale, min_scale, length.out = length(birth_cohorts))
	pdata  = tibble(yob = birth_cohorts, scale = scalev, skew = ref$skew, shape = ref$shape)
}

pdata %<>% mutate(median = qskewlogis(.5, scale, shape, skew))

# Generate biases
bias_f <- switch(params$bias,
    "none" = function(age) 0,
    "norm" = function(age, sd = 0.3) rnorm(1, 0, sd),
    "logis" = function(age, max = 4, r = .5, mid = 23) max / (1 + exp(-r * (age - mid))) - max / 2,
    "women" = function(age) .5 - bias_lgt(age, , 2.5),
    "men" = function(age) .5 + bias_lgt(age, r = .25)
)

# Generate pooled and survey data
## Generate true afs
afsd = my_pop %>%
    left_join(pdata, "yob") %>%
    arrange(yob) %>%
    group_by(yob) %>%
    mutate(afs = rskewlogis(n(), scale, shape, skew)) %>%
    select(id, yob, afs)

afsd <- afsd %>%
    uncount(params$nsv, .id = "svy") %>%
    mutate(
        svy = chosen_svy[svy],
        age = svy - yob
    ) %>%
    filter(age %in% 15:49) %>%
        mutate(
            bias = bias_f(age) + rnorm(n(), 0, 0.3),
            biased_afs = afs + bias,
            sampling_weight = age_weight(age),
            event = if_else(biased_afs <= age, 1, 0)
        )

# Fit model
source("get_posterior.R")
 
# parallel within this
post = get_posterior(
    data = afsd, sample_size = params$sample_size, K = params$theK, S = 50,
    check = T, 
    yob_order = 2
)

attributes(post)$ref_par = ref 

saveRDS(post, mypath)
