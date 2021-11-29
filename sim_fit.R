---
title: "DB Method"
author: "Kinh Nguyen"
date: "`r Sys.Date()`"
params:
  nsv: 3
  sample_size: 1000
  bias: "none"
  trend: "none"
---

```{r include=FALSE, message=FALSE, warning=FALSE}
# params = list(nsv=2, sample_size=1000, bias='none', trend='none')

library(abind)
library(TMB)
library(data.table)
library(tidyverse)
# remotes::install_github('kklot/ktools')
library(ktools)
library(loo)                # 

# Check if input from bash is correct
myname = unlist(params) %>% paste0(names(.), .) %>% paste0(collapse="_") %>% paste0(".rds")
if (file.exists(myname))
	stop("did it")
message(myname)

set.seed(123)
options(mc.cores = parallel::detectCores()-2)
```

## Pool population

```{r}
N             <- 10^6
birth_cohorts <- 1900:2021
elig_age      <- 15:49 # age eligible for including in the year of surveys

# take age distribution of SZ
sz = data.frame(
	age=c(15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32,
				33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49),
	n=c(525, 566, 514, 469, 443, 464, 351, 385, 321, 316, 262, 301, 273, 241,
			222, 245, 211, 188, 184, 185, 151, 184, 164, 235, 123, 151, 136, 144,
			109, 116, 108, 136, 124, 136, 81))

sz %>% {approxfun(.[,1], .[,2]/sum(.[,2]), yleft=0, yright=0)} -> age_weight

my_pop <- data.table(id  = 1:N, yob = sample(birth_cohorts, N, T), afs=0)
```

# Load dll

```{r}
openmp(1) # we do multiple fit on parallel so avoiding using it in TMB
compile('model.cpp')
dyn.load(dynlib('model'))
invisible(config(tape.parallel=FALSE, DLL='model'))

source("tmb_sampling.R")
```

# AFS parameters and sampling

Following skew log-logistic distribution, at the begining, year 1900

```{r}
ref = list(scale = 0.06204, shape = 10, skew  = 1.5)
qskewlogis(.5, ref$scale, ref$shape, ref$skew) 
min_scale = 0.0703  # 15
qskewlogis(.5, min_scale, ref$shape, ref$skew) 
max_scale = 0.05859 # 18
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

## Generate true afs
afsd = my_pop %>% 
	left_join(pdata, "yob") %>% 
  arrange(yob) %>% 
  group_by(yob) %>% 
  mutate(afs = rskewlogis(n(), scale, shape, skew)) %>% 
  select(id, yob, afs)
```

# Generate biases

```{r}
bias_none = function(age) 0
bias_norm = function(age, sd = 0.3) rnorm(1, 0, sd)
bias_logis = function(age, max = 4, r = .5, mid = 23) max / (1 + exp(-r*(age - mid))) - max/2
# bias_logis(15:49, r = .5) %>% plott(15:49, .)

if (params$bias == "none") bias_f = bias_none
if (params$bias == "norm") bias_f = bias_norm
if (params$bias == "logis") bias_f = bias_logis

plott(bias_f(15:49))
```

# Generate survey data

```{r}
afsd = afsd %>%
  uncount(params$nsv, .id = "svy") %>% 
  mutate(svy = 2000 + 5 * (svy - 1))  %>% 
  mutate(age = svy - yob) %>%
	filter(age %in% 15:49) %>% 
  mutate(bias = bias_f(age), biased_afs = afs + bias) %>% 
  mutate(weight = age_weight(age)) %>% 
  mutate(event = if_else(biased_afs <= age, 1, 0))
```

# Fit model

```{r get_post}
source("get_posterior.R")
 
# parallel
post = get_posterior(data=afsd, sample_size=params$sample_size, K=100)
attributes(post)$ref_par = ref 
```

# Get some stats

```{r get_stat}
X = attributes
ekld = function(x) x %>% aperm(c(1,3,2)) %>% rowMeans(dims = 2) %>% colSums %>% mean

# apply to one scenario with 1000 posterior samples

x = post[[1]]

stats = function(x, agref = 1) {
	realskew = x[, 'log_skew'] %>% exp
	# skewe = realskew %>% quantile95 %>% t %>% as.data.table(1) %>% rename(skew_lo=lo, skew_up=up, skew_med=med)
	realshape = x[, 'log_shape'] %>% exp
	# shapee = realshape %>% quantile95 %>% t %>% as.data.table(1) %>% rename(shape_lo=lo, shape_up=up, shape_med=med)

	ageref  = tibble(age=X(x)$age, id=X(x)$age_id) %>% mutate(agr = findInterval(age, seq(15, 50, 5))) 
	agecoef = x[, ageref %>% filter(agr==agref) %>% pull(id)] %>% rowMeans

	realscale = exp(sweep(x[, attributes(x)$yob_id], 1, x[, 'intercept'], '+') + agecoef) 
	qskewlogis(.5, realscale[1, ], realshape[1], realskew[1]) %>% histt


	# qskewlogis(.5, 1/20, 9, 1.1)
	kld = lapply(1:1000, getQld)
	list(kld = kld, cvr = cvr)
}

getstats = function(aref) 
{
	tst = parallel::mclapply(post, stats, agref = aref) 
	kldi = sapply(seq_along(tst), function(x) 
								sapply(1:1000, function(y) 
											 tst[[x]]$kld[[y]]) %>% colSums %>% mean) %>% mean

	lapply(seq_along(tst), function(x) tst[[x]]$cvr) %>% rbindlist(idcol = "sim") %>% 
		group_by(sim) %>% 
		summarise(ishape = sum(inside_shape)/n(),
							iscale = sum(inside_scale)/n(),
							iskew = sum(inside_skew)/n(), 
							shape_b = sum(shape_med - shape)/n(), 
							skew_b = sum(skew_med - skew)/n(), 
							scale_b = sum(scale_med - scale)/n(), 
							)  %>% 
		select(-sim)  %>% 
		colMeans  %>% 
		c(kldi, nsv = params$nsv, size = params$sample_size, trend = params$trend, bias = params$bias) 
}
	
out1 = getstats(1)
out2 = getstats(2)
out3 = getstats(3)
out4 = getstats(4)

out = list(out1, out2, out3, out4)

getwd()
saveRDS(out, myname)
# saveRDS(out, paste0("output/", myname))
```

