library(ktools)
library(tidyverse)
library(data.table)
birth_cohorts <- 1900:2021
elig_age      <- 15:49 # age eligible for including in the year of surveys
ref = list(scale = 0.06204, shape = 10, skew  = 1.5)
qskewlogis(.5, ref$scale, ref$shape, ref$skew) 
min_scale = 0.0703  # 15
qskewlogis(.5, min_scale, ref$shape, ref$skew) 
max_scale = 0.05859 # 18
qskewlogis(.5, max_scale, ref$shape, ref$skew) 

refnone = tibble(yob = birth_cohorts, scale = ref$scale, skew = ref$skew, shape = ref$shape)
scalev = seq(min_scale, max_scale, length.out = length(birth_cohorts))
refincrease = tibble(yob = birth_cohorts, scale = scalev, skew = ref$skew, shape = ref$shape)
scalev = seq(max_scale, min_scale, length.out = length(birth_cohorts))
refdecrease  = tibble(yob = birth_cohorts, scale = scalev, skew = ref$skew, shape = ref$shape)

X = attributes
ekld = function(x) x %>% aperm(c(1,3,2)) %>% rowMeans(dims = 2) %>% colSums %>% mean

toatable = function(x) 
{
	readRDS(x) %>% 
		lapply(t) %>% 
		lapply(as.data.table, keep.rownames=1) %>% 
		list(.id = "age_ref") %>% 
		do.call("bind_rows", .)
}

os = list.files("fuchs", '.*rds$', full.names=T) %>%
	lapply(toatable) %>% 
	do.call("bind_rows", .) %>% 
	mutate(across(age_ref:size, as.double)) %>% 
	mutate(age_ref = factor(age_ref, 1:4, c('15-19', '20-24', '25-29', '30-34')))

osw = os %>% 
	rename(shape_cov = ishape, scale_cov = iscale, skew_cov = iskew) %>% 
	pivot_longer(2:7, names_to = char(par, stat), names_sep = "_") %>% 
	mutate(par = str_to_title(par))  

newpal = c('#EC768C', '#F6D3D9', '#DCF4EA', '#DAF1F9', '#B8D6B4', '#95DEC1',
					 '#8DD5EB', '#FAE29B', '#CCCAA1', '#DACF9B', '#E3DEB7', '#E8EFCF')

library(ktools)
options(ggplot2.discrete.fill =  hcl.colors(8), 
				ggplot2.discrete.color = hcl.colors(8))


g = os %>% 
	rename(KLD = V7) %>%
	mutate(nsv = paste0('Number of surveys = ', nsv)) %>% 
	mutate(size = paste0('Sample size/survey = ', size)) %>% 
	ggplot() +
	geom_line(aes(age_ref, KLD, color = trend, linetype = bias, group = interaction(trend, bias))) +
	geom_point(aes(age_ref, KLD, color = trend, linetype = bias, group = interaction(trend, bias))) +
	facet_wrap(char(nsv, size)) + labs(x = "Age's effect reference group") +
	scale_y_continuous(trans = 'log', labels = function(x) sprintf("%.2f",x)) +
	theme(legend.position = 'bottom', legend.spacing.x = unit(.01, "cm"),
				strip.background = element_rect(fill='grey95', color = NA)) +
	coord_cartesian(xlim=c(0.9, 4.2), expand=FALSE, clip = "off")

quartz_off(g, 'fig/KLD', 7, 7)

library(ggpubr)

bap = osw %>% 
	filter(stat == 'b')  %>% 
	filter(nsv == 2 & size == 1000) %>% 
	rename(Bias = value) %>% 
	ggplot() +
	geom_line(aes(age_ref, Bias, color = trend, linetype = bias, group=interaction(trend, bias)), size=.71) + 
	geom_point(aes(age_ref, Bias, color = trend, linetype = bias, group=interaction(trend, bias)), size=.71) + 
	facet_wrap('par', scale = 'free_y') +
	labs(x = "Reference age-group", color = "Trend", linetype = "Bias")
covp = osw %>% 
	filter(stat == 'cov')  %>% 
	filter(nsv == 2 & size == 1000) %>% 
	rename(Coverage = value) %>% 
	ggplot() +
	geom_line(aes(age_ref, Coverage, color = trend, linetype = bias, group=interaction(trend, bias)), size=.71) + 
	geom_point(aes(age_ref, Coverage, color = trend, linetype = bias, group=interaction(trend, bias)), size=.71) + 
	facet_wrap('par') +
	labs(x = "", color = "Trend", linetype = "Bias")

g = ggarrange(covp, bap, nrow = 2, ncol = 1, common.legend = T, legend = "bottom")
quartz_off(g, "fig/Bias_Coverage", 7, 7) 
# resize(7,5)
