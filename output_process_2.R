library(ktools)
library(tidyverse)
library(data.table)

newpal = c("#EC768C", "#F6D3D9", "#DCF4EA", "#DAF1F9", "#B8D6B4", "#95DEC1", "#8DD5EB", "#FAE29B", "#CCCAA1", "#DACF9B", "#E3DEB7", "#E8EFCF")

options(
    ggplot2.discrete.fill = hcl.colors(8),
    ggplot2.discrete.color = hcl.colors(8)
)

setwd("/Users/knguyen/GitHub/db_paper 2/simulation_study/")

real_diff_1970_2005 <- tibble(
    trend = char(none, increase, decrease),
    real_diff = c(0, 2.38565, -1.917808)
)

processed <- lapply(list.files("fuchs7", full.names = T), readRDS) %>%
    bind_rows(.id = "scenario") %>%
    mutate(
        ageref = factor(ageref, 1:4, c("15-19", "20-24", "25-29", "30-34")),
        nsv = factor(nsv)
    ) %>%
    left_join(real_diff_1970_2005, "trend") %>%
    mutate(
        trend = paste0("Imposed trend: ", trend),
        bias = paste0("Imposed bias: ", bias)
    )
    
processed %>%
    ggplot() +
    geom_boxplot(aes(ageref, ave_trend, fill = nsv)) +
    geom_hline(aes(yintercept = real_diff), linetype = "dashed") +
    facet_grid(vars(bias), vars(trend),scales = 'free') +
    theme_bw() +
    theme(legend.position = "bottom", strip.background = element_rect(fill = "grey95", color = "grey")) +
    guides(fill = guide_legend(nrow = 1)) +
    scale_fill_viridis_d() +
    # scale_y_continuous(breaks = seq(-2, 5.5, 1)) +
    labs(
        x = "Referenced age-group", y = "Differences", fill = "Number of surveys",
        title = "Differences in median AFS between 1970-2005 birth cohorts"
    ) -> g

quartz_off(g, "trend_diff", 7, 7, open = 1)

g_last_diff <- processed %>%
    ggplot() +
    geom_boxplot(aes(ageref, ave_last, fill = nsv)) +
    geom_hline(aes(yintercept = 0), linetype = "dashed") +
    facet_grid(vars(bias), vars(trend)) +
    theme_bw() +
    theme(legend.position = "bottom", strip.background = element_rect(fill = "grey95", color = "grey")) +
    guides(fill = guide_legend(nrow = 1)) +
    scale_y_continuous(breaks = seq(-2, 5.5, 1)) +
    scale_fill_viridis_d() +
    labs(
        x = "Referenced age-group", y = "Differences",
        title = "Differences in median AFS of the 2005 birth cohort", 
		fill = "Number of surveys"
    )

quartz_off(g_last_diff, "last_diff", 7, 7, open = 1)

processed %>%
    group_by(ageref, trend, bias, nsv) %>%
    summarise(cov_95 = mean(cov_95), cov_iqr = mean(cov_iqr)) %>%
    ggplot() +
    geom_line(aes(ageref, cov_95, color = nsv, group = nsv)) +
    facet_grid(vars(bias), vars(trend)) +
    theme_bw() +
    theme(legend.position = "bottom", strip.background = element_rect(fill = "grey95", color = "grey")) +
    guides(color = guide_legend(nrow = 1)) +
    scale_fill_viridis_d() +
	scale_y_continuous(labels = scales::percent) +
    labs(
        x = "Referenced age-group", y = "Coverage",
        subtitle = "Percentage of simulation included the true trend difference in 95% CrI", 
		color = "Number of surveys"
    ) -> g_cov
quartz_off(g_cov, "cov_diff", 7, 5, open = 1)

processed %>%
    group_by(ageref, trend, bias, nsv) %>%
    summarise(cov_95 = mean(cov_95), cov_iqr = mean(cov_iqr)) %>%
    ggplot() +
    geom_line(aes(ageref, cov_iqr, color = nsv, group = nsv)) +
    facet_grid(vars(bias), vars(trend)) +
    theme_bw() +
    theme(legend.position = "bottom", strip.background = element_rect(fill = "grey95", color = "grey")) +
    guides(color = guide_legend(nrow = 1)) +
    scale_fill_viridis_d() +
	scale_y_continuous(labels = scales::percent) +
    labs(
        x = "Referenced age-group", y = "Coverage",
        subtitle = "Percentage of simulation included the true trend difference in IQR estimate", 
		color = "Number of surveys"
    ) -> g_cov
quartz_off(g_cov, "cov_diff_iqr", 7, 5, open = 1)

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
