## plotting initial yeast growth curves!


# load packages -----------------------------------------------------------

library(tidyverse)
library(janitor)
library(broom)
library(lubridate)
library(plotrix)
library(purrr)
library(stringr)


# read in data ------------------------------------------------------------

yeast_raw <- read_csv("data-raw/bioscreen/30C_march8_2017.csv")
yeast_raw2 <- read_csv("data-raw/bioscreen/30C_march9_2017.csv")
yeast_raw3 <- read_csv("data-raw/bioscreen/30C_march10_2017.csv")
yeast_raw4 <- read_csv("data-raw/bioscreen/30C_march11_2017.csv")
yeast_raw5 <- read_csv("data-raw/bioscreen/30C_march12_2017.csv")
id <- read_csv("data-raw/bioscreen/30C_bioscreen_08062017.csv")



yeast_raw_all <- bind_rows(yeast_raw, yeast_raw2, yeast_raw3, yeast_raw4, yeast_raw5)
yeast_raw_all$sample <- rownames(yeast_raw_all)
# 
# yeast_raw_all %>% 
# 	select(sample, everything()) %>% View
# yeast_raw <- yeast_raw %>% 
# 	select(sample, everything()) %>% 
# 	mutate(sample = as.numeric(sample))
# 	
# yeast_raw$sample <- rownames(yeast_raw)	
# yeast_raw2$sample <- rownames(yeast_raw2)	
# yeast_raw3$sample <- rownames(yeast_raw3)
# yeast_raw4$sample <- rownames(yeast_raw4)	
# yeast_raw5$sample <- rownames(yeast_raw5)	
# 
# 
# 
# 
# yeast_raw2 <- yeast_raw2 %>% 
# 	mutate(sample = as.numeric(sample)) %>% 
# 	mutate(sample = sample + 94) %>% 
# 	select(sample, everything()) %>%
# 	filter(sample != 95)
# 
# 
# yeast_raw3 <- yeast_raw3 %>% 
# 	mutate(sample = as.numeric(sample)) %>% 
# 	mutate(sample = sample + 194) %>% 
# 	select(sample, everything())
# 
# yeast_raw_all <- bind_rows(yeast_raw, yeast_raw2)

yeast <- yeast_raw_all %>% 
	clean_names() %>% 
	gather(key = "well", value = "OD", starts_with("well")) %>% 
	separate(well, into = c("well_number", "well_id")) %>% 
	select(-well_number) %>% 
	mutate(well_id = as.numeric(well_id))

yeast2 <- left_join(yeast, id)

yeast3 <- yeast2 %>% 
	mutate(sample = as.numeric(sample)) %>% 
	mutate(time_since_innoc_days = (((sample *15)/60/24))) %>% 
	select(time_since_innoc_days, everything()) %>% 
	filter(sample != "95") %>% 
	filter(sample != "1") %>% 
	filter(sample != "195") %>% 
	filter(sample != "287") %>% 
	filter(sample != "362")
 
yeast3 %>% 
	# filter(time > 10) %>% 
	filter(well_id < 251, well_id > 221) %>% 
	filter(OD>0.12) %>% 
	filter(!grepl("blank", treatment)) %>% 
	ggplot(aes(x = sample, y = OD, color = treatment)) + geom_point() + 
	# geom_line() +
	facet_wrap( ~ well_id) +
	theme_bw() + 
	theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("figures/30C_ypg_growth_day1.png")


str(yeast2)
yeast2$start.time <- hms("00:00:00")
yeast2$time_since_innoc <- interval(yeast2$start.time, yeast2$time)
glimpse(yeast2)

# yeast3 <- yeast2 %>% 
# 	mutate(time_since_innoc_days = time/ddays(1)) %>% 
# 	mutate(time_since_innoc_hours = time/dhours(1))



### estimate r
yeast3 %>% 
	group_by(treatment, well_id) %>% 
	do(tidy(nls(OD ~ 0.15 * (1+a)^(time_since_innoc),
							data= .,  start=list(a=0.01),
							control = nls.control(maxiter=100, minFactor=1/204800000)))) %>% 
	ungroup() %>%
	group_by(treatment) %>% 
	summarise_each(funs(mean, std.error), estimate) %>%
	ggplot(aes(x = treatment, y = mean)) + geom_point(size = 4) + 
	geom_errorbar(aes(ymin = mean - std.error, ymax = mean + std.error), width = 0.2) + 
	theme_bw() + xlab("treatment") + ylab("intrinsic growth rate (r)") + 
	theme(text = element_text(size=18))

### log OD vs time
yeast3 %>% 
	filter(time_since_innoc_days > 0.5, time_since_innoc_days < 2) %>% 
	mutate(hours = time_since_innoc_days*24) %>% 
	# filter(!grepl("blank", treatment)) %>% 
	# filter(well_id == 131) %>%
	# ggplot(aes(x = time_since_innoc, y = log(OD), color = treatment_b)) + geom_point() +
	# facet_wrap( ~ well_setup) + geom_smooth(method = "lm")
	group_by(treatment, well_id) %>% 
	do(tidy(lm(log(OD) ~ hours, data = .), conf.int = TRUE)) %>% 
	ungroup() %>%
	filter(term != "(Intercept)") %>%
	filter(!grepl("blank", treatment)) %>% 
	group_by(treatment) %>% 
	summarise_each(funs(mean, std.error), estimate) %>%
	ggplot(aes(x = treatment, y = mean)) + geom_point(size = 4) + 
	geom_errorbar(aes(ymin = mean - std.error, ymax = mean + std.error), width = 0.2) + 
	theme_bw() + xlab("treatment") + ylab("intrinsic growth rate (r)") + 
	theme(text = element_text(size=18))

	yeast3 %>% 
		# filter(well_setup == 54) %>% 
		ggplot(aes(x = time_since_innoc, y = OD, color = treatment)) + geom_point() +
		facet_wrap( ~ treatment)


yeast %>% 
	mutate(treatment = ifelse(well_id %in% c("1", "10", "91", "100"), "blank", "yeast")) %>%
	# filter(well_id %in% c(1:10)) %>% 
	ggplot(aes(x = time, y = OD, color = treatment)) + geom_point() + 
	geom_line() +
	facet_wrap( ~ well_id, scales = "free_y") + theme_bw() + 
	theme(axis.text.x = element_blank(), axis.text.y = element_blank()) + 
	theme(strip.background = element_blank(),
				strip.text.x = element_blank()) + 
	theme(text = element_text(size = 18))

ggsave("figures/16C_growth_all.png", width = 20, height = 15)


### trying Nathaniel's code

nderiv<-function(fit, x, eps=1e-5){(predict(fit, x + eps) - predict(fit, x - eps))/(2 * eps)}

spline.slope.single<-function(x, y,  n=101, eps=1e-5, span=0.5){
	max(nderiv(loess(log(y) ~ x, degree=1, span=span), seq(min(df$x), max(df$x), length=n)), na.rm=TRUE)
}

spline.slope<-function(df, n=101, eps=1e-5, span=0.2){
	x <- df$hours
	y <- df$OD
	max(nderiv(loess(log(y) ~ x, degree=1, span=span), seq(min(x), max(x), length=n)), na.rm=TRUE)
}


str(yeast3)


snippet <- yeast3 %>% 
	# filter(time_since_innoc > 1, time_since_innoc < 2) %>% 
	# filter(!grepl("blank", treatment)) %>% 
	filter(well_id == 131) %>% 
	mutate(hours = time_since_innoc) %>% 
	select(hours, OD) 


spline.slope.single(snippet$hours, snippet$OD)




### now for the whole dataset


yeast_split <- yeast3 %>% 
	mutate(hours = time_since_innoc_days*24) %>% 
	filter(time_since_innoc_days > 0.5, time_since_innoc_days < 2) %>% 
	filter(!grepl("blank", treatment)) %>% 
	select(hours, OD, well_id) %>% 
	split(.$well_id)

str(yeast_split)

growth_rates <- yeast_split %>% 
	map_df(spline.slope) %>% 
	gather(key = "well_id", value = "growth_rate", 1:189) %>% 
	mutate(well_id = as.integer(well_id)) %>% 
	mutate(temperature = 30)


yeast4 <- left_join(growth_rates, id)

str(yeast4)
yeast4 %>% 
	group_by(treatment) %>% 
	summarise_each(funs(mean, std.error), growth_rate) %>%
	ggplot(aes(x = treatment, y = mean)) + geom_point(size = 4) + 
	geom_errorbar(aes(ymin = mean - std.error, ymax = mean + std.error), width = 0.2) + 
	theme_bw() + xlab("treatment") + ylab("growth rate") + 
	theme(text = element_text(size=18))


# 16C on YAPD -------------------------------------------------------------

march6 <- read_csv("data-raw/bioscreen/16C_test_yeast_march6.csv")
march7 <- read_csv("data-raw/bioscreen/16C_test_yeast_march7.csv")
march8 <- read_csv("data-raw/bioscreen/16C_test_yeast_march8.csv")


march6$sample_time <- rownames(march6)
march7$sample_time <- rownames(march7)
march8$sample_time <- rownames(march8)


march6b <- march6 %>%
	mutate(sample_time = as.numeric(sample_time)) 

march7b <- march7 %>%
	mutate(sample_time = as.numeric(sample_time)) %>% 
	mutate(sample_time = sample_time + 88) %>% 
	select(sample_time, everything())

march8b <- march8 %>% 
	mutate(sample_time = as.numeric(sample_time)) %>% 
	mutate(sample_time = sample_time + 186) %>% 
	select(sample_time, everything())

yeast_raw_all_16 <- bind_rows(march6b, march7b, march8b)

yeast16 <- yeast_raw_all_16 %>% 
	clean_names() %>% 
	gather(key = "well", value = "OD", starts_with("well")) %>% 
	separate(well, into = c("well_number", "well_id")) %>% 
	select(-well_number) %>% 
	mutate(well_id = as.numeric(well_id))

yeast16 %>% 
	mutate(sample_time = (sample_time*15)/60) %>% 
	mutate(well_id = as.numeric(well_id) - 100) %>% 
	mutate(treatment = ifelse(well_id %in% c("1", "10", "91", "100"), "blank", "yeast")) %>%
	filter(treatment == "yeast") %>% 
	ggplot(aes(x = sample_time, y = OD, color = treatment)) + geom_point() +
	# facet_wrap( ~ well_id, scales = "free_y") +
	theme_bw() + 
	theme(axis.text.x = element_blank(), axis.text.y = element_blank()) + 
	theme(strip.background = element_blank(),
				strip.text.x = element_blank()) + 
	theme(text = element_text(size = 18))

### now with Nathaniel's code

yeast16 <- yeast_raw_all_16 %>% 
	clean_names() %>% 
	gather(key = "well", value = "OD", starts_with("well")) %>% 
	separate(well, into = c("well_number", "well_id")) %>% 
	select(-well_number) %>% 
	mutate(well_id = as.numeric(well_id)) %>% 
	mutate(treatment = ifelse(well_id %in% c("101", "110", "191", "200"), "blank", "yeast")) %>%
	filter(treatment == "yeast")

length(unique(yeast16$well_id))
yeast_split_16 <- yeast16 %>% 
	mutate(sample_time = as.numeric(sample_time)) %>% 
	mutate(time_since_innoc_days = (((sample_time *15)/60/24))) %>% 
	mutate(hours = time_since_innoc_days*24) %>% 
	# filter(time_since_innoc_days > 0.5, time_since_innoc_days < 2) %>% 
	select(hours, OD, well_id) %>% 
	split(.$well_id)

str(yeast_split_16)

growth_rates_16 <- yeast_split_16 %>% 
	map_df(spline.slope) %>% 
	gather(key = "well_id", value = "growth_rate", 1:96) %>% 
	mutate(well_id = as.integer(well_id)) %>% 
	mutate(temperature = "16") %>% 
	mutate(temperature = as.numeric(temperature)) %>% 
	mutate(treatment = "16C YAPD")


## now join all the growth rates
str(growth_rates_16)
growth_rates_all <- bind_rows(growth_rates_16, yeast4)


growth_rates_all %>% 
	mutate(treatment = str_replace(treatment, "glyc_1", "30C, YPG, 1% glycerol")) %>% 
	mutate(treatment = str_replace(treatment, "glyc_3", "30C, YPG, 3% glycerol")) %>% 
	mutate(treatment = str_replace(treatment, "YAPD", "16C YAPD")) %>% 
	group_by(treatment) %>% 
	summarise_each(funs(mean, std.error), growth_rate) %>%
	ggplot(aes(x = treatment, y = mean)) + geom_point(size = 4) + 
	geom_errorbar(aes(ymin = mean - std.error, ymax = mean + std.error), width = 0.2) + 
	theme_bw() + xlab("") + ylab("max growth rate (/hr)") + 
	theme(text = element_text(size=18)) +
	theme(axis.text.x = element_text(angle = 45, hjust = 1))




### now try with Nathaniel's example data

## read in data
example_data_raw <- read_csv("data-raw/bioscreen/example_bioscreen_30C_YPAD.csv")

example_data_raw$sample_time <- rownames(example_data_raw) ## cheating a bit here with the times, to better deal with importing multiple csvs from the same bioscreen run

## clean up data
example_data <- example_data_raw %>% 
	clean_names() %>% 
	gather(key = "well", value = "OD", starts_with("well")) %>% 
	separate(well, into = c("well_number", "well_id")) %>% 
	select(-well_number) %>% 
	mutate(well_id = as.numeric(well_id)) %>% 
	mutate(sample_time = as.numeric(sample_time)) %>% 
	mutate(time_since_innoc_days = (((sample_time *15)/60/24))) %>% 
	mutate(hours = time_since_innoc_days*24) %>% 
	filter(hours < 21, hours > 0.5) %>% 
	select(hours, OD, well_id) 

## split example df by well_id
example_split <-  example_data %>% 
	split(.$well_id)


### define Nathaniel's fitting function
nderiv<-function(fit, x, eps=1e-5){(predict(fit, x + eps) - predict(fit, x - eps))/(2 * eps)}

spline.slope <- function(df, n=101, eps=1e-5, span=0.5){
	x <- df$hours
	y <- df$OD
	time_point <- seq(min(x), max(x), length=n)
	growth <- nderiv(loess(log(y) ~ x, degree=1, span=span), time_point)
	output <- top_n(data.frame(growth, time_point), n = 1, wt = growth)
return(output)
}


## fit each well
growth_rates_example <- example_split %>% 
	map_df(spline.slope, .id = "well_id")

## join the growth rate results back to initial df
example_results <- left_join(example_data, growth_rates_example)

## plot it, to visualize where the max growth rate was found along the time series
example_results %>% 
	filter(well_id < 131) %>% 
	mutate(time = ifelse(hours >= (time_point -0.5) & hours <= (time_point + 0.5), "max_growth_range", "time")) %>%
	ggplot(aes(x = hours, y = OD, color = time)) + geom_point(size = 0.5) +
	facet_wrap( ~ well_id) + theme_bw()



### now trying with my 30C data
yeast_trim <- yeast3 %>% 
	mutate(hours = time_since_innoc_days*24) %>% 
	filter(time_since_innoc_days > 0.05) %>% 
	filter(!grepl("blank", treatment))

yeast_split <- yeast_trim %>% 
	select(hours, OD, well_id) %>% 
	split(.$well_id)

str(yeast_split)

growth_rates <- yeast_split %>% 
	map_df(spline.slope, .id = "well_id") %>% 
	mutate(well_id = as.integer(well_id))


yeast4 <- left_join(growth_rates, id)

yeast_30_all <- left_join(yeast4, yeast_trim, by = "well_id")

yeast_30_all %>% 
	filter(well_id < 132) %>% 
	mutate(time = ifelse(hours >= (time_point -1) & hours <= (time_point + 1), "max_growth_range", "time")) %>%
	ggplot(aes(x = hours, y = OD, color = time)) + geom_point(size = 0.5) +
	facet_wrap( ~ well_id) + theme_bw()

yeast_30_all %>% 
	filter(time_point < 50) %>% 
	group_by(treatment.x) %>% 
	summarise_each(funs(mean, std.error), time_point) %>% 
	ggplot(aes(x = treatment.x, y = mean)) + geom_point(size = 4) + 
	geom_errorbar(aes(ymin = mean - std.error, ymax = mean + std.error), width = 0.2)


yeast_30_all %>% 
	filter(time_point < 50) %>% 
	group_by(treatment.x) %>% 
	summarise_each(funs(mean, std.error), growth) %>%
	ggplot(aes(x = treatment.x, y = mean)) + geom_point(size = 4) + 
	geom_errorbar(aes(ymin = mean - std.error, ymax = mean + std.error), width = 0.2) + 
	theme_bw() + xlab("") + ylab("max growth rate (/hr)") + 
	theme(text = element_text(size=18)) +
	theme(axis.text.x = element_text(angle = 45, hjust = 1))
	


growth_rates_all <- bind_rows(growth_rates_16, yeast4, growth_rates_example)


growth_rates_all %>% 
	mutate(treatment = str_replace(treatment, "glyc_1", "30C, YPG, 1% glycerol")) %>% 
	mutate(treatment = str_replace(treatment, "glyc_3", "30C, YPG, 3% glycerol")) %>% 
	group_by(treatment) %>% 
	summarise_each(funs(mean, std.error), growth_rate) %>%
	ggplot(aes(x = treatment, y = mean)) + geom_point(size = 4) + 
	geom_errorbar(aes(ymin = mean - std.error, ymax = mean + std.error), width = 0.2) + 
	theme_bw() + xlab("") + ylab("max growth rate (/hr)") + 
	theme(text = element_text(size=18)) +
	theme(axis.text.x = element_text(angle = 45, hjust = 1))


### plot them all together


example_trim <- example_data %>% 
	mutate(sample_time = as.numeric(sample_time)) %>% 
	mutate(time_since_innoc_days = (((sample_time *15)/60/24))) %>% 
	mutate(hours = time_since_innoc_days*24) %>% 
	filter(hours < 21, hours > 0.5) %>% 
	select(hours, OD, well_id) %>% 
	mutate(treatment = "30C YPAD")

trim_16 <-  yeast16 %>% 
	mutate(sample_time = as.numeric(sample_time)) %>% 
	mutate(time_since_innoc_days = (((sample_time *15)/60/24))) %>% 
	mutate(hours = time_since_innoc_days*24) %>% 
	filter(hours > 0.5) %>% 
	select(hours, OD, well_id) %>% 
	mutate(treatment = "16C YPAD")
	

all_YPAD <- bind_rows(example_trim, trim_16)


all_YPAD %>% 
	ggplot(aes(x = hours, y = OD, color = treatment, group = well_id)) + geom_point(size = 0.25) + 
	theme_bw()
