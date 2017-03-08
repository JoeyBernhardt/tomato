## plotting initial yeast growth curves!


# load packages -----------------------------------------------------------

library(tidyverse)
library(janitor)


# read in data ------------------------------------------------------------

yeast_raw <- read_csv("data-raw/bioscreen/16C_test_yeast_march6.csv")


yeast_raw %>% 
	clean_names() %>% 
	filter(time > 905) %>%
	gather(key = "well", value = "OD", starts_with("well")) %>% 
	separate(well, into = c("well_number", "well_id")) %>% 
	select(-well_number) %>% 
	mutate(well_id = as.numeric(well_id)) %>% 
	ggplot(aes(x = time, y = OD)) + geom_point() + 
	geom_line() +
	facet_wrap( ~ well_id, scales = "free_y") + theme_bw() + 
	theme(axis.text.x = element_text(angle = 45, hjust = 1))


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
