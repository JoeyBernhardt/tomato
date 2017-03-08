### looking at bioscreen temperature data


# load packages -----------------------------------------------------------
library(tidyverse)
library(lubridate)
library(janitor)


# read in data ------------------------------------------------------------

ibutton12_6 <- read_csv("data-raw/12-6.csv")
ibutton13_5 <- read_csv("data-raw/13-5.csv")
ibutton27_a <- read_csv("data-raw/27-a.csv")
ibutton8_1 <- read_csv("data-raw/8-1.csv")

ibutton12_6  <- clean_names(ibutton12_6)
ibutton13_5  <- clean_names(ibutton13_5)
ibutton27_a  <- clean_names(ibutton27_a) %>% 
	mutate(time = hour(date_time)) %>% View
ibutton8_1  <- clean_names(ibutton8_1)

# separate out temperature chunks -----------------------------------------
?lubridate

ibutton27_a %>% 
	filter(date_time <-  )

