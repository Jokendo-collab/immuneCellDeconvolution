
library(tidyverse)      # data manipulation and visualization
library(lubridate)      # easily work with dates and times
library(fpp2)           # working with time series data
library(zoo)

require(smooth)
require(Mcomp)

smoothScatter(dm2$EGFR, h=18, silent=FALSE)
#Load the data

dm <- read.table("C:/Users/Javan_Okendo/Desktop/Kagisho_Data/analysis_data.txt",header = T,sep = '\t')
dim(dm)

dm2=na.omit(dm) #remove the mising values


head(dm2) #view the first few rows

dim(dm2) #get to know the dataframe dimensions

savings <- dm2 %>%
  select(dm2$Year_since_started, srate = EGFR) %>%
  mutate(srate_ma01 = rollmean(srate, k >= 90, fill = NA),
         srate_ma02 = rollmean(srate, k < 60, fill = NA),
         srate_ma03 = rollmean(srate, k < 89, fill = NA))
