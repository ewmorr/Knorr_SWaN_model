###############################################################################
###############################################################################
###############################################################################
#This script runs empirical model of temperature-respiration, and then uses model coefficients to
#predict fluxes using hourly measurements of soil temperature

#code prepared by Alix Contosta 06/21/2018
#updated 01/12/2020 to include years 2018 and 2019
#updated 01/26/2021 to include year 2020
#updated 4/21/2023 to include year 2021 and 2022

#updated 03/19/2021 to do the following:
#1. filter out temperatures from warmed plots that do 
#not have a delta of > 2.5 or < 7.5
#2. model respiration as a function of temperature
#by treatment instead of using treatment as a model
#covariate. 
#3. estimate annual fluxes using temperature data that may contain
#missing values and scale to 365 days 
#4. calculate cumulative differences from control 
#to compare with changes in soil carbon

#updated 04/21/2023 by EM for the following
#1. Split the loops that model respiration 
# (i.e., the field respiration measurements)
# and the loops that predict respiration based on
# hourly temperature values. These are now separate loops.
#2. As part of the update the indexing was switched from
# vector-based to dplyr filtering.
#3. Model outputs are now saved in a list and output as an RDS.
#4. We also removed N plots 11,16, and 23 from the calculation of statistics, 
# plot generation, and the output model predictions. These plots are not monitored for
# temperature and we previously took an average of the 3 monitored plots
# but this artificially increases the number of measurements. These plots are
# included in the modeling.

###############################################################################
#There are four stages to the process

#Stage 1: data set up: load libraries, read data, and do initial data processing to prepare for modeling and annual flux estimation. The temperature data are quality controlled by removing any values in the heated plots where the difference between the heated plot and it's reference control plot is less than 2.5C or greater than 7.5C. These values indicate an error with the heating system or the data logger system. Individual plots with missing values are then back-filled by taking an average of the other plots in that treatment for individual timepoints. Larger blocks of missing data where all plots have missing values, indicating a system shut down, are estimated by linear interpolation within plots.

#Stage 2: run model of temperature-respiration

#Stage 3: predict respiration based on soil temperature, treatment, and year

#Stage 4: run statistical analyses

###############################################################################
###############################################################################
###############################################################################

#####Stage 1: data set up: load libraries, read data, and do initial data processing to prepare for modeling and annual flux estimation

#####load libraries
library(data.table)
library(dplyr)
library(caTools)
library(multcomp)
library(nlme)
library(MASS)
library(reshape2)
library(stringr)
library(zoo)
#library(ggplot2)

#####set working directory and import data

#setwd("C:\\Users\\Alix\\Documents\\GitHub\\swan")
#setwd("~/Documents/R/swan_2020")
#setwd("~/repo/swan")

resp = read.table("data/SWN_annualflux_2006-22_cc.csv", head = TRUE, sep = ",")
temps = read.table("data/daily_avg_2006_2022.csv", head = TRUE, sep = ",") #produced using avg_daily_temps.R
#Note that the daily temp file is not used for the final predictions and statistics (in favor of the estimates based on the hourly data)
#However, it is handy for resting some of the longer running parts of the script (e.g., QA/QC steps)
hourly = read.table("data/cleaned_hourly_swntemp_2006_2022.csv", head = TRUE, sep = ",")
plot_pairs = read.table("data/heating_pairs.txt", header = T)

###################################################
#####initial processing for field respiration data
###################################################

#force doy to be day bc it's not truly an ordinal date
resp$day = resp$doy

#create date object
resp$Date = paste(resp$month, resp$day, resp$year, sep = "/")

#make date posix compliant
resp$DATETIME = strptime(resp$Date, "%m/%d/%Y",tz="EST")

#extract ordinal day from DATETIME
resp$doy = as.numeric(strftime(resp$DATETIME, "%j"))

#make Year_doy
resp$ydoy = as.integer(ifelse(resp$doy < 10, paste(resp$year, 0, 0, resp$doy, sep = ""),
            ifelse(resp$doy < 100, paste(resp$year, 0, resp$doy, sep = ""),
            paste(resp$year, resp$doy, sep = ""))))

#log transform daily fluxes for fitting exponential models
resp$LnDFlux = log(resp$dflux)
resp$LnFlux = log(resp$corrflux)

#End field respiration processing
#################################

##################################################
#####initial processing for soil temperature data
##################################################

#use Sys.setenv to set system environmental variables, in this case, GMT time
Sys.setenv(TZ='GMT')

#add columns for plots 11, 16, and 23 by averaging over other N temperatures from plots 9 and 11

temps$N_11 = (temps$N_9 + temps$N_21) / 2
temps$N_16 = (temps$N_9 + temps$N_21) / 2
temps$N_23 = (temps$N_9 + temps$N_21) / 2

hourly$N_11 = (hourly$N_6 + hourly$N_9 + hourly$N_21) / 3
hourly$N_16 = (hourly$N_6 + hourly$N_9 + hourly$N_21) / 3
hourly$N_23 = (hourly$N_6 + hourly$N_9 + hourly$N_21) / 3

#create new columns for Day and Month that add zeroes onto days or months that are <10 (important for merging tables in a later step)
temps$Day.1 = ifelse(temps$Day < 10, paste(0, temps$Day, sep = ""), temps$Day)
temps$Month.1 = ifelse(temps$Month < 10, paste(0, temps$Month, sep = ""), temps$Month)

hourly$Day.1 = ifelse(hourly$DAY < 10, paste(0, hourly$DAY, sep = ""), hourly$DAY)
hourly$Month.1 = ifelse(hourly$MONTH < 10, paste(0, hourly$MONTH, sep = ""), hourly$MONTH)
hourly$TIME.1 = ifelse(hourly$TIME < 1000, paste(0, hourly$TIME, sep = ""), hourly$TIME)

sst <- data.frame(do.call(rbind, stringr::str_split(hourly$TIME.1, "")))
sst$hh = as.numeric(paste(sst$X1, sst$X2, sep = ""))

#Set hour 24 to hour 0 for POSIX compliance
hourly$hour = ifelse(sst$hh == 24, 0, sst$hh)

#make a Date column for merging in a later step
temps$Date = as.factor(paste(temps$Year, temps$Month.1, temps$Day.1, sep = "-"))
hourly$Date = as.factor(paste(hourly$YEAR, hourly$Month.1, hourly$Day.1, hourly$hour, sep = "-"))

#create a vector of all the possible days in the time series from Jan 1, 2007 to Dec 31, 2020 including leap years
doys = as.character(seq(from = as.Date("2007-01-01"), to = as.Date("2022-12-31"), by = 1, tz = "UTC"))
doyhs = seq(from = as.POSIXct("2007-01-01 0:00:00"), to = as.POSIXct("2022-12-31 23:00:00"), by = "hour", tz = "GMT")

#make doys into a data.frame
doys = data.frame(doys)
doyhs = data.frame(doyhs)

#make date posix compliant
doys$DATETIME = strptime(doys$doys, "%Y-%m-%d", tz="UTC")

#extract year, month, and day from DATETIME
doys$Year = as.numeric(strftime(doys$DATETIME, "%Y"))
doys$Month = as.numeric(strftime(doys$DATETIME, "%m"))
doys$Day = as.numeric(strftime(doys$DATETIME, "%d"))
doys$doy = as.numeric(strftime(doys$DATETIME, "%j"))

doyhs$Year = as.numeric(strftime(doyhs$doyhs, tz = "UTC", "%Y"))
doyhs$Month = as.numeric(strftime(doyhs$doyhs, tz = "UTC", "%m"))
doyhs$Day = as.numeric(strftime(doyhs$doyhs, tz = "UTC", "%d"))
doyhs$doy = as.numeric(strftime(doyhs$doyhs, tz = "UTC", "%j"))
doyhs$hour = as.numeric(strftime(doyhs$doyhs, tz = "UTC", "%H"))

doyhs$Day.1 = ifelse(doyhs$Day < 10, paste(0, doyhs$Day, sep = ""), doyhs$Day)
doyhs$Month.1 = ifelse(doyhs$Month < 10, paste(0, doyhs$Month, sep = ""), doyhs$Month)

doyhs$Date = as.factor(paste(doyhs$Year, doyhs$Month.1, doyhs$Day.1, doyhs$hour, sep = "-"))

#merge doys wth temps to enable a complete time series (there are some days totally missing from the temps data set)
temps.1 = merge(doys, temps, by.x = "doys", by.y = "Date", all.x = T, all.y = F)
hourly.1 = merge(doyhs, hourly, by.x = "Date", by.y = "Date", all.x = T, all.y = F)

#extract cols 11:34 and 14:37(the soil temperature columns in temps.1 and hourly.1)
plots = temps.1[ , 11:33]
hplots = hourly.1[ , 14:37]

##################################################
#####temperature data QA

#create empty dataframes for holding the results of heated plots QA
H_plots.QA = data.frame(matrix(nrow = nrow(plots), ncol = 11))
H_hplots.QA = data.frame(matrix(nrow = nrow(hplots), ncol = 11))

#calculate delta between heated and control plots.
#omit heated plot temperatures where delta is < 2.5 or > 7.5

#for daily average temperatures
for(i in 1:nrow(plot_pairs)) {
  C_plot = plots[ , as.character(plot_pairs$Control[i])]
  H_plot = plots[ , as.character(plot_pairs$Heated[i])]

  H_plot.QA = ifelse(H_plot - C_plot < 2.5 | H_plot - C_plot > 7.5,
                     NA, H_plot)
  
  H_plots.QA[[i]] = H_plot.QA 
  names(H_plots.QA)[i] = as.character(plot_pairs$Heated[i])
}

#These indices must correctly match the column numbers of the control and N plots
plots.QA = cbind(plots[ , c(1:6, 7:9, 21:23)], H_plots.QA)

#hourly temperatures
for(i in 1:nrow(plot_pairs)) {
  C_hplot = hplots[ , as.character(plot_pairs$Control[i])]
  H_hplot = hplots[ , as.character(plot_pairs$Heated[i])]
  
  H_hplot.QA = ifelse(H_hplot - C_hplot < 2.5 | H_hplot - C_hplot > 7.5,
                     NA, H_hplot)
  
  H_hplots.QA[[i]] = H_hplot.QA 
  names(H_hplots.QA)[i] = as.character(plot_pairs$Heated[i])
}

#These indices must correctly match the column numbers of the control and N plots
hplots.QA = cbind(hplots[ , c(1:6, 7:9, 22:24)], H_hplots.QA)

####################################
#After removing temp values based on the QA we back fill values for plots with missing temps where possible, by using the average temp value in a given treatment. This fills individual time points where some of the plots had acceptable temperature values.

C_plots = c("C_1", "C_12", "C_14", "C_19", "C_20", "C_24")
N_plots = c("N_6", "N_9", "N_21", "N_11", "N_16", "N_23")

H_plots = c("H_2", "H_5", "H_8", "H_13", "H_18") #one heated plot is excluded because the heating system is off
HN_plots = c("HN_4", "HN_7", "HN_10", "HN_3", "HN_17", "HN_15")

#daily data gap filling based on remaining treatment values

plots.QA.C = plots.QA[,C_plots]
plots.QA.N = plots.QA[,N_plots]
plots.QA.H = plots.QA[,H_plots]
plots.QA.HN = plots.QA[,HN_plots]


for(i in 1:nrow(plots.QA.C)){
    if(sum( is.na(plots.QA.C[i,]) ) > 0 & sum( is.na(plots.QA.C[i,]) ) < ncol(plots.QA.C) ){ #if the number of NAs is gt 0 but lt the number of columns
        plots.QA.C[i, is.na( plots.QA.C[i,] ) ] = sum(plots.QA.C[i, !is.na( plots.QA.C[i,] ) ]) / length(plots.QA.C[i, !is.na( plots.QA.C[i,] ) ])
    }else{
        plots.QA.C[i,] = plots.QA.C[i,]
    }
}

for(i in 1:nrow(plots.QA.N)){
    if(sum( is.na(plots.QA.N[i,]) ) > 0 & sum( is.na(plots.QA.N[i,]) ) < ncol(plots.QA.N) ){
        plots.QA.N[i, is.na( plots.QA.N[i,] ) ] = sum(plots.QA.N[i, !is.na( plots.QA.N[i,] ) ]) / length(plots.QA.N[i, !is.na( plots.QA.N[i,] ) ])
    }else{
        plots.QA.N[i,] = plots.QA.N[i,]
    }

}

for(i in 1:nrow(plots.QA.H)){
    if(sum( is.na(plots.QA.H[i,]) ) > 0 & sum( is.na(plots.QA.H[i,]) ) < ncol(plots.QA.H) ){
        plots.QA.H[i, is.na( plots.QA.H[i,] ) ] = sum(plots.QA.H[i, !is.na( plots.QA.H[i,] ) ]) / length(plots.QA.H[i, !is.na( plots.QA.H[i,] ) ])
    }else{
        plots.QA.H[i,] = plots.QA.H[i,]
    }

}

for(i in 1:nrow(plots.QA.HN)){
    if(sum( is.na(plots.QA.HN[i,]) ) > 0 & sum( is.na(plots.QA.HN[i,]) ) < ncol(plots.QA.HN) ){
        plots.QA.HN[i, is.na( plots.QA.HN[i,] ) ] = sum(plots.QA.HN[i, !is.na( plots.QA.HN[i,] ) ]) / length(plots.QA.HN[i, !is.na( plots.QA.HN[i,] ) ])
    }else{
        plots.QA.HN[i,] = plots.QA.HN[i,]
    }

}

plots.QA = cbind(plots.QA.C, plots.QA.N, plots.QA.H, plots.QA.HN)

########################
#hourly data gap filling
hplots.QA.C = hplots.QA[,C_plots]
hplots.QA.N = hplots.QA[,N_plots]
hplots.QA.H = hplots.QA[,H_plots]
hplots.QA.HN = hplots.QA[,HN_plots]

for(i in 1:nrow(hplots.QA.C)){
    if(sum( is.na(hplots.QA.C[i,]) ) > 0 & sum( is.na(hplots.QA.C[i,]) ) < ncol(hplots.QA.C) ){ #if the number of NAs is gt 0 but lt the number of columns
        hplots.QA.C[i, is.na( hplots.QA.C[i,] ) ] = sum(hplots.QA.C[i, !is.na( hplots.QA.C[i,] ) ]) / length(hplots.QA.C[i, !is.na( hplots.QA.C[i,] ) ])
    }else{
        hplots.QA.C[i,] = hplots.QA.C[i,]
    }
}

for(i in 1:nrow(hplots.QA.N)){
    if(sum( is.na(hplots.QA.N[i,]) ) > 0 & sum( is.na(hplots.QA.N[i,]) ) < ncol(hplots.QA.N) ){
        hplots.QA.N[i, is.na( hplots.QA.N[i,] ) ] = sum(hplots.QA.N[i, !is.na( hplots.QA.N[i,] ) ]) / length(hplots.QA.N[i, !is.na( hplots.QA.N[i,] ) ])
    }else{
        hplots.QA.N[i,] = hplots.QA.N[i,]
    }
}

for(i in 1:nrow(hplots.QA.H)){
    if(sum( is.na(hplots.QA.H[i,]) ) > 0 & sum( is.na(hplots.QA.H[i,]) ) < ncol(hplots.QA.H) ){
        hplots.QA.H[i, is.na( hplots.QA.H[i,] ) ] = sum(hplots.QA.H[i, !is.na( hplots.QA.H[i,] ) ]) / length(hplots.QA.H[i, !is.na( hplots.QA.H[i,] ) ])
    }else{
        hplots.QA.H[i,] = hplots.QA.H[i,]
    }
}

for(i in 1:nrow(hplots.QA.HN)){
    if(sum( is.na(hplots.QA.HN[i,]) ) > 0 & sum( is.na(hplots.QA.HN[i,]) ) < ncol(hplots.QA.HN) ){
        hplots.QA.HN[i, is.na( hplots.QA.HN[i,] ) ] = sum(hplots.QA.HN[i, !is.na( hplots.QA.HN[i,] ) ]) / length(hplots.QA.HN[i, !is.na( hplots.QA.HN[i,] ) ])
    }else{
        hplots.QA.HN[i,] = hplots.QA.HN[i,]
    }
}

hplots.QA = cbind(hplots.QA.C, hplots.QA.N, hplots.QA.H, hplots.QA.HN)

#Data formatting
#transform the plots file so that all of the soil temperature columns become rows
plots.1 <- reshape2::melt(plots.QA)
names(plots.1) <- c("Plot_Trt", "temp.1")

hplots.1 <- reshape2::melt(hplots.QA)
names(hplots.1) <- c("Plot_Trt", "temp.1")

#add year, Month, Day, and doy to plots.1 and year, Month, Day, doy, and hour to hplots.1
#note that the date/time cols are wrapped along all of the data which is now stacked by plot

plots.1$year = temps.1$Year.x
plots.1$Month = temps.1$Month.x
plots.1$Day = temps.1$Day.x
plots.1$doy = temps.1$doy

hplots.1$year = hourly.1$Year
hplots.1$Month = hourly.1$Month
hplots.1$Day = hourly.1$Day
hplots.1$doy = hourly.1$doy
hplots.1$hour = hourly.1$hour.x

#remove temperature data from 6/24/2013 1900 to 7/9/2013 1000 and 
#from 5/30/2016 0500 to 5/30/2016 1300 due to system malfunction
plots.1$temp.1 = ifelse(plots.1$year == 2013 & plots.1$Month == 6 & plots.1$Day >= 24, NA,
                         ifelse(plots.1$year == 2013 & plots.1$Month == 7 & plots.1$Day <= 9, NA,
                                ifelse(plots.1$year == 2016 & plots.1$Month == 5 & plots.1$Day == 30, NA,
                                       plots.1$temp.1)))

hplots.1$temp.1 = ifelse(hplots.1$year == 2013 & hplots.1$Month == 6 & hplots.1$Day >= 24, NA,
                  ifelse(hplots.1$year == 2013 & hplots.1$Month == 7 & hplots.1$Day <= 9, NA,
                  ifelse(hplots.1$year == 2016 & hplots.1$Month == 5 & hplots.1$Day == 30, NA,
                  hplots.1$temp.1)))

#End QA and gap filling
########################

#split the Plot_Trt column into Plot and Trt
sp = as.character(plots.1$Plot_Trt)
sps = read.table(textConnection(sp), sep = "_")
names(sps) <- c("Trt", "Plot")

#add Plot and Trt onto plots.1 data.frame, making Plot an integer (important for using the plots.1 data.frame when predicting fluxes)
plots.1$Plot = sps$Plot
plots.1$Trt = sps$Trt

#split the Plot_Trt column into Plot and Trt
sp = as.character(hplots.1$Plot_Trt)
sps = read.table(textConnection(sp), sep = "_")
names(sps) <- c("Trt", "Plot")

#add Plot and Trt onto hplots.1 data.frame, making Plot an integer (important for using the hplots.1 data.frame when predicting fluxes)
hplots.1$Plot = sps$Plot
hplots.1$Trt = sps$Trt

######################################
######################################
#Temp interpolation block
#gap fill missing data with linear interpolation

#make sure data are in plot by date order
plots.1 = plots.1[order(plots.1$Plot, plots.1$year, plots.1$Month, plots.1$Day), ]
hplots.1 = hplots.1[order(hplots.1$Plot, hplots.1$year, hplots.1$Month, hplots.1$Day, hplots.1$hour), ]

#use na.approx inside mutate to interpolate within plots
plots.2 = plots.1 %>%
group_by(Plot) %>%
mutate(temp =  na.approx(temp.1, na.rm=FALSE))

hplots.2 = hplots.1 %>%
group_by(Plot) %>%
mutate(temp =  na.approx(temp.1, na.rm=FALSE))

#The interpolated column is "temp", not "temp.1"
plots.1 = plots.2
hplots.1 = hplots.2
#End temp interpolation
###############################################

#make ydoy for downstream processing
plots.1$ydoy = as.integer(ifelse(plots.1$doy < 10, paste(plots.1$year, 0, 0, plots.1$doy, sep = ""),
            ifelse(plots.1$doy < 100, paste(plots.1$year, 0, plots.1$doy, sep = ""),
            paste(plots.1$year, plots.1$doy, sep = ""))))

hplots.1$ydoy = as.integer(ifelse(hplots.1$doy < 10, paste(hplots.1$year, 0, 0, hplots.1$doy, sep = ""),
            ifelse(hplots.1$doy < 100, paste(hplots.1$year, 0, hplots.1$doy, sep = ""),
            paste(hplots.1$year, hplots.1$doy, sep = ""))))

#remove leap years to facilitate iterative model runs
plots.1 = plots.1[plots.1$doy != 366, ]
hplots.1 = hplots.1[hplots.1$doy != 366, ]

#remove plot 22 (has not heated since 2009) to avoid poor model fits
plots.1 = plots.1[plots.1$Plot != 22, ]
hplots.1 =  hplots.1[hplots.1$Plot != 22, ]
resp = resp[resp$Plot != 22, ]

#End Stage 1
############


##############################################################################
##############################################################################
##############################################################################
#Stage 2: run models of temperature-respiration using two year moving windows

#order dataframes by year (just in case :) )
resp = resp[order(resp$year), ]
plots.1 = plots.1[order(plots.1$year, plots.1$Trt, plots.1$Plot, plots.1$Month, plots.1$Day), ]
hplots.1 = hplots.1[order(hplots.1$year, hplots.1$Trt, hplots.1$Plot, hplots.1$Month, hplots.1$Day, hplots.1$hour), ]

#create a dataframe with possible start and end years based on a sliding two-year window
start.year = seq(2006, 2021)
end.year = seq(2007, 2022)
wins = data.frame(cbind(start.year, end.year))

############################
#Respiration ~ temp model loop

#list for models
resp_mod_list = list()

#loop through 2-yr windows
for(i in 1:nrow(wins)){
    
    #2 years of data
    temp.resp = resp %>%
        filter(year == wins[i,1] | year == wins[i,2])
    
    #loop for each treatment
    #filter data to trt
    #run model and store in list
    trts = unique(resp$Trt)
    for(j in 1:length(trts)){
        temp.resp.trt = temp.resp %>%
            filter(Trt == trts[j])
        #note that we are using LnFlux NOT LnDFlux in the lm. The latter is 
        # Ln(Flux*24). The spot measurements are accurate on the short term basis.
        # Scaling these to the day makes sense if using the daily average temp
        # for prediction. However, the hourly temps are what end up used
        # in the final calcs. If necessary can multiply predictions based on the
        # daily average by 24.
        resp_mod_list[[ as.character(wins[i,2]) ]][[ trts[j] ]] = lm(
            LnFlux ~ temp,
            data = temp.resp.trt
        )
    }
}
#End respiration modeling
resp_mod_list

saveRDS(resp_mod_list, "data/model_output/2022respiration_lm_list.RDS")

#End stage 2
########################


##############################################################################
##############################################################################
##############################################################################

#####Stage 3: Predict respiration based on soil temperature, Trt, and year


################################################################
#Flux prediction using either daily avg temp or hourly
#Note that the daily ave predictions should be multiplied by 24
#bc we modeled based on the spot measurements

###################
#Daily preds

#list for each data subset
daily_preds = list()

#We subset the data based on the end point year in wins and use the 2-year
# window resp model with the same end year for prediction
for(i in 1:nrow(wins)){
    #1 year window
    temp.plots.1 = plots.1 %>%
        filter(year == wins[i,2])
    
    #list for holding treatment-level preds for simpler bind_rows()
    trt_preds = list()
    #loop treatments
    trts = unique(temp.plots.1$Trt)
    for(j in 1:length(trts)){
        temp.plots.1.trt = temp.plots.1 %>%
            filter(Trt == trts[j])
        temp.preds = data.frame(
            year = temp.plots.1.trt$year,
            doy = temp.plots.1.trt$doy,
            Plot = temp.plots.1.trt$Plot,
            Trt = temp.plots.1.trt$Trt,
            predflux = 24 * exp( #multiply the preds by 24 for daily
                predict(
                    resp_mod_list[[ as.character(wins[i,2]) ]][[ trts[j] ]],
                    newdata = temp.plots.1.trt
                )
            )
        )
        trt_preds[[ trts[j] ]] = temp.preds
    }
    daily_preds[[ as.character(wins[i, 2]) ]] = bind_rows(trt_preds)
}
all_pred = bind_rows(daily_preds)

###############
#Hourly

#list for each data subset
hourly_preds = list()

#We subset the data based on the end point year in wins and use the 2-year
# window resp model with the same end year for prediction
for(i in 1:nrow(wins)){
    #1 year window
    temp.hplots.1 = hplots.1 %>%
        filter(year == wins[i,2])
    
    #list for holding treatment-level preds for simpler bind_rows()
    trt_preds = list()
    #loop treatments
    trts = unique(temp.hplots.1$Trt)
    for(j in 1:length(trts)){
        temp.hplots.1.trt = temp.hplots.1 %>%
            filter(Trt == trts[j])
        temp.preds = data.frame(
            year = temp.hplots.1.trt$year,
            doy = temp.hplots.1.trt$doy,
            hour = temp.hplots.1.trt$hour,
            Plot = temp.hplots.1.trt$Plot,
            Trt = temp.hplots.1.trt$Trt,
            predflux = exp(
                predict(
                    resp_mod_list[[ as.character(wins[i,2]) ]][[ trts[j] ]],
                    newdata = temp.hplots.1.trt
                )
            )
        )
        trt_preds[[ trts[j] ]] = temp.preds
    }
    hourly_preds[[ as.character(wins[i, 2]) ]] = bind_rows(trt_preds)
}
all_pred.h = bind_rows(hourly_preds)

#End prediction loops
######################

#at some point the treatment designations in the data file were switched from numeric to letters.
#we substitute back to numeric here so the remaining code works as written

all_pred$Trt[all_pred$Trt == "C"] = 1
all_pred$Trt[all_pred$Trt == "H"] = 2
all_pred$Trt[all_pred$Trt == "HN"] = 3
all_pred$Trt[all_pred$Trt == "N"] = 4

all_pred.h$Trt[all_pred.h$Trt == "C"] = 1
all_pred.h$Trt[all_pred.h$Trt == "H"] = 2
all_pred.h$Trt[all_pred.h$Trt == "HN"] = 3
all_pred.h$Trt[all_pred.h$Trt == "N"] = 4


########################################
#Calculate annual sums of respiration

##############
#daily fluxes
#add together daily predicted respiration within each plot x year combination
plot.sum = aggregate(all_pred$predflux, by = list(all_pred$Trt, all_pred$Plot, all_pred$year), sum, na.rm = T)
names(plot.sum) = c("Trt", "Plot", "Year", "cum_flux")

#determine the number of days within each Year / Plot that 
#went into annual sum calculation

#flag non-NA values
all_pred$fluxNA = ifelse(is.na(all_pred$predflux) == T, 0, 1)
#add instances of non NA values together. This is the number of days measured
plot.comp = aggregate(all_pred$fluxNA, by = list(all_pred$Trt, all_pred$Plot, all_pred$year), sum, na.rm = T)
names(plot.comp) = c("Trt", "Plot", "Year", "comp_flux")

#add comp_flux to plot.sum
plot.sum$comp = plot.comp$comp_flux

#scale estimated fluxes to be 365 days
plot.sum$cum_flux_365 = (plot.sum$cum_flux / plot.sum$comp) * 365

#convert units from mg to g
plot.sum$cum_g = plot.sum$cum_flux_365 / 1000

#order by year
plot.sum = plot.sum[order(plot.sum$Year), ]

################
#hourly fluxes
#add together predicted respiration within each plot x year combination
plot.sum.h = aggregate(all_pred.h$predflux, by = list(all_pred.h$Trt, all_pred.h$Plot, all_pred.h$year), sum, na.rm = T)
names(plot.sum.h) = c("Trt", "Plot", "Year", "cum_flux")

#determine the number of days within each Year / Plot that 
#went into annual sum calculation

#flag non-NA values
all_pred.h$fluxNA = ifelse(is.na(all_pred.h$predflux) == T, 0, 1)
#add instances of non NA values together. This is the number of hours measured
plot.comp = aggregate(all_pred.h$fluxNA, by = list(all_pred.h$Trt, all_pred.h$Plot, all_pred.h$year), sum, na.rm = T)
names(plot.comp) = c("Trt", "Plot", "Year", "comp_flux")

#add comp_flux to plot.sum.h
plot.sum.h$comp = plot.comp$comp_flux

#scale estimated fluxes to be 365 days
plot.sum.h$cum_flux_365 = (plot.sum.h$cum_flux / plot.sum.h$comp) * (365 * 24)

#convert units from mg to g
plot.sum.h$cum_g = plot.sum.h$cum_flux_365 / 1000

#order by year
plot.sum.h = plot.sum.h[order(plot.sum.h$Year), ]

#End annual sums
##################

###################
#Data formatting

#REMOVING THE N PLOTS WITH NO TEMP MEASUREMENTS
plot.sum.h = plot.sum.h[!plot.sum.h$Plot %in% c(11, 16, 23),]

#assign replicates to plots
# We are removing the N plots 11, 16, and 23
# These had an average temperature assigned based on the other plots
# However, including these artificially reduces the error estimate for N
# because it essentially produces three "replicates" with the same value

plot.sum.rep = plot.sum.h
plot.sum.rep$Rep = NA

plot.sum.rep.1 = filter(plot.sum.rep, Plot == 24 | Plot == 13 | Plot == 10 | Plot == 6)
plot.sum.rep.1$Rep = 1

plot.sum.rep.2 = filter(plot.sum.rep, Plot == 8 | Plot == 12 | Plot == 15) #| Plot == 16
plot.sum.rep.2$Rep = 2

plot.sum.rep.3 = filter(plot.sum.rep, Plot == 2 | Plot == 4 | Plot == 14 | Plot == 21)
plot.sum.rep.3$Rep = 3

plot.sum.rep.4 = filter(plot.sum.rep, Plot == 1 | Plot == 17 | Plot == 18) # | Plot == 11 
plot.sum.rep.4$Rep = 4

plot.sum.rep.5 = filter(plot.sum.rep, Plot == 5 | Plot == 7 | Plot == 9 | Plot == 19)
plot.sum.rep.5$Rep = 5

plot.sum.rep.6 = filter(plot.sum.rep, Plot == 3 | Plot == 20 | Plot == 22) #| Plot == 23
plot.sum.rep.6$Rep = 6

plot.sum.rep = rbind(plot.sum.rep.1, plot.sum.rep.2, plot.sum.rep.3, plot.sum.rep.4, plot.sum.rep.5, plot.sum.rep.6)

#make year x rep column
plot.sum.rep$yr.rep = paste(plot.sum.rep$Year, plot.sum.rep$Rep, sep = " ")

###################################
#Calculate differences from control

#subset data to only include control plots
rep.c = plot.sum.rep[plot.sum.rep$Trt == 1, ]

#select columns needed for merging
rep.sub = rep.c[ , c("yr.rep", "cum_g")]
names(rep.sub) = c("yr.rep", "cum_c")

#merge rep.sub and plot.sum.rep
rep.all = merge(rep.sub, plot.sum.rep, by.x = "yr.rep", by.y = "yr.rep", all.x = T, all.y = T)

#calculate difference from control
rep.all$cdif = rep.all$cum_g - rep.all$cum_c

############################################
#calculate cumulative difference from control
rep.all$csum = ave(rep.all$cdif, rep.all$Plot, FUN = cumsum)

#make into data table for data reduction
repa = data.table(rep.all)

#create Year by Trt combination
repa$ID = paste(repa$Year, repa$Trt, sep = " ")

#calculate median and IQR
repsum = repa[ , list(Year = unique(Year), Trt = unique(Trt),
                      mcsum = median(csum), maxcsum = quantile(csum, 0.75), mincsum = quantile(csum, 0.25)), 
               by = ID]

repsum.mean_se = repa[ , list(Year = unique(Year), Trt = unique(Trt),
                mean.csum = mean(csum), SE.sum = sd(csum)/sqrt(length(csum))),
                by = ID]

#End stage 3
###############

###############################################################################
###############################################################################
###############################################################################

#####Stage 4: Run statistical analyses

########################################
#mixed effects models to examine differences among treatments, years, and treatment x years

plot.sum.h$Trt = factor(plot.sum.h$Trt)
plot.sum.h$Year = factor(plot.sum.h$Year)


gls.af = gls(cum_g ~ Trt * Year, data = plot.sum.h)
lme.af = lme(fixed = cum_g ~ Trt * Year, random = ~1|Plot, data = plot.sum.h)
anova(gls.af, lme.af)
#model fit improved with fixed and random effects, check which model gives the lower AIC

ac.af = update(lme.af, correlation=corAR1(form = ~ 1 | Plot))
AIC(lme.af, ac.af)
#model fit improved with correlation structure
var.af.1 = update(lme.af, weights = varIdent( ~ Plot | Trt))
var.af.2 = update(lme.af, weights = varIdent( ~ Plot | Year))
var.af.3 = update(lme.af, weights = varIdent( ~ Plot | Trt * Year))
AIC(lme.af, ac.af, var.af.1, var.af.2, var.af.3)


#model residuals and qq plots for model verification
plot(ac.af)
qqnorm(ac.af)

#pairwise comparisons

#create models for all single factor and interaction terms

#make a Trt x Year column
plot.sum.h$TrtYear = as.factor(paste(plot.sum.h$Trt, plot.sum.h$Year, sep = " "))

MCTrt <- lme(fixed = cum_g ~ Trt, correlation=corAR1(form = ~ 1 | Plot), random = ~1|Plot, data = plot.sum.h)
MCYear <- lme(fixed = cum_g ~ Year, correlation=corAR1(form = ~ 1 | Plot), random = ~1|Plot, data = plot.sum.h)
MCTrtYear <- lme(fixed = cum_g ~ TrtYear, correlation=corAR1(form = ~ 1 | Plot), random = ~1|Plot, data = plot.sum.h)

#use glht function and acquire summary of multiple comparisons
#set test adjustment as appropriate.

summary((glht(MCTrt, linfct = mcp(Trt = "Tukey"))))
summary((glht(MCYear, linfct = mcp(Year = "Tukey"))))

########################################

#t-tests to examine significant differences from zero

#reconfigure dataframe

rep.H = rep.all[rep.all$Trt == 2, ]
H.2007 = rep.H[rep.H$Year == 2007, ]
H.2008 = rep.H[rep.H$Year == 2008, ]
H.2009 = rep.H[rep.H$Year == 2009, ]
H.2010 = rep.H[rep.H$Year == 2010, ]
H.2011 = rep.H[rep.H$Year == 2011, ]
H.2012 = rep.H[rep.H$Year == 2012, ]
H.2013 = rep.H[rep.H$Year == 2013, ]
H.2014 = rep.H[rep.H$Year == 2014, ]
H.2015 = rep.H[rep.H$Year == 2015, ]
H.2016 = rep.H[rep.H$Year == 2016, ]
H.2017 = rep.H[rep.H$Year == 2017, ]
H.2018 = rep.H[rep.H$Year == 2018, ]
H.2019 = rep.H[rep.H$Year == 2019, ]
H.2020 = rep.H[rep.H$Year == 2020, ]
H.2021 = rep.H[rep.H$Year == 2021, ]
H.2022 = rep.H[rep.H$Year == 2022, ]

rep.HN = rep.all[rep.all$Trt == 3, ]
HN.2007 = rep.HN[rep.HN$Year == 2007, ]
HN.2008 = rep.HN[rep.HN$Year == 2008, ]
HN.2009 = rep.HN[rep.HN$Year == 2009, ]
HN.2010 = rep.HN[rep.HN$Year == 2010, ]
HN.2011 = rep.HN[rep.HN$Year == 2011, ]
HN.2012 = rep.HN[rep.HN$Year == 2012, ]
HN.2013 = rep.HN[rep.HN$Year == 2013, ]
HN.2014 = rep.HN[rep.HN$Year == 2014, ]
HN.2015 = rep.HN[rep.HN$Year == 2015, ]
HN.2016 = rep.HN[rep.HN$Year == 2016, ]
HN.2017 = rep.HN[rep.HN$Year == 2017, ]
HN.2018 = rep.HN[rep.HN$Year == 2018, ]
HN.2019 = rep.HN[rep.HN$Year == 2019, ]
HN.2020 = rep.HN[rep.HN$Year == 2020, ]
HN.2021 = rep.HN[rep.HN$Year == 2021, ]
HN.2022 = rep.HN[rep.HN$Year == 2022, ]

rep.N = rep.all[rep.all$Trt == 4, ]
N.2007 = rep.N[rep.N$Year == 2007, ]
N.2008 = rep.N[rep.N$Year == 2008, ]
N.2009 = rep.N[rep.N$Year == 2009, ]
N.2010 = rep.N[rep.N$Year == 2010, ]
N.2011 = rep.N[rep.N$Year == 2011, ]
N.2012 = rep.N[rep.N$Year == 2012, ]
N.2013 = rep.N[rep.N$Year == 2013, ]
N.2014 = rep.N[rep.N$Year == 2014, ]
N.2015 = rep.N[rep.N$Year == 2015, ]
N.2016 = rep.N[rep.N$Year == 2016, ]
N.2017 = rep.N[rep.N$Year == 2017, ]
N.2018 = rep.N[rep.N$Year == 2018, ]
N.2019 = rep.N[rep.N$Year == 2019, ]
N.2020 = rep.N[rep.N$Year == 2020, ]
N.2021 = rep.N[rep.N$Year == 2021, ]
N.2022 = rep.N[rep.N$Year == 2022,]


tran.H = data.frame(cbind("Plot" = H.2007$Plot, "delt.2007" = H.2007$cdif, "delt.2008" = H.2008$cdif, "delt.2009" = H.2009$cdif,
         "delt.2010" = H.2010$cdif, "delt.2011" = H.2011$cdif, "delt.2012" = H.2012$cdif, "delt.2013" = H.2013$cdif,
         "delt.2014" = H.2014$cdif, "delt.2015" = H.2015$cdif, "delt.2016" = H.2016$cdif, "delt.2017" = H.2017$cdif, "delt.2018" = H.2018$cdif, "delt.2019" = H.2019$cdif, "delt.2020" = H.2020$cdif, "delt.2021" = H.2021$cdif, "delt.2022" = H.2022$cdif))

tran.HN = data.frame(cbind("Plot" = HN.2007$Plot, "delt.2007" = HN.2007$cdif, "delt.2008" = HN.2008$cdif, "delt.2009" = HN.2009$cdif,
         "delt.2010" = HN.2010$cdif, "delt.2011" = HN.2011$cdif, "delt.2012" = HN.2012$cdif, "delt.2013" = HN.2013$cdif,
         "delt.2014" = HN.2014$cdif, "delt.2015" = HN.2015$cdif, "delt.2016" = HN.2016$cdif, "delt.2017" = HN.2017$cdif, "delt.2018" = HN.2018$cdif, "delt.2019" = HN.2019$cdif, "delt.2020" = HN.2020$cdif, "delt.2021" = HN.2021$cdif, "delt.2022" = HN.2022$cdif))

tran.N = data.frame(cbind("Plot" = N.2007$Plot, "delt.2007" = N.2007$cdif, "delt.2008" = N.2008$cdif, "delt.2009" = N.2009$cdif,
         "delt.2010" = N.2010$cdif, "delt.2011" = N.2011$cdif, "delt.2012" = N.2012$cdif, "delt.2013" = N.2013$cdif,
         "delt.2014" = N.2014$cdif, "delt.2015" = N.2015$cdif, "delt.2016" = N.2016$cdif, "delt.2017" = N.2017$cdif, "delt.2018" = N.2018$cdif, "delt.2019" = N.2019$cdif, "delt.2020" = N.2020$cdif, "delt.2021" = N.2021$cdif, "delt.2022" = N.2022$cdif))


#apply t-tests across columns

#These indices need to be updated when new years of data are added i.e, tran.H[2:X]
H.test <- lapply (tran.H[2:ncol(tran.H)], function(x)
              t.test(x, y = NULL, alternative = "two.sided", mu = 0, data = tran.H))

HN.test <- lapply (tran.HN[2:ncol(tran.HN)], function(x)
              t.test(x, y = NULL, alternative = "two.sided", mu = 0, data = tran.HN))

N.test <- lapply (tran.N[2:ncol(tran.N)], function(x)
              t.test(x, y = NULL, alternative = "two.sided", mu = 0, data = tran.N))

#extract p-values
H.test.p <- sapply(H.test , '[', 'p.value')
H.test.res <- as.data.frame(unlist(H.test.p))
names(H.test.res) <- "p_val"

HN.test.p <- sapply(HN.test , '[', 'p.value')
HN.test.res <- as.data.frame(unlist(HN.test.p))
names(HN.test.res) <- "p_val"

N.test.p <- sapply(N.test , '[', 'p.value')
N.test.res <- as.data.frame(unlist(N.test.p))
names(N.test.res) <- "p_val"

H.test.p
HN.test.p
N.test.p

###########################################################
###########################################################
###########################################################

##### Output results

#quality controled temperature values
write.csv(plots.1, "data/model_output/2022_daily_temps.csv", quote = F, row.names = F)
write.csv(hplots.1, "data/model_output/2022_hourly_temps.csv", quote = F, row.names = F)

#annual sums of flux and deltas (hourly data)
write.csv(plot.sum.h, "data/model_output/2022_plot.sum.h.csv", quote = F, row.names = F)
write.csv(rep.all, "data/model_output/2022_rep.all.csv", quote = F, row.names = F)

#treatment median and IQR of deltas
write.csv(repsum, "data/model_output/2022_repsum.csv", quote = F, row.names = F)
#treatment mean and SE of deltas
write.csv(repsum.mean_se, "data/model_output/2022_repsum.mean_se.csv", quote = F, row.names = F)

