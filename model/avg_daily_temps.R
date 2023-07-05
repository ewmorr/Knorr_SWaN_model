library(dplyr)

htemps = read.csv("2021_data/cleaned_hourly_swntemp_2006_2022.csv")
head(htemps)

dtemp = htemps %>% 
  group_by(YEAR, MONTH, DAY) %>%
  summarize(
    Count = n(),
    across(-TIME, mean)
  ) 

colnames(dtemp)[1:3] = c("Year", "Month", "Day")
write.csv(dtemp, "daily_avg_swntemp.csv", row.names = F)
