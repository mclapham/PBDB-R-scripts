#read stage match
time_match <- read.csv("https://docs.google.com/spreadsheets/d/1DPWsRF3DmGEo0eJXjvotuDJLphFaSYqjfCpYXWSgDvY/export?gid=0&format=csv")

single.stage <- function(occurrence_data, reso = "stage") {
  
  if (reso == "stage") {
    occurrence_data$max_stage <- time_match$matched_stage[match(occurrence_data$early_interval, time_match$int_name)]
    occurrence_data$min_stage <- time_match$matched_stage[match(occurrence_data$late_interval, time_match$int_name)]
  } else if (reso == "series") {
    occurrence_data$max_stage <- time_match$matched_series[match(occurrence_data$early_interval, time_match$int_name)]
    occurrence_data$min_stage <- time_match$matched_series[match(occurrence_data$late_interval, time_match$int_name)]
  }

  
  occurrence_data <- subset(occurrence_data, is.na(occurrence_data$max_stage)==F) #must remove ones where max_stage doesn't match
  
  single_int <- subset(occurrence_data, occurrence_data$late_interval=="")
  
  single_int$min_stage <- single_int$max_stage #must fill in min stage for single interval
  
  mult_int <- subset(occurrence_data, occurrence_data$late_interval!="")
  
  mult_int <- subset(mult_int, mult_int$max_stage==mult_int$min_stage)
  
  rbind(single_int, mult_int)
}
