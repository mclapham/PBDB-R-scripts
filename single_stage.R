#read stage match
time_match <- read.csv("https://docs.google.com/spreadsheets/d/1DPWsRF3DmGEo0eJXjvotuDJLphFaSYqjfCpYXWSgDvY/export?gid=0&format=csv")

single.stage <- function(occurrence_data) {
  
  #match PBDB interval to international stage
  occurrence_data$max_stage <- time_match$matched_stage[match(occurrence_data$early_interval, time_match$int_name)]
  occurrence_data$min_stage <- time_match$matched_stage[match(occurrence_data$late_interval, time_match$int_name)]
  
  #must remove ones where max_stage doesn't match
  occurrence_data <- subset(occurrence_data, is.na(occurrence_data$max_stage)==F) 
  
  #ones where already belongs to one stage
  single_int <- subset(occurrence_data, occurrence_data$late_interval=="")
  
  #must fill in min stage for single interval
  single_int$min_stage <- single_int$max_stage 
  
  #ones where PBDB has both max and min interval
  mult_int <- subset(occurrence_data, occurrence_data$late_interval!="")
  
  #finds where max international interval is same as min international interval 
  mult_int <- subset(mult_int, mult_int$max_stage==mult_int$min_stage)
  
  rbind(single_int, mult_int)
}
