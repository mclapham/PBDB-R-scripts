read.pbdb <- function(include_taxon,maxinterval="Phanerozoic",mininterval="Phanerozoic",temp_res="stage",environ="all") {
  
  #DATA ACQUISITION
  #reads names and ages of time intervals
  time_int <- read.csv("http://paleobiodb.org/data1.1/intervals/list.txt?scale=all&limit=all")
  
  #finds maximum age of oldest user-specified interval
  max_interval_ma <- subset(time_int$early_age,time_int$interval_name==maxinterval)
  
  #finds minimum age of youngest user-specified interval
  min_interval_ma <- subset(time_int$late_age,time_int$interval_name==mininterval)
  
  #reads occurrences based on specified taxa and age range
  occurrences <- read.csv(paste("http://paleobiodb.org/data1.1/occs/list.txt?base_name=",include_taxon,"&min_ma=",min_interval_ma,"&max_ma=",max_interval_ma,"&show=time,phylo,geo,ident,paleoloc,&limit=all",sep=""))
  
  #finds only those collections resolved to a single time interval
  if (temp_res=="stage") {
    resolved_occs <- subset(occurrences,occurrences$cx_int_no %in% subset(time_int$interval_no,time_int$level==5))
    
    #creates variable in data frame with stage number
    resolved_occs$time_int <- resolved_occs$cx_int_no
  } else {resolved_occs <- occurrences}
  
  #list of marine environments
  reef_env <- c("reef, buildup or bioherm","perireef or subreef","intrashelf/intraplatform reef","platform/shelf-margin reef","slope/ramp reef","basin reef")
  carbonate_env <- c("carbonate indet.","peritidal","shallow subtidal indet.","open shallow subtidal","lagoonal/restricted shallow subtidal","sand shoal","reef, buildup or bioherm","perireef or subreef","intrashelf/intraplatform reef","platform/shelf-margin reef","slope/ramp reef","basin reef","deep subtidal ramp","deep subtidal shelf","deep subtidal indet.","offshore ramp","offshore shelf","offshore indet.","slope","basinal (carbonate)","basinal (siliceous)")
  siliciclastic_env <- c("marine indet.","marginal marine indet.","coastal indet.","estuary/bay","lagoonal","paralic indet.","delta plain","interdistributary bay","delta front","prodelta","deltaic indet.","foreshore","shoreface","transition zone/lower shoreface","offshore","submarine fan","basinal (siliciclastic)","deep-water indet.")
  marine_env <- c(carbonate_env,siliciclastic_env)
  
  #list of terrestrial environments
  terrestrial_env <- c("terrestrial indet.","fluvial indet.","alluvial fan","channel lag","coarse channel fill","fine channel fill","channel","wet floodplain","dry floodplain","floodplain","crevasse splay","levee","mire/swamp","fluvial-lacustrine indet.","delta plain","fluvial-deltaic indet.","lacustrine-large","lacustrine-small","pond","crater lake","lacustrine delta plain","lacustrine interdistributary bay","lacustrine delta front","lacustrine prodelta","lacustrine deltaic indet.","lacustrine indet.","dune","interdune","loess","eolian indet.","cave","fissure fill","sinkhole","karst indet.","tar","spring","glacial")
  
  #finds collections with appropriate environments
  if(environ!="all") {
    resolved_occs$environment <- gsub("\"","",resolved_occs$environment) #some terrestrial environments have quotations in the name, which must be removed
    resolved_occs <- subset(resolved_occs,resolved_occs$environment %in% get(noquote(paste(environ,"_env",sep=""))))
  }
  
  resolved_occs
  
}
