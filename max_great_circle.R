#GEOGRAPHIC RANGE CALCULATION - MAXIMUM GREAT CIRCLE DISTANCE
#Uses the haversine formula to calculate the maximum great circle distance between any two occurrences of a taxon in a time interval
#Version 1.0 (Matthew Clapham, mclapham@ucsc.edu)

#GETTING STARTED
#To use this, copy the entire text from this document and paste it all at the prompt
#Wait for it to finish and then type the geog_range function with your parameters (a taxon name and time interval is the minimum requirement)
#Parameters are described in more detail below
#Here is an example: geog_range("Rhynchonelliformea","Wuchiapingian","Olenekian","genus","marine")

#REQUIRED PARAMETERS
#Specify a taxon or taxa to include in the analysis in quotations, separated by commas (if multiple taxa)
#Taxa do not have to be at the same rank. E.g. geog_range("Theropoda,Equidae,Apatosaurus",...)
#To analyze multiple time intervals (separate geographic ranges calculated will be calculated for each), specify a name in maxinterval and mininterval
#E.g., geog_range("Theropoda,Equidae,Apatosaurus","Early Albian","Barstovian",...)
#Any PBDB time interval can be used, including regional scales, but analyses will be split by official ICS stages

#OPTIONAL PARAMETERS
#The analysis can calculate the geographic range at the species, genus, or family level by specifying one of tax_level="species", "genus", or "family" (genus is default)
#Environment can be specified as one of environ="reef", "siliciclastic", "carbonate", "marine", or "terrestrial" (the default is to include all occurrences)
#Analysis can be restricted to species that have formal opinions classifying them, by specifying formal_id="yes" (warning: some intervals have very few formally classified species, especially for marine inverts)

#THE OUTPUT
#Results are stored as max_gc_dist, which you can view in R by typing that name at the text prompt
#To write the results to a file viewable in a spreadsheet, paste this command at the text prompt: write.csv(max_gc_dist,"geog_range_results.csv")
#The file, named geog_range_results, will be saved to the current working directory, which is most often your Documents folder. To find that folder, type getwd() at the R prompt


geog_range<-function(include_taxon,maxinterval,mininterval=maxinterval,tax_level="genus",environ="all",formal_id="no") {

#DATA ACQUISITION
  #reads names and ages of time intervals
  time_int<-read.csv("http://paleobiodb.org/data1.1/intervals/list.txt?scale=all&limit=all")
  
  #finds maximum age of oldest user-specified interval
  max_interval_ma<-subset(time_int$early_age,time_int$interval_name==maxinterval)
  
  #finds minimum age of youngest user-specified interval
  min_interval_ma<-subset(time_int$late_age,time_int$interval_name==mininterval)
  
  #reads occurrences based on specified taxa and age range
  occurrences<-read.csv(paste("http://paleobiodb.org/data1.1/occs/list.txt?base_name=",include_taxon,"&min_ma=",min_interval_ma,"&max_ma=",max_interval_ma,"&show=time,paleoloc,phylo,geo,ident&limit=all",sep=""))
  
  #finds only those collections resolved to a stage
  resolved_occs<-subset(occurrences,occurrences$cx_int_no %in% subset(time_int$interval_no,time_int$level==5))
 
  #list of marine environments
  reef_env<-c("reef, buildup or bioherm","perireef or subreef","intrashelf/intraplatform reef","platform/shelf-margin reef","slope/ramp reef","basin reef")
  carbonate_env<-c("carbonate indet.","peritidal","shallow subtidal indet.","open shallow subtidal","lagoonal/restricted shallow subtidal","sand shoal","reef, buildup or bioherm","perireef or subreef","intrashelf/intraplatform reef","platform/shelf-margin reef","slope/ramp reef","basin reef","deep subtidal ramp","deep subtidal shelf","deep subtidal indet.","offshore ramp","offshore shelf","offshore indet.","slope","basinal (carbonate)","basinal (siliceous)")
  siliciclastic_env<-c("marine indet.","marginal marine indet.","coastal indet.","estuary/bay","lagoonal","paralic indet.","delta plain","interdistributary bay","delta front","prodelta","deltaic indet.","foreshore","shoreface","transition zone/lower shoreface","offshore","submarine fan","basinal (siliciclastic)","deep-water indet.")
  marine_env<-c(carbonate_env,siliciclastic_env)
  
  #list of terrestrial environments
  terrestrial_env<-c("terrestrial indet.","fluvial indet.","alluvial fan","channel lag","coarse channel fill","fine channel fill","channel","wet floodplain","dry floodplain","floodplain","crevasse splay","levee","mire/swamp","fluvial-lacustrine indet.","delta plain","fluvial-deltaic indet.","lacustrine-large","lacustrine-small","pond","crater lake","lacustrine delta plain","lacustrine interdistributary bay","lacustrine delta front","lacustrine prodelta","lacustrine deltaic indet.","lacustrine indet.","dune","interdune","loess","eolian indet.","cave","fissure fill","sinkhole","karst indet.","tar","spring","glacial")
  
  #finds collections with appropriate environments
  if(environ!="all") {
    resolved_occs$environment<-gsub("\"","",resolved_occs$environment) #some terrestrial environments have quotations in the name, which must be removed
    resolved_occs<-subset(resolved_occs,resolved_occs$environment %in% get(noquote(paste(environ,"_env",sep=""))))
  }
      
  #ASSORTED SUBFUNCTIONS FOR CALCULATIONS
  #1 and 2: functions to convert radians to degrees or degrees to radians
  as.radians<-function(degrees) {degrees*pi/180}
  as.degrees<-function(radians) {radians*180/pi}
  
  #3: Haversine formula to calculate great circle distance
  great.circle<-function(lat1,lat2,long1,long2) {6371*(acos(sin(as.radians(lat1))*sin(as.radians(lat2))+cos(as.radians(lat1))*cos(as.radians(lat2))*cos(as.radians(long2)-as.radians(long1))))}
  
  #4: Function to calculate maximum great circle distance for a taxon, using Haversine formula above
  max.dist<-function(latdata,longdata) {
    
    if (length(unique(latdata))>1 | length(unique(longdata))>1) {
      
      #creates all possible pairs of latitudes and longitudes
      lat.pairs<-combn(latdata,2)
      long.pairs<-combn(longdata,2)
      
      #calculates great circle distance for all possible lat-long pairs
      gc.dist<-great.circle(lat.pairs[1,],lat.pairs[2,],long.pairs[1,],long.pairs[2,])
      
      #finds maximum great circle distance from list, or sets distance as 1 km (in case of all occurrences being at a single coordinate)
      maxdist<<-max(gc.dist,na.rm=T) } else if (length(unique(latdata))==1 & length(unique(longdata))==1) {
      
      maxdist<<-1 } else maxdist<<-NA
  }
  #End 4: maximum great circle distance function

#FILE PREPARATION AND DATA CLEANING
  #file preparation and data cleaning for species level analysis
  if(tax_level=="species")
    {
    #removes rows where species qualified by cf. or aff., question mark, ex gr. or quotation mark
    resolved_occs<-subset(resolved_occs,resolved_occs$species_reso=="" | resolved_occs$species_reso=="n. sp.")
    
    if (formal_id=="yes") {
      #deletes occurrences not resolved to species level using formally-classified species
      resolved_occs<-subset(resolved_occs,resolved_occs$matched_rank==3)
      
      #prepares a data frame with only the necessary columns
      cleaned_occs<-data.frame(taxon_name=resolved_occs$matched_name,interval_name=resolved_occs$early_interval,interval_no=resolved_occs$cx_int_no,paleolng=resolved_occs$paleolng,paleolat=resolved_occs$paleolat)
   
    } else {
      #deletes occurrences not resolved to species level using classified and unclassified species
      resolved_occs<-subset(resolved_occs,resolved_occs$taxon_rank=="species")
        
      #strips all extraneous qualifiers from taxon name
      resolved_occs$taxon_name<-gsub("\" ","",resolved_occs$taxon_name)
      resolved_occs$taxon_name<-gsub("cf. ","",resolved_occs$taxon_name)
      resolved_occs$taxon_name<-gsub("aff. ","",resolved_occs$taxon_name)
      resolved_occs$taxon_name<-gsub("\\? ","",resolved_occs$taxon_name)
      resolved_occs$taxon_name<-gsub("n. gen. ","",resolved_occs$taxon_name)
      resolved_occs$taxon_name<-gsub(" n. sp.","",resolved_occs$taxon_name)
      resolved_occs$taxon_name<-gsub(" sensu lato","",resolved_occs$taxon_name)
      resolved_occs$taxon_name<-gsub(" informal","",resolved_occs$taxon_name)
      
      #separates classified and unclassified occurrences
      classified_occs<-subset(resolved_occs,resolved_occs$matched_rank==3)
      unclassified_occs<-subset(resolved_occs,resolved_occs$matched_rank!=3)
      
      #prepares a data frame with only the necessary columns for classified and unclassified occurrences
      cleaned_occs_class<-data.frame(taxon_name=classified_occs$matched_name,interval_name=classified_occs$early_interval,interval_no=classified_occs$cx_int_no,paleolat=classified_occs$paleolat,paleolng=classified_occs$paleolng)
      cleaned_occs_unclass<-data.frame(taxon_name=unclassified_occs$taxon_name,interval_name=unclassified_occs$early_interval,interval_no=unclassified_occs$cx_int_no,paleolat=unclassified_occs$paleolat,paleolng=unclassified_occs$paleolng)
      
      #combines classified and unclassified occurrences
      cleaned_occs<-rbind(cleaned_occs_class,cleaned_occs_unclass)
    }
  
    #or file preparation and data cleaning for family level analysis
  } else if (tax_level=="family") {
 
    #deletes occurrences not resolved to at least family level using classified and unclassified species
    resolved_occs<-subset(resolved_occs,resolved_occs$matched_rank<=9)
    
    #deletes occurrences where genus or family are qualified with question mark, quotations, cf. or aff.
    if(length(grep("\\? ",resolved_occs$taxon_name))>0) resolved_occs<-resolved_occs[-grep("\\? ",resolved_occs$taxon_name),]
    if(length(grep("\" ",resolved_occs$taxon_name))>0) resolved_occs<-resolved_occs[-grep("\" ",resolved_occs$taxon_name),]
    if(length(grep("cf. ",resolved_occs$taxon_name))>0) resolved_occs<-resolved_occs[-grep("cf. ",resolved_occs$taxon_name),]
    if(length(grep("aff. ",resolved_occs$taxon_name))>0) resolved_occs<-resolved_occs[-grep("aff. ",resolved_occs$taxon_name),]
    
    #prepares a data frame with only the necessary columns
    cleaned_occs<-data.frame(taxon_name=resolved_occs$family,interval_name=resolved_occs$early_interval,interval_no=resolved_occs$cx_int_no,paleolat=resolved_occs$paleolat,paleolng=resolved_occs$paleolng)
    
    #or file preparation and data cleaning for genus level analysis
  } else {
    
    #deletes occurrences not resolved to at least genus level using classified and unclassified species
    resolved_occs<-subset(resolved_occs,resolved_occs$matched_rank<=5)
    
    #deletes occurrences where genus is qualified with question mark, quotations, cf. or aff.
    resolved_occs<-subset(resolved_occs,resolved_occs$genus_reso=="" | resolved_occs$genus_reso=="n. sp.")
      
    #prepares a data frame with only the necessary columns
    cleaned_occs<-data.frame(taxon_name=resolved_occs$genus_name,interval_name=resolved_occs$early_interval,interval_no=resolved_occs$cx_int_no,paleolat=resolved_occs$paleolat,paleolng=resolved_occs$paleolng)
    
  }

#GREAT CIRCLE DISTANCE ANALYSIS
  #calculates maximum great circle distance for occurrences of each taxon in each stage
  max_gc_dist<<-sapply(split(cleaned_occs,cleaned_occs$interval_no),function(x) sapply(split(x,x$taxon_name),function(y) max.dist(y$paleolat,y$paleolng)))

#RESULTS ORDERING
  #the intervals are in ascending numerical order by interval_no, but not necessarily chronological (this finds the chronological order)
  interval_order<-order(sapply(as.numeric(colnames(max_gc_dist)),function(x) which(time_int$interval_no==x)),decreasing=T)
  
  #moves columns into chronological order
  max_gc_dist<<-max_gc_dist[,interval_order]
   
  #renames columns with interval name
  colnames(max_gc_dist)<<-time_int$interval_name[sapply(colnames(max_gc_dist),function(x) which(time_int$interval_no==x))] 
  #removes empty rows
  max_gc_dist<<-subset(max_gc_dist,apply(max_gc_dist,1,function(x) sum(x,na.rm=T))>0)
}
