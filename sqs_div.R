#SHAREHOLDER QUORUM SUBSAMPLING
#Based on exact SQS method described in Alroy (2014)
#No guarantee that this works correctly. Fossilworks (www.fossilworks.org) offers SQS diversity curve generation

#GETTING STARTED
#To use this, copy the entire text from this document and paste it all at the prompt
#Wait for it to finish and then type the extinction.calc function with your parameters
#a taxon name, sampling threshold (quorum, between 0 and 1), and number of trials are the minimum requirements
#Parameters are described in more detail below
#Here is an example:
#sqs.calc("Bivalvia","Roadian","Olenekian",q=0.7,trials=10,"genus","marine")

#REQUIRED PARAMETERS
#Specify a taxon or taxa to include in the analysis in quotations, separated by commas (if multiple taxa)
#Taxa do not have to be at the same rank. E.g. extinction.calc("Theropoda,Equidae,Apatosaurus",...)

#To analyze a specific time range, specify a name in maxinterval and mininterval
#E.g., sqs.calc("Theropoda,Equidae,Apatosaurus","Early Albian","Barstovian",...)
#Any PBDB time interval can be used, including regional scales, but rates will be calculated at stage or 10 Myr bin

#Specify a quorum, between 0 and 1 (see Alroy 2010 for more details)
#Specify the number of resampling trials desired, e.g. trials=10

#OPTIONAL PARAMETERS
#Diversity can be calculated at stage (default) or "10 Myr bin" (sensu Alroy et al., 2008) resolution by specifying temp_res="myrbin"
#The analysis can calculate diversity at the species, genus, or family level by specifying one of tax_level="species", "genus", or "family" (genus is default)
#Environment can be specified as one of environ="reef", "siliciclastic", "carbonate", "marine", or "terrestrial" (the default is to include all occurrences)
#Analysis can be restricted to species that have formal opinions classifying them, by specifying formal_id="yes" (warning: some intervals have very few formally classified species, especially for marine inverts)

#THE OUTPUT
#Results are stored as diversity_results, which you can view in R by typing that name at the text prompt
#To write the results to a file viewable in a spreadsheet, paste this command at the text prompt: write.csv(diversity_results,"sqs_diversity.csv")
#The file, named sqs_diversity.csv, will be saved to the current working directory, which is most often your Documents folder. To find that folder, type getwd() at the R prompt

sqs.calc<-function(include_taxon,maxinterval="Phanerozoic",mininterval="Phanerozoic",q,trials,temp_res="stage",tax_level="genus",environ="all",formal_id="no") {
  
  #DATA ACQUISITION
  print("Requesting data from server")
  
  #initializes timing
  ptm<-proc.time()
  
  #reads names and ages of time intervals
  time_int<-read.csv("http://paleobiodb.org/data1.1/intervals/list.txt?scale=all&limit=all")
  
  #finds maximum age of oldest user-specified interval
  max_interval_ma<-subset(time_int$early_age,time_int$interval_name==maxinterval)
  
  #finds minimum age of youngest user-specified interval
  min_interval_ma<-subset(time_int$late_age,time_int$interval_name==mininterval)
  
  #reads occurrences based on specified taxa and age range
  occurrences<-read.csv(paste("http://paleobiodb.org/data1.1/occs/list.txt?base_name=",include_taxon,"&min_ma=",min_interval_ma,"&max_ma=",max_interval_ma,"&show=time,phylo,geo,ident&limit=all",sep=""))
  
  #reports time
  print(paste("Finished data acquisition: elapsed time",round(as.numeric((proc.time()-ptm)[3]),1),"seconds"))
  
  #finds only those collections resolved to a single time interval
  if (temp_res=="myrbin") {
    #loads file with matches between PBDB intervals and "10 myr bins"
    time_conv<-read.csv("https://github.com/mclapham/PBDB-R-scripts/raw/master/time_convers.csv")
    
    #deletes occurrences where at least one of the time intervals is not contained within a single 10-myr-bin
    occurrences<-subset(occurrences,occurrences$early_interval %in% time_conv$interval_name & occurrences$late_interval %in% time_conv$interval_name)
    
    #matches occurrence name (for max and min interval) and enters 10-myr-bin number in new field
    occurrences$myr_bin_max<-time_conv$myr_bin_no[match(occurrences$early_interval,time_conv$interval_name)]
    occurrences$myr_bin_min<-time_conv$myr_bin_no[match(occurrences$late_interval,time_conv$interval_name)]
    
    #finds occurrences confined to single 10-myr-bin
    resolved_occs<-subset(occurrences,occurrences$myr_bin_max==occurrences$myr_bin_min)
    
    #creates variable in data frame with 10-myr-bin number
    resolved_occs$time_int<-resolved_occs$myr_bin_max
    
  } else {
    resolved_occs<-subset(occurrences,occurrences$cx_int_no %in% subset(time_int$interval_no,time_int$level==5))
    
    #creates variable in data frame with stage number
    resolved_occs$time_int<-resolved_occs$cx_int_no
  }
  
  
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
  
  #FILE PREPARATION AND DATA CLEANING
  #file preparation and data cleaning for species level analysis
  if(tax_level=="species")
  {
    #removes rows where species qualified by cf. or aff., question mark, ex gr. or quotation mark
    resolved_occs<-subset(resolved_occs,resolved_occs$species_reso=="" | resolved_occs$species_reso=="n. sp.")
    
    if (formal_id=="yes") {
      #deletes occurrences not resolved to species level using formally-classified species
      cleaned_occs<-subset(resolved_occs,resolved_occs$matched_rank==3)
      
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
      
      #combines classified and unclassified occurrences
      cleaned_occs<-rbind(classified_occs,unclassified_occs)
    }
    
    #or file preparation and data cleaning for family level analysis
  } else if (tax_level=="family") {
    
    #deletes occurrences where genus or family are qualified with question mark, quotations, cf. or aff.
    if(length(grep("\\? ",resolved_occs$taxon_name))>0) resolved_occs<-resolved_occs[-grep("\\? ",resolved_occs$taxon_name),]
    if(length(grep("\" ",resolved_occs$taxon_name))>0) resolved_occs<-resolved_occs[-grep("\" ",resolved_occs$taxon_name),]
    if(length(grep("cf. ",resolved_occs$taxon_name))>0) resolved_occs<-resolved_occs[-grep("cf. ",resolved_occs$taxon_name),]
    if(length(grep("aff. ",resolved_occs$taxon_name))>0) resolved_occs<-resolved_occs[-grep("aff. ",resolved_occs$taxon_name),]
    
    #deletes occurrences not resolved to at least family level
    cleaned_occs<-subset(resolved_occs,resolved_occs$matched_rank<=9)
    
    cleaned_occs$matched_name<-cleaned_occs$family
    
    #or file preparation and data cleaning for genus level analysis
  } else {
    
    #deletes occurrences not resolved to at least genus level using classified and unclassified species
    cleaned_occs<-subset(resolved_occs,resolved_occs$matched_rank<=5)
    
    #deletes occurrences where genus is qualified with question mark, quotations, cf. or aff.
    cleaned_occs<-subset(cleaned_occs,cleaned_occs$genus_reso=="" | cleaned_occs$genus_reso=="n. gen.")
    
    #extracts genus name from matched_name string
    cleaned_occs$matched_name<-gsub(" .*","",cleaned_occs$matched_name)
  }
    
  goodsu<-function(taxon_list) {
    #GOOD'S U CALCULATION
    #Single-occurrence taxa
    p1<-length(which(table(taxon_list)==1))
    
    #Most abundant taxon
    o1<-max(table(taxon_list))
    
    O<-length(taxon_list) #Total number of occurrences
    
    #Good's u (coverage of data)
    1-p1/(O-o1)
    
  }
  
  #19 to 32
  occurrence_data<-subset(cleaned_occs,cleaned_occs$time_int==24)
  
  #EXACT SQS METHOD
  
  sqs.method<-function(occurrence_data) {
    
    if(nrow(occurrence_data)>25) { #an arbitrary cutoff below which the subsampling is unlikely to work (or be meaningful)
      
      #finds three random collections per reference (if fewer than three, samples all)
      rand_coll_list<-sapply(split(occurrence_data$collection_no,occurrence_data$reference_no),function(x) unique(x)[sample(seq(1:length(unique(x))))][1:3])
      
      #orders references in random order
      rand_coll_list<-rand_coll_list[,sample(seq(1,ncol(rand_coll_list)))]
      
      #vector of reference_no-collection_no combinations from random list
      ref_coll_comb<-paste((rep(colnames(rand_coll_list),apply(rand_coll_list,2,function(x) length(subset(x,x>0))))),na.omit(as.vector(rand_coll_list)))
      
      #subsets occurrences that match both reference number and collection number
      sub_occs<-occurrence_data[which(paste(occurrence_data$reference_no,occurrence_data$collection_no) %in% ref_coll_comb),]
      
      #orders occurrences by random reference order (and then by previously-chosen random collection order within reference)
      sub_occs<-sub_occs[order(match(paste(sub_occs$reference_no,sub_occs$collection_no),ref_coll_comb)),]
      
      #finds starting position of each reference/collection combination
      collection_start<-sort(match(ref_coll_comb,paste(sub_occs$reference_no,sub_occs$collection_no)))
      
      u<-numeric(0)
      div_ct<-numeric(0)
      n=1
      #calculates Good's u after adding each collection
      for (i in collection_start) {
        data_sub<-sub_occs$matched_name[1:i]
        u[n]<-goodsu(data_sub)
      
        if (n>2) {
          if(is.nan(u[n])==F & is.nan(u[n-1])==F) {
            #finds diversity at points where Good's U crosses above and below target q
            if (u[n]>q & u[n-1]<q) {
              div_ct<-c(div_ct,length(unique(data_sub)))
            }
          
            if (u[n]<q & u[n-1]>q) {
              div_ct<-c(div_ct,length(unique(data_sub)))
            }
          }
        }
        n=n+1
      }
      
      #finds median diversity 
      median(div_ct)
    
    } else {return(NA)}
  }
  
  sqs.apply<-function(occs_data) {
    #presence-absence matrix for sampled taxa in each time interval
    sapply(split(occs_data,occs_data$time_int),sqs.method)
  }
  
  #repeates subsampling and diversity calculation 
  sqs_diversity<<-replicate(trials,sqs.apply(cleaned_occs))
  
  #calculates average per time interval
  div_avg<-apply(sqs_diversity,1,function(x) mean(x,na.rm=T))
  
  if (temp_res=="stage") {
    int_names<-time_int$interval_name[match(as.numeric(rownames(sqs_diversity)),time_int$interval_no)]
    int_min<-rowMeans(cbind(time_int$late_age[match(as.numeric(rownames(sqs_diversity)),time_int$interval_no)],time_int$early_age[match(as.numeric(rownames(sqs_diversity)),time_int$interval_no)]))
  } else if (temp_res=="myrbin") {
    int_names<-time_conv$myr_bin_name[match(as.numeric(rownames(sqs_diversity)),time_conv$myr_bin_no)]
    myr_stages<-c("Fortunian","Early Cambrian","Middle Cambrian","Furongian","Tremadocian","Arenig","Llanvirn","Caradoc","Hirnantian","Telychian","Pridoli","Pragian","Emsian","Givetian","Frasnian","Famennian","Tournaisian","Asbian","Serpukhovian","Moscovian","Gzhelian","Sakmarian","Kungurian","Capitanian","Changhsingian","Olenekian","Ladinian","Carnian","Rhaetian","Sinemurian","Pliensbachian","Aalenian","Bathonian","Kimmeridgian","Tithonian","Valanginian","Barremian","Aptian","Albian","Cenomanian","Santonian","Campanian","Maastrichtian","Thanetian","Lutetian","Priabonian","Chattian","Serravallian","Late Pleistocene")
    myr_endpt<-time_int$late_age[match(myr_stages,time_int$interval_name)]-0.5*diff(c(542,time_int$late_age[match(myr_stages,time_int$interval_name)]))
    int_min<-myr_endpt[as.numeric(rownames(sqs_diversity))]
  }
  
  #results table with interval names and average extinction rate
  diversity_results<<-data.frame(interval=int_names,age_ma=int_min,div_avg)
  diversity_results<<-diversity_results[order(diversity_results$age_ma),]
}


