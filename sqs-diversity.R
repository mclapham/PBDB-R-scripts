#SHAREHOLDER QUORUM SUBSAMPLING DIVERSITY ANALYSIS
#Based on the method described by Alroy (2010)
#(I think this should work correctly!)
#Enter the taxonomic groups to be included (comma-separated list accepted) as include_taxon
#Enter the oldest (maxinterval) and youngest (mininterval) time intervals, using any PaleoDB time interval name
#Enter the quorum (q) you wish to sample (a number between 0 and 1)
  #TO ADD: an error check in case q is lower than the coverage of a time interval
#Enter the number of replicates (n_rep); default is 10
#Optional: select taxonomic level (species, genus, family)
#Optional: environment (terrestrial, marine, siliciclastic, carbonate, reef)
#Optional: restrict analysis to formally-classified species (only relevant if tax_level="species")

sqs.div<-function(include_taxon,maxinterval,mininterval=maxinterval,q,n_rep=10,tax_level="genus",environ="all",formal_id="no") {

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
  
  #reads references corresponding to occurrences
  references<-read.csv(paste("http://paleobiodb.org/data1.1/occs/refs.txt?base_name=",include_taxon,"&min_ma=",min_interval_ma,"&max_ma=",max_interval_ma,"&limit=all",sep=""))
    
  #reports time
  print(paste("Finished data acquisition: elapsed time",round(as.numeric((proc.time()-ptm)[3]),1),"seconds"))
  
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
      cleaned_occs<-data.frame(collection_no=resolved_occs$collection_no,taxon_name=resolved_occs$matched_name,interval_no=resolved_occs$cx_int_no,reference_no=resolved_occs$reference_no)
   
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

      cleaned_occs<-data.frame(collection_no=resolved_occs$collection_no,taxon_name=resolved_occs$matched_name,interval_no=resolved_occs$cx_int_no,reference_no=resolved_occs$reference_no)
      
      #prepares a data frame with only the necessary columns for classified and unclassified occurrences
      cleaned_occs_class<-data.frame(collection_no=classified_occs$collection_no,taxon_name=classified_occs$matched_name,interval_no=classified_occs$cx_int_no,reference_no=classified_occs$reference_no)
      cleaned_occs_unclass<-data.frame(collection_no=unclassified_occs$collection_no,taxon_name=unclassified_occs$taxon_name,interval_no=unclassified_occs$cx_int_no,reference_no=unclassified_occs$reference_no)
      
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
    cleaned_occs<-data.frame(collection_no=resolved_occs$collection_no,taxon_name=resolved_occs$family,interval_no=resolved_occs$cx_int_no,reference_no=resolved_occs$reference_no)
    
    #or file preparation and data cleaning for genus level analysis
  } else {
    
    #deletes occurrences not resolved to at least genus level using classified and unclassified species
    resolved_occs<-subset(resolved_occs,resolved_occs$matched_rank<=5)
    
    #deletes occurrences where genus is qualified with question mark, quotations, cf. or aff.
    resolved_occs<-subset(resolved_occs,resolved_occs$genus_reso=="" | resolved_occs$genus_reso=="n. gen.")
    
    resolved_occs$matched_name<-gsub(" .*","",resolved_occs$matched_name)
    
    #prepares a data frame with only the necessary columns
    cleaned_occs<-data.frame(collection_no=resolved_occs$collection_no,taxon_name=resolved_occs$matched_name,interval_no=resolved_occs$cx_int_no,reference_no=resolved_occs$reference_no)
    
  }
  
  #function to replace occurrence reference numbers with earliest dated collection reference
  collection.ref<-function(coll_refs) {
    pubyr<-subset(references$pubyr,references$reference_no %in% sort(coll_refs))
    which(pubyr==min(pubyr))  
  }
  
  #finds oldest publication per collection (chooses smallest reference_no if there are two or more publications from same year)
  coll_references<-lapply(sapply(split(cleaned_occs$reference_no,cleaned_occs$collection_no),unique),function(x) min(sort(x)[collection.ref(x)]))
  
  #creates collection_ref column with the oldest publication for that collection
  cleaned_occs$collection_ref<-apply(cleaned_occs,1,function(x) as.numeric(coll_references[which(names(coll_references)==as.numeric(x[1]))]))
  
  
  #GOOD'S U CALCULATION ALGORITHM
  goods.u<-function(data) {
    #Occurrences of taxa in single references
    data1p<-names(subset(sapply(split(data$collection_ref,data$taxon_name),function(x) length(unique(x))),sapply(split(data$collection_ref,data$taxon_name),function(x) length(unique(x)))==1))
    
    #Taxa in largest collection
    largest<-subset(data$taxon_name,data$collection_no==names(sort(sapply(split(data$taxon_name,data$collection_no),length),decr=T)[1]))
    
    #Single reference taxa excluding largest collection
    data1p<-subset(data1p,!data1p %in% largest)
    
    o1<-as.numeric(sort(sapply(split(data$collection_no,data$taxon_name),length),decr=T)[1]) #Occurrence count of most common taxon
    O<-nrow(data) #Total number of occurrences
    p1<-length(data1p) #Count of single reference taxa excluding largest collection
    
    #Good's u (coverage of data)
    1-p1/(O-o1)
  }
  #END GOOD'S U CALCULATION ALGORITHM
  
  
  #REFERENCE THROWBACK ALGORITHM
  #Counts of collections per reference
  ref.throw<-function(data) {
    refcount<-data.frame("Ref"=names(split(data$collection_no,data$collection_ref)),"Count"=as.numeric(sapply(split(data$collection_no,data$collection_ref),function (x) length(unique(x)))))
    
    #Assigns count of collections per reference to each collection
    dataloc_raw<-data.frame("Collection"=data$collection_no,"Refcount"=apply(data,1,function(x) subset(refcount$Count,refcount$Ref==as.numeric(x[which(names(cleaned_occs)=="collection_ref")]))))
    
    #Creates table with collection numbers and the inverse of the number of collections in the collection reference
    data.frame("collection_no"=names(sapply(split(dataloc_raw$Refcount,dataloc_raw$Collection),unique)),"inv_count"=1/sapply(split(dataloc_raw$Refcount,dataloc_raw$Collection),unique))
    #In SQS collections will be randomly sampled with a weight equal to the inverse of the number of collections described in the reference 
  }
  #END THROWBACK ALGORITHM
  
  
  #runs reference throwback to creates table of collection numbers and 1/number of references
  ref.ct.table<-lapply(split(cleaned_occs,cleaned_occs$interval_no),ref.throw)
  
  
  #BEGIN SQS CALCULATION
  sqs.calc<-function(occ_sub) {
        
    #obtains collection reference count information for pertinent time interval
    dataloc<-ref.ct.table[[paste(unique(occ_sub$interval_no))]]
    
    #randomly samples collections based on weight equal to inverse of the number of collections described in reference
    temp_loc<-sample(dataloc$collection_no,prob=dataloc$inv_count) 
    
    #calculates the proportion of occurrences belonging to each genus
    taxon_count<-data.frame(taxon_name=levels(occ_sub$taxon_name),occ_prop=sapply(split(occ_sub$taxon_name,occ_sub$taxon_name),length)/nrow(occ_sub))
    
    #normalizes occurrence frequencies to exclude the most common taxon
    taxon_count$occ_prop<-taxon_count$occ_prop/(1-max(taxon_count$occ_prop))
    
    #adjusts user-defined quorum q by Good's u value
    q_adj<-q/goods.u(occ_sub)
    
    #creates list of taxa in each collection
    coll_taxon_lists<-sapply(temp_loc,function(x) unique(subset(occ_sub$taxon_name,occ_sub$collection_no==x)))
    
    #adds first collection
    sampled_taxa<-unlist(coll_taxon_lists[1])
    
    #increments through list of collections
      for (i in seq(2,length(coll_taxon_lists))) {
      
      #creates a list of sampled taxa including the next collection
      temp_sampled_taxa<-unique(c(sampled_taxa,coll_taxon_lists[[i]]))
      
      if (max_added==0) { #quicker to do the %in% search only once and keep a counter
        
        #checks if largest taxon has been added to the list; if so, removes its occurrence frequency from the quorum counter
        if (which(taxon_count$occ_prop==max(taxon_count$occ_prop)) %in% sampled_taxa) {
           temp_running_q<-sum(taxon_count$occ_prop[temp_sampled_taxa])-max(taxon_count$occ_prop)
            max_added<-1
          } else { temp_running_q<-sum(taxon_count$occ_prop[temp_sampled_taxa]) }
        
        #if the largest taxon has been added (max_added=1), remove its occurrence frequency from the quorum counter
        } else {temp_running_q<-sum(taxon_count$occ_prop[temp_sampled_taxa])-max(taxon_count$occ_prop)
      }
      
      #if adding the collection does not overshoot the coverage target, do it
      if (temp_running_q<=q_adj) {
        
        if (length(temp_sampled_taxa)>length(sampled_taxa)) {last_coll_id<-i}
        
        sampled_taxa<-temp_sampled_taxa
        running_q<-temp_running_q
      } else if (temp_running_q-q_adj<q_adj-running_q) { #if collection overshoots but overshoot < undershoot by not adding, add it
        sampled_taxa<-temp_sampled_taxa
        running_q<-temp_running_q
      }
    }
    length(sampled_taxa)
  }

  subsamp_div<-sapply(split(cleaned_occs,cleaned_occs$interval_no),function(x) replicate(n_rep,sqs.calc(x)))

#RESULTS ORDERING
  #the intervals are in ascending numerical order by interval_no, but not necessarily chronological (this finds the chronological order)
  interval_order<-order(sapply(as.numeric(colnames(subsamp_div)),function(x) which(time_int$interval_no==x)),decreasing=T)
  
  #moves columns into chronological order
  subsamp_div<-subsamp_div[,interval_order]
   
  #renames columns with interval name
  colnames(subsamp_div)<-time_int$interval_name[sapply(colnames(subsamp_div),function(x) which(time_int$interval_no==x))]
  
  sqsdiv<<-apply(subsamp_div,2,mean)
  sqsdiv
}
