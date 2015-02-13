taxonClean <- function(occurrences, tax_level="genus", formal_id="no") {
  
  #FILE PREPARATION AND DATA CLEANING
  #file preparation and data cleaning for species level analysis
  if(tax_level=="species")
  {
    #removes rows where species qualified by cf. or aff., question mark, ex gr. or quotation mark
    resolved_occs <- subset(occurrences, occurrences$species_reso=="" | resolved_occs$species_reso=="n. sp.")
    
    if (formal_id=="yes") {
      #deletes occurrences not resolved to species level using formally-classified species
      cleaned_occs <- subset(resolved_occs, resolved_occs$matched_rank==3)
      cleaned_occs$final_taxon <- cleaned_occs$matched_name
      
    } else {
      #deletes occurrences not resolved to species level using classified and unclassified species
      resolved_occs <- subset(resolved_occs, resolved_occs$taxon_rank=="species")
      
      #strips all extraneous qualifiers from taxon name
      resolved_occs$taxon_name <- gsub("\" ", "", resolved_occs$taxon_name)
      resolved_occs$taxon_name <- gsub("cf. ", "", resolved_occs$taxon_name)
      resolved_occs$taxon_name <- gsub("aff. ", "", resolved_occs$taxon_name)
      resolved_occs$taxon_name <- gsub("\\? ", "", resolved_occs$taxon_name)
      resolved_occs$taxon_name <- gsub("n. gen. ", "", resolved_occs$taxon_name)
      resolved_occs$taxon_name <- gsub(" n. sp.", "", resolved_occs$taxon_name)
      resolved_occs$taxon_name <- gsub(" sensu lato", "", resolved_occs$taxon_name)
      resolved_occs$taxon_name <- gsub(" informal", "", resolved_occs$taxon_name)
      
      #separates classified and unclassified occurrences
      classified_occs <- subset(resolved_occs, resolved_occs$matched_rank==3)
      unclassified_occs <- subset(resolved_occs, resolved_occs$matched_rank!=3)
      
      classified_occs$final_taxon <- classified_occs$matched_name
      unclassified_occs$final_taxon <- unclassified_occs$taxon_name
      
      #combines classified and unclassified occurrences
      cleaned_occs <- rbind(classified_occs, unclassified_occs)
    }
    
    #or file preparation and data cleaning for family level analysis
  } else if (tax_level=="family") {
    
    #deletes occurrences not resolved to at least family level
    resolved_occs <- subset(occurrences, occurrences$matched_rank <= 9)    
    
    #deletes occurrences where genus or family are qualified with question mark, quotations, cf. or aff.
    if(length(grep("\\? ", resolved_occs$taxon_name)) > 0) resolved_occs <- resolved_occs[-grep("\\? ", resolved_occs$taxon_name),]
    if(length(grep("\" ", resolved_occs$taxon_name)) > 0) resolved_occs <- resolved_occs[-grep("\" ", resolved_occs$taxon_name),]
    if(length(grep("cf. ", resolved_occs$taxon_name)) > 0) resolved_occs <- resolved_occs[-grep("cf. ", resolved_occs$taxon_name),]
    if(length(grep("aff. ", resolved_occs$taxon_name)) > 0) resolved_occs <- resolved_occs[-grep("aff. ", resolved_occs$taxon_name),]
    
    cleaned_occs <- resolved_occs
    
    cleaned_occs$final_taxon <- cleaned_occs$family
    
    #or file preparation and data cleaning for genus level analysis
  } else {
    
    #deletes occurrences not resolved to at least genus level using classified and unclassified species
    resolved_occs <- subset(occurrences, occurrences$matched_rank <= 5)
    
    #deletes occurrences where genus is qualified with question mark, quotations, cf. or aff.
    cleaned_occs <- subset(resolved_occs, resolved_occs$genus_reso=="" | resolved_occs$genus_reso=="n. gen.")
    
    #extracts genus name from matched_name string
    cleaned_occs$final_taxon <- gsub(" .*", "", cleaned_occs$matched_name)
    
  }
  
  cleaned_occs

}
