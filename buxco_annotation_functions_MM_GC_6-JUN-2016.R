library(plethy)
library(plyr)
library(flux)

# function to combine data
combine_data = function(data=buxco_annot_penh, FUN=log, variables="Penh", db_dir){
  for(i in 1:length(dir(db_dir))){
    #print(i)
    data[[i]] <- get_db(i, FUN=FUN, variables=variables, db_dir)
    if(nrow(data[[i]]) == 0){
      
    }else{
      data[[i]]$Date = paste(sapply(strsplit(dir(db_dir)[i],"_"),"[",1),sapply(strsplit(dir(db_dir)[i],"_"),"[",2),sep="_")
      
    }
  }
  return(do.call(rbind,data))
}

# function to read in data base and format data
get_db = function(x=1, FUN=log, variables="Penh", db_dir){ 
  db <- makeBuxcoDB(db.name=file.path(db_dir, dir(db_dir)[x]))
  #print(x)
  dat <- retrieveData(db, variables=variables)
  format_data(dat, FUN=FUN) -> dat
  return(dat)
}

# retrieve data from each buxco batch and get means
format_data = function(data, FUN=log){
  
  # summarize over Break type labels == EXP
  if(nrow(data) == 0){
    return(data)
  }else{
    data[which(data[,"Break_type_label"] == "EXP"),] -> data_v2
    
    # add Days_PI column (calculated for each sample name, Days, and Rec_Exp_date)
    data_v2$Days_PI = NA
    
    # Set day annotation column
    #timepoint_name = names(data_v2)[which(sapply(1:ncol(data_v2),function(x)sum(grep("D|day",data_v2[,x])) > 0))]
    timepoint_name = "Rec_Exp_date"
    
    for (s in unique(data_v2[, "Sample_Name"])) {
      #print(s)
      min_tp = min(as.numeric(as.character(gsub("D|day ","",data_v2[data_v2$Sample_Name == s, timepoint_name]))))
      data_v2[data_v2$Sample_Name == s,"Days_PI"] <- data_v2[data_v2$Sample_Name == s, "Days"] + min_tp
    }
    
    # Format sample names
    data_v2[,"Sample_Name"] <- gsub("  "," ",data_v2[,"Sample_Name"])
    data_v2[,"Sample_Name"] <- gsub(" ","_",data_v2[,"Sample_Name"])
    
    data_v3 <- data_v2
    
    ## Set blank values to NA
    data_v3[which(data_v3[,"Value"] == ""),"Value"] <- NA
    
    # Get values in heatmap (mean of transformed values for each animal)
    # infinite values are ignored (e.g. log(0))
    data_v4 <- ddply(.data=data_v3, .variables=c("Days_PI", "Sample_Name"), .fun=function(x) data.frame(mean_per_day=mean(FUN(as.numeric(as.character(x$Value)))[is.finite(FUN(as.numeric(as.character(x$Value))))], na.rm=T)))
    
    data_v3[which(data_v3[,1] == "3015x5306_f91_SARS"),]
    ## ** please use column names here **
    #data_v3[, c(1,4,6,8,9,10,12)] -> data_v3_full
    cols_to_keep = c("Sample_Name", "Variable_Name", "Rec_Exp_date", "Virus", "Days", 
                     "Break_type_label", "Days_PI")
    data_v3[, cols_to_keep] -> data_v3_full
    
    # Remove duplicates, so there is only one row per sample/day
    data_v3_full[!duplicated(data_v3_full),] -> data_v3_full
    
    merge(data_v3_full, data_v4, by=c("Days_PI","Sample_Name"), all.x=T) -> data_v4
    
    # median values in heatmap (median of mean of log values for each animal)
    data_v4_med_heat <- ddply(.data=data_v4, .variables=c("Sample_Name"), .fun=function(x) data.frame(median_mean_over_all_days=median(x$mean_per_day, na.rm=T)))
      
    # heatmap: median of transformed values for each day per virus
    virus_day_median <- ddply(.data=data_v3, .variables=c("Days_PI",'Virus'), .fun=function(x) data.frame(virus_median_per_day=median(FUN(as.numeric(as.character(x$Value)))[is.finite(FUN(as.numeric(as.character(x$Value))))],na.rm=T)))
          
    merge(data_v4, data_v4_med_heat[,c("Sample_Name","median_mean_over_all_days")], by="Sample_Name") -> data_v5

    data_v5$Mating = sapply(strsplit(data_v5$Sample_Name,"_"),"[",1)
    
    data_v5$RIX_ID = gsub("f|NA","",sapply(strsplit(data_v5$Sample_Name,"_"),"[",2))
    
    data_v5$ID = paste(data_v5$Mating, data_v5$RIX_ID, sep='_')
    
    data_v5$Sex = 'F'
    
    data_v5$Sample_Name = paste0(data_v5$Mating, '_', tolower(data_v5$Sex), data_v5$RIX_ID, '_', data_v5$Virus)

    # add virus median per day
    merge(data_v5, virus_day_median, by=c("Days_PI","Virus")) -> data_v6
    
    # order
    data_v6[,c("Sample_Name","ID","Mating","Sex","RIX_ID","Virus","Rec_Exp_date","Days","Days_PI","Variable_Name","Break_type_label","mean_per_day","median_mean_over_all_days","virus_median_per_day")] -> data_v7

    data_v7[order(data_v7[,"ID"], data_v7[,"Days"]),] -> data_v8
      
    ## Annotate break type labels for each sample
    unique(data[,c("Sample_Name","Break_type_label")]) -> label_id
  
    label_id$Mating = sapply(strsplit(label_id$Sample_Name," "),"[",1)
    
    label_id$RIX_ID = gsub("f|NA","",sapply(strsplit(label_id$Sample_Name," "),"[",2))
        
    label_id$Sex = 'F'
    
    label_id$Virus = sapply(strsplit(label_id$Sample_Name," "),"[",3)
    
    label_id$Sample_Name = paste0(label_id$Mating, '_', tolower(label_id$Sex), label_id$RIX_ID, '_', label_id$Virus)
    
    label_id$Break_type_label_all = NA
    
    for(i in 1:length(unique(label_id[,1]))){
      label_id[label_id[,1]%in%unique(label_id[,1])[i],"Break_type_label_all"] <- paste(label_id[label_id[,1]%in%unique(label_id[,1])[i],2],sep=",",collapse=",")
    }
    
    unique(label_id[,c("Sample_Name","Break_type_label_all")]) -> label_id_all
        
    data_v8$Break_type_label_all = NA
    
    for(i in 1:dim(label_id_all)[1]){
      #print(i)
      data_v8[which(data_v8[,"Sample_Name"] %in% label_id_all[i,1]),"Break_type_label_all"] <- as.character(label_id_all[i,2])
    }
    
    stopifnot(sum(is.na(data_v8[,"Break_type_label_all"]))==0)
    
    return(data_v8)
  }
}

# get mock averages and AUC
mock_mean = function(data){
  
  mock_average <- ddply(.data=data, .variables=c("Mating","Virus","Days_PI","Date"), .fun=function(x) data.frame(mock_mean_per_day=mean(x$mean_per_day)))
  
  mock_sd <- ddply(.data=data, .variables=c("Mating","Virus","Days_PI","Date"), .fun=function(x) data.frame(mock_sd_per_day=sd(x$mean_per_day)))
  
  mock_average[which(mock_average[,2] == "Mock"),] -> mock_average_v2
  mock_sd[which(mock_sd[,2] == "Mock"),] -> mock_sd_v2
  
  cbind(mock_average_v2, mock_sd_per_day=mock_sd_v2[,5]) -> mock_dist
  
  merge(data, mock_dist[,-2], by=c("Mating","Days_PI","Date"), all.x=T) -> data_v2
  
  names(mock_average)[5] <- "line_virus_mean_per_day"
  
  merge(data_v2, mock_average, by=c("Mating","Virus","Days_PI","Date"), all.x=T) -> data_v2
  
  ## Calculate difference between each animal's value and the same-day mock average
  data_v2$mean_diff_infected_mock_per_day = data_v2[,"mean_per_day"] - data_v2[,"mock_mean_per_day"]
  
  # add auc columns
  data_v2$AUC = NA
  
  
  # added in if statement, values have to be non-empty
  for(i in 1:length(unique(data_v2[,"ID"]))){
    #print(i)
    if(length(na.omit(data_v2[which(data_v2[,"ID"] == unique(data_v2[,"ID"])[i]),"mean_per_day"])) > 1){

      data_v2[which(data_v2[,"ID"] == unique(data_v2[,"ID"])[i]),"AUC"] <- auc(data_v2[which(data_v2[,"ID"] == unique(data_v2[,"ID"])[i]),"Days_PI"],data_v2[which(data_v2[,"ID"] == unique(data_v2[,"ID"])[i]),"mean_per_day"])
      
    }else{
      data_v2[which(data_v2[,"ID"] == unique(data_v2[,"ID"])[i]),"AUC"] <- NA
    }
  }
    
  # get mean AUC for mock (updated code MM)
  mock_auc <- ddply(.data=data_v2, .variables=c("Mating","Virus","Date"), .fun=function(x) data.frame(AUC_mock_mean=mean(as.numeric(as.character(x$AUC[!duplicated(x$ID)])), na.rm=T)))
  
  mock_auc[which(mock_auc[,"Virus"] == "Mock"),] -> mock_auc_v2
  
  merge(data_v2, mock_auc_v2[,-2], by=c("Mating","Date"),all.x=T) -> data_v3
  
  # get diff of AUC for infected vs mock
  data_v3$AUC_diff_infected_mock = data_v3$AUC - data_v3$AUC_mock_mean
  
  data_v3[order(data_v3[,"ID"],data_v3[,"Days_PI"]),] -> data_v3
  
  # add AUC per day per animal
  data_v3$AUC_per_day = NA
  
  for(j in 1:length(unique(data_v3[,"ID"]))){
    #print(unique(data_v3[,"ID"])[j])
    for(k in 1:length(unique(data_v3[,"Variable_Name"]))){
      #print(unique(data_v3[,"Variable_Name"])[k])
      if(sum(data_v3[,"ID"] %in% unique(data_v3[,"ID"])[j] & data_v3[,"Variable_Name"] %in% unique(data_v3[,"Variable_Name"])[k] & !is.na(data_v3[,"mean_per_day"])) >1){
        data_v3[data_v3[,"ID"] %in% unique(data_v3[,"ID"])[j] & data_v3[,"Variable_Name"] %in% unique(data_v3[,"Variable_Name"])[k] & !is.na(data_v3[,"mean_per_day"]),c("Days_PI","mean_per_day")] -> test_data
        
        auc_val = vector("list",nrow(test_data)-1)
        for(i in 1:(nrow(test_data)-1)){
          #print(i)
          auc_val[i] = auc(test_data[seq(i,i+1),1],test_data[seq(i,i+1),2])
        }
        
        data_v3[data_v3[,"ID"] %in% unique(data_v3[,"ID"])[j] & data_v3[,"Variable_Name"] %in% unique(data_v3[,"Variable_Name"])[k] & !is.na(data_v3[,"mean_per_day"]),c("AUC_per_day")] <- c(NA,unlist(auc_val))
      }
    }
  }
  
  return(data_v3)
}

inverse = function(x) 1/x
