# Parse Buxco Data: Example August 2013

################################################

# Source functions for getting plot data frame

################################################

source('utilities_MM.R')

################################################

# Load Buxco Data and clean

################################################

# Set file path of data
aug_2013 = '2013_iAugust - buxco.txt'

# Read in data
aug_2013 = read.delim(aug_2013, sep=',', header=T, as.is=T, skip=1)

# Save number of columns
length(names(aug_2013)) -> output_cols

# Set data header and footer
table_header='Table WBPth'
table_footer='Table Metabolism'
table_footer_cols = as.character(tail(aug_2013, 1))[!is.na(tail(aug_2013, 1))]

# Subset to exclude last two rows 
aug_2013 = aug_2013[1:(dim(aug_2013)[1]-2), ]

# Annotate Mating, RIX, Virus, and Sex
aug_2013$Mating = sapply(strsplit(aug_2013$Subject, ' +'), '[', 1)
aug_2013$RIX_ID = sapply(strsplit(aug_2013$Subject, ' +'), '[', 2)
aug_2013$Virus = sapply(strsplit(aug_2013$Subject, ' +'), '[', 3)
aug_2013$Sex = NA
aug_2013$Sex[grepl('f|F', aug_2013$RIX_ID)] = 'F'
aug_2013$Sex[grepl('m|M', aug_2013$RIX_ID)] = 'M' 
aug_2013$RIX_ID = gsub('f|m', '', aug_2013$RIX_ID)

# Check if Virus is annotated correctly
unique(aug_2013[grep("x",aug_2013[,"Subject"]),"Virus"])

# Code chunk to annotate virus based on weight sheets
# aug_2013$Virus2 = NA
# for(i in 1:dim(virus_annot)[1]){
#   print(i)
#   aug_2013[paste(aug_2013[,"Mating"],aug_2013[,"RIX_ID"],sep="_") %in% virus_annot[i,1],"Virus2"] <- as.character(virus_annot[i,2])
# }
# 
# aug_2013[which(is.na(aug_2013[,"Virus"])),"Virus"] <- aug_2013[which(is.na(aug_2013[,"Virus"])),"Virus2"]

# Standardize Virus column: SARS, FLU, Mock
aug_2013[which(aug_2013[,"Virus"] == "sars"),"Virus"] <- "SARS"
aug_2013[which(aug_2013[,"Virus"] == "flu"),"Virus"] <- "FLU"
aug_2013[which(aug_2013[,"Virus"] == "mock"),"Virus"] <- "Mock"

# Clean sample names
aug_2013$Subject_cleaned = aug_2013$Subject
aug_2013[grep("x",aug_2013[,"Subject"]),"Subject_cleaned"] = paste(aug_2013[grep("x",aug_2013[,"Subject"]),"Mating"], paste0(tolower(aug_2013[grep("x",aug_2013[,"Subject"]),"Sex"]), aug_2013[grep("x",aug_2013[,"Subject"]),"RIX_ID"]), aug_2013[grep("x",aug_2013[,"Subject"]),"Virus"], sep=' ')
aug_2013$Subject = aug_2013$Subject_cleaned

# Check Unique Sample Names
# Default rows that are removed in parsing process: burn.in.lines=c("Measurement", "Create measurement", "Waiting for","Site Acknowledgement Changed")
# "Subject" and blanks are also removed.
# Note any others that do not have the form Mating RIX Virus (i.e. "Responding to")
unique(aug_2013$Subject)

# remove variable RH, EXP, value >= 100
# remove variabe Tc, EXP, negative values and -50
if(length(which(aug_2013[,"RH"] >=100)) > 0){
  aug_2013[-which(aug_2013[,"RH"] >= 100),] -> aug_2013
}

if(length(which(aug_2013[,"Tc"] < 0)) > 0){
  aug_2013[-which(aug_2013[,"Tc"] < 0),] -> aug_2013
}

# Output in Buxco format
colnames(aug_2013)[output_cols] = ''
cleaned_file = '2013_iAugust - buxco_cleaned.txt'
cat(table_header, file=cleaned_file, sep='\n')
write.table(aug_2013[,1:output_cols], file=cleaned_file, append=T, sep=',', col.names=T, row.names=F, quote=F, na='')
cat(table_footer, file=cleaned_file, append=T, sep='\n')
cat(table_footer_cols, file=cleaned_file, append=T, sep=',')
cat(',\n', file=cleaned_file, append=T)

################################################

# Create Buxco Database

################################################

# Set the file path to the cleaned buxco data
aug_2013_path <- "2013_iAugust - buxco_cleaned.txt"

# Set the file size to the number of rows in the file
chunk.size <- dim(aug_2013)[1]

# Set the file path of where to save the data base
db.name <- file.path("August2013_database.db")

# Parse the data, add "Responding to" in the burn.in.lines if necessary (Note: This takes a few minutes to run)
parse.buxco(file.name=aug_2013_path, chunk.size=chunk.size, db.name=db.name, verbose=FALSE)

# Note any parsing warnings that get printed (Sample Name and Break number), none in this case

################################################

# Add Annotation

################################################

# Read in the data base that was created
August2013_database.db <- makeBuxcoDB(db.name=file.path("August2013_database.db"))

# Add the Virus, Day and Break type level (EXP, ACC, ERR, or UNK)
addAnnotation(August2013_database.db, query=virus.query, index=TRUE)  
addAnnotation(August2013_database.db, query=day.infer.query, index=TRUE)  
addAnnotation(August2013_database.db, query=break.type.query, index=TRUE)

################################################

# Save parsing warnings, error, and unknown rows

################################################

# Check Break type levels
annoLevels(August2013_database.db)

# Subset the parsing warning and error rows in each file
# Example:
acc.aug2013 <- retrieveData(August2013_database.db, variables=variables(August2013_database.db), Break_type_label = 'ACC')

# acc.feb2013[which(acc.feb2013[,1] == "16513x16188 f105 FLU" & acc.feb2013[,"Break_number"] == 158),] -> acc.feb2013_warning

################################################

# Buxco Annotation Workflow

################################################

# Source Functions
source("buxco_annotation_functions_MM_GC_6-JUN-2016.R")

# Set working directory to file that contains Buxco databases
db_dir = "Database"

# Initialize empty list to the number of databases in the directory
buxco_annot_penh_log = vector("list",length(dir(db_dir)))

# Retrieve data from all batches and combine together with means, use specific transformation for each variable (i.e. Penh: log transform)
combine_data(data=buxco_annot_penh_log, FUN=log, variables="Penh",db_dir=db_dir) -> buxco_annot_penh_data_log

# Get mock means and AUC
mock_mean(data=buxco_annot_penh_data_log) -> buxco_annot_penh_data_log_mock

# Save annotation file
# write.table(file="./buxco_annotation_penh_log.txt",x=buxco_annot_penh_data_log_mock, sep="\t", row.names=F, quote=F, na="")

# View first few columns of annotation file
head(buxco_annot_penh_data_log_mock)
