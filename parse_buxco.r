# Parse Buxco Data: Example August 2013

# Load Buxco Data
read.delim(file="./Buxco_Data/2013_iAugust - buxco.txt", sep=",", skip=1, header=T, colClasses=c(rep(NA,21),"NULL")) -> aug_2013_data

# Check Unique Sample Names
# default burn.in.lines=c("Measurement", "Create measurement", "Waiting for","Site Acknowledgement Changed")
# "Subject" and blanks are also removed.
# Note any others that do not have the form Mating RIX Virus (i.e. "Responding to")
unique(aug_2013_data[,"Subject"]) 

# Create Buxco Database
# Set the file path to the buxco data
aug_2013 <- "./Buxco_Data/2013_iAugust - buxco.txt"
# Set the file size the number of rows in the file
chunk.size <- dim(aug_2013_data)[1]
# Set the file path of where to save the data base
db.name <- file.path("./Buxco_Data/Database/August2013_database.db")
# Parse the data, add "Responding to" in the burn.in.lines if necessary
parse.buxco(file.name=aug_2013, chunk.size=chunk.size, db.name=db.name, verbose=FALSE)
# Note any parsing warnings that get printed (Sample Name and Break number), none in this case

# Add Annotation
# Read in the data base that was created
August2013_database.db <- makeBuxcoDB(db.name=file.path("./Buxco_Data/Database/August2013_database.db"))
# Add the Day and Break type level (EXP, ACC, ERR, or UNK)
addAnnotation(August2013_database.db, query=day.infer.query, index=TRUE)  
addAnnotation(August2013_database.db, query=break.type.query, index=TRUE)

# Check Break type levels
annoLevels(August2013_database.db)
# $Days
# [1]  0  4  5  6  7  9 11 14 19 26
# 
# $Break_type_label
# [1] "ACC" "ERR" "EXP"

# Save parsing warnings, error, and unknown rows
acc.aug2013 <- retrieveData(August2013_database.db, variables=variables(August2013_database.db), Break_type_label = 'ACC')

exp.aug2013 <- retrieveData(August2013_database.db, variables=variables(August2013_database.db), Break_type_label = 'EXP')

err.aug2013 <- retrieveData(August2013_database.db, variables=variables(August2013_database.db), Break_type_label = 'ERR')

# save error rows
#write.table(file="./Buxco_parsing_warnings/aug_2013_error_rows.txt", x=err.aug2013, sep="\t", quote=F, row.names=F)

# Subset the parsing warning rows in each file
# Example:
# acc.feb2013[which(acc.feb2013[,1] == "16513x16188 f105 FLU" & acc.feb2013[,"Break_number"] == 158),] -> acc.feb2013_warning
