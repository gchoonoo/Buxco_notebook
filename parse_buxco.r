# Parse Buxco Data: Example August 2013

################################################

# Source functions for plotting

################################################

source('utilities_MM.R')

################################################

# Load Buxco Data

################################################

read.delim(file="./Buxco_Data/2013_iAugust - buxco.txt", sep=",", skip=1, header=T, colClasses=c(rep(NA,21),"NULL")) -> aug_2013_data

################################################

# Check Unique Sample Names

################################################

# default burn.in.lines=c("Measurement", "Create measurement", "Waiting for","Site Acknowledgement Changed")
# "Subject" and blanks are also removed.
# Note any others that do not have the form Mating RIX Virus (i.e. "Responding to")
unique(aug_2013_data[,"Subject"]) 

################################################

# Create Buxco Database

################################################

# Set the file path to the buxco data
aug_2013 <- "./Buxco_Data/2013_iAugust - buxco.txt"

# Set the file size to the number of rows in the file
chunk.size <- dim(aug_2013_data)[1]

# Set the file path of where to save the data base
db.name <- file.path("./Buxco_Data/Database/August2013_database.db")

# Parse the data, add "Responding to" in the burn.in.lines if necessary (Note: This takes a few minutes to run)
parse.buxco(file.name=aug_2013, chunk.size=chunk.size, db.name=db.name, verbose=FALSE)

# Note any parsing warnings that get printed (Sample Name and Break number), none in this case

################################################

# Add Annotation

################################################

# Read in the data base that was created
August2013_database.db <- makeBuxcoDB(db.name=file.path("August2013_database.db"))
# Add the Day and Break type level (EXP, ACC, ERR, or UNK)
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

# Create Boxplot

################################################

# Observe variables
variables(August2013_database.db)
# [1] "f"     "TVb"   "MVb"   "Penh"  "PAU"   "Rpef"  "Comp"  "PIFb" 
# [9] "PEFb"  "Ti"    "Te"    "EF50"  "Tr"    "Tbody" "Tc"    "RH"   
# [17] "Rinx" 

# Choose the variable and category to visualize
exp.penh <- retrieveData(August2013_database.db, variables="Penh", Break_type_label = 'EXP')

# Get table of categories
with(exp.penh, table(Days, Break_type_label))

# Create Boxplot
boxplot(Value~Sample_Name, data=exp.penh)

################################################

# Create Heatmap

################################################

# Output results for experiment rows only (EXP) and stratify by virus
mvtsplot(August2013_database.db, outer.group.name='Virus', break.types=c('EXP'))
