library(tidyr)
library(dplyr)
library(synapser)
library(jsonlite)

synLogin()

#################################
# Load data
#################################

# load the file view
# I manually added "name" to the fileview in the web UI
queryResults <- synTableQuery("select * from syn72030052")
data <- as.data.frame(queryResults)
colnames(data)[28] = 'filename'

# load nuc hashing metadata
file.metadata.path = synGet('syn53446426')$path
file.metadata = read.csv(file.metadata.path)
head(file.metadata)

# load the individual metadata file
ind.metadata.path = synGet('syn55251012')$path
ind.metadata = read.csv(ind.metadata.path)
head(ind.metadata)

#################################
# IndividualtoFileMap table
#################################

# Pivot longer to stack all the fastq filename columns
df_long <- file.metadata %>%
  pivot_longer(
    cols = starts_with("pool"),  # Select all columns starting with "pool"
    names_to = "pool_info",       # Column name for the original column names
    values_to = "fastq_filename"  # Column name for the fastq filenames
  ) %>%
  select(IndividualID, fastq_filename) %>%  # Keep only the columns you want
  filter(!is.na(fastq_filename) & fastq_filename != "")  # Remove empty values if any

# View the result
head(df_long)
colnames(df_long)[2] = 'filename'

# Adding synID gives us the individualToFileMap (id = synid, individualId)
indivtofile  <- merge(df_long, data[, c("filename", "id")], by = "filename", all.x = TRUE)
indivtofile2 = indivtofile[c('id', 'IndividualID')]
write.csv(indivtofile2,'~/cohort_builder/nps_ad/driver_tables/output/npsad_individualtoFileMap_nometadata.csv')

#manually downloaded and uploaded to synapse

# TODO: Add the metadata files to this as well - going to exclude for now as they are excluded in ELITE implementation

missingfromindtofile = data[which(!(data$id %in% indivtofile2$id )),] # all metadata, plus 1 raw seq file, 1 processed file

#################################
# Backend Individuals table
#################################
# need to think about what things folks would want to search on for cohort builder

backend_individuals = ind.metadata
which(is.na(backend_individuals$individualID)) #none
which(is.na(backend_individuals$cohort)) #none
which(is.na(backend_individuals$species)) #none
which(is.na(backend_individuals$sex)) #none
which(is.na(backend_individuals$race)) #none
which(is.na(backend_individuals$apoeGenotype)) #none
length(which(is.na(backend_individuals$ageDeath))) # 157 samples have no ageDeath
which(is.na(backend_individuals$Braak)) #none
which(is.na(backend_individuals$apoe4Status)) #none
which(is.na(backend_individuals$dataContributionGroup)) #none
which(is.na(backend_individuals$study)) #none


backend_individuals = ind.metadata[c('individualID', 'cohort',
                                     'species', 'sex', 'race',
                                     'apoeGenotype', 'ageDeath',
                                     'Braak', 'apoe4Status',
                                     'dataContributionGroup','study')]

# missing individuals
missing.ind = backend_individuals[which(!(backend_individuals$individualID %in% indivtofile2$IndividualID)),] #are there ind in backend_ind not in the indivtofilemap, yes 22
which(!(missing.ind$individualID %in% file.metadata$IndividualID)) # all of them

backend_individuals = backend_individuals %>% filter(!(backend_individuals$individualID %in% missing.ind$individualID))
intersect(backend_individuals$individualID,missing.ind$individualID)

write.csv(backend_individuals,'~/cohort_builder/nps_ad/driver_tables/output/npsad_backend_individuals.csv')

#manually downloaded and uploaded to synapse


###############################################

# This will help you update the Fileview with the individualIDs, need to do the same for sex

indivtofile2$IndividualID = gsub( " ", "", indivtofile2$IndividualID) 

df_grouped <- indivtofile2 %>%
  group_by(id) %>%
  summarise(IndividualID = list(IndividualID)) %>%
  ungroup()

# View the result
head(df_grouped)
df_grouped$IndividualID <- sapply(df_grouped$IndividualID, function(x) rjson::toJSON(x))

colnames(df_grouped)[2] = 'individualID'

data2  <- merge(data, df_grouped[, c("individualID", "id")], by = "id", all.x = TRUE)
data2$individualID.x = NULL
colnames(data2)[28] = "individualID"

data2 <- data2 %>% 
  select("ROW_ID","ROW_VERSION","ROW_ETAG", "id","individualID", everything())

colnames(data2)[28] = 'name'

synStore(Table("syn72030052", data2)) # this worked, but I needed to change the char size limit for individualID to 150 manually in the web UI

