library(tidyverse)


ViCAT_reads_mapped = "combined_readscount.tsv"
mydata <- read_tsv(ViCAT_reads_mapped, col_names = FALSE)
mydata <- mydata %>% separate(X1, into=c("sample_name", "RefID"), sep=".REF_")
mydata <- mydata %>% pivot_wider(names_from = sample_name, values_from=X2, values_fill = 0)
write.table(mydata, "ViCAT_otutable.csv", sep=",", col.names=TRUE, row.names=FALSE)

###Writes samples to file for metadata in order they are listed in mydata
samples <- colnames(mydata[,-1])
write.table(samples, "metadata.csv", sep=",", col.names=c("PlantID"), row.names=FALSE)
metaformat <- read_csv("metadata.csv")
metaformat <- cbind(metaformat, "SampleID" = 1:nrow(metaformat), "Genotype" = 1:nrow(metaformat), "Year_Sampled"= 1:nrow(metaformat), "Location"= 1:nrow(metaformat), "Other"= 1:nrow(metaformat))
write.table(metaformat, "metadata.csv", sep=",", col.names=TRUE, row.names=FALSE)