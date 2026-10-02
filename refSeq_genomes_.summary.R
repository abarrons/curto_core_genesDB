## pool together summary info for out group genomes downloded from NCBI ##
library(readr)
library(tidyverse)
setwd("~/SyncFolder/CA_curto/curto_core_genesDB/")

root_folder <- "~/SyncFolder/CA_curto/curto_core_genesDB/og_refseq_genomes/"

# Recursively search for files with .tsv extension in subfolders
file_list <- list.files(path = root_folder, pattern = "\\.tsv$", 
                        full.names = TRUE, recursive = TRUE)

# Print the list of file paths
print(file_list)



col_types <- c("character", rep("numeric", 22))

#f <- read_tsv("S.11/transposed_report.tsv", comment = "", col_names = TRUE)
col_types <- cols(
  Assembly = col_character()
)


for (file in file_list){
  
  # if the merged dataset doesn't exist, create it
  if (!exists("all.data")){
    all.data <- read_tsv(file, col_names = TRUE, comment = "")
  }
  
  # if the merged dataset does exist, append to it
  if (exists("all.data")){
    temp_dataset <- read_tsv(file, col_names = TRUE, comment = "")
    all.data<-unique(plyr::rbind.fill(all.data, temp_dataset))
    rm(temp_dataset)
  }
}

all.data <- all.data %>%
  separate(`Organism Scientific Name`, into = c("Genus", "Species"), sep = " ", extra = "merge", remove = F) %>% 
  mutate(file_name = paste(all.data$`Assembly Accession`, "_", all.data$`Assembly Name`, "_genomic.fna", sep = "")) %>%
  mutate(seq_id = gsub(".fna$", "", file_name))

mapping_data <- all.data %>%
  select(Genus, file_name, seq_id) %>%
  rename(taxon = Genus)


write.table(mapping_data, "refSeq_genomes_summary.txt", quote = F, sep = "\t", row.names = F)
