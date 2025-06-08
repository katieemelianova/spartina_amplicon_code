library(dada2)
library(dplyr)
library(magrittr)
library(readr)
library(stringr)
library(tidyverse)

convert_to_filtered_file_name <- function(x, reads_dir) {
  x %>%
    basename() %>%
    str_replace("\\.fastq$", ".filtered.fastq") %>% 
    file.path(reads_dir, .)
}

# provide fastq path, read in fastqs and set filtered filenames
reads_directory <- ("/lisc/scratch/spartina/JR_amplicons")
read_1_fastq_files <- list.files(reads_directory, pattern = "_1\\.fastq$", full.names = TRUE)
read_2_fastq_files <- list.files(reads_directory, pattern = "_2\\.fastq$", full.names = TRUE)
read_1_filtered_fastq_files <- read_1_fastq_files %>% convert_to_filtered_file_name(reads_directory)
read_2_filtered_fastq_files <- read_2_fastq_files %>% convert_to_filtered_file_name(reads_directory)



out <- filterAndTrim(read_1_fastq_files, read_1_filtered_fastq_files, read_2_fastq_files, read_2_filtered_fastq_files, truncLen=c(220,150),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=FALSE) 


# read them in again to get existing files only, I guess some files do not pass filtering and dont get created
read_1_filtered_fastq_files <- list.files(reads_directory, pattern = "_1\\.filtered.fastq$", full.names = TRUE)
read_2_filtered_fastq_files <- list.files(reads_directory, pattern = "_2\\.filtered.fastq$", full.names = TRUE)

# dereplicate
derep1s <- derepFastq(read_1_filtered_fastq_files, verbose=TRUE)
derep2s <- derepFastq(read_2_filtered_fastq_files, verbose=TRUE)

# error estimation
err1 <- learnErrors(derep1s, multithread=TRUE)
err2 <- learnErrors(derep2s, multithread=TRUE)
plotErrors(err1, nominalQ=TRUE)

# run dada inference
dada1 <- dada(derep1s, err=err1, multithread=TRUE)
dada2 <- dada(derep2s, err=err1, multithread=TRUE)

# merge paired reads
merged <- mergePairs(
  dadaF = dada1,
  dadaR = dada2,
  derepF = derep1s,
  derepR = derep2s,
  justConcatenate = TRUE)


# make an ASV table
seqtab <- makeSequenceTable(merged)

# remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)


# way to subset the seqtable
#getSequences(seqtab.nochim)[1:5]

taxonomy_levels <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")

taxa <- assignTaxonomy(
  seqs = seqtab.nochim,
  refFasta = "/lisc/scratch/spartina/spartina_amplicon_code/silva_nr99_v138.2_toGenus_trainset.fa.gz",
  taxLevels = taxonomy_levels,
  multithread=12
)

taxa <- addSpecies(taxa, "/lisc/scratch/spartina/spartina_amplicon_code/silva_v138.2_assignSpecies.fa.gz")

write.table(taxa, "taxonomy_species_america", quote=FALSE)


taxonomy_species_america



