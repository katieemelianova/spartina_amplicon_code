library(dada2)
library(dplyr)
library(magrittr)
library(readr)
library(stringr)
library(tidyverse)

convert_to_filtered_file_name <- function(x, reads_dir) {
  x %>%
    basename() %>%
    str_replace("\\.fastq\\.gz$", ".filtered.fastq.gz") %>% 
    file.path(reads_dir, .)
}




#reads_directory <- ("/Users/katieemelianova/Desktop/Spartina/JMF_results/JMF-2503-07__all__rRNA_SSU_515_806__JMFR_MSRI_LWRV6")
reads_directory <- ("/lisc/scratch/spartina/Results_JMF-2503-07/JMF-2503-07__all__rRNA_SSU_515_806__JMFR_MSRI_LWRV6")


read_1_fastq_files <- list.files(reads_directory, pattern = "\\.1\\.fastq.gz$", full.names = TRUE)
read_2_fastq_files <- list.files(reads_directory, pattern = "\\.2\\.fastq.gz$", full.names = TRUE)

read_1_filtered_fastq_files <- read_1_fastq_files %>% convert_to_filtered_file_name(reads_directory)
read_2_filtered_fastq_files <- read_2_fastq_files %>% convert_to_filtered_file_name(reads_directory)



out <- filterAndTrim(read_1_fastq_files, read_1_filtered_fastq_files, read_2_fastq_files, read_2_filtered_fastq_files, truncLen=c(220,150),
              maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
              compress=TRUE, multithread=TRUE) 

# read them in again to get existing files only, I guess some files do not pass filtering and dont get created
read_1_filtered_fastq_files <- list.files(reads_directory, pattern = "\\.1\\.filtered.fastq.gz$", full.names = TRUE)
read_2_filtered_fastq_files <- list.files(reads_directory, pattern = "\\.2\\.filtered.fastq.gz$", full.names = TRUE)

# dereplicate
derep1s <- derepFastq(read_1_filtered_fastq_files, verbose=TRUE)
derep2s <- derepFastq(read_2_filtered_fastq_files, verbose=TRUE)

# split out the different files depending on which primers (?) were used
F_read_1_dereplicated <- derep1s[names(derep1s) %>% str_detect("^[^_].+\\.F\\.1\\.filtered\\.fastq\\.gz$")]
R_read_1_dereplicated <- derep1s[names(derep1s) %>% str_detect("^[^_].+\\.R\\.1\\.filtered\\.fastq\\.gz$")]
F_read_2_dereplicated <- derep2s[names(derep2s) %>% str_detect("^[^_].+\\.F\\.2\\.filtered\\.fastq\\.gz$")]
R_read_2_dereplicated <- derep2s[names(derep2s) %>% str_detect("^[^_].+\\.R\\.2\\.filtered\\.fastq\\.gz$")]

# infer sample names from filename and set sample names
names(F_read_1_dereplicated) <- F_read_1_dereplicated %>% names() %>% str_split_i("\\.", 1) %>% str_split_i("-", 4)
names(R_read_1_dereplicated) <- R_read_1_dereplicated %>% names() %>% str_split_i("\\.", 1) %>% str_split_i("-", 4)
names(F_read_2_dereplicated) <- F_read_2_dereplicated %>% names() %>% str_split_i("\\.", 1) %>% str_split_i("-", 4)
names(R_read_2_dereplicated) <- R_read_2_dereplicated %>% names() %>% str_split_i("\\.", 1) %>% str_split_i("-", 4)

# error estimation (on the unsplit reads, 1 and 2, as per JMF code also)
err1 <- learnErrors(derep1s, multithread=TRUE)
err2 <- learnErrors(derep2s, multithread=TRUE)
plotErrors(err1, nominalQ=TRUE)

# run dada inference on split reads using the error estimates from unsplit reads - check if this is correct?
dadaF1 <- dada(F_read_1_dereplicated, err=err1, multithread=TRUE)
dadaR1 <- dada(R_read_1_dereplicated, err=err1, multithread=TRUE)
dadaF2 <- dada(F_read_2_dereplicated, err=err2, multithread=TRUE)
dadaR2 <- dada(R_read_2_dereplicated, err=err2, multithread=TRUE)

# merge paired reads
merged <- mergePairs(
  dadaF = c(dadaF1, dadaR1),
  dadaR = c(dadaF2, dadaR2),
  derepF = c(F_read_1_dereplicated, R_read_1_dereplicated),
  derepR = c(F_read_2_dereplicated, R_read_2_dereplicated))


length(merged)
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

write.table(taxa, "taxonomy_species", quote=FALSE)



