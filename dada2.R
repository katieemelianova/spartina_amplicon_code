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


errF <- learnErrors(read_1_filtered_fastq_files, multithread=TRUE)
errR <- learnErrors(read_2_filtered_fastq_files, multithread=TRUE)
plotErrors(errF, nominalQ=TRUE)

# read them in again to get existing files only, I guess some files do not pass filtering and dont get created
read_1_filtered_fastq_files <- list.files(reads_directory, pattern = "\\.1\\.filtered.fastq.gz$", full.names = TRUE)
read_2_filtered_fastq_files <- list.files(reads_directory, pattern = "\\.2\\.filtered.fastq.gz$", full.names = TRUE)

# infer sample names from filename
sample_names <- read_1_filtered_fastq_files %>% basename() %>% str_split_i("\\.", 1) %>% str_split_i("-", 4)

# dereplicate
derepFs <- derepFastq(read_1_filtered_fastq_files, verbose=TRUE)
derepRs <- derepFastq(read_2_filtered_fastq_files, verbose=TRUE)


  
  
# set sample names
names(derepFs) <- sample_names
names(derepRs) <- sample_names

# run dada
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)




# merge paired reads
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)

# make an ASV table
seqtab <- makeSequenceTable(mergers)

# remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)


getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
head(track)

nrow(track)

length(sample_names)

getSequences(seqtab.nochim)[1:5]

taxonomy_levels <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")

seqtab.nochim

start.time <- Sys.time()
taxonomy_raw <- assignTaxonomy(
  seqs = seqtab.nochim,
  refFasta = "/lisc/scratch/spartina/spartina_amplicon_code/silva_nr99_v138.2_toGenus_trainset.fa.gz",
  taxLevels = taxonomy_levels,
  multithread=12
)
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken

dim(seqtab.nochim)


write.table(taxonomy_raw, "taxonomy_raw", quote=FALSE)

