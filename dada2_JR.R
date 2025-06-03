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


