library(dada2)
library(dplyr)
library(magrittr)
library(readr)
library(stringr)
library(tidyverse)
library(phyloseq)



convert_to_filtered_file_name <- function(x, reads_dir) {
  x %>%
    basename() %>%
    str_replace("\\.fastq.gz$", ".filtered.fastq.gz") %>% 
    file.path(reads_dir, .)
}

# provide fastq path, read in fastqs and set filtered filenames
reads_directory <- ("/Users/katieemelianova/Desktop/Spartina/spartina_amplicon_code/jr_amplicons")
#reads_directory <- ("/lisc/scratch/spartina/JR_amplicons")
read_1_fastq_files <- list.files(reads_directory, pattern = "_1\\.fastq.gz$", full.names = TRUE)
read_2_fastq_files <- list.files(reads_directory, pattern = "_2\\.fastq.gz$", full.names = TRUE)
read_1_filtered_fastq_files <- read_1_fastq_files %>% convert_to_filtered_file_name(reads_directory)
read_2_filtered_fastq_files <- read_2_fastq_files %>% convert_to_filtered_file_name(reads_directory)

out <- filterAndTrim(read_1_fastq_files, read_1_filtered_fastq_files, read_2_fastq_files, read_2_filtered_fastq_files, truncLen=c(220,150),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=FALSE) 


# read them in again to get existing files only, I guess some files do not pass filtering and dont get created
read_1_filtered_fastq_files <- list.files(reads_directory, pattern = "_1\\.filtered.fastq.gz$", full.names = TRUE)
read_2_filtered_fastq_files <- list.files(reads_directory, pattern = "_2\\.filtered.fastq.gz$", full.names = TRUE)

# dereplicate
derep1s <- derepFastq(read_1_filtered_fastq_files, verbose=TRUE)
derep2s <- derepFastq(read_2_filtered_fastq_files, verbose=TRUE)

# error estimation
err1 <- learnErrors(derep1s, multithread=TRUE)
err2 <- learnErrors(derep2s, multithread=TRUE)
plotErrors(err1, nominalQ=TRUE)

# run dada inference
dada1 <- dada(derep1s, err=err1, multithread=TRUE)
dada2 <- dada(derep2s, err=err2, multithread=TRUE)

# merge paired reads
merged_jr <- mergePairs(
  dadaF = dada1,
  dadaR = dada2,
  derepF = derep1s,
  derepR = derep2s,
  justConcatenate = TRUE)


# make an ASV table
seqtab_jr <- makeSequenceTable(merged_jr)

# remove chimeras
seqtab.nochim_jr <- removeBimeraDenovo(seqtab_jr, method="consensus", multithread=TRUE, verbose=TRUE)


# way to subset the seqtable
#getSequences(seqtab.nochim)[1:5]

taxonomy_levels <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")

taxa_jr <- assignTaxonomy(
  seqs = seqtab.nochim_jr,
  refFasta = "/lisc/scratch/spartina/spartina_amplicon_code/silva_nr99_v138.2_toGenus_trainset.fa.gz",
  taxLevels = taxonomy_levels,
  multithread=12
)

taxa_jr <- addSpecies(taxa_jr, "/lisc/scratch/spartina/spartina_amplicon_code/silva_v138.2_assignSpecies.fa.gz")

write.table(taxa_jr, "taxonomy_species_america", quote=FALSE)






taxonomy_species_america <- read_delim("taxonomy_species_america") %>% set_colnames(c("sequence", "Domain", "Phylum", "Class", "Order", "Family")) %>% column_to_rownames("sequence") %>% separate(Family, sep=" ", c("Family", "Genus"))
samplesheet <- read_delim("/Users/katieemelianova/Desktop/Spartina/spartina_amplicon_code/JR_metadata.csv")

metadata_jr <- samplesheet %>% dplyr::select(Run, HOST, isolation_source) %>% 
  mutate(depth=case_when(grepl("5-10", isolation_source) ~ "5-10",
                         grepl("0-5", isolation_source) ~ "0-5")) %>%
  mutate(height=case_when(grepl("Tall", isolation_source) ~ "Tall",
                          grepl("Short", isolation_source) ~ "Short",
                          grepl("Medium", isolation_source) ~ "Medium")) %>%
  mutate(compartment=case_when(grepl("Endosphere", isolation_source) ~ "endosphere",
                               grepl("Rhizosphere", isolation_source) ~ "rhizosphere",
                               grepl("sediment", isolation_source) ~ "sediment")) %>%
  mutate(sample_id=Run) %>%
  dplyr::select(-isolation_source) %>%
  column_to_rownames("Run")
  



OTU_JR <- otu_table(seqtab.nochim_jr, taxa_are_rows=FALSE)
TAX_JR <- tax_table(taxonomy_species_america) %>% set_rownames(rownames(taxonomy_species_america)) %>% set_colnames(c("Domain", "Phylum", "Class", "Order", "Family", "Genus"))
sample_names(OTU_JR) <- str_remove(sample_names(OTU_JR), "_1.filtered.fastq.gz")

ps_jr <- phyloseq(OTU_JR, 
               sample_data(metadata_jr), 
               TAX_JR)

ps_jr <- subset_samples(ps_jr, height != "Medium")



plot_richness(ps_jr, x="height", measures=c("Shannon", "Simpson"), color="compartment") + geom_point(size = 3.5)



ps.prop_jr <- transform_sample_counts(ps_jr, function(otu) otu/sum(otu))
ord.nmds.bray_jr <- ordinate(ps.prop_jr, method="NMDS", distance="bray")
plot_ordination(ps.prop_jr, ord.nmds.bray_jr, color="height", title="Bray NMDS", shape="compartment") + 
  geom_point(size = 7) +
  theme(strip.text.x = element_text(size=25),
        axis.text.x = element_text(size=25),
        axis.text.y = element_text(size=20),
        axis.title = element_text(size=25),
        legend.title = element_text(size=20),
        legend.text = element_text(size=20))


ps_rel_jr <- transform_sample_counts(ps_jr, function(OTU) OTU/sum(OTU) )
ps_rel_jr <- subset_taxa(ps_rel_jr, Genus == "Sedimenticola") %>% subset_taxa(!(is.na(Genus))) %>% subset_taxa(!(Genus == "NA"))
## did this to remove black slices as per https://github.com/joey711/phyloseq/issues/721
phylumGlommed = tax_glom(ps_rel_jr, "Genus")
plot_bar(phylumGlommed, fill="Genus") + facet_wrap(~height+compartment, scales="free_x", ncol=2) + 
  theme(strip.text.x = element_text(size=25),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size=20),
        axis.title = element_text(size=25),
        legend.title = element_blank(),
        legend.text = element_text(size=20),
        axis.ticks.x = element_blank()) +
  ylab("Relative Abundance")




