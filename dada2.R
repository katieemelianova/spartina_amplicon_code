library(dada2)
library(dplyr)
library(magrittr)
library(readr)
library(stringr)
library(tidyverse)
library(msa)
library(seqinr)

convert_to_filtered_file_name <- function(x, reads_dir) {
  x %>%
    basename() %>%
    str_replace("\\.fastq\\.gz$", ".filtered.fastq.gz") %>% 
    file.path(reads_dir, .)
}




reads_directory <- ("/Users/katieemelianova/Desktop/Spartina/JMF_results/JMF-2503-07__all__rRNA_SSU_515_806__JMFR_MSRI_LWRV6")
#reads_directory <- ("/lisc/scratch/spartina/Results_JMF-2503-07/JMF-2503-07__all__rRNA_SSU_515_806__JMFR_MSRI_LWRV6")


read_1_fastq_files <- list.files(reads_directory, pattern = "\\.1\\.fastq.gz$", full.names = TRUE)
read_2_fastq_files <- list.files(reads_directory, pattern = "\\.2\\.fastq.gz$", full.names = TRUE)

read_1_filtered_fastq_files <- read_1_fastq_files %>% convert_to_filtered_file_name(reads_directory)
read_2_filtered_fastq_files <- read_2_fastq_files %>% convert_to_filtered_file_name(reads_directory)



out <- filterAndTrim(read_1_fastq_files, read_1_filtered_fastq_files, read_2_fastq_files, read_2_filtered_fastq_files, truncLen=c(220,150),
              maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
              compress=TRUE, multithread=FALSE) 

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
err1 <- learnErrors(derep1s, multithread=FALSE)
err2 <- learnErrors(derep2s, multithread=FALSE)
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
  derepR = c(F_read_2_dereplicated, R_read_2_dereplicated),
  justConcatenate = TRUE)



matrix_to_data_table <- function(x, names) {
  x <- as.data.table(x, keep.rownames = TRUE)
  x <- melt(x, id.vars = "rn", variable.factor = FALSE)
  setnames(x, names)
  x
}

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

write.table(taxa, "taxonomy_species_europe", quote=FALSE)





library(phyloseq)
taxonomy_species_europe <- read_delim("taxonomy_species_europe") %>% set_colnames(c("sequence", "Domain", "Phylum", "Class", "Order", "Family", "Genus")) %>% column_to_rownames("sequence")

samplesheet <- readxl::read_xlsx("/Users/katieemelianova/Desktop/Spartina/JMF_results/Results_2025-04-30/JMF-2503-07.xlsx")[1:80,]
samplesheet_marshinfo <- readxl::read_xlsx("/Users/katieemelianova/Desktop/Spartina/Sequencing/March2025_LeFaouAmplicon/LeFaou_samplesheet.xlsx")

# !!! I tried to set the sample ID as the sample number e.e.g 0001 etc and this led me on a wild goose chase when the plot_bar failed to work
# turns out phyloseq doesnt like sample names which are integers, writing this down for posterity
sample_metadata <- inner_join(samplesheet_marshinfo, samplesheet, by="User sample ID") %>% 
  #mutate(sample_id=str_split_i(`JMF sample ID.y`, "-", 4), 
         mutate(sample_id=`JMF sample ID.y` %>% str_replace_all("-", "."), 
         compartment=case_when(`Sample description.y` == "soil from saltmarsch" ~ "Sediment", 
                               `Sample description.y` == "spartina roots" ~ "Root"),
         elevation = case_when(`Sample description.x` == "Sediment marsh" ~ "Low Elevation",
                               `Sample description.x` == "Sediment dry" ~ "High Elevation",
                               `Sample description.x` == "Root marsh" ~ "Low Elevation",
                               `Sample description.x` == "Root dry" ~ "High Elevation",
                               `Sample description.x` %in% c("Sedment unknown", "Root unknown") ~ "remove"),
         compartment = case_when(compartment == "Sediment" ~"Rhizosphere",
                                 compartment == "Root" ~"Endosphere")) %>% 
  dplyr::select(elevation, `concentration (ng/µL) (NanoDrop).x`, compartment, `User sample ID`, sample_id) %>% 
  set_colnames(c("elevation", "concentration", "compartment", "user_id", "sample_id"))
sample_metadata %<>% data.frame() %>% set_rownames(sample_metadata$sample_id)

rownames(seqtab.nochim) <- rownames(seqtab.nochim) %>% str_remove("[a-zA-Z]$")



# make a phyloseq object 

ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=TRUE), 
               sample_data(sample_metadata), 
               tax_table(taxonomy_species_europe))
               
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=TRUE), 
               sample_data(sample_metadata), 
               tax_table(taxonomy_species_europe) %>% set_rownames(rownames(taxonomy_species_europe)) %>% set_colnames(c("Domain", "Phylum", "Class", "Order", "Family", "Genus")))


test <- tax_table(taxonomy_species_europe) %>% set_rownames(rownames(taxonomy_species_europe)) %>% set_colnames(c("Domain", "Phylum", "Class", "Order", "Family", "Genus"))

taxa_names(test)
taxonomy_species_europe %>% set_rownames(rownames(taxonomy_species_europe)) %>% set_colnames(c("Domain", "Phylum", "Class", "Order", "Family", "Genus")) %>% dplyr::select(Genus)

# remove samples with unknown marsh elevation
ps <- subset_samples(ps, elevation != "remove")

# convert concentration tonumeric from character
ps@sam_data$concentration <- as.numeric(ps@sam_data$concentration)

# I get a warning about singletons with this command but see https://github.com/benjjneb/dada2/issues/214
# here I plot diversity vs concentration, coloured by compartment and eleveation
plot_richness(ps, x="concentration", measures=c("Shannon", "Simpson"), color="elevation", shape="compartment") + geom_point(size = 3.5)

ps.prop <- transform_sample_counts(ps, function(otu) otu/sum(otu))
ord.nmds.bray <- ordinate(ps.prop, method="NMDS", distance="bray")
plot_ordination(ps.prop, ord.nmds.bray, color="elevation", title="Bray NMDS", shape="compartment") + geom_point(size = 7)




##### giving up on raw data plottingbecause I think the samples need further library F/R merging 
##### need to ask JMF how to merge without using sequence IDs (I want to keep the ASV seqs as seqIDs)
##### running analysis plots using the final datasets provided by JMF for now

otumat <- read.table("/Users/katieemelianova/Desktop/Spartina/JMF_results/Results_2025-04-30/JMF-2503-07__all__rRNA_SSU_515_806__JMFR_MSRI_LWRV6/DADA2_counts_as_matrix.tsv", header=TRUE) %>%
  column_to_rownames(var="Sequence_ID")
taxmat <- read_delim("/Users/katieemelianova/Desktop/Spartina/JMF_results/Results_2025-04-30/JMF-2503-07__all__rRNA_SSU_515_806__JMFR_MSRI_LWRV6/DADA2_ASVs.rRNA_SSU.SILVA_reference.DADA2_classified.tsv") %>%
  column_to_rownames(var="Sequence_ID")

OTU = otu_table(otumat, taxa_are_rows = TRUE)
sample_names(OTU) %<>% str_remove("A")
TAX = tax_table(taxmat)
rownames(TAX) <- rownames(taxmat)
colnames(TAX) <- colnames(taxmat)

ps <- phyloseq(OTU, 
               sample_data(sample_metadata), 
               TAX)



# remove samples with unknown marsh elevation
ps <- subset_samples(ps, elevation != "remove")

# convert concentration tonumeric from character
ps@sam_data$concentration <- as.numeric(ps@sam_data$concentration)


# I get a warning about singletons with this command but see https://github.com/benjjneb/dada2/issues/214
# here I plot diversity vs concentration, coloured by compartment and eleveation
plot_richness(ps, x="concentration", measures=c("Shannon", "Simpson"), color="elevation", shape="compartment") + geom_point(size = 3.5)

ps.prop <- transform_sample_counts(ps, function(otu) otu/sum(otu))
ord.nmds.bray <- ordinate(ps.prop, method="NMDS", distance="bray")
plot_ordination(ps.prop, ord.nmds.bray, color="elevation", title="Bray NMDS", shape="compartment") + 
  geom_point(size = 7) +
  theme(strip.text.x = element_text(size=25),
        axis.text.x = element_text(size=25),
        axis.text.y = element_text(size=20),
        axis.title = element_text(size=25),
        legend.title = element_text(size=20),
        legend.text = element_text(size=20))




ps_rel <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU) )


barplot_thiodiazotropha <- subset_taxa(ps_rel, Genus == "Candidatus Thiodiazotropha") %>% subset_taxa(!(is.na(Genus))) %>% 
  tax_glom("Genus") %>%
  plot_bar(fill="Genus") + facet_wrap(~elevation+compartment, scales="free_x", ncol=2) + 
  theme(strip.text.x = element_text(size=20),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size=20),
        axis.title = element_text(size=25),
        legend.title = element_blank(),
        legend.text = element_text(size=20),
        axis.ticks.x = element_blank()) +
  ylab("Relative Abundance")

barplot_sedimenticola <- subset_taxa(ps_rel, Genus == "Sedimenticola") %>% subset_taxa(!(is.na(Genus))) %>% 
  tax_glom("Genus") %>%
  plot_bar(fill="Genus") + facet_wrap(~elevation+compartment, scales="free_x", ncol=2) + 
  theme(strip.text.x = element_text(size=20),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size=20),
        axis.title = element_text(size=25),
        legend.title = element_blank(),
        legend.text = element_text(size=20),
        axis.ticks.x = element_blank()) +
  ylab("Relative Abundance")
  








#################################################################
#               get sequences of Sedimenticola                  #
#################################################################

asv_mapping <- read_tsv("/Users/katieemelianova/Desktop/Spartina/JMF_results/Results_2025-04-30/JMF-2503-07__all__rRNA_SSU_515_806__JMFR_MSRI_LWRV6/DADA2_results.tsv")
sedimenticola_asv <- subset_taxa(ps, Genus == "Sedimenticola")
thiodiazotropha_asv <- subset_taxa(ps, Genus == "Candidatus Thiodiazotropha")
sedimenticola_asv_seqs <- asv_mapping %>% filter(Sequence_ID %in% rownames(sedimenticola_asv@otu_table)) %>% pull(Sequence) %>% unique()
thiodiazotropha_asv_seqs <- asv_mapping %>% filter(Sequence_ID %in% rownames(thiodiazotropha_asv@otu_table)) %>% pull(Sequence) %>% unique()


write.fasta(as.list(sedimenticola_asv_seqs), paste0("france_", 1:length(sedimenticola_asv_seqs)), "sedimenticola_asv_seqs_france.fasta")
write.fasta(as.list(thiodiazotropha_asv_seqs), paste0("france_", 1:length(sedimenticola_asv_seqs)), "thiodiazotropha_asv_seqs_france.fasta")



thiodiazotropha_asv_seqs_all_tree.newick






##############################################
#               sedimenticola                #
##############################################

# aligned with mafft and made a tree using geneious tree builder
tree_all_sedi <- ggtree::read.tree("sedimenticola_asv_seqs_all tree.newick")
tree_all_sedi$tip.label <- str_remove_all(tree_all_sedi$tip.label, "'")
tree_all_sedi$range <- case_when(startsWith(tree_all_sedi$tip.label, "usa") ~ "USA",
                            startsWith(tree_all_sedi$tip.label, "france") ~ "France")

ggtree(tree_all_sedi)  + geom_tiplab()

dd <- data.frame(taxa=tree_all_sedi$tip.label,
                 Sedimenticola=tree_all_sedi$range)
p<-ggtree(tree_all_sedi, size=0.3, colour="gray60")

pdf("sedimenticola_asv_seqs_all_sourcetree.pdf", height=15, width=10)
p %<+% dd + geom_tippoint(aes(color=Sedimenticola, size=Sedimenticola)) + 
  theme(legend.title = element_text(size=25),
        legend.text = element_text(size=20)) +
  scale_colour_manual(values=c("deeppink2", "darkolivegreen4")) +
  scale_size_manual(values=c(4, 1.2)) +
  guides(size = guide_legend(override.aes = list(size = 7)))
dev.off()


##############################################
#             thiodiazotropha                #
##############################################

tree_all_thio <- ggtree::read.tree("thiodiazotropha_asv_seqs_all_tree.newick")
tree_all_thio$tip.label <- str_remove_all(tree_all_thio$tip.label, "'")
tree_all_thio$range <- case_when(startsWith(tree_all_thio$tip.label, "usa") ~ "USA",
                                 startsWith(tree_all_thio$tip.label, "france") ~ "France")

ggtree(tree_all_thio)  + geom_tiplab()

dd <- data.frame(taxa=tree_all_thio$tip.label,
                 Thiodiazotropha=tree_all_thio$range)
p<-ggtree(tree_all_thio, size=0.3, colour="gray60")

pdf("thiodiazotropha_asv_seqs_all_sourcetree.pdf", height=15, width=10)
p %<+% dd + geom_tippoint(aes(color=Thiodiazotropha, size=Thiodiazotropha)) + 
  theme(legend.title = element_text(size=25),
        legend.text = element_text(size=20)) +
  scale_colour_manual(values=c("deeppink2", "darkolivegreen4")) +
  scale_size_manual(values=c(4, 1.2)) +
  guides(size = guide_legend(override.aes = list(size = 7)))
dev.off()



####################################################
#          functional annotation of taxa           #
####################################################


functional <- readxl::read_xlsx("/Users/katieemelianova/Desktop/Spartina/JR_functionally_annotated_genera.xlsx", skip=2) %>% 
  dplyr::select("S/Sulfate reducing genera", "Sulfur oxidizing genera", "Iron oxidizing genera", "Nitrifying genera") %>%
  set_colnames(c("sulfate_reducers", "sulfur_oxidisers", "iron_oxidisers", "nitrifiers"))

sulfate_reducers <- functional %>% dplyr::select(sulfate_reducers) %>% drop_na() %>% pull() %>% str_replace_all("_", " ")
sulfur_oxidisers <- functional %>% dplyr::select(sulfur_oxidisers) %>% drop_na() %>% pull() %>% str_replace_all("_", " ")
iron_oxidisers <- functional %>% dplyr::select(iron_oxidisers) %>% drop_na() %>% pull() %>% str_replace_all("_", " ")
nitrifiers <- functional %>% dplyr::select(nitrifiers) %>% drop_na() %>% pull() %>% str_replace_all("_", " ")




############## sulfur oxidisers ###########




transform_and_glom <- function(ps_object, taxon_list, filter_value){
  ps_object <<- ps_object
  taxon_list <<- taxon_list
  filter_value <<- filter_value
  transformed_and_glommed <- subset_taxa(ps_object, Genus %in% taxon_list) %>%
    subset_samples(compartment == filter_value) %>%
    tax_glom("Genus", NArm = F) %>% 
    transform_sample_counts(function(x) {x * 100/sum(x)})
  return(transformed_and_glommed)
}


get_mean_abundance <- function(ps_object, taxon_list, filter_value, transpose=FALSE){
  transformed_and_glommed <- transform_and_glom(ps_object, taxon_list, filter_value)
  if (transpose == TRUE){
    mean_abundance <- data.frame(Genus = tax_table(transformed_and_glommed)[,"Genus"], Mean = rowMeans(otu_table(transformed_and_glommed) %>% t()), row.names = NULL)
  } else {
    mean_abundance <- data.frame(Genus = tax_table(transformed_and_glommed)[,"Genus"], Mean = rowMeans(otu_table(transformed_and_glommed)), row.names = NULL)
  }
  return(mean_abundance)
}



taxon_list <- sulfur_oxidisers
filter_value <- "Rhizosphere"     # I have to declare these variables first because see https://github.com/joey711/phyloseq/issues/1695
france_sulfur_oxidisers_rhizospere_mean <- get_mean_abundance(ps_rel, taxon_list, filter_value)
usa_sulfur_oxidisers_rhizospere_mean <- get_mean_abundance(ps_rel_jr, taxon_list, filter_value, transpose=TRUE)
filter_value <- "Endosphere"
france_sulfur_oxidisers_endosphere_mean <- get_mean_abundance(ps_rel, taxon_list, filter_value)
usa_sulfur_oxidisers_endosphere_mean <- get_mean_abundance(ps_rel_jr, taxon_list, filter_value, transpose=TRUE)
sulfur_oxidisers_all <- inner_join(france_sulfur_oxidisers_rhizospere_mean,
                                   france_sulfur_oxidisers_endosphere_mean, 
                                   usa_sulfur_oxidisers_rhizospere_mean,
                                   usa_sulfur_oxidisers_endosphere_mean, by="Genus") %>% 
  set_colnames(c("Genus", "france_mean", "usa_mean")) %>%
  reshape2::melt() %>%
  mutate("role" = "Sulfur Oxidiser")

inner_join(inner_join(france_sulfur_oxidisers_rhizospere_mean,
           france_sulfur_oxidisers_endosphere_mean,
           by="Genus"),
           inner_join(usa_sulfur_oxidisers_rhizospere_mean,
           usa_sulfur_oxidisers_endosphere_mean,
           by="Genus"),
           by="Genus")





france_sulfur_oxidisers_rhizosphere <- subset_taxa(ps_rel, Genus %in% sulfur_oxidisers) %>%
  subset_samples(compartment == "Rhizosphere") %>%
  tax_glom("Genus", NArm = F) %>% 
  transform_sample_counts(function(x) {x * 100/sum(x)})

france_sulfur_oxidisers_endosphere <- subset_taxa(ps_rel, Genus %in% sulfur_oxidisers) %>%
  subset_samples(compartment == "Endosphere") %>%
  tax_glom("Genus", NArm = F) %>% 
  transform_sample_counts(function(x) {x * 100/sum(x)})

usa_sulfur_oxidisers_rhizosphere <- subset_taxa(ps_rel_jr, Genus %in% sulfur_oxidisers) %>%
  subset_samples(compartment == "Rhizosphere") %>%
  tax_glom("Genus", NArm = F) %>% 
  transform_sample_counts(function(x) {x * 100/sum(x)})

usa_sulfur_oxidisers_endosphere <- subset_taxa(ps_rel_jr, Genus %in% sulfur_oxidisers) %>%
  subset_samples(compartment == "Endosphere") %>%
  tax_glom("Genus", NArm = F) %>% 
  transform_sample_counts(function(x) {x * 100/sum(x)})


france_sulfur_oxidisers_rhizospere_mean <- data.frame(Genus = tax_table(france_sulfur_oxidisers_rhizosphere)[,"Genus"], Mean = rowMeans(otu_table(france_sulfur_oxidisers_rhizosphere)), row.names = NULL) %>% mutate()
france_sulfur_oxidisers_endosphere_mean <- data.frame(Genus = tax_table(france_sulfur_oxidisers_endosphere)[,"Genus"], Mean = rowMeans(otu_table(france_sulfur_oxidisers_endosphere)), row.names = NULL)

usa_sulfur_oxidisers_rhizospere_mean <- data.frame(Genus = tax_table(usa_sulfur_oxidisers_rhizosphere)[,"Genus"], Mean = rowMeans(otu_table(usa_sulfur_oxidisers_rhizosphere) %>% t()), row.names = NULL)
usa_sulfur_oxidisers_endosphere_mean <- data.frame(Genus = tax_table(usa_sulfur_oxidisers_endosphere)[,"Genus"], Mean = rowMeans(otu_table(usa_sulfur_oxidisers_endosphere) %>% t()), row.names = NULL)




sulfur_oxidisers_all <- inner_join(france_sulfur_oxidisers_rhizospere_mean,
                                   france_sulfur_oxidisers_endosphere_mean, 
                                   usa_sulfur_oxidisers_rhizospere_mean,
                                   usa_sulfur_oxidisers_endosphere_mean, by="Genus") %>% 
  set_colnames(c("Genus", "france_mean", "usa_mean")) %>%
  reshape2::melt() %>%
  mutate("role" = "Sulfur Oxidiser")

############## sulfate reducers ###########

france_sulfate_reducers <- subset_taxa(ps_rel, Genus %in% sulfate_reducers) %>%
  tax_glom("Genus", NArm = F) %>% 
  transform_sample_counts(function(x) {x * 100/sum(x)})

usa_sulfate_reducers <- subset_taxa(ps_rel_jr, Genus %in% sulfate_reducers) %>%
  tax_glom("Genus", NArm = F) %>% 
  transform_sample_counts(function(x) {x * 100/sum(x)})

france_sulfate_reducers_mean <- data.frame(Genus = tax_table(france_sulfate_reducers)[,"Genus"], Mean = rowMeans(otu_table(france_sulfate_reducers)), row.names = NULL)
usa_sulfate_reducers_mean <- data.frame(Genus = tax_table(usa_sulfate_reducers)[,"Genus"], Mean = rowMeans(otu_table(usa_sulfate_reducers) %>% t()), row.names = NULL)

sulfate_reducers_all <- inner_join(france_sulfate_reducers_mean, usa_sulfate_reducers_mean, by="Genus") %>% 
  set_colnames(c("Genus", "france_mean", "usa_mean")) %>%
  reshape2::melt() %>% 
  mutate("role"="Sulfur Reducer")


############## iron oxidisers ###########

france_iron_oxidisers <- subset_taxa(ps_rel, Genus %in% iron_oxidisers) %>%
  tax_glom("Genus", NArm = F) %>% 
  transform_sample_counts(function(x) {x * 100/sum(x)})

france_iron_oxidisers@otu_table[is.na(france_iron_oxidisers@otu_table)] <- 0

usa_iron_oxidisers <- subset_taxa(ps_rel_jr, Genus %in% iron_oxidisers) %>%
  tax_glom("Genus", NArm = F) %>% 
  transform_sample_counts(function(x) {x * 100/sum(x)})

usa_iron_oxidisers@otu_table[is.na(usa_iron_oxidisers@otu_table)] <- 0

france_iron_oxidisers_mean <- data.frame(Genus = tax_table(france_iron_oxidisers)[,"Genus"], Mean = rowMeans(otu_table(france_iron_oxidisers)), row.names = NULL)
usa_iron_oxidisers_mean <- data.frame(Genus = tax_table(usa_iron_oxidisers)[,"Genus"], Mean = rowMeans(otu_table(usa_iron_oxidisers) %>% t()), row.names = NULL)

iron_oxidisers_all <- inner_join(france_iron_oxidisers_mean, usa_iron_oxidisers_mean, by="Genus") %>% 
  set_colnames(c("Genus", "france_mean", "usa_mean")) %>%
  reshape2::melt() %>% 
  mutate("role"="Iron Oxidiser")


############## nitrifiers ###########

france_nitrifiers <- subset_taxa(ps_rel, Genus %in% nitrifiers) %>%
  tax_glom("Genus", NArm = F) %>% 
  transform_sample_counts(function(x) {x * 100/sum(x)})

france_nitrifiers@otu_table[is.na(france_nitrifiers@otu_table)] <- 0

usa_nitrifiers <- subset_taxa(ps_rel_jr, Genus %in% nitrifiers) %>%
  tax_glom("Genus", NArm = F) %>% 
  transform_sample_counts(function(x) {x * 100/sum(x)})

usa_nitrifiers@otu_table[is.na(usa_nitrifiers@otu_table)] <- 0

france_nitrifiers_mean <- data.frame(Genus = tax_table(france_nitrifiers)[,"Genus"], Mean = rowMeans(otu_table(france_nitrifiers)), row.names = NULL)
usa_nitrifiers_mean <- data.frame(Genus = tax_table(usa_nitrifiers)[,"Genus"], Mean = rowMeans(otu_table(usa_nitrifiers) %>% t()), row.names = NULL)

nitrifiers_all <- inner_join(france_nitrifiers_mean, usa_nitrifiers_mean, by="Genus") %>% 
  set_colnames(c("Genus", "france_mean", "usa_mean")) %>%
  reshape2::melt() %>% 
  mutate("role"="Nitrifier")




rbind(sulfate_reducers_all, sulfur_oxidisers_all, iron_oxidisers_all, nitrifiers_all) %>%
  mutate(Genus = fct_reorder(Genus, as.numeric(as.factor(role)))) %>%
  ggplot(aes(y = Genus, x = variable)) + 
  geom_point(aes(size = sqrt(value), fill = role), alpha = 0.75, shape = 21)


