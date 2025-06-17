

path <- "/Users/katieemelianova/Desktop/Spartina/spartina_amplicon_code/MiSeq_SOP" # CHANGE ME to the directory containing the fastq files after unzipping.
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- basename(fnFs) %>% str_split_i("_", 1)


filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(240,160),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=FALSE)

errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)

# Naming the derep objects may be a point where I go wrong
names(derepFs) <- sample.names
names(derepRs) <- sample.names

dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)

mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)

seqtab_tutorial <- makeSequenceTable(mergers)
seqtab.nochim_tutorial <- removeBimeraDenovo(seqtab_tutorial, method="consensus", multithread=TRUE, verbose=TRUE)
taxa_tutorial <- assignTaxonomy(seqtab.nochim_tutorial, "/Users/katieemelianova/Desktop/Spartina/spartina_amplicon_code/silva_nr99_v138.2_toGenus_trainset.fa.gz", multithread=TRUE)
taxa_tutorial <- addSpecies(taxa_tutorial, "/Users/katieemelianova/Desktop/Spartina/spartina_amplicon_code/silva_v138.2_assignSpecies.fa.gz")


samples.out_tutorial <- rownames(seqtab.nochim_tutorial)
subject <- sapply(strsplit(samples.out_tutorial, "D"), `[`, 1)
gender <- substr(subject,1,1)
subject <- substr(subject,2,999)
day <- as.integer(sapply(strsplit(samples.out_tutorial, "D"), `[`, 2))
samdf <- data.frame(Subject=subject, Gender=gender, Day=day)
samdf$When <- "Early"
samdf$When[samdf$Day>100] <- "Late"
rownames(samdf) <- samples.out_tutorial

otu_table(seqtab.nochim_tutorial, taxa_are_rows=FALSE)
sample_data(samdf)


ps_tutorial<- phyloseq(otu_table(seqtab.nochim_tutorial, taxa_are_rows=FALSE), 
               sample_data(samdf), 
               tax_table(taxa_tutorial))

ps_tutorial

sample_names(ps_tutorial) %>% is.matrix()
rank_names(ps_tutorial) %>% is.matrix()
sample_variables(ps_tutorial) %>% is.matrix()
otu_table(ps_tutorial) %>% is.matrix()
tax_table(ps_tutorial) %>% is.matrix()
taxa_names(ps_tutorial) %>% is.matrix()
plot_bar(ps_tutorial)


top20_tutorial <- names(sort(taxa_sums(ps_tutorial), decreasing=TRUE))[1:20]
ps.top20_tutorial <- transform_sample_counts(ps_tutorial, function(OTU) OTU/sum(OTU))
ps.top20_tutorial <- prune_taxa(top20_tutorial, ps.top20_tutorial)
plot_bar(ps.top20_tutorial, x="Day", fill="Family") + facet_wrap(~When, scales="free_x")





