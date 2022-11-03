library(dada2); packageVersion("dada2")

library(readxl)
install.packages("unmarked")
Sample_information <- read_excel("/mnt/nfs/bioinfdata/ngs/ME/kuramae_group/marcio/project-amazon-saf/Brazil_samples/Sample_information.xls")

RcppParallel::setThreadOptions(numThreads = 16)
setDadaOpt(OMEGA_A = 1e-20, OMEGA_C = 1e-20, BAND_SIZE = 32)

path <- "/mnt/nfs/bioinfdata/ngs/ME/kuramae_group/marcio/project-amazon-saf/Brazil_samples/18s" # CHANGE ME to the directory containing the fastq files after unzipping.

files1 <- list.files(path, pattern="_R1.fastq.gz", full.names = TRUE)
files2 <- list.files(path, pattern="_R2.fastq.gz", full.names = TRUE)

# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
#16S data
#Selecting first the bacteria
fnFs <- sort(files1)
fnRs <- sort(files2)

# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

plotQualityProfile(fnFs[1:2])
plotQualityProfile(fnRs[1:2])

# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered_18S", paste0(sample.names, "_18S_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered_18S", paste0(sample.names, "_18S_R_filt.fastq.gz"))

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(240,160),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE
head(out)

errF <- learnErrors(filtFs, multithread=16, nbases=1e+14, verbose=TRUE, MAX_CONSIST=20)
errR <- learnErrors(filtRs, multithread=16, nbases=1e+14, verbose=TRUE, MAX_CONSIST=20)

plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)

dada2:::checkConvergence(errF[[1]])

#Dereplication
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names

#Sample Inference
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)

dadaFs[[1]]

mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)

seqtab <- makeSequenceTable(mergers)
dim(seqtab)

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))

#remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)
#[1] 0.8178506
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)

#Assign taxonomy

taxa <- assignTaxonomy(seqtab.nochim, "/data/db/dada2/silva/132/silva_nr_v132_train_set.fa.gz", multithread=TRUE)
taxa <- addSpecies(taxa, "/data/db/dada2/silva/132/silva_species_assignment_v132.fa.gz")

#Assign taxonomy Fungi
unite.ref <- "/data/db/unite/20181118/sh_general_release_dynamic_s_02.02.2019.fasta"
taxa.fungi <- assignTaxonomy(seqtab.nochim, unite.ref, multithread = 16, tryRC = TRUE)
taxa.print.fungi <- taxa.fungi  # Removing sequence rownames for display only
rownames(taxa.print.fungi) <- NULL
head(taxa.print.fungi)

#Construct phylogenetic tree
seqs.fungi <- getSequences(seqtab.nochim)
names(seqs.fungi) <- seqs.fungi # This propagates to the tip labels of the tree
alignment.fungi <- AlignSeqs(DNAStringSet(seqs.fungi), anchor=NA)

#The phangorn R package
library(phangorn)
BiocManager::install("phangorn")
phang.align.fungi <- phyDat(as(alignment.fungi, "matrix"), type="DNA")
dm.fungi <- dist.ml(phang.align.fungi)
treeNJ.fungi <- NJ(dm.fungi) # Note, tip order != sequence order. also takes really long!
fit.fungi = pml(treeNJ.fungi, data=phang.align.fungi)

## negative edges length changed to 0!

fitGTR.fungi <- update(fit.fungi, k=4, inv=0.2)
fitGTR.fungi <- optim.pml(fitGTR.fungi, model="GTR", optInv=TRUE, optGamma=TRUE,
                          rearrangement = "stochastic", control = pml.control(trace = 0))
detach("package:phangorn", unload=TRUE)


# #Alternative DECIPHER
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("DECIPHER")
# library(DECIPHER); packageVersion("DECIPHER")
# dna <- DNAStringSet(getSequences(seqtab.nochim)) # Create a DNAStringSet from the ASVs
# load("~/tax/IDTaxa/SILVA_SSU_r132_March2018.RData") # CHANGE TO THE PATH OF YOUR TRAINING SET
# ids <- IdTaxa(dna, trainingSet, strand="top", processors=NULL, verbose=FALSE) # use all processors
# ranks <- c("domain", "phylum", "class", "order", "family", "genus", "species") # ranks of interest
# # Convert the output object of class "Taxa" to a matrix analogous to the output from assignTaxonomy
# taxid <- t(sapply(ids, function(x) {
#   m <- match(ranks, x$rank)
#   taxa <- x$taxon[m]
#   taxa[startsWith(taxa, "unclassified_")] <- NA
#   taxa
# }))
# colnames(taxid) <- ranks; rownames(taxid) <- getSequences(seqtab.nochim)

taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)

# #Construct phylogenetic tree
# seqs <- getSequences(seqtab.nochim)
# names(seqs) <- seqs # This propagates to the tip labels of the tree
# alignment <- AlignSeqs(DNAStringSet(seqs), anchor=NA)
# 
# #The phangorn R package
# BiocManager::install("phangorn")
# library(phangorn)
# aphang.align <- phyDat(as(alignment, "matrix"), type="DNA")
# dm <- dist.ml(aphang.align)
# treeNJ <- NJ(dm) # Note, tip order != sequence order. also takes really long!
# fit = pml(treeNJ, data=aphang.align)
# 
# ## negative edges length changed to 0!
# 
# fitGTR <- update(fit, k=4, inv=0.2)
# fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
#                     rearrangement = "stochastic", control = pml.control(trace = 0))
# detach("package:phangorn", unload=TRUE)

#Combine data into a phyloseq object
summary(mutate_all(data.frame(taxa.print), as.factor))
#Handoff to phyloseq
# BiocManager::install("phyloseq")
library(phyloseq); packageVersion("phyloseq")
library(ggplot2); packageVersion("ggplot2")
theme_set(theme_bw())

#Construct sample dataframe
# samples.out <- rownames(seqtab.nochim)
# subject <- sapply(strsplit(samples.out, "D"), `[`, 1)
# gender <- substr(subject,1,1)
# subject <- substr(subject,2,999)
# day <- as.integer(sapply(strsplit(samples.out, "D"), `[`, 2))
# samdf <- data.frame(Subject=subject, Gender=gender, Day=day)
# samdf$When <- "Early"
# samdf$When[samdf$Day>100] <- "Late"
# rownames(samdf) <- samples.out

#contruct phyloseq object
# new_map_wood <- read.delim("/mnt/nfs/bioinfdata/ngs/ME/kuramae_group/marcio/BacteriaUSADADA2/new_map_wood.txt")
# samdf <- data.frame(new_map_wood)
# rownames(samdf) <- samdf$Description
#samdf <- samdf[samples.out,]
library(phyloseq)

ps <- phyloseq(tax_table(taxtab), sample_data(samdf),
               otu_table(seqtab, taxa_are_rows = FALSE)) #End of the bioinformatics

ps.bac <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
                   sample_data(samdf), 
                   tax_table(taxa))
ps.bac
phy_tree(ps.bac)
plot_richness(ps.bac, x="Treatment", measures=c("Shannon", "Simpson"), color="Category")

#Combine data into a phyloseq object
seqtab.final <- seqtab.nochim
colnames(seqtab.final) <- paste0("sp", 1:ncol(seqtab.final))
mimarks_path <- "/mnt/nfs/bioinfdata/ngs/ME/kuramae_group/dilution_project/DADA2_analysis/analysis/metada.txt"
samdf <- Sample_information %>% 
  mutate(Plot = paste0("P",`Sample information`)) %>% 
  select(-...2) %>% 
  rename(SampleID = sample)
samdf$SampleID
samdf <- samdf[samdf$SampleID %in% rownames(seqtab),]
#samdf <- samdf[!duplicated(samdf$SampleID),] # Remove dupicate entries for reverse reads
rownames(seqtab) # Check labels discrepancy
all(rownames(seqtab) %in% samdf$SampleID) # TRUE

taxInfo <- taxa.print
barplot(unlist(lapply(apply(taxInfo[,-1], 2, unique), length)))

samdf.final <- sample_data(samdf)
sample_names(samdf.final) <- samdf$SampleID
ps.fungi <- phyloseq(tax_table(taxInfo), sample_data(samdf.final),
               otu_table(seqtab.final, taxa_are_rows = FALSE)#,phy_tree(fitGTR$tree) #In case of phylogenetic info
               ) #End of the bioinformatics

####-----------------------------------------------#
#phyloseq
library("phylos12eq")
library("gridExtra")
ps = readRDS("data/ps.rds")
ps

#Taxonomic Filtering

# Show available ranks in the dataset
rank_names(ps.bac)
rank_names(ps.fungi)

get_taxa_unique(ps.fungi, "Kingdom")
get_taxa_unique(ps.fungi, "Phylum")
get_taxa_unique(ps.fungi, "Genus")

# Create table, number of features for each phyla
table(tax_table(ps.fungi)[, "Phylum"], exclude = NULL)
subset_taxa(ps.fungi, is.na(Phylum) & !Phylum %in% c("", "uncharacterized"))
ps0.fungi <- subset_taxa(ps.fungi, !is.na(Phylum) & Phylum %in% c("Ascomycota", "Basidiomycota", "Blastocladiomycota", "Chytridiomycota",
                                                                  "Cryptomycota", "LKM15", "Mucoromycota", "Neocallimastigomycota", "Peronosporomycetes",
                                                                  "Zoopagomycota"))
View(table(tax_table(ps0.fungi)[, "Phylum"], exclude = NULL))

View(table(tax_table(ps0.fungi)[, "Genus"], exclude = NULL))

# # Compute prevalence of each feature, store as data.frame
# prevdf.bac = apply(X = otu_table(ps0.fungi),
#                    MARGIN = ifelse(taxa_are_rows(ps0.fungi), yes = 1, no = 2),
#                    FUN = function(x){sum(x > 0)})
# # Are there phyla that are comprised of mostly low-prevalence features? Compute the total and average prevalences of
# # the features in each phylum 
# plyr::ddply(prevdf.bac, "Phylum", function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))})
# 
# # Define phyla to filter
# filterPhyla = c("Fusobacteria", "Deinococcus-Thermus")
# # Filter entries with unidentified Phylum.
# ps1 = subset_taxa(ps0, !Phylum %in% filterPhyla)
# ps1
# 
# #
# dna.bac <- Biostrings::DNAStringSet(taxa_names(ps0.bac))
# names(dna.bac) <- taxa_names(ps0.bac)
# ps0.bac <- merge_phyloseq(ps0.bac, dna.bac)
# taxa_names(ps0.bac) <- paste0("ASV", seq(ntaxa(ps0.bac)))
# ps0.bac
samDF <- data.frame(sample_data(ps0.fungi))
samDF$`Site-code` <- substr(samDF$Plot, start = 1, stop = 4) %>% 
  gsub("A", "",.) %>%
  gsub("B", "",.) %>% 
  gsub("C", "",.) %>% 
  gsub("D", "",.) %>% 
  gsub("E", "",.) %>% 
  gsub("F", "",.)
DescriptionSAF <- read_excel("/mnt/nfs/bioinfdata/ngs/ME/kuramae_group/marcio/project-amazon-saf/Brazil_samples/DescriptionSAF.xlsx")
View(DescriptionSAF)

SitesInfo <- unique(DescriptionSAF)

samDF <- samDF %>% 
  left_join(SitesInfo, by = "Site-code")
rownames(samDF) <- samDF$SampleID
ps1.fungi <- ps0.fungi
sample_data(ps1.fungi) <- samDF
plot_richness(ps1.fungi, x="Description", measures=c("Shannon", "Simpson"), color="State")

# Transform data to proportions as appropriate for Bray-Curtis distances
ps.prop.bac <- transform_sample_counts(ps1.fungi, function(otu) otu/sum(otu))
ord.nmds.bray.bac <- ordinate(ps.prop.bac, method="NMDS", distance="bray")

SAFdata.F <- ps1.fungi
colnames(tax_table(SAFdata.F))
# colnames(tax_table(SAFdata.F)) = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
get_taxa_unique(SAFdata.F, "Kingdom")
get_taxa_unique(SAFdata.F, "Phylum")
get_taxa_unique(SAFdata.F, "Class")
get_taxa_unique(SAFdata.F, "Family")
get_taxa_unique(SAFdata.F, "Genus")
get_taxa_unique(SAFdata.F, "Species")



save(SAFdata.F,ps1.fungi,DescriptionSAF,SitesInfo, samDF,ps0.fungi,ps.fungi,taxInfo,seqtab.final,samdf.final,samdf,
           file = "SAF18Sfinal.RData")

# Define the categorical variable
get_variable(SAFdata.F)
vartotal.tab <- get_variable(SAFdata.F) %>% 
  mutate(farm_type = ifelse(land_use == "forest", "forest", farm_type))

#Plotting

p = plot_bar(SAFdata.F, "Kingdom", fill="Phylum", facet_grid=~Description)
p + geom_bar(aes(color=Phylum, fill=Phylum), stat="identity", position="stack")

p = plot_bar(bactdata, "Phylum", fill="Class", facet_grid=~salt)
p + geom_bar(aes(color=Phylum, fill=Class), stat="identity", position="stack")

p = plot_bar(archdata, "Phylum", fill="Phylum", facet_grid=~Treatments)
p + geom_bar(aes(color=Phylum, fill=Phylum), stat="identity", position="stack")


#Fungi
fungidata <- SAFdata.F
# fungidata <- subset_taxa(fungidata, Phylum!="NA")
# fungidata <- subset_taxa(fungidata, Class!="unclassified")
# fungidata = prune_taxa(taxa_sums(fungidata) > 15, fungidata) # keep OTUs with >15 counts across all the samples
any(is.na(otu_table(fungidata))) # check for NAs, should be false
any(otu_table(fungidata) < 0) # check for weird OTUs, should be false
sample_names(fungidata)

#
fdata.phyloseq <- tax_glom(fungidata, taxrank="Genus", NArm=FALSE)  # agglomerate samples up to Genus level

fdata <- tax_glom(fungidata, taxrank="Genus", NArm=FALSE)  # agglomerate samples up to Genus level
# fdata <- prune_taxa(taxa_sums(fdata) > 10, fdata)
ftemp <- as.matrix(cbind(t(otu_table(fdata)),tax_table(fdata)))
row.names(ftemp) <- as.vector(paste(ftemp[,"Phylum"], ftemp[,"Genus"], sep = "_")) # check that these are the same first
ftemp <- subset(ftemp, select=-c(Kingdom, Phylum, Class, Order, Family, Genus, Species))
class(ftemp) <- "numeric"
fdata <- as.data.frame(t(ftemp))
unique(colnames(fdata))
SampleDataFungi <- get_variable(fungidata)
rownames(fdata) <- SampleDataFungi$Plot

# new pdata at genus level
pdata2.phyloseq <- tax_glom(procdata, taxrank="Genus", NArm=FALSE)  # agglomerate samples up to Class level

pdata2 <- tax_glom(procdata, taxrank="Genus", NArm=FALSE)  # agglomerate samples up to Class level
# pdata2 <- prune_taxa(taxa_sums(pdata) > 10, pdata2)
ptemp <- as.matrix(cbind(otu_table(pdata2),tax_table(pdata2)))
vartotal <- get_variable(procdata, "Description")
row.names(ptemp) <- as.vector(paste(ptemp[,"Phylum"], ptemp[,"Genus"], sep = "_")) # check that these are the same first
ptemp <- subset(ptemp, select=-c(Domain, Phylum, Class, Order, Family, Genus, Species))
class(ptemp) <- "numeric"
pdata2 <- as.data.frame(t(ptemp))
unique(colnames(pdata2))
rownames(pdata2)

SampleDataBac <- get_variable(procdata)%>% 
  mutate(Plot = gsub("\\.", "",X.SampleID))

SampleTot <-  SampleDataFungi %>% 
  left_join(SampleDataBac, by = "Plot")
View(SampleTot)
write.table(SampleTot, file = "SampleTot.txt", sep = "\t", quote = FALSE)
write.table(SampleDataBac, file = "SampleDataBac.txt", sep = "\t", quote = FALSE)

#Solve pseudo-replicates


repSAF <- subset_samples(fungidata, Site.code %in% c("P117", "P68", "P70", "P119", "P120", "P121"))
sample_data(repSAF)$Sample.information
p = plot_bar(repSAF, "Kingdom", fill="Phylum", facet_grid=~Sample.information)
p = plot_bar(repSAF, x= "Sample.information", y="Kingdom", fill="Phylum")
p + ggplot2::geom_bar(ggplot2::aes(color=Phylum, fill=Phylum), stat="identity", position="stack")
repTab <- data.frame(TotAbund = colSums(t(otu_table(repSAF))))
rownames(repTab) <- sample_data(repSAF)$Sample.information
barplot(t(repTab),las=2)


SampleTot2 <- SampleDataFungi %>% 
  mutate(Sample.information = gsub("117C1", "117C", Sample.information)) %>% #View
  mutate(Sample.information = gsub("68D2", "68D", Sample.information)) %>% 
  mutate(Sample.information = gsub("68E2", "68E", Sample.information)) %>% 
  mutate(Sample.information = gsub("70E1", "70E", Sample.information)) %>% 
  mutate(Sample.information = gsub("119C2", "119C", Sample.information)) %>% 
  mutate(Sample.information = gsub("120B3", "120B", Sample.information)) %>% 
  mutate(Sample.information = gsub("120C3", "120C", Sample.information)) %>% 
  mutate(Sample.information = gsub("120E1", "120E", Sample.information)) %>% 
  mutate(Sample.information = gsub("121D1", "121D", Sample.information)) %>% 
  mutate(Plot = paste0("P", Sample.information)) %>% 
  left_join(SampleDataBac, by = "Plot") #%>% View

envdataT <- phy_chm_all %>% 
  rename(Site.code = Code)

#Final info table
AllInfo <- SampleTot2 %>% 
  left_join(envdataT, by = "Site.code")# %>% View

#Combining bacteria and fungi data
y <- pdata2[,-1] %>% 
  separate(taxonomy, c("Kingdom","Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = "; ") %>%
  filter(Kingdom %in% c("k__Archaea", "k__Bacteria")) %>% 
  group_by(Kingdom, Phylum,  Class,   Order,   Family,  Genus,   Species) %>% 
  summarise_all(sum) %>% 
  ungroup
dim(y)

dim(fdata)
yF <- fdata[,-1] %>% 
  separate(taxonomy, c("Kingdom","Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = "; ") %>% 
  filter(Kingdom == "k__Fungi") %>% 
  filter(Phylum %in% c("p__Ascomycota","p__Basidiomycota","p__Chytridiomycota","p__Glomeromycota",    
                       "p__Mortierellomycota","p__Mucoromycota","p__Olpidiomycota","p__Rozellomycota")) %>% 
  group_by(Kingdom, Phylum,  Class,   Order,   Family,  Genus,   Species) %>% 
  summarise_all(sum) %>% 
  ungroup

dim(yF)
apply(y[,1:7],2,unique)
apply(yF[,1:7],2,unique)
View(yF)

taxInfo.B.gen <- data.frame(tax_table(pdata2.phyloseq)@.Data)
dim(taxInfo.B.gen)
taxInfo.F.gen <- data.frame(tax_table(fdata.phyloseq)@.Data) %>% 
  rename(Domain = Kingdom)
dim(taxInfo.F.gen)

TaxInfo.T <- rbind(taxInfo.B.gen, taxInfo.F.gen) %>% 
  select(Domain, Phylum,  Class,   Order,   Family,  Genus,   Species) %>% 
  mutate(TaxLab = paste(Phylum,Genus, sep = "_")) %>% 
  mutate(TaxAbbr = paste(abbreviate(Phylum,6),abbreviate(Genus,6))) %>% 
  mutate(TaxAbbr = ifelse(TaxAbbr %in% TaxAbbr[duplicated(TaxAbbr)], paste(abbreviate(Phylum,6),abbreviate(Family,4),abbreviate(Genus,6)),TaxAbbr)) %>% 
  mutate(TaxAbbr = ifelse(TaxAbbr %in% TaxAbbr[duplicated(TaxAbbr)], paste(abbreviate(Phylum,6),abbreviate(Order,4),abbreviate(Family,4),abbreviate(Genus,6)),TaxAbbr)) %>% 
  mutate(TaxAbbr = ifelse(TaxAbbr %in% TaxAbbr[duplicated(TaxAbbr)], paste(abbreviate(Phylum,6),abbreviate(Class,4),abbreviate(Order,4),abbreviate(Family,4),abbreviate(Genus,6)),TaxAbbr))
View(TaxInfo.T)
View(TaxInfo.T[duplicated(TaxInfo.T$TaxAbbr),])
TaxInfo.T$TaxAbbr[duplicated(TaxInfo.T$TaxAbbr)]

View(colnames(pdata2)[duplicated(colnames(pdata2))])
View(colnames(fdata)[duplicated(colnames(fdata))])

View(data.frame(TaxInfo.T, Colnamestb = c(colnames(pdata2),colnames(fdata))))
class(pdata2) <- "numeric"
ydataP <- pdata2  
  
colnames(ydataP) <- TaxInfo.T[TaxInfo.T$TaxLab %in% colnames(ydataP), "TaxAbbr"]
ydataP <- mutate_all(ydataP, function(x) as.numeric(as.character(x)))
View(ydataP)
barplot(rowSums(ydataP))
View(rowSums(ydataP))
summary(rowSums(ydataP))
ydataP2 <- ydataP %>% 
  mutate(SampleID = rownames(pdata2), TotBac = rowSums(ydataP)) %>% 
  filter(TotBac>2000)
View(ydataP2)
ydataP2$SampleID

ydataF <- fdata
barplot(rowSums(ydataF))
summary(rowSums(ydataF))
rownames(ydataF) <- TaxInfoF$TaxAbbr
ydataF <- t(ydataF)
View(ydataF)
colnames(ydataF) <- TaxInfo.T[TaxInfo.T$TaxLab %in% colnames(ydataF), "TaxAbbr"]

#Covariates Data
colnames(AllInfo)
View(AllInfo)
cov.DT <- AllInfo %>%
  mutate(Description.x = ifelse(Site.code == "P121", "CapM", Description.x))%>% 
  mutate(Region = ifelse(Site.code == "P121", "GUR", Region)) %>% 
  mutate(State = ifelse(Site.code == "P121", "MA", State)) %>% 
  data.frame()
# rownames(cov.DT) <- cov.DT$Plate_location
View(cov.DT)

#GJAM
library(gjam)
colnames(cov.DT)
levels(factor(cov.DT$Description.x))
xdata  <- cov.DT %>% 
  mutate(LandUse = factor(Description.x, levels = c("FLO", "CapA", "CapM", "CapB", "CapEn", "SAFC", "SAFQ","SAFAb"))) %>% 
  mutate(Region = factor(Region), State = factor(State))
summary(xdata)
plot(xdata[,-1])
rownames(xdata) <- rownames(cov.DT)

# otu    <- ydata[rownames(cov.DT),]
# dim(otu)
# 
# par(mfrow=c(1,3), bty='n', mar=c(1,1,1,1), oma = c(0,0,0,0), 
#     mar = c(3,2,2,1), tcl = -0.5, mgp = c(3,1,0), family='')
# hist(otu, nclass=100, ylab = 'Reads', main='each observation')
# nobs <- gjamTrimY(otu, minObs = 1, OTHER = F)$nobs
# hist(nobs, nclass=100, ylab = 'Total reads per OTU', main='Full sample')
# dev.off()

tmpP <- gjamTrimY(ydataP2[,1:941], minObs = 150)$y %>% 
  data.frame() %>% 
  rename(otherP = other)
dim(ydataP)               # all OTUs
dim(tmpP)                       # trimmed data
tail(colnames(tmpP))            # 'other' class added
colSums(tmpP)
View(tmpP)
rownames(tmpP) <- gsub("\\.", "", ydataP2$SampleID)
rownames(tmpP) <- ifelse(rownames(tmpP) == "P32C1", "P32C", rownames(tmpP))
rownames(tmpP) <- ifelse(rownames(tmpP) == "P113B1", "P113B", rownames(tmpP))
rowInfo

tmpF <- gjamTrimY(ydataF, minObs = 40)$y %>% 
  data.frame() %>% 
  rename(otherF = other)
dim(ydataF)               # all OTUs
dim(tmpF)                       # trimmed data
tail(colnames(tmpF))            # 'other' class added
colSums(tmpF)
View(rownames(tmpF))
rownames(tmpP)
rownames(tmpF)

tmpP$Site.code <- rownames(tmpP)
tmpF$Site.code <- rownames(tmpF) %>% 
  gsub("P117C1", "P117C", .) %>% #View
  gsub("68D2", "68D",.) %>% 
  gsub("68E2", "68E",.) %>% 
  gsub("70E1", "70E",.) %>% 
  gsub("119C2", "119C",.) %>% 
  gsub("120B3", "120B", .) %>% 
  gsub("120C3", "120C", .) %>% 
  gsub("120E1", "120E",.) %>% 
  gsub("121D1", "121D", .) 
  
ydataT <- tmpF %>% 
  left_join(tmpP, by = "Site.code")
View(ydataT)
ydataT$Site.code
summary(ydataT[,1:5])
View(colnames(ydataT))
heatmap(as.matrix(vegan::decostand(as.matrix(ydataT[,-97]), MARGIN = 2, method = "standardize", na.rm= TRUE)))
dim(ydataT)
dim(ydataT[complete.cases(ydataT),])
ydataT$Site.code
cov.DT$Sample.code <- paste0("P", cov.DT$Sample.information)




ydataT.cov.DT <- ydataT %>% 
  rename(Sample.code = Site.code) %>% 
  left_join(cov.DT, by = "Sample.code")
View(ydataT.cov.DT[,colnames(cov.DT)])

dim(cov.DT)

#GJAM model
S     <- ncol(ydata)
colnames(ydata)
typeNames    <- c(rep('CC',S-3),"DA","CA","CA")   # composition count data
typeNames


#Standardize continuous variables
xdata
colnames(xdata)
xdata$Plot <- factor(xdata$Plot)
ydataFinal <- ydataT.cov.DT %>% 
  select(Ascmyc.Asprgl:otherF,Vrrcmc.DA10.NA:otherP,C_big_10cm.:C_total_abv_bel) %>% 
  mutate(Silte.argila = 100-Gross_sand-areia_fina) %>% 
  mutate(Gross_sand = Gross_sand/100, areia_fina = areia_fina/100, Silte.argila = Silte.argila/100) %>% 
  select(-Silte, -Argila,-Silte.argila.1)
colnames(ydataFinal)

#shortdata
colnames(ydataT.cov.DT)
ydataFinal <- ydataT.cov.DT %>% 
  select(Ascmyc.Asprgl:otherF,Vrrcmc.DA10.NA:otherP,C_big_10cm.:C_total_abv_bel) %>% 
  mutate(Silte.argila = 100-Gross_sand-areia_fina) %>% 
  mutate(Gross_sand = Gross_sand/100, areia_fina = areia_fina/100, Silte.argila = Silte.argila/100) %>% 
  select(-Silte, -Argila,-Silte.argila.1)
colnames(ydataFinal)

typeNames    <- c(rep('CC',ncol(select(ydataFinal,Ascmyc.Asprgl:otherF))),   # composition count data
                  rep('CC',ncol(select(ydataFinal,Vrrcmc.DA10.NA:otherP))),# composition count data
                  rep('CA',ncol(select(ydataFinal,C_big_10cm.:Saturation))),
                  rep('FC',ncol(select(ydataFinal,Gross_sand:Silte.argila))),
                  rep('CA',ncol(select(ydataFinal,M.O:C_total_abv_bel)))
)

CCFBgroups <- c(rep(1,ncol(select(ydataFinal,Ascmyc.Asprgl:otherF))),   # composition count data
                rep(2,ncol(select(ydataFinal,Vrrcmc.DA10.NA:otherP))),# composition count data
                rep(0,ncol(select(ydataFinal,C_big_10cm.:Saturation))),
                rep(0,ncol(select(ydataFinal,Gross_sand:Silte.argila))),
                rep(0,ncol(select(ydataFinal,M.O:C_total_abv_bel)))
)

rl <- list(r = 15, N = 30)
ml <- list(ng = 2000, burnin = 500, typeNames = typeNames, reductList = rl#, #random = "Plot",
           #CCgroups = CCFBgroups
           )

output <- gjam(~ LandUse, xdata, ydataFinal, modelList = ml)
ydataFinalC <- ydataFinal[complete.cases(ydataFinal),]
randvalues <- matrix(rexp(2460), 15)
cov.DT2
ydataFinalC[,244:269] <- ydataFinalC[,244:269]+(randvalues/10)
ydataFinalC[,255:257] <- ydataFinal[complete.cases(ydataFinal),255:257]
ydataFinalC <- ydataFinalC[-c(99,101,102),]
dim(ydataFinalC)
xdataC <- xdata[rownames(ydataFinalC),]
dim(xdataC)
xdataC$LandUse <- factor(xdataC$LandUse, levels = c("FLO", "CapA", "CapM", "CapB", "CapEn", "SAFC", "SAFQ"))
View(xdataC$LandUse)
#short
outputC <- gjam(~ LandUse, xdataC, ydataFinalC, modelList = ml)
outputC
ml.r <- list(ng = 4000, burnin = 1000, typeNames = typeNames, reductList = rl, random = "Region",
             CCgroups = CCFBgroups
)
outputC.r <- gjam(~ LandUse, xdataC, ydataFinalC, modelList = ml.r)
outputC.r
View(outputC.r$parameters$betaStandXTable)

#CLR transformed data
d.czm <- zCompositions::cmultRepl(ydataFinalCLR[,1:243],  label=0, method="CZM")
ydataFinalCLR <- ydataFinalC
d.clr <- apply(d.czm, 2, function(x){log(x) - mean(log(x))})
ydataFinalCLR[,1:243] <- d.clr
typeNamesCLR    <- c(rep('CON',ncol(select(ydataFinalCLR,Ascmyc.Asprgl:otherF))),   # composition count data
                              rep('CON',ncol(select(ydataFinalCLR,Vrrcmc.DA10.NA:otherP))),# composition count data
                              rep('CA',ncol(select(ydataFinalCLR,C_big_10cm.:Saturation))),
                              rep('FC',ncol(select(ydataFinalCLR,Gross_sand:Silte.argila))),
                              rep('CA',ncol(select(ydataFinalCLR,M.O:C_total_abv_bel)))
)
ml.clr <- list(ng = 4000, burnin = 1000, typeNames = typeNamesCLR, reductList = rl, random = "Region"
)
outputCLR <- gjam(~ LandUse, xdataC, ydataFinalCLR, modelList = ml.clr)
outputCLR
ml.r <- list(ng = 4000, burnin = 1000, typeNames = typeNames, reductList = rl, random = "Region",
             CCgroups = CCFBgroups
)
outputC.r <- gjam(~ LandUse, xdataC, ydataFinalC, modelList = ml.r)
outputC.r
View(outputC.r$parameters$betaStandXTable)
#RandomEffects
ydataFinal.R <- ydataT.cov.DT %>% 
  select(Ascmyc.Asprgl:otherF,Vrrcmc.DA10.NA:otherP)

typeNames.R    <- c(rep('CC',ncol(select(ydataFinal,Ascmyc.Asprgl:otherF))),   # composition count data
                  rep('CC',ncol(select(ydataFinal,Vrrcmc.DA10.NA:otherP)))# composition count data
                  # rep('CA',ncol(select(ydataFinal,C_big_10cm.:Saturation))),
                  # rep('FC',ncol(select(ydataFinal,Gross_sand:Silte.argila))),
                  # rep('CA',ncol(select(ydataFinal,M.O:C_total_abv_bel)))
)

CCFBgroups.R <- c(rep(1,ncol(select(ydataFinal,Ascmyc.Asprgl:otherF))),   # composition count data
                rep(2,ncol(select(ydataFinal,Vrrcmc.DA10.NA:otherP)))# composition count data
                # rep(0,ncol(select(ydataFinal,C_big_10cm.:Saturation))),
                # rep(0,ncol(select(ydataFinal,Gross_sand:Silte.argila))),
                # rep(0,ncol(select(ydataFinal,M.O:C_total_abv_bel)))
)


ml.r <- list(ng = 5000, burnin = 2000, typeNames = typeNames.R, reductList = rl, random = "Plot",
             CCgroups = CCFBgroups.R)
ml.r <- list(ng = 5000, burnin = 2000, typeNames = typeNames.R, reductList = rl, random = "Region",
             CCgroups = CCFBgroups.R)
# ml.r <- list(ng = 5000, burnin = 2000, typeNames = typeNames, reductList = rl, random = c("Region","Plot"),
#              CCgroups = CCFBgroups.R)
output.r <- gjam(~ LandUse, xdata, ydataFinal.R, modelList = ml.r)


#Graphs
#Sensitivity
selected.model <- output
selected.model <- outputCLR
selected.model$inputs$y
ynames <- colnames(selected.model$inputs$y)
ynames
# group  <- 
bact  <- gsub("_","", colnames(select(ydataFinalC,Vrrcmc.DA10.NA:otherP)))
fungi  <- gsub("_","", colnames(select(ydataFinalC,Ascmyc.Asprgl:otherF)))
soil  <- gsub("_","", colnames(select(ydataFinalC,Soil_Dens:C_stock_solo)))
biom  <- gsub("_","", c(colnames(select(ydataFinalC,C_big_10cm.:C_total_aboveground)),colnames(select(ydataFinalC,C_total_abv_bel))))
ynames %in% c(colnames(select(ydataFinalC,C_big_10cm.:C_total_aboveground)),colnames(select(ydataFinalC,C_total_abv_bel)))
# full <- gjamSensitivity(selected.model)
# cc   <- gjamSensitivity(selected.model, group)
# nt <- ncol(full)
# 
# ylim <- range( full )
# boxplot( full, boxwex = 0.25, # at = 1:nt, 
#          col='gray', log='y',
#          xaxt = 'n', ylim = ylim, 
#          xlab = 'Predictors', ylab='Sensitivity')
# axis(1,at=1:nt,
#      labels=colnames(full),las=1,cex = 0.8)
# 
# # ylim <- range( rbind(cc,bc) )
# ylim <- range( rbind(cc,full) )
# Labels <- gsub("farmtype.x","", colnames(full)) %>% gsub("landuse.x","", .)
# boxplot( full, boxwex = 0.25,  at = 1:nt - .21, col='blue', log='y',
#          xaxt = 'n', ylim = ylim, 
#          xlab = 'Predictors', ylab='Sensitivity')
# boxplot(cc, boxwex = 0.25, at = 1:nt + .2, col='forestgreen', add=T,
#         xaxt = 'n')
# axis(1,at=1:nt,labels=abbreviate(Labels,6),las=2,cex = 0.8)
# legend('bottomleft',c('full','microbes'),
#        text.col=c('blue','forestgreen'))

bS   <- gjamSensitivity(selected.model, bact)
fS   <- gjamSensitivity(selected.model, fungi)
sS   <- gjamSensitivity(selected.model, soil)
bioS   <- gjamSensitivity(selected.model, biom)
# 
# colnames(bS)
# 
# ylim <- range( rbind(bS,fS,sS,pS,CWMS) )
# boxplot(bS, boxwex = 0.25,  at = 1:nt - .5, col='blue', log='y',ylim = ylim, xlim = c(0,nt+1),
#         xaxt = 'n', xlab = 'Predictors', ylab='Sensitivity')
# boxplot(fS, boxwex = 0.25, at = 1:nt -.25, col='forestgreen', add=T,
#         xaxt = 'n')
# boxplot(sS, boxwex = 0.25, at = 1:nt, col='brown', add=T,
#         xaxt = 'n')
# boxplot(pS, boxwex = 0.25, at = 1:nt + .25, col='red', add=T,
#         xaxt = 'n')
# boxplot(CWMS, boxwex = 0.25, at = 1:nt + .5, col='orange', add=T,
#         xaxt = 'n')
# axis(1,at=1:nt,labels=Labels,las=2)
# legend('bottomleft',c('bac','fun', 'soil', 'prod', 'CWM'),
#        text.col=c('blue','forestgreen', 'brown', 'red', 'orange'))

SensTab <- data.frame(GroupVar = c(rep("Bact", nrow(bS)),rep("Fung", nrow(fS)),rep("Soil", nrow(sS)),rep("Biom", nrow(bioS))),
                      rbind(bS,fS, sS, bioS))
View(SensTab)
SensTab %>% 
  gather(key = "TreatVar", value = "Sensitivity", -GroupVar) %>% 
  mutate(GroupVar = factor(GroupVar, levels = c("Bact", "Fung", "Soil", "Biom"))) %>% 
  mutate(TreatVar = gsub("LandUse", "",TreatVar)) %>% 
  mutate(TreatVar = factor(TreatVar, levels = c("FLO", "CapA", "CapM", "CapB", "CapEn", "SAFC", "SAFQ"))) %>% 
  # mutate(TreatVar = gsub("farmtype.x", "FT",TreatVar)) %>% 
  # mutate(TreatVar = gsub("landuse.x", "LU",TreatVar)) %>% 
  # mutate(TreatVar = factor(TreatVar, levels = c("FTagroecological", "FTconventional", "FTlarge.scale", "LUcoffee", "LUpasture",
  #                                               "FTagroecological.LUcoffee", "FTlarge.scale.LUcoffee"))) %>% #summary
  ggplot(aes(x = TreatVar, y = Sensitivity, fill = GroupVar))+
  geom_boxplot()+
  scale_y_log10()+
    theme_cowplot()
colnames(SensTab)
SensTab %>% 
  gather(key = "TreatVar", value = "Sensitivity", -GroupVar) %>% 
  mutate(GroupVar = factor(GroupVar, levels = c("Bact", "Fung", "Soil", "Biom"))) %>% #View
  mutate(TreatVar = gsub("LandUse", "",TreatVar)) %>% 
  mutate(TreatVar = factor(TreatVar, levels = c("FLO", "CapA", "CapM", "CapB", "CapEn", "SAFC", "SAFQ"))) %>% 
  # mutate(TreatVar = gsub("farmtype.x", "FT",TreatVar)) %>% 
  # mutate(TreatVar = gsub("landuse.x", "LU",TreatVar)) %>% 
  mutate(TreatVar = factor(TreatVar, levels = colnames(SensTab[,-1]))) %>% #summary
  ggplot(aes(x = TreatVar, y = Sensitivity, fill = GroupVar))+
  geom_boxplot()+
  scale_y_log10()


boxplot(output$ematrix)
boxplot(output.r$ematrix)

#Coefficients
heatmap(selected.model$parameters$betaStandXWmu)
library("gplots")
netTab <- selected.model$parameters$betaStandXWmu
colnames(netTab)
# netTab[,1:241] <- netTab[,1:241]*2
rownames(netTab) <- gsub("LandUse","", rownames(netTab))
# netTab[,colnames(outputrTnew.MiNULL$parameters$betaStandXWmu)] <- outputrTnew.MiNULL$parameters$betaStandXWmu
heatmap.2(t(netTab), scale = "none", col = rev(bluered(100)), cexRow =0.8,#lhei=lhei, lwid=lwid,lmat =lmat,
          trace = "none", density.info = "none")
coeffTabfinal <- selected.model$parameters$betaStandXWTable
# coeffTabfinal[rownames(outputrTnew.MiNULL$parameters$betaStandXWTable),] <- outputrTnew.MiNULL$parameters$betaStandXWTable
NullSigEff <- coeffTabfinal %>% 
  mutate(VarEffect = rownames(selected.model$parameters$betaStandXWTable)) %>% 
  tidyr::separate(VarEffect, into = c("DepVar", "Treat"), sep = "_") %>%
  #tidyr::separate(DepVar, into = c("Phylum", "Class")) %>% 
  # mutate(DepVar = ifelse(DepVar == "HAl", "H_Al", DepVar)) %>% #View
  mutate(DepVar = factor(DepVar, levels = unique(DepVar))) %>% 
  filter(sig95 == "*")



unique(NullSigEff$DepVar)
labelsENG <- colnames(netTab[,unique(NullSigEff$DepVar)]) %>% 
  gsub(".NA", "",.)
labelsENG[207:226] <- c("Plants > 10cm (dbh)" ,          "Plants < 10cm (dbh)" ,       "Dead logs" ,  "Leaf litter"  ,           "Twigs litter",
                        "TAGB"  ,  "Water (%)"  ,      "Porosity"  ,        "Saturation" ,         "Org. Matter" ,                "pH",
                        "P"  ,                 "K" ,                  "Ca" ,                 "Mg."   ,              "H+Al" ,                "Na",
                        "Al" ,  "C stock soil"   ,       "LAGB" )

TaxInfo.TGJAM <- TaxInfo.T %>% 
  mutate(TaxAbbrGJAM = gsub(" ", ".", TaxAbbr)) %>% 
  mutate(TaxAbbrGJAM = gsub("_", "", TaxAbbrGJAM)) %>% 
  mutate(TaxAbbrGJAM = gsub("-", ".", TaxAbbrGJAM)) %>%
  mutate(TaxAbbrGJAM = gsub("\\(", ".", TaxAbbrGJAM))

TaxInfo.T.Sig <- TaxInfo.T %>% 
  mutate(TaxAbbrGJAM = gsub(" ", ".", TaxAbbr)) %>% 
  mutate(TaxAbbrGJAM = gsub("_", "", TaxAbbrGJAM)) %>% 
  mutate(TaxAbbrGJAM = gsub("-", ".", TaxAbbrGJAM)) %>%
  mutate(TaxAbbrGJAM = gsub("\\(", ".", TaxAbbrGJAM)) %>%
  filter(TaxAbbrGJAM %in% unique(NullSigEff$DepVar)) %>% 
  mutate(Color = Domain)
TaxInfo.T.Sig$TaxAbbrGJAM
unique(NullSigEff$DepVar)
View(TaxInfo.T.Sig[match(TaxInfo.T.Sig$TaxAbbrGJAM, unique(NullSigEff$DepVar)[1:206]),])
length(TaxInfo.T.Sig$Color)
levels(TaxInfo.T.Sig$Color) <- c("#e41a1c", "#377eb8", "#4daf4a")

ColortabGJAM <- data.frame(SigVars = unique(NullSigEff$DepVar)[1:206], TaxAbbrGJAM = unique(NullSigEff$DepVar)[1:206]) %>% 
  left_join(TaxInfo.T.Sig, by = "TaxAbbrGJAM")

colorVec <- c(as.character(ColortabGJAM$Color), rep("#984ea3",6),rep("#ff7f00",12),rep("#984ea3",2))

heatmap.2(netTab[,unique(NullSigEff$DepVar)], scale = "none", col = rev(bluered(100)),cexCol = 0.4, cexRow =0.6,#lhei=lhei, lwid=lwid,lmat =lmat,
          trace = "none", density.info = "none", labCol = labelsENG, colCol = colorVec)

TaxInfo.T.Sig
labelsTab <- data.frame(TaxAbbrGJAM = unique(NullSigEff$DepVar)) %>% 
  left_join(TaxInfo.T.Sig, by = "TaxAbbrGJAM") %>% 
  mutate(finalLabel = ifelse(!is.na(Genus), as.character(Genus), 
                             ifelse(!is.na(Family), as.character(Family),
                                    ifelse(!is.na(Order), as.character(Order),
                                           ifelse(!is.na(Class), as.character(Class),as.character(Phylum))))))
labelsTab$finalLabel[207:226] <- c("Plants > 10cm (dbh)" ,          "Plants < 10cm (dbh)" ,       "Dead logs" ,  "Leaf litter"  ,           "Twigs litter",
                               "TAGB"  ,  "Water (%)"  ,      "Porosity"  ,        "Saturation" ,         "Org. Matter" ,                "pH",
                               "P"  ,                 "K" ,                  "Ca" ,                 "Mg."   ,              "H+Al" ,                "Na",
                               "Al" ,  "C stock soil"   ,       "LAGB" )
View(labelsTab)
# distance & hierarchical clustering

SelectedVars <- as.vector(unique(NullSigEff$DepVar)[-c(203,223,224)])
SelectedVars <- c(SelectedVars[1:208],"X.CDESTRUTIVA",SelectedVars[209:223])
labelsTabFig <- labelsTab[-c(203,223,224),]
CoeffTabCluster <- data.frame(t(netTab[,SelectedVars]))
CoeffTabCluster$Labels <- c(labelsTabFig$finalLabel[1:208],"Shrubs",labelsTabFig$finalLabel[209:223])
View(CoeffTabCluster)

write.table(CoeffTabCluster, file = "exportCoefftabTree.txt", sep = "\t", row.names = TRUE, quote = FALSE)

distanceSig= dist(netTab[,SelectedVars], method ="euclidean")    
hclusterSig = hclust(distanceSig, method ="ward.D")
plot(hclust(distanceSig, method ="ward.D2"))
distanceVar= dist(t(netTab[,SelectedVars]), method ="euclidean")    
hclusterVar = hclust(distanceVar, method ="ward.D")
summary(hclust(distanceVar, method ="ward.D"))
plot(hclusterVar)
g <- data.frame(cutree(hclusterVar, k = c(2,4,8)))
plot(hclusterVar, labels =labelsTab$finalLabel, type = "fan", tip.color = colorVec)
hclusterVar$labels <- CoeffTabCluster$Labels
plot(ape::as.phylo(hclusterVar), type = "fan", tip.color = colorVec,
     label.offset = 1, cex = 0.7)
plot(ape::as.phylo(hclusterVar),    "u", lab4ut = "axial", tip.color = colorVec,
     label.offset = 1, cex = 0.7)
legend("topright",      # location of the legend on the heatmap plot
       legend = c("Fungi", "Bacteria", "Archea",  "Plant biomass", "Soil"), # category labels
       col = unique(colorVec),  # color key
       lty= 1,             # line style
       lwd = 10            # line width
)
g$TaxAbbrGJAM <- rownames(g)
View(g)
TaxInfo.T.Sig %>% 
  left_join(g, by = "TaxAbbrGJAM") %>% write.table("ClusterTaxInfo.txt", sep = "\t", row.names = FALSE, quote = FALSE)

library(ape)
my_tree <- as.phylo(hclusterVar) 
write.tree(phy=my_tree, file="exported_tree_new.newick") 

table(grp2 = g[,"2"], grp4 = g[,"4"], grp8 = g[,"8"])
heatmap.2(netTab[,SelectedVars],  
          # main = paste( "test"),  
          trace="none",          
          # margins =c(5,7),      
          col=rev(bluered(100)),        
          # breaks=col_breaks,     
          dendrogram="both",      
          Rowv = as.dendrogram(hclusterSig),  
          Colv = as.dendrogram(hclusterVar), 
          key.title = "Regression Coefficient",
          density.info=c("density"),
          labRow = c("MF", "OSF", "MSF", "YSF", "EFA", "CPA", "HA"),
          labCol = CoeffTabCluster$Labels, 
          colCol = c(colorVec[1:208],"#984ea3",colorVec[209:223]),
          ColSideColors =c(colorVec[1:208],"#984ea3",colorVec[209:223]),
          # legend = c("Fungi", "Archea", "Bacteria", "Plant biomass", "Soil"), 
          # cexRow =0.6,
          cexCol = 0.4,
          na.rm = TRUE ) 
par(lend = 1)           # square line ends for the color legend
legend("topright",      # location of the legend on the heatmap plot
       legend = c("Fungi", "Bacteria", "Archea",  "Plant biomass", "Soil"), # category labels
       col = unique(colorVec),  # color key
       lty= 1,             # line style
       lwd = 10            # line width
)
#Barplot for summary of effects
NullSigEff %>% 
  mutate(TaxAbbrGJAM = DepVar, Effect = ifelse(Estimate >0, "pos", "neg")) %>% 
  left_join(TaxInfo.T.Sig, by = "TaxAbbrGJAM") %>% 
  ggplot()+
  geom_bar(aes(x = Phylum, y = ..count.., fill = Effect))+
  geom_bar(aes(x = Phylum, y = -..count.., fill = Effect))+
  # geom_bar(subset = .(Effect == "pos"), aes(y = ..count.., fill = Effect), stat = "identity") + 
  # geom_bar(subset = .(Effect == "neg"), aes(y = -..count.., fill = Effect), stat = "identity")+
  facet_wrap(~Treat)

NullSigEff %>% 
  mutate(TaxAbbrGJAM = DepVar, Effect = ifelse(Estimate >0, "pos", "neg")) %>% 
  left_join(TaxInfo.T.Sig, by = "TaxAbbrGJAM") %>%
  ggplot(aes(Phylum)) +
  geom_bar(subset = .(Effect == "pos"), aes(y = ..count.., fill = Effect), stat = "identity") + 
  geom_bar(subset = .(Effect == "neg"), aes(y = -..count.., fill = Effect), stat = "identity") + 
  xlab("")

df.m <- NullSigEff %>% 
  mutate(TaxAbbrGJAM = DepVar, Effect = ifelse(Estimate >0, "pos", "neg")) %>% 
  left_join(TaxInfo.T.Sig, by = "TaxAbbrGJAM") %>% 
  mutate(Treat = gsub("LandUse", "", Treat)) %>% 
  mutate(Treat = factor(Treat, levels = c("CapB", "CapM", "CapA", "FLO", "CapEn", "SAFQ", "SAFC"))) %>% 
  mutate(Effect = factor(Effect, levels = c("pos", "neg"))) %>% 
  mutate(Phylum = factor(Phylum, levels = unique(Phylum)))
levels(df.m$Phylum)
View(unique(df.m[,10:11]))
LU.names <- as_labeller(c(`CapB` = "YSF",`CapM` = "MSF",`CapA` = "OSF",`FLO` = "MF",`CapEn` = "EFA", `SAFQ` = "HA", `SAFC` = "CPA"))
df.m %>% 
  mutate(Phylum = factor(Phylum, levels = c("Acidobacteria","Actinobacteria","Bacteroidetes", "Chloroflexi",
                                            "Cyanobacteria","Elusimicrobia","Firmicutes","Gemmatimonadetes","Nitrospirae",
                                            "Planctomycetes","Proteobacteria","Verrucomicrobia","Thaumarchaeota","Ascomycota",
                                            "Basidiomycota","Mucoromycota","Chytridiomycota","Cryptomycota","LKM15"))) %>% 
  ggplot(aes(x = Phylum)) + 
  geom_bar(data = subset(df.m, Effect == "pos"), 
           aes(y = ..count.., fill = Effect)) +
  geom_bar(data = subset(df.m, Effect == "neg"), 
           aes(y = -..count.., fill = Effect)) + 
  geom_hline(yintercept = 0,colour = "grey90")+
  scale_x_discrete(limits = c("Acidobacteria","Actinobacteria","Bacteroidetes", "Chloroflexi",
                              "Cyanobacteria","Elusimicrobia","Firmicutes","Gemmatimonadetes","Nitrospirae",
                              "Planctomycetes","Proteobacteria","Verrucomicrobia","Thaumarchaeota","Ascomycota",
                              "Basidiomycota","Mucoromycota","Chytridiomycota","Cryptomycota","LKM15"))+
  facet_wrap(~Treat, ncol = 4,labeller=LU.names)+
  geom_text(stat='count',data = subset(df.m, Effect == "pos"), 
            aes(label = ..count..), vjust=-0.1)+
  geom_text(stat='count', data = subset(df.m, Effect == "neg"), position = "stack", 
            aes(label = -..count..))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        panel.background = element_rect(fill = "white", colour = "grey50", linetype = "solid"))

df.m %>% 
  mutate(Phylum = factor(Phylum, levels = c("Acidobacteria","Actinobacteria","Bacteroidetes", "Chloroflexi",
                                            "Cyanobacteria","Elusimicrobia","Firmicutes","Gemmatimonadetes","Nitrospirae",
                                            "Planctomycetes","Proteobacteria","Verrucomicrobia","Thaumarchaeota","Ascomycota",
                                            "Basidiomycota","Mucoromycota","Chytridiomycota","Cryptomycota","LKM15"))) %>% 
  ggplot(aes(x = Phylum, y = Estimate, color = Effect)) + 
  geom_jitter()+
  scale_x_discrete(limits = c("Acidobacteria","Actinobacteria","Bacteroidetes", "Chloroflexi",
                              "Cyanobacteria","Elusimicrobia","Firmicutes","Gemmatimonadetes","Nitrospirae",
                              "Planctomycetes","Proteobacteria","Verrucomicrobia","Thaumarchaeota","Ascomycota",
                              "Basidiomycota","Mucoromycota","Chytridiomycota","Cryptomycota","LKM15"))+
  geom_hline(yintercept = 0,colour = "grey90")+
  facet_wrap(~Treat, ncol = 4,labeller=LU.names)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        panel.background = element_rect(fill = "white", colour = "grey50", linetype = "solid"))


library(cowplot)
title = unique(NullSigEff$Treat)
varSelected = title[1]
labelsTab <- mutate_all(labelsTab, list(as.character))

NullSigEff %>% 
  mutate(sig95 = factor(sig95), Treat = factor(Treat)) %>% 
  mutate(TaxAbbrGJAM = as.character(DepVar)) %>% 
  left_join(labelsTab, by = "TaxAbbrGJAM") %>% 
  mutate(Groups = ifelse(DepVar %in% bact, "bact",
                         ifelse(DepVar %in% fungi, "fungi",
                                ifelse(DepVar %in% soil, "soil",
                                       ifelse(DepVar %in% biom, "biom", "Error"))))) %>% #View
  mutate(Groups = ifelse(Domain == "Archaea", "archaea", as.character(Groups))) %>% #View
  mutate(Groups = #ifelse(DepVar %in% bact, "bact",
                         ifelse(DepVar %in% fungi, "fungi",
                                ifelse(DepVar %in% soil, "soil",
                                       ifelse(DepVar %in% biom, "biom", as.character(Groups))))) %>% #View
  mutate(Groups = factor(Groups, levels = c("bact", "fungi","archaea", "biom", "soil"))) %>% 
  mutate(Treat = gsub("LandUse", "", Treat)) %>% 
  mutate(Treat = factor(Treat, levels = c("FLO", "CapA", "CapM", "CapB", "CapEn", "SAFC", "SAFQ"))) %>% 
  #filter(Treat == varSelected) %>% 
  mutate(rank = rank(Estimate, ties.method = "first")) %>%
  ungroup() %>% 
  mutate(rank2 = 1:nrow(.)) %>% 
  mutate(DepVarAbb = factor(finalLabel, levels = unique(finalLabel))) %>% 
  # mutate(DepVarAbb = factor(abbreviate(DepVar, 10), levels = unique(abbreviate(DepVar, 10)))) %>% 
  ggplot(aes(x = reorder(DepVarAbb,Phylum), y = Estimate, color = Groups))+
  geom_point()+
  geom_errorbar(aes(ymin=CI_025, ymax=CI_975), width=.1)+
  ylab(varSelected)+
  geom_hline(yintercept=0, linetype="dashed")+
  coord_flip()+
  scale_color_manual(values = unique(colorVec))+
  facet_grid(Groups+Phylum~Treat,scales = "free",labeller=LU.names, 
             space = "free")+
  theme_cowplot(font_size = 8)

#environmental Response
envRes <- selected.model$parameters$ematrix
# envRes[rownames(outputrTnew.MiNULL$parameters$ematrix),colnames(outputrTnew.MiNULL$parameters$ematrix)] <- outputrTnew.MiNULL$parameters$ematrix
heatmap(selected.model$parameters$ematrix)
# heatmap(output$parameters$ematrix)
# heatmap(output.r$parameters$ematrix)

heatmap(envRes)

#Network graph
selected.model$parameters$betaStandXWTable
betarTnew <- selected.model$parameters$betaStandXWTable

beta.graphviz <-betarTnew %>% 
  mutate(VarName = rownames(betarTnew)) %>% 
  # mutate(VarName = gsub('\\b\\.','',VarName)) %>% mutate(VarName = gsub('\\b\\:','',VarName)) %>% mutate(VarName = gsub('\\b\\-','',VarName)) %>% 
  filter(sig95 == "*") %>% 
  separate(VarName, into = c("DepVar","TreatEff"), sep = "_") %>% 
  select(DepVar, TreatEff, Estimate) 
graph_from_data_frame(beta.graphviz,directed=TRUE, vertices=TreatEff)
print(graph_from_data_frame(beta.graphviz), e=TRUE, v=TRUE)

beta.graphviz %>% 
  mutate(color = ifelse(Estimate >0 , "blue", "red")) %>% 
  mutate(digraph = paste(TreatEff, "->", DepVar, "[color=", color, ", penwidth=3.24, label =", Estimate, " ];")) %>% 
  select(digraph) %>% write.table(file = "digraphP1.txt", quote = FALSE, row.names = FALSE)

data.frame(VarCode = c(unique(beta.graphviz$TreatEff), unique(beta.graphviz$DepVar))) %>% 
  mutate(VarName = paste(VarCode, "[fontsize = 28, fontname = Helvetica];")) %>% 
  select(VarName)%>% write.table(file = "digraphP2.txt", quote = FALSE, row.names = FALSE)

beta.network <- betarTnew %>% 
  mutate(VarName = rownames(betarTnew)) %>% 
  filter(sig95 == "*") %>% 
  separate(VarName, into = c("DepVar","TreatEff"), sep = "_") %>% 
  select(DepVar, TreatEff, Estimate) %>% 
  spread(key = DepVar, value = Estimate, fill = 0)

#Network with igraph
library(igraph) # Load the igraph package
net <- graph(envRes,  weighted=TRUE, add.rownames="code") 

biggroups <- data.frame(VarNames = ynames, 
                        GroupVars = c(rep("Bact", length(bact)),rep("Fung", length(fungi)),rep("Soil", length(soil)),rep("Biom", length(biom)))) #%>% 

install.packages("qgraph")
library(qgraph)
graphtab <- plyr::rbind.fill(data.frame(selected.model$ematrix),beta.network[,-1])
diag(envRes) <- 0
qgraph(envRes,#Check that
       posCol = "blue", negCol = "red", minimum = 0.5, cut = 0.8,
       vsize = 3, 
       groups = biggroups$GroupVars, 
       legend = TRUE, borders = FALSE)
title("Big 5 correlations", line = 2.5)

qgraph(envRes, layout = "circle",  
       groups = biggroups$GroupVars,# rotation = "promax", 
       minimum = 0.5, cut = 0.9, borders = FALSE, vTrans = 200)

emptyDF <- data.frame(matrix(0, 6, 6), row.names = rownames(netTab))
colnames(emptyDF) <- rownames(netTab)
TestTab <- data.frame(rbind(envRes[-103,-103],netTab)) %>% #View()
  cbind(data.frame(rbind(t(netTab),emptyDF)))
qgraph(TestTab , minimum = 0.25, cut = 0.4,posCol = "blue", negCol = "red",
       #vsize = 1.5, 
       groups = factor(c(as.vector(biggroups$GroupVars),rep("Treat",6))), 
       legend = TRUE, borders = FALSE)
qgraph(TestTab , minimum = 0.25, cut = 0.4,layout = "circle",
       vsize = 1.5, groups = factor(c(as.vector(biggroups$GroupVars),rep("Treat",6))), 
       legend = TRUE, borders = FALSE)

PlotCorr <- envRes
colnames(PlotCorr) <- abbreviate(colnames(PlotCorr), 10)
rownames(PlotCorr) <- abbreviate(rownames(PlotCorr), 10)
install.packages("network")
cluteredPlot <- corrplot::corrplot(PlotCorr,tl.cex = 0.5,#tl.col=network::as.color(biggroups$GroupVars, opacity = 1), 
                                   method = "color",#type = "lower", 
                                   order = "hclust", 
                                   addrect = 4)
rownames(biggroups) <- abbreviate(biggroups$VarNames, 10)
cluterlabels <- biggroups[rownames(cluteredPlot),]
corrplot::corrplot(PlotCorr,tl.cex = 0.5,tl.col=network::as.color(cluterlabels$GroupVars, opacity = 1), method = "color",#type = "lower", 
                   #order = "hclust", 
                   addrect = 4)

#Model predictions
#What should I change to make all the systems have the same microbiome as found in the Mature forest?
beta <- selected.model$parameters$betaTable
ws   <- grep('FLO',rownames(beta))  # find coefficients for status
beta[ws,]

output$parameters$sigMu['FLO',]

#What is the microbiome of a MF?
gjamPredict(selected.model, y2plot = colnames(selected.model$inputs$y)) #predict the data in-sample
xdataMF     <- selected.model$inputs$xdata
xdataMF$LandUse <- rep("FLO", length(xdataMF$LandUse))     # mean for x[,3]
newdata   <- list(xdata = xdataMF[1:2,c(1,ncol(xdataMF))], nsim = 50 )
p1 <- gjamPredict(selected.model, newdata = newdata)
View(p1$sdList$yMu)
summary(p1$sdList$yMu)
origdata   <- list(xdata = xdataC[,c(1,ncol(xdataC))], nsim = 50 )
porig <- gjamPredict(selected.model, newdata = origdata)
heatmap(porig$sdList$yMu[xdataC$LandUse == "FLO",])
#What should I change?
# colnames(ydataEST)
# nvar <- 123 #hoursmowning
# nvar <- 104:112 #soilFertility
# nvar <- 106 #Psoil
p1DF <- porig$sdList$yMu[xdataC$LandUse == "FLO",]
summary(p1DF)
yMF <-apply(p1DF[,c(bact, fungi)],2,median)
yMFsd <-apply(p1DF[,c(bact, fungi)],2,sd)
summary(yMFsd)

AllLists <-apply(as.data.frame(yMF),1,function(x) cbind(x,yMF))
View(AllLists)
dim(AllLists)
dim(xdataC)
newSB    <- list(ydataCond = AllLists[1:161,], nsim=50)
mature <- gjamPredict(selected.model, newdata = newSB) 
heatmap(mature$sdList$yMu[,c(bact, fungi)])

colnames(mature$sdList$yMu)
heatmap(mature$sdList$yMu[,244:269])
maturemicobe <- mature$sdList$yMu[,244:269]
porigPS <- porig$sdList$yMu[,244:269]

xdataC$LandUse != "FLO"
change <- maturemicobe-porigPS
SiteSystem <- xdataC$LandUse[xdataC$LandUse != "FLO"]
heatmap.2(change[xdataC$LandUse != "FLO",], scale = "none", col = rev(bluered(100)),cexCol = 0.4, cexRow =0.6,#lhei=lhei, lwid=lwid,lmat =lmat,
          trace = "none", density.info = "none", labRow = SiteSystem)

change2 <- xdataC
colnames(change2) <- gsub("_","",colnames(change2))
dim(change2[,colnames(maturemicobe)])
change2 <- change2[,colnames(maturemicobe)]
change3 <- change2 %>% 
  mutate(Mg. = ifelse(Mg. <0, -1*Mg., Mg.))
View(change3)
changeorig <- maturemicobe-change3
dim(changeorig)
dim(maturemicobe)
heatmap.2(as.matrix(changeorig[xdataC$LandUse != "FLO",]), scale = "none", col = rev(bluered(100)),cexCol = 0.4, cexRow =0.6,#lhei=lhei, lwid=lwid,lmat =lmat,
          trace = "none", density.info = "none", labRow = SiteSystem)
SimData <- data.frame(LU =SiteSystem, change[xdataC$LandUse != "FLO",]) 
SimData <- data.frame(LU =SiteSystem, changeorig[xdataC$LandUse != "FLO",]) 
# infor <- SimData[xdataC$LandUse %in% c("CapEn", "SAFQ", "SAFC"),-1]*2
# SimData[xdataC$LandUse %in% c("CapEn", "SAFQ", "SAFC"),-1] <- infor
colnames(SimData)
SimData[,-c(13:15)] %>% 
  gather(key = "TreatVar", value = "Sensitivity", -LU) %>% 
  mutate(Sensitivity = ifelse(LU %in% c("CapEn", "SAFQ", "SAFC"), Sensitivity*2,Sensitivity)) %>% 
  # mutate(GroupVar = factor(GroupVar, levels = c("Bact", "Fung", "Soil", "Biom"))) %>% #View
  # mutate(TreatVar = gsub("LandUse", "",TreatVar)) %>% 
  mutate(LU = factor(LU, levels = c("FLO", "CapA", "CapM", "CapB", "CapEn", "SAFC", "SAFQ"))) %>% 
  # mutate(TreatVar = gsub("farmtype.x", "FT",TreatVar)) %>% 
  # mutate(TreatVar = gsub("landuse.x", "LU",TreatVar)) %>% 
  # mutate(TreatVar = factor(TreatVar, levels = colnames(SensTab[,-1]))) %>% #summary
  ggplot(aes(x = TreatVar, y = Sensitivity, fill = LU))+
  geom_boxplot()+
  facet_wrap(~LU)+
  #ylim(-150,100)+
  ylim(-4,25)+ # For the ratio of changes
  theme(axis.text.x = element_text(angle=90, hjust=1))

#ratio of changes
changeratio <- maturemicobe/change3
SimData.ratio <- data.frame(LU =SiteSystem, changeratio[xdataC$LandUse != "FLO",]) 
colnames(SimData.ratio)
SimData.ratio[,-c(13:15)] %>% 
  select(-CsoilWB,-Cstocksolo,-Ctotalabvbel,-SoilDens,-Na,-Al,-Ctotalaboveground) %>% 
  gather(key = "TreatVar", value = "Sensitivity", -LU) %>% 
  # mutate(Sensitivity = ifelse(LU %in% c("CapEn", "SAFQ", "SAFC"), Sensitivity*2,Sensitivity)) %>% 
  # mutate(GroupVar = factor(GroupVar, levels = c("Bact", "Fung", "Soil", "Biom"))) %>% #View
  # mutate(TreatVar = gsub("LandUse", "",TreatVar)) %>% 
  mutate(LU = factor(LU, levels = c("FLO", "CapA", "CapM", "CapB", "CapEn", "SAFC", "SAFQ"))) %>% 
  mutate(TreatVar = factor(TreatVar, levels = unique(TreatVar))) %>% 
  # mutate(TreatVar = gsub("farmtype.x", "FT",TreatVar)) %>% 
  # mutate(TreatVar = gsub("landuse.x", "LU",TreatVar)) %>% 
  # mutate(TreatVar = factor(TreatVar, levels = colnames(SensTab[,-1]))) %>% #summary
  ggplot(aes(x = TreatVar, y = Sensitivity, fill = LU))+
  geom_boxplot()+
  facet_wrap(~LU,labeller=LU.names2)+
  #ylim(-150,100)+
  geom_hline(yintercept=1, linetype="dashed", color = "black")+
  ylim(-0.5,17)+ # For the ratio of changes
  ylab("Ratio of change")+
  xlab("")+
  scale_x_discrete(labels= labels.ratio[labels.ratio!="TAGB"])+
  # scale_y_continuous(limits=c(0,NA),expand = c(0.05,0.05), breaks=labelsR, labels=labelsR)+
  # scale_y_continuous(breaks=labelsR, labels=labelsR, expand=c(0.075,0))+
  theme(axis.text.x = element_text(angle=90, hjust=1),
        panel.background = element_blank(),legend.position = "none",
        axis.line = element_line(colour = "black"))
labels.ratio <- c(labelsTab$finalLabel[207:211], "Shrubs", labelsTab$finalLabel[212:222])
LU.names2 <- as_labeller(c(`CapB` = "Young Secondary Forest (YSF)",
                          `CapM` = "Mid-age Secondary Forest (MSF)",
                          `CapA` = "Old Secondary Forest (OSF)",`FLO` = "Mature Forest (MF)",
                          `CapEn` = "Enriched Fallow (EFA)", `SAFQ` = "Homegarden (HA)", 
                          `SAFC` = "Commercial Plantations (CPA)"))
tablvalues <- na.omit(SimData.ratio[,-c(1,13:15)])
tablvalues[tablvalues < Inf]
step <- 5
scale <- floor(35/ 45) - 1
breaksR <- seq(0, max(tablvalues[tablvalues < Inf]), step)
labelsR <- seq(0, 15+step, step)
labelsR <- c(0,5,10,15,20,35,40)
labelsR <- append(labelsR, scale * seq(from=ceiling((-0.5 + step) / step) * step, 
                                       length.out=length(breaksR) - length(labelsR), by=step))

SimData.ratio[,-c(13:15)] %>% 
  select(-CsoilWB,-Cstocksolo,-Ctotalabvbel,-SoilDens,-Na,-Al,-Ctotalaboveground) %>% 
  gather(key = "TreatVar", value = "Sensitivity", -LU) %>% 
  # mutate(Sensitivity = ifelse(LU %in% c("CapEn", "SAFQ", "SAFC"), Sensitivity*2,Sensitivity)) %>% 
  # mutate(GroupVar = factor(GroupVar, levels = c("Bact", "Fung", "Soil", "Biom"))) %>% #View
  # mutate(TreatVar = gsub("LandUse", "",TreatVar)) %>% 
  mutate(LU = factor(LU, levels = c("FLO", "CapA", "CapM", "CapB", "CapEn", "SAFC", "SAFQ"))) %>% 
  mutate(TreatVar = factor(TreatVar, levels = unique(TreatVar))) %>%
  filter(Sensitivity<17) %>% 
  group_by(LU, TreatVar) %>% 
  summarise(Median = median(Sensitivity), Q1 = quantile(Sensitivity,0.25), Q3 = quantile(Sensitivity,0.75)) %>% 
  ungroup() %>% 
  mutate(Result = paste0(round(Median,2)," [", round(Q1,2), "-", round(Q3,2), "]")) %>% 
  select(LU, TreatVar,Result) %>% 
  spread(key = "LU", value = "Result") %>% View
  

#Standardize the data
SimData.std <- data.frame(LU = SimData$LU,vegan::decostand(SimData[,-c(1,13:15)], method = "standardize"))
SimData.std %>% 
  gather(key = "TreatVar", value = "Sensitivity", -LU) %>% 
  mutate(Sensitivity = ifelse(LU %in% c("CapEn", "SAFQ", "SAFC"), Sensitivity*2,Sensitivity)) %>% 
  # mutate(GroupVar = factor(GroupVar, levels = c("Bact", "Fung", "Soil", "Biom"))) %>% #View
  # mutate(TreatVar = gsub("LandUse", "",TreatVar)) %>% 
  mutate(LU = factor(LU, levels = c("FLO", "CapA", "CapM", "CapB", "CapEn", "SAFC", "SAFQ"))) %>% 
  # mutate(TreatVar = gsub("farmtype.x", "FT",TreatVar)) %>% 
  # mutate(TreatVar = gsub("landuse.x", "LU",TreatVar)) %>% 
  # mutate(TreatVar = factor(TreatVar, levels = colnames(SensTab[,-1]))) %>% #summary
  ggplot(aes(x = TreatVar, y = Sensitivity, fill = LU))+
  geom_boxplot()+
  facet_wrap(~LU)+
  #ylim(-150,100)+
  theme(axis.text.x = element_text(angle=90, hjust=1))


SimDataFLO <- data.frame(LU =SiteSystem, change[xdataC$LandUse == "FLO",]) 
# infor <- SimData[xdataC$LandUse %in% c("CapEn", "SAFQ", "SAFC"),-1]*2
# SimData[xdataC$LandUse %in% c("CapEn", "SAFQ", "SAFC"),-1] <- infor
SimDataFLO %>% 
  gather(key = "TreatVar", value = "Sensitivity", -LU) %>% 
  mutate(Sensitivity = ifelse(LU %in% c("CapEn", "SAFQ", "SAFC"), Sensitivity*2,Sensitivity)) %>% 
  # mutate(GroupVar = factor(GroupVar, levels = c("Bact", "Fung", "Soil", "Biom"))) %>% #View
  # mutate(TreatVar = gsub("LandUse", "",TreatVar)) %>% 
  mutate(LU = factor(LU, levels = c("FLO", "CapA", "CapM", "CapB", "CapEn", "SAFC", "SAFQ"))) %>% 
  # mutate(TreatVar = gsub("farmtype.x", "FT",TreatVar)) %>% 
  # mutate(TreatVar = gsub("landuse.x", "LU",TreatVar)) %>% 
  # mutate(TreatVar = factor(TreatVar, levels = colnames(SensTab[,-1]))) %>% #summary
  ggplot(aes(x = TreatVar, y = Sensitivity, fill = LU))+
  geom_boxplot()+
  facet_wrap(~LU)+
  theme(axis.text.x = element_text(angle=90, hjust=1))+
  theme_cowplot()
