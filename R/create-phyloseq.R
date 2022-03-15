### CREATE A PHYLOSEQ OBJECT FROM DADA2 DATA

# EVS 3/2022

# load packages that we need
library(phyloseq)
install.packages("tidyverse") # if it needs to be installed
library(tidyverse)

### ---- make phyloseq ----

# load ASV table
load("dada2-data/metforminasv-table.RData")

# load tax table
load("dada2-data/metformintax-table.RData")

# matrices in R
dim(seqtab.nochim)
dim(tax)

# create phyloseq object
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows = FALSE),
               tax_table(tax))

# look at the phyloseq object
ps
# get number of taxa
ntaxa(ps)
# get taxa ranks
rank_names(ps)

# access the data "slots" with @
head(ps@tax_table)
head(ps@otu_table)

# fix ASV names
### from dada2 tutorial: fix ASV names
dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
ps

# get taxa names
head(taxa_names(ps))
head(sample_names(ps))

## ---- metadata ----

# read in metadata file
dat <- readxl::read_xlsx(path = "data/sampleIDs-with-sequencing-controls.txt",
                         sheet = 1)
# read in text file
dat <- read.table("data/sampleIDs-with-sequencing-controls.txt", sep = "\t", header = TRUE)

# look at data
head(dat)
str(dat)

# fix sample names to get ONLY the sample ID
names <- str_extract(sample_names(ps), "M(\\d){1,3}")

# change the sample names to NAMES
sample_names(ps) <- names

# did it work?
head(sample_names(ps))

# format our data to add to phyloseq
sampdf <- dat %>% 
  column_to_rownames(var = "NovaGeneID")

# add to phyloseq
sample_data(ps) <- sampdf

# did it work?
ps

## save our output

# re-name our phyloseq
psraw <- ps

# save as RImage
save(psraw, file = "data/phyloseq-raw.RData")
