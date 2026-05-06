# Script R for 16S data analysis

# Installing dada2 and cutadapt
# https://benjjneb.github.io/dada2/tutorial.html
# dada2

#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("dada2")

library("dada2")
packageVersion("dada2") # version used 1.24.0

#cutadpat
# https://benjjneb.github.io/dada2/ITS_workflow.html
cutadapt <- "cutadapt.exe"
system2(cutadapt, args = "--version") # version used 4.0

#Other packages
library(ShortRead)
packageVersion("ShortRead")
library(Biostrings)
packageVersion("Biostrings")

### Filtering N and taking out primers

### I am considering ITS pipeline to use cutadapt
#https://benjjneb.github.io/dada2/ITS_workflow.html
# setting files path
getwd()
path <- "DENG/"
list.files(path)

#Setting forward and reverse files
fnFs <- sort(list.files(path, pattern = "_1.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern = "_2.fastq", full.names = TRUE))

#Veryfing the samples
samples_sequenced <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
samples_sequenced

#Identifying primers
FWD <- "ACTCCTACGGGAGGCAGCAG"  ## primer 16S F 
REV <- "GGACTACHVGGGTWTCTAAT"   ## primer 16S R 

# primers used in DENG et al. 2017
# 338F (5'-ACTCCTACGGGAGGCAGCAG-3') and 806R (5'-GGACTACHVGGGTWTCTAAT-3')

#In theory if you understand your amplicon sequencing setup, this is sufficient to continue. 
#However, to ensure we have the right primers, and the correct orientation of the primers on the reads,
#we will verify the presence and orientation of these primers in the data.

allOrients <- function(primer) {
  # Create all orientations of the input sequence
  require(Biostrings)
  dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
  orients <- c(Forward = dna, Complement = Biostrings::complement(dna), Reverse = Biostrings::reverse(dna),
               RevComp = Biostrings::reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector
}
FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)
FWD.orients
REV.orients

#The presence of ambiguous bases (Ns) in the sequencing reads makes accurate mapping of short primer 
#sequences difficult. Next we are going to "pre-filter" the sequences just to remove those with Ns, 
#but perform no other filtering.

fnFs.filtN <- file.path(path, "filtN", basename(fnFs)) # Put N-filterd files in filtN/ subdirectory
fnRs.filtN <- file.path(path, "filtN", basename(fnRs))
filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 0, multithread = FALSE)

#We are now ready to count the number of times the primers appear in the forward and reverse read, 
#while considering all possible primer orientations. Identifying and counting the primers on one set of paired 
#end FASTQ files is sufficient, assuming all the files were created using the same library preparation, 
#so we'll just process the first sample.

primerHits <- function(primer, fn) {
  # Counts number of reads in which the primer is found
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}

rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.filtN[[1]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.filtN[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.filtN[[1]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.filtN[[1]]))

#We now create output filenames for the cutadapt-ed files, and define the parameters we are going to give 
#the cutadapt command. The critical parameters are the primers, and they need to be in the right orientation, 
#i.e. the FWD primer should have been matching the forward-reads in its forward orientation, and the REV primer 
#should have been matching the reverse-reads in its forward orientation. Warning: A lot of output will be written
#to the screen by cutadapt!

path.cut <- file.path(path, "cutadapt")
if(!dir.exists(path.cut)) dir.create(path.cut)
fnFs.cut <- file.path(path.cut, basename(fnFs))
fnRs.cut <- file.path(path.cut, basename(fnRs))

FWD.RC <- dada2:::rc(FWD)
REV.RC <- dada2:::rc(REV)

# Trim FWD and the reverse-complement of REV off of R1 (forward reads)
R1.flags <- paste("-g", FWD, "-a", REV.RC) 
# Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
R2.flags <- paste("-G", REV, "-A", FWD.RC) 
# Run Cutadapt
for(i in seq_along(fnFs)) {
  system2(cutadapt, args = c(R1.flags, R2.flags, "-n", 2, "-m", 1, # -n 2 required to remove FWD and REV from reads
                             "-o", fnFs.cut[i], "-p", fnRs.cut[i], # output files
                             fnFs.filtN[i], fnRs.filtN[i])) # input files
}

# Sanity check
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.cut[[1]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.cut[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.cut[[1]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.cut[[1]]))

#If there are still a few of primers in the wrong orientation left, I can just ignore them
# https://github.com/marcelm/cutadapt/issues/588

#The primer-free sequence files are now ready to be analyzed through the DADA2 pipeline. 
#Similar to the earlier steps of reading in FASTQ files, we read in the names of the cutadapt-ed FASTQ files 
#and applying some string manipulation to get the matched lists of forward and reverse fastq files.

# Forward and reverse fastq filenames have the format:
cutFs <- sort(list.files(path.cut, pattern = "_1.fastq", full.names = TRUE))
cutRs <- sort(list.files(path.cut, pattern = "_2.fastq", full.names = TRUE))

### Inspecting samples and generating ASVs

#Inspecting read quality profile

plotQualityProfile(cutFs[1:2])
plotQualityProfile(cutRs[1:2])

#Filtering and trimming
#The samples of Bruna had a low quality, that is why I tested several different parameters
# and now the parameters that I use is: maxEE (8,8), and truncLen = c(260,230)
#Assigning the filenames for the output of the filtered reads to be stored as fastq.gz files.

filtFs <- file.path(path.cut, "filtered", basename(cutFs))
filtRs <- file.path(path.cut, "filtered", basename(cutRs))

# Name the filt objects by the sample names
names(filtFs) <- samples_sequenced
names(filtRs) <- samples_sequenced

out <-  filterAndTrim(cutFs, filtFs, cutRs, filtRs, truncLen = c(270,230),
                      maxN = 0, maxEE = c(4, 6), truncQ = 2, rm.phix = TRUE, 
                      compress = TRUE, multithread = FALSE)  # on windows, set multithread = FALSE

head(out)
out

# From this point on, I will switch to the 16S protocol because it has more details and it is the same
#for the both pipelines (https://benjjneb.github.io/dada2/tutorial.html)
#Learning the error rates
#The DADA2 algorithm makes use of a parametric error model (err) and every amplicon dataset has a different 
#set of error rates. The learnErrors method learns this error model from the data, by alternating estimation 
#of the error rates and inference of sample composition until they converge on a jointly consistent solution.

errF <- learnErrors(filtFs, multithread=TRUE)

errR <- learnErrors(filtRs, multithread=TRUE)

plotErrors(errF, nominalQ=TRUE)


#Sample inference
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)

dadaRs <- dada(filtRs, err=errR, multithread=TRUE)

dadaFs[[1]]


#Merge paired reads
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])

#Construct sequence table
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))


#Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)

sum(seqtab.nochim)/sum(seqtab)


##Track reads through pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- samples_sequenced
head(track)
track


#Assign taxonomy (I will use the sample 2 - better results)
#Using the database created by Denny. I should update it.

silva_16S <- "silva_nr99_v138.1_train_set.fa.gz"  # CHANGE ME to location on your machine
taxa_16S <- assignTaxonomy(seqtab.nochim, silva_16S, multithread = TRUE, tryRC = TRUE)

#inspecting
taxa.print <- taxa_16S  # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)


# Exporting tables
getwd()
write.csv(seqtab.nochim, file = "Deng_seqtab.nochim_16S.csv")
write.csv(taxa_16S, file = "Deng_taxa_16S.csv")

