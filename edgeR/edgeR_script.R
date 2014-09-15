# Quantitative Evolution Workshop Sep 14-22 2014, Sirince
# Script for processing and doing differential expression analysis on read counts from RNAseq
# To be used for the first 104 libraries for Protea repens from common garden
# update September 2014

# Let's install the necessary pacakages

source("http://bioconductor.org/biocLite.R")
biocLite("edgeR")
biocLite("limma")

library(limma); library(edgeR)

# Now change your working directory to the QEvolution folder

getwd() # This will show you the active working directory
setwd("~/QEvolution/edgeR") #This will set your new working directory

#### Prepare the count data for analysis 

# load the raw count data
reads = read.csv("merged_counts.csv") # csv should have the contig names in the first column,
dim(reads)
# Let's have a look at the overall picture
counts = reads[,c(2:115)]
barplot(colSums(counts)*1e-6, names =1:114, ylab="Library size (millions)")
# You will realize that there are some libraries with very low counts
hist(colSums(counts))
#Let's remove the 6 libraries with counts< 500K, these are:
# morning: POT_1_, RND_44
# afternoon: CER_1_, GAR_41, KLM_20, POT_1_
to_remove = rep(FALSE, 115)
to_remove[c(20, 25, 29, 41, 98, 106)]=TRUE
reads = reads[,!to_remove]
dim(reads)

# Let's see if it worked out fine
counts = reads[,c(2:109)]
barplot(colSums(counts)*1e-6, names =1:108, ylab="Library size (millions)")
# Now there is only one library with ~600K reads. This is a personal choice you can always be more stringent

#Let's give the names for the contigs
contigs=reads[,1]
head(contigs) # does it look ok?
head(counts) # does it look ok?

# Indicates which population each individual is from
pops = as.vector(sapply(names(counts), substr, start=1, stop=3)) # population ids in order of the columns of the reads matrix
indiv = as.vector(sapply(names(counts), substr, start=1, stop=6)) # pull out unique individual identifier
# And if all data are combined, which time course 
# Alternatively we can look at differential expression across pops at the two times separately, and see how much overlap there is between the two timepoints. 
time = c(rep("afternoon", 53), rep("morning", 55))# morning or afternoon indicator for each column

# NOTE some of the libraries in the "mornings" sequencing run were actually afternoon libraries. 
#These were:BAN_3, BRD_22, KAR_1, KSW_4, MGU_57, RND_20,SWA_45, VAN_2,
extra_pm = c(62, 67, 78, 86, 92, 99, 102, 108) #if low count samples not removed (above + POT_1)
time[extra_pm] = "afternoon"
 

#### edgeR analysis

# edgeR works with a DGEList object composed of lists(contig names, counts) and explanotory factors (we will add these as we go through the analysis). Let's create the DEGList object

dge = DGEList(counts, genes=contigs) # done! now we can do our processing

# remove contigs with low counts (Keep genes with least 1 count-per-million reads (cpm) in at least half of samples samples)
# it makes more sense for this to be done seperately for afternoons and mornings (you will see what I mean in a bit)
expr = rowSums(cpm(dge)>1)>=57
dge = dge[expr, keep.lib.sizes=FALSE]
dim(dge)
# How many contigs are removed?

# Let's normalize our libraries for different library sizes

reads.norm = calcNormFactors(dge)
reads.norm$samples # This will give the sample library sizes and corresponding normalization factors

# Let's calculate the common dispersion for all libraries. Read the edgeR user's guide for more information

reads.norm = estimateCommonDisp(reads.norm)
reads.norm$common.dispersion

#mds plot
mds <- plotMDS(reads.norm, top=500, col=(time=="afternoon")+1,  cex=0.5, ndim=3, labels=indiv, xlab="MDS axis 1 (MDS1)", ylab="MDS axis 2 (MDS2)") # this is the default way, pairwise gene selection and top 500 genes

mds <- plotMDS(reads.norm, top=15000, col=(time=="afternoon")+1,  cex=0.5, ndim=3, labels=indiv, xlab="MDS axis 1 (MDS1)", ylab="MDS axis 2 (MDS2)") # gene selection pariwise

mds <- plotMDS(reads.norm, gene.selection="common", col=(time=="afternoon")+1,  cex=0.5, ndim=3, labels=indiv, xlab="MDS axis 1 (MDS1)", ylab="MDS axis 2 (MDS2)") # gene selection common

# What do these results imply? Shall we analyze the 2 time-points seperately? Maybe a better idea. Let's analyze the morning data only.

### Analyzing morning data only ###

# We will load this new dataset now the same way we did for the previous dataset

reads.M = read.csv("morning_counts.csv")
names(reads.M)
counts = reads.M[,c(2:58)]
barplot(colSums(counts)*1e-6, names =1:57, ylab="Library size (millions)")

# Note that 24, 40 and 48 have low counts
# NOTE also that some of the libraries in the "mornings" sequencing run were actually afternoon libraries. 
# These, to be removed, are:BAN_3, BRD_22, KAR_1, KSW_4, MGU_57, RND_20,SWA_45, VAN_2 (10,15,26,34,40,48,52,58)
# Let's remove both afternoons and low count reads
to_remove = rep(FALSE, 58)
to_remove[c(10,15,26,34,40,48,52,58,25,41,49)]=TRUE
reads.M = reads.M[,!to_remove]
ncol(reads.M)
counts.M = reads.M[,c(2:47)]
barplot(colSums(counts.M)*1e-6, names =1:46, ylab="Library size (millions)") #does it look ok?

dim(counts.M)
contigs=reads.M[,1]


#create DGEList object

dge.M = DGEList(counts.M, genes=contigs)
# remove contigs with low counts(Keep genes with least 1 count-per-million reads (cpm) in at least 23 samples)
expr = rowSums(cpm(dge.M)>1)>=23
dge.M = dge.M[expr, keep.lib.sizes=FALSE]
dim(dge.M)
# How many contigs are removed?


# Indicates which population each individual is from
pops.M = as.vector(sapply(names(counts.M), substr, start=1, stop=3)) # population ids in order of the columns of the reads matrix
indiv.M = as.vector(sapply(names(counts.M), substr, start=1, stop=6)) # pull out unique individual identifier

# Normalize libraries
reads.norm.M = calcNormFactors(dge.M)
reads.norm.M$samples

# Calculate common dispersion
reads.norm.M = estimateCommonDisp(reads.norm.M)
reads.norm.M$common.dispersion

# How does the MDS plot look like?

mds.M <- plotMDS(reads.norm.M, top=500, gene.selection="common", cex=0.5, ndim=3, labels=pops.M, xlab="MDS axis 1 (MDS1)", ylab="MDS axis 2 (MDS2)") # this is the default way, pairwise gene selection and top 500 genes
# Do some populations show more variation then others? Which axis could explain the population effects better?

# Could there be any correlations between the environmental and trait variation and these MDS axes? Let's have a look.

# read in the environment and trait data, individuals in the same order as the count data 

envtdata = read.csv("envtandtraitdata_mornings_only.csv")

envtdata$mds1 = mds.M$x
envtdata$mds2 = mds.M$y

# Compare to environmental variables and traits

# Add the variable names to see if MDS axes correlate with and environmental variable (regressions)
names(envtdata)

summary(lm(mds1~raincon, envtdata))
summary(lm(mds2~raincon, envtdata))
summary(lm(mds1~map, envtdata))
summary(lm(mds2~map, envtdata))



# Which other environmental variables would be explanatory in general expression patterns? Make a summary of MDS and trait correlations


# Let's see them visually

par(mfrow=c(2,2))
plot(mds1~frost, envtdata)
plot(mds2~frost, envtdata)
plot(mds1~julmint, envtdata)
plot(mds2~julmint, envtdata)

summary(lm(mds1~Wet_leaf_Mass_g_2012, envtdata))
summary(lm(mds2~Wet_leaf_Mass_g_2012, envtdata))
plot(mds1~Wet_leaf_Mass_g_2012, envtdata)
plot(mds2~Wet_leaf_Mass_g_2012, envtdata)

# It seems like gene regulation is effected by soruce population environmental variables? How could this be explained? How about the trait values? ARe there differences between years?

### Gene by gene regressions ###
library(MASS)

We can now look at the correlations between gene expression and trait/envr variables gene by gene. We will use counts per million of reads for this analysis.

reads.cpm = cpm(reads.norm.M)

# Function to fit a negative binomial glm relating one row of expression data (counts for one contig for each population) to an environmental or trait variable measured for those same populations).
fit.glm <- function(readcount, envtvar) {
    m1 = glm.nb(round(readcount)~envtvar) # fit negative binomial regression using MASS library
    #m1 = glm(readsnorm~envtvar, family="poisson")
    return(c(coef(m1), summary(m1)$coef[2,4], 1-m1$deviance/m1$null.deviance))
}

# We can use this function to look at simple regressions between gene expression and environmental/trait data. Just replace the variable name

# Example: for map. This will take some time to run all these models, don't be alarmed!

envtvar.map = scale(envtdata$map) # center and scale the explantory variable
n.contigs = dim(reads.cpm)[1] # how many contigs there are in the normalized data set
glm.tests.map = {} # create blank variable to hold the results
for (i in 1:n.contigs) {
    glm.tests.map = rbind(glm.tests.map, fit.glm(reads.cpm[i,], envtvar.map))
}
glm.tests.map = as.data.frame(glm.tests.map)
names(glm.tests.map) = c("intercept", "slope", "p.value", "pct.dev.explained")
row.names(glm.tests.map) = reads.norm.M$genes
write.table(glm.tests.map, "glm_tests_am_map.csv", sep=",", row.names=TRUE, col.names=TRUE)
# export the contig numbers for the "significantly" associated loci
contigs.pm.map = reads.norm$genes$genes[which(glm.tests.map$p.value < 0.001)]
write(contigs.pm.map, "glmcontigs_am_map.txt")
#Find a way to print contigs names here!!

# Explore the results; insert ".variablename" after glm.tests
hist(glm.tests.map$slope)
hist(glm.tests.map$p.value)
sum(glm.tests.map$p.value < 0.001)# how many signif at a given level?
hist(glm.tests.map$pct.dev.explained)
sum(glm.tests.map$p.value < (0.001) & glm.tests.map$pct.dev.explained >= 0.2) # how many significant with pretty good R2?
z = which(glm.tests.map$p.value < 0.001) # theshold by p-value only
z = which(glm.tests.map$p.value < 0.001 & glm.tests.map$pct.dev.explained >= 0.15 & abs(glm.tests.map$slope) >= 0.1) # multidim criterion

# Now what can you do with these genes that correlate with any given environmental variable? BLAST, GO analysis etc.
