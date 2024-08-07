#---------------------#
####     Setup     ####
#---------------------#
###### Setup options ######
## Create the options object
opt<-list()
opt$wd <- path.expand(file.path("~","Bio321G_DiffExp"))

## Setup package options
opt$pckDir       <- "/stor/scratch/Bio321G_NB_Spring2023/R" #This is my package library made for this assignment with (largely) compatible software versions.
opt$createPckDir <- T # This toggles the if statement that adjusts the package directory
opt$biocPackages <- c("ggplot2","ggrepel","TCGAbiolinks","DESeq2") # These are the minimum packages you need to run this script file


## For downloading GDC data
opt$gdcPath      <- "/stor/scratch/Bio321G_NB_Spring2023/GDCdata" #We can all download GDC files to this file. Please don't put other files here.
opt$createGdcDir <- T # This toggles the if statement that adjusts the GDC data directory. Set this to TRUE if you are using a non-EduPod system

## Set query filtering options
opt$minPerGroup <- 12           # Used to control warning about small sample size
opt$maxPerGroup <- 50           # Used to control how large before warning about size is produced & sampling size
opt$sampleToMaxLogic <- T       # Used to control whether query is reduced to a smaller size when large
opt$sampleToMax.SeedLogic <- T  # Used to control whether seed is set prior to sampling
opt$sampleToMax.Seed <- 3242023 # Used to specify the seed if set


## Set cutoffs for determining significance
opt$fdr.cut.off   <- 0.01 ## Maximum false discovery rate adjusted p-value cutoff:
opt$lfc.cut.off   <- 1.50 ## Minimum Log fold DIFFERENCE (absolute value of change) cutoff to be considered

## For visualization
opt$visScriptFile <- "./volcanoPlotFun.R" # You won't have this file as it is part of your homework
opt$nPerSide      <- 10                  # How many significant gene names to print per x-axis side. Used in the visScriptFile.


###### Setup location to store packages ######
## Setup the working directory
dir.create(opt$wd)
setwd(opt$wd)

###### Install packages ######
## Adjust package directory if needed
if(!dir.exists(opt$pckDir)){
    warning("Changing package directory...")
    opt$pckDir    <- file.path(opt$wd,"R_GDC")
    dir.create(opt$pckDir)
}

## Tell R where to find packages
.libPaths(opt$pckDir)

## You should only need to install once
if(!all(opt$biocPackages%in%installed.packages())){
    ## Install an installer tool
    # Installation goes into custom library path due to file permission issues with some packages
    if(!"BiocManager"%in%installed.packages()){
        install.packages("BiocManager", lib = opt$pckDir)
    }
    
    ## Update the installed packages
    update.packages(instlib = opt$pckDir,ask = F)
    
    ## Install on modern non-EduPod version of R
    # If on Windows and you do not have Rtools installed, say no to the pop-up dialogue box (or install Rtools)
    if(!R.Version()$version.string=="R version 4.0.3 (2020-10-10)"){
        BiocManager::install(
            opt$biocPackages, lib = opt$pckDir,
            ask = F,update = T
        )
    }
}
    

###### Load packages ######
## Don't forget to change your package dir:  
.libPaths(opt$pckDir)
## Just check the following for warnings like the following:
# Ex: Error in library(gibberish) : there is no package called ‘gibberish’
# Most other warning messages mean the package was installed, but it gave feedback.
library(ggplot2)
library(ggrepel)
library(TCGAbiolinks) #Takes a bit to load
library(DESeq2)       #Produces lots of "warnings"

#---------------------------#
#### Download TCGA Files ####
#---------------------------#
###### Search for TCGA Files ######
#Read more about this using browseVignettes("TCGAbiolinks")
# For this analysis of adrenal gland cancer samples, the project was specified
query1 <- GDCquery(
    project = c("TCGA-PCPG",'TCGA-ACC'),
    data.category = "Transcriptome Profiling",
    experimental.strategy = "RNA-Seq",
    data.type = "Gene Expression Quantification",
    sample.type = "Primary Tumor" #Remove this one if searching for non
)

###### Pull out IDs from initial search   ######
## Pull sample specific data out of complicated query object
samDf <- query1$results[[1]]
table(samDf$sample_type) # Primary tumor: 258

###### Load the IDs identified on the GDC ######

# dowloading data from computer: 
group1 <- read.table('muc16_nonMut.tsv', fill=T,header = 1)
group2 <- read.table('muc16_mutants.tsv',fill=T,header = 1)
library(dplyr)

## Check group size and if the exclusion analysis worked
if(nrow(group1)<opt$minPerGroup){warning("Group1: You should have a larger number of samples for this kind of analysis")}
if(nrow(group2)<opt$minPerGroup){warning("Group2: You should have a larger number of samples for this kind of analysis")}

if(nrow(group1)>opt$maxPerGroup){warning("Group1: You should have a smaller number of samples if working on the edupod")}
# 87 samples originally in group without a MUC16 mutation (group 1)
if(nrow(group2)>opt$maxPerGroup){warning("Group2: You should have a smaller number of samples if working on the edupod")}

if(!all(!group1$ID.1%in%group2$ID.1)&all(!group2$ID.1%in%group1$ID.1)){
    stop("The group IDs were not exclusive...")
}

###### Compare the searched-for & GDC IDs ######
## Subset samDf to just those rows with "sample.submitter_id" values present in groupIds
groupIds<-c(group1$ID.1,group2$ID.1)
samDf <- samDf[samDf$sample.submitter_id%in%groupIds,]

###### Sample if query too large for EduPod ###### 
if(opt$sampleToMaxLogic){
    samDf.1 <- samDf[samDf$sample.submitter_id%in%group1$ID.1,]
    samDf.2 <- samDf[samDf$sample.submitter_id%in%group2$ID.1,]
    if(nrow(samDf.1)>opt$maxPerGroup){
        if(opt$sampleToMax.SeedLogic){set.seed(opt$sampleToMax.Seed)}
        samDf.1 <- samDf.1[sample(1:nrow(samDf.1),size = opt$maxPerGroup),]
    }
    if(nrow(samDf.2)>opt$maxPerGroup){
        if(opt$sampleToMax.SeedLogic){set.seed(opt$sampleToMax.Seed*2+1)} # Just avoiding using numbers I am likely to change opt$sampleToMax.Seed to
        samDf.2 <- samDf.2[sample(1:nrow(samDf.2),size = opt$maxPerGroup),]
    }
    samDf<-rbind(samDf.1,samDf.2)
}

###### Determine corresponding barcodes   ######
desiredBarcodes <- samDf$cases
class(desiredBarcodes)
barCode

###### Search for just these barcodes     ######
query2 <- GDCquery(
    project = c("TCGA-PCPG",'TCGA-ACC'),
    data.category = "Transcriptome Profiling",
    experimental.strategy = "RNA-Seq",
    data.type = "Gene Expression Quantification",
    barcode = desiredBarcodes,
    sample.type = "Primary Tumor"
)

##### Download and prepare the data    ######
## Adjust the storage directory if needed 
if(!file.exists(opt$gdcPath)){
    if(opt$createGdcDir){
        warning("Creating directory to store GDC data in the wd.")
        opt$gdcPath <- "./GDCdata/"
    }else{
        stop("Oh, no! The file you are trying to store GDCdata in doesn't exist.")
    }
}

## Download data
GDCdownload(query2,method = "api",
            files.per.chunk = 10,
            directory = opt$gdcPath)

###### Prepare the downloaded data        ######
dds <- GDCprepare(query = query2,directory = opt$gdcPath)

## Check if it worked
View(as.data.frame(colData(dds)))

# saving barcodes analyzed
barcodeData <- sampleData[,c(1, 175)]
class(barcodeData)

which(colnames(barcodeCorrected) == 'paper_ADS')
which(colnames(barcodeCorrected) == 'ajcc_pathologic_stage')
barcodeLess <- barcodeCorrected[,!(55:127)]

write.csv(resOutput, "Table_GeneData.csv", row.names=FALSE)

write.csv(barcodeData, 'Table_BarcodeData.csv', row.names = F)
#---------------------------------------------#
#### Filter the samples and loci analyzed  ####
#---------------------------------------------#
###### Notes on exploring the summarized experiment object
## Visualize some of the clinical data
# View(as.data.frame(colData(dds)))

## Visualize the genomic data
View(as.data.frame(rowRanges(dds)))

## See what data has been added to the observation section
names(assays(dds))

## See what the unnormalized counts (what we want) look like.
View(as.data.frame(assays(dds)[["counts"]]))

## See what the tpm-transformed values look like (we won't use these).
View(as.data.frame(assays(dds)[["tpm_unstrand"]]))

###### Remove multiple observations / undesired tissues, etc #####
if(sum(duplicated(dds$sample_submitter_id))>0){
    warning("Some IDs present more than once. Consider subsetting if this was unintentional.")
}
if(length(unique(dds$sample_type))>1){
    warning("More than one type of sample tissue type present. Consider subsetting if this was unintentional.")
}

###### Identify the groupings of samples        ######
## Create new "columns" in the sample data based on grouping
dds$group1 <- dds$sample_submitter_id%in%group1$ID.1
dds$group2 <- dds$sample_submitter_id%in%group2$ID.1

## These should be the opposite of each other if barcode based subsetting was successful
if(!all(dds$group1==!dds$group2)){
    stop("Your groupings are not mutually exclusive")
}

# The way you communicate groupings to DESeq is the sample descriptor's
#    column name corresponding to a factor vector.
# Deseq2 uses the order of factor levels to set "control" vs. "manipulated".
# The first level is set as the control.
# Since the log fold change is manipulated - control, this effects visualizations.

# in this siutation it would be MUC16 mut - non MUC 16 mut
dds$comp <- factor(dds$group2,levels = c("FALSE","TRUE"))

# DESeq doesn't like spaces in levels, so you can remove them at the "level" level.
## This isn't needed if not using levels with spaces (e.g. FALSE and TRUE)
## It is included here only because it's a tricky to debug error
levels(dds$comp)<-gsub(" ","_",levels(dds$comp))

###### Convert to DESeq object                  ######
dds <- DESeqDataSet(dds, design = ~ comp)


###### Normalize the read depth values          ######
dds <- estimateSizeFactors(dds)

###### Filter loci with extremely low RD ######
## Isolate the raw counts
rawCounts  <- as.data.frame(counts(dds, normalized = F))
View(as.data.frame(rawCounts[,dds$group1]==0))

## Filter based on % of each grouping with 0 RD
grp1PctWith0<-rowMeans((rawCounts[,dds$group1]==0))
grp2PctWith0<-rowMeans((rawCounts[,dds$group2]==0))
maxPct0Cutoff <- 0.9

## Visualize low rd loci
hist(c(grp1PctWith0,grp2PctWith0),1000)
abline(v = maxPct0Cutoff,col="red")

##Do the subset
pctWith0Logic <- grp1PctWith0 < maxPct0Cutoff & grp2PctWith0 < maxPct0Cutoff

dds <- dds[pctWith0Logic,]

#onesRemoved <- dds[!pctWith0Logic,]

##### Filter samples based on unnormalized RD ##### 
## Omit samples with extremely different RD 
hist(colMeans(rawCounts),1000,border="blue")
hist(colMeans(rawCounts[,dds$group1]),1000,border="red",add=T)

# badSamples <- c(
#     "TCGA-AH-6903-01A-11R-1928-07"# From: 
#which.max(colMeans(rawCounts)) # based on histogram)

##### Remove redundant count data ####
rm("rawCounts")

##### Filter samples based on normalized RD ####
## Isolate the normalized read depths
normCounts <- as.data.frame(counts(dds, normalized = T))
## Determining the most variable rows
perLocusVar   <- apply(normCounts,1,var)
## Calculating the 25% quantile of this data
pcaLocusLogic <- perLocusVar>quantile(perLocusVar,probs=0.25)
## Use a principal component analysis to visualize the distribution
pca <- prcomp(t(normCounts[pcaLocusLogic,]), scale.  = T)
names(pca)
## Calculate means from the standard normal
pca1Cutoffs <- mean(pca$x[,1])+c(-1,1)*4*sd(pca$x[,1])
# initial plot used to visualize, not used for final report
plot(pca$x[,1],pca$x[,2],asp = 1)
abline(v = pca1Cutoffs,col="red")

# Use logic to isolate and name "bad" samples 
pcaLogic <- pca$x[,1]>pca1Cutoffs[1]&pca$x[,1]<pca1Cutoffs[2]
names(which(!pcaLogic))

##### Utilize the bad sample vector to subset dds ####
dds <- dds[,!dds$barcode%in%badSamples]

##### Remove redundant count data ####
rm("normCounts")

#-------------------------------------------#
#### Do differential expression analysis ####
#-------------------------------------------#
###### Convert to DESeq object ######
dds <- DESeqDataSet(dds, design = ~ comp)

## Because we have significantly changed the loci used, recalculate normalization
dds <- estimateSizeFactors(dds)

## Run the Deseq2 package analysis
dds <- DESeq(dds)
res <- results(dds) # Organizes the results

sampleData <- as.data.frame(colData(dds))

saveRDS(sampleData,file = "colData_nmc2355.rds")

str(sampleData)
noLists <- data.frame()

for (col in names(sampleData)) {
  if (!is.list(sampleData[[col]])) {
    noLists[[col]] <- sampleData[[col]]
  }
}


write.csv(sampleData, "Table_SampleData.csv", row.names=FALSE)

##### Accumulate data into a data.frame #####
## Add in gene data from the rowRanges section of the SummarizedExperiment object
## Adding columns with the same name causes problems with ggplot2, so
## It only adds other columns
colsAdded <- !colnames(as.data.frame(rowRanges(dds)))%in%colnames(as.data.frame(res))
resOutput <- cbind(
    as.data.frame(rowRanges(dds))[,colsAdded],
    as.data.frame(res)
)
saveRDS(resOutput,file = "resOutput_nmc2355.rds")

#saving resOutput as a csv file
write.csv(resOutput, "Table_GeneData.csv", row.names=FALSE)

# Finding which rows in the data meet the conditions to be considered significantly
# differentially expressed:
conditions <- 
  -log10(resOutput$padj) > -log10(opt$fdr.cut.off) & 
  abs(resOutput$log2FoldChange) > opt$lfc.cut.off

resSignificant <- resOutput[conditions, ]

#removing redunant NA rows 
resSignificant <- resSignificant[!is.na(resSignificant$seqnames), ]  

nrow(resSignificant) #1225 rows meet the conditions of being considered significant

# ##### Check cutoffs and plot ####
# # The following loads a script file that you will make as part of a question
# # This will not work for your system
# source(opt$visScriptFile)
# volc1<-volcPlot(resOutput)
# volc1
# 
# #---------------------------------------------#
# #### Save R environment for later analysis ####
# #---------------------------------------------#
# ##### Save the plot ####
# ggsave(filename = "apcPlot.png",volc1,
#        width = 9,height = 9,dpi = 600)
# 
# 
# #You should be able to download my version of the image from the website
# save.image(paste0(
#     "afterDeseq_",
#     Sys.Date(),
#     "_",
#     as.numeric(Sys.time()),
#     ".rimage")
# )

# removing rows with NA values for every column
resOutput_Plots <- resOutput[!is.na(resOutput$seqnames), ]

#finding the ones that meet the conditions of significance 
resOutput_Plots$`Significantly Differentially Expressed?` <- ifelse((-log10(resOutput$padj) > -log10(opt$fdr.cut.off) & 
                                                                       abs(resOutput$log2FoldChange) > opt$lfc.cut.off),
                                                                    'Yes', 'No')


# 1225 observations that are significant and have greater lfc than the preprescribed 
# number 
table(resOutput_Plots$`Significantly Differentially Expressed?`)

####PLOTS####
#####Volcano Plot#####
volcanoPlot <- ggplot(data = resOutput_Plots, aes(x = log2FoldChange,
                                                  y = -log10(padj))) +
  geom_point(aes(shape = `Significantly Differentially Expressed?`, col = `Significantly Differentially Expressed?`)) +
  geom_vline(xintercept = c(-1.5, 1.5), linetype = 2, col = 'red') +
  geom_hline(yintercept =  -log10(.01), linetype = 2, col = 'red') + 
  geom_text_repel(data = subset(resOutput_Plots, -log(padj) > 15 
                                & abs(resOutput_Plots$log2FoldChange) > opt$lfc.cut.off ), 
                  aes(label = gene_name)) +
  scale_color_viridis(discrete = T, option = 'turbo') + 
  theme_bw() + 
  labs(
    y = '-log10(FDR-Adjusted p-values)',
    x = 'log2(Fold Change)', title = 'Differential Gene Expressions among TCGA-PCPG and TCGA-ACC \nwith MUC16 Mutations (12) - Without MUC16 Mutations (50)',
  ) + 
  theme(plot.title = element_text(hjust = 0.5)) 

# visualizing the plot: 
volcanoPlot

ggsave(filename = "apcPlot.png",volcanoPlot,
       width = 9,height = 9,dpi = 600)


#####PCA Plot####

# testing if the right points are being mapped 
all(rownames(pca$x) == colnames(dds))

#making new df to be able to color the points based on the comparison grouping
mutGroup = data.frame(
  mutation_group = ifelse(dds$comp == TRUE, 'TRUE', 'FALSE'),
  x = pca$x[,1],
  y = pca$x[,2] )

length(dds$comp)
length(pca$x[,1])
# ggfortify package function to plot pca 
library(ggfortify)
pcaPlot <- autoplot(pca, data = mutGroup, col = 'mutation_group') + 
  labs(title = 'Normalized Sample Expression Counts plotted onto PC1 and PC2', 
       col = 'MUC16 Mutation Present') + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme_bw()

ggsave(filename = "pcaPlot.png",pcaPlot,
       width = 9,height = 9,dpi = 600)


#####Boxplot####
colnames(normCounts) == sampleData$barcode
all(rownames(resOutput)==rownames(normCounts))

orderVec <- order(resOutput$padj, -abs(resOutput$log2FoldChange), resOutput$gene_name)
resOutput_ordered <- resOutput[orderVec, ]
#View(resOutput_ordered)
gene1 <- rownames(resOutput_ordered)[1]
wnt1Gene1 <- normCounts[gene1, sampleData$barcode]

sampleData$wnt1Gene1 <- unlist(wnt1Gene1)
View(sampleData)
sampleData$wnt1Gene1
colnames(sampleData)
sampleData$comp

wnt1_DistributionPlot <- ggplot(sampleData, aes(comp, log10(wnt1Gene1))) + 
  geom_boxplot() + 
  geom_point() +
  labs(x = 'MUC16 Mutation Present',
       y = 'log10(Normalized counts at ENSG00000125084.12)', 
       title = 'Distribution of WNT1 gene expression counts between samples with
       and without a MUC16 mutation') + 
  theme(plot.title = element_text(hjust = 0.5))


ggsave(filename = "wnt1DistributionPlot.png",wnt1DistributionPlot,
       width = 9,height = 9,dpi = 600)

ncol(sampleData)
#Subsetting sampleData so that only relevant, non-list, non-NA data is included #
barcodeData <- sampleData[,c(1, 175)]


which(colnames(sampleData) == 'comp')
# read here for is.list use: https://stackoverflow.com/questions/23880630/why-are-my-data-frame-columns-lists
list_cols <- sapply(smallerSampleData, is.list)
noListSampleData <- smallerSampleData[,!list_cols]
write.csv(noListSampleData, "Table_SampleData.csv", row.names=FALSE)
