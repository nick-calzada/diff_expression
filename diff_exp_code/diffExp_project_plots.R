
##### Filter samples based on normalized RD ####
## Isolate the normalized read depths
normCounts <- as.data.frame(counts(dds, normalized = T))
## Determining the most variable rows
perLocusVar   <- apply(normCounts,1,var)
## Calculating the 25% quantile of this data
pcaLocusLogic <- perLocusVar>quantile(perLocusVar,probs=0.25)
## Use a principal component analysis to visualize the distribution
pca <- prcomp(t(normCounts[pcaLocusLogic,]),scale.  = T)
#given code in assignment 
pca1Cutoffs <- mean(pca$x[,1])+c(-1,1)*4*sd(pca$x[,1])
plot(pca$x[,1],pca$x[,2],asp = 1, col = )
abline(v = pca1Cutoffs,col="red")

# subsetting in order to be able to color points
mutGroup <- as.data.frame(ifelse(dds$comp == TRUE, 'TRUE', 'FALSE'))
# renaming the column name 
colnames(mutGroup)[1] = "mutation_group"

#
autoplot(pca, data = mutGroup, col = 'mutation_group') + 
  labs(title = 'Normalized sample expression counts plotted onto PC1 and PC2', 
       col = 'MUC16 Mutation Present')
