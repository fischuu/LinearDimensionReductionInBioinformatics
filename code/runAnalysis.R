# Run the Analysis for the manuscript
# Fischer, D., Nordhausen, K., Oja, H.: On Linear Dimension Reduction Based on
# Diagonalization of Scatter Matrices for Bioinformatics Downstream Analyses
###############################################################################

# Set the project parameter and load required parameters
  set.seed(1)
# Change this path to the current location of the cloned git repository
  projFolder <- c("/home/fischuu/git/LinearDimensionReductionInBioinformatics")

# Change this, if Figures should be exported to HDD
  exportFigures <- FALSE
  
# Import the libraries and scripts
  library(ICS)
  library(ICSNP)
  library(ICtest)
  library(SpatialNP)
  source(file.path(projFolder,"code","SIR.R"))
  
# First, import the data
  PrCa <- read.table(file.path(projFolder, "data", "qn.probes.dataset.csv"), sep=";")
  PrCa.pheno <- read.table(file.path(projFolder, "data", "persons.csv"), sep=";")

# Run the PCA
################################################################################
  
  PrCa.pca <- prcomp(PrCa)
# Store the rotated data
  PrCa.pca.rot <- PrCa.pca$x
# Store the explained variance per component
  PrCa.pca.var <- PrCa.pca$sdev^2
  PrCa.pca.cumsum <- cumsum(PrCa.pca.var)/sum(PrCa.pca.var)
  PrCa.pca.80 <- sum(PrCa.pca.cumsum<0.80)+1
  PrCa.pca.99 <- sum(PrCa.pca.cumsum<0.99)+1
  
  # Just as a reminder, to get the origial data, one need to multiply the rotated data (=x here)
  # with the transposed rotation matrix (= the matrix that contains the eigenvectors). As prcomp
  # centers the data by default, this center has to be added then after he multiplication
  #  origData <- t(t(PrCa.pca$x %*% t(PrCa.pca$rotation))+ PrCa.pca$center)
  
  # If the data is scaled in addition prcomp(x, scale=TRUE), then the matrix has to be multiplied
  # in addition to it like this
  #  origData <- t(t(PrCa.pca$x %*% t(PrCa.pca$rotation))*PrCa.pca$scale + PrCa.pca$center)
  
# Table of first 10 squared singular values of SVD:
  tableOut <- data.frame(value=PrCa.pca.var,
                         cumsum=PrCa.pca.cumsum)
  
  tableOut[1:10,]
  
# Figure of explained variance
  if(exportFigures) png(file=file.path(projFolder,"Results","PrCa-PCA-Screeplot-SVD.png"), width=1000, height=1000)
    barplot(PrCa.pca.var[1:98])
  if(exportFigures)  dev.off()

  if(exportFigures) png(file=file.path(projFolder,"Results","PrCa-PCA-Scatterplot.png"), width=1000, height=1000)
    pairs(PrCa.pca.rot[,1:4])
  if(exportFigures)  dev.off()
  
  # Run the ICS 
###############################################################################
  
# Ths ICS-step
# First run ics with the default methods, using the first 10 principle components as rather arbitrary amount of components.
# Further we use 80 percent of the explained variance (=First four components) and also 99%, to cover the whole dataset (=98 Components).
    
  FOBI <- ics(PrCa.pca.rot[,1:10])
  FOBI.80 <- ics(PrCa.pca.rot[,1:PrCa.pca.80])
  FOBI.99 <- ics(PrCa.pca.rot[,1:PrCa.pca.99])
    
# Robust ICS with different scatter functionals
  FOBI.80.robust <- ics(PrCa.pca.rot[,1:PrCa.pca.80],duembgen.shape, symm.huber)

    
# Who is deviating in the first component?
  which(ics.components(FOBI.80)[,1]>2)
  
# It sems that all persons from slides 24,25,26 are in that group...
  if(exportFigures)  png(file=file.path("Results","PrCa-ICS-Scatterplot-80.png"), width=1000, height=1000)
    plot(FOBI.80, col=as.numeric(PrCa.pheno$slideNumber>23)+1,pch=20, index=1:4, main="Based on 4 components from PCA")
  if(exportFigures)dev.off()
  
# Robust ICA
  if(exportFigures)  png(file=file.path("Results","PrCa-robustICS-Scatterplot-80.png"), width=1000, height=1000)
    plot(FOBI.80.robust, col=as.numeric(PrCa.pheno$slideNumber>23)+1, index=1:min(PrCa.pca.80,10), 
         main="Based on 80% Variancefrom PCA - Duembgen/Symm Huber")
  if(exportFigures)dev.off()
    
# Automatic detection of components with ICtest. The k=4 comes from FOBIasym, with the smallest k that does not reject the test.
  FOBIasymp(PrCa.pca.rot[,1:10], k=3)
  if(exportFigures) png(file=file.path("Results","PrCa-ICS-FOBIasymp-Scatterplot-5.png"), width=1000, height=1000)
    plot(FOBIasymp(PrCa.pca.rot[,1:5], k=2), col=as.numeric(PrCa.pheno$slideNumber>23)+1, which="k", main="Based on 10 components from PCA, dim=first non-sig.")
  if(exportFigures) dev.off()
  if(exportFigures)png(file=file.path("Results","PrCa-ICS-FOBIasymp-Scatterplot-10.png"), width=1000, height=1000)
    plot(FOBIasymp(PrCa.pca.rot[,1:10], k=4), col=as.numeric(PrCa.pheno$slideNumber>23)+1, which="k", main="Based on 10 components from PCA, dim=first non-sig.")
  if(exportFigures)dev.off()