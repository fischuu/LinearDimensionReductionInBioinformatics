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
  PrCa.pca.99 <- sum(PrCa.pca.cumsum<0.9999)+1
  
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

# Data for Table 1  
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

# Data for table 2:
  apply(cbind(FOBI.80@gKurt, FOBI.80.robust@gKurt),2,round,4)

# Data for table 3  
    round(c(FOBIasymp(PrCa.pca.rot[,1:PrCa.pca.80], k=0, model="ICA")$p.value,
          FOBIasymp(PrCa.pca.rot[,1:PrCa.pca.80], k=1, model="ICA")$p.value,
          FOBIasymp(PrCa.pca.rot[,1:PrCa.pca.80], k=2, model="ICA")$p.value),4)

  
  
  # SIR for cancer class  
################################################################################    
set.seed(1234)
# Getting the phenotype class
  y <- PrCa.pheno[,7]  
    
# Performing the SIR (80 refers to 4 components and 99 to 98, depending on the variance explained by these components)
  RES.SIR.10 <- SIR(PrCa.pca.rot[,1:10], y, natH=TRUE)
  RES.SIR.80 <- SIR(PrCa.pca.rot[,1:PrCa.pca.80], y, natH=TRUE)
  RES.SIR.99 <- SIR(PrCa.pca.rot[,1:PrCa.pca.99], y, natH=TRUE)

# Plotting the results
    if(exportFigures)  png(file=file.path("Results","PrCa-SIR-Scatterplot-99.png"), width=1000, height=1000)
      pairs(cbind(y, RES.SIR.99$S[,1:min(5,PrCa.pca.99)]), col=as.numeric(PrCa.pheno[,7])+1)
    if(exportFigures)dev.off()

    if(exportFigures)  png(file=file.path("Results","PrCa-SIR-Scatterplot-80.png"), width=1000, height=1000)
      pairs(cbind(y, RES.SIR.80$S[,1:min(5,PrCa.pca.80)]), col=as.numeric(PrCa.pheno[,7])+1)
    if(exportFigures)dev.off()
      
# Data for table 4
      