##################################################################################################################
# Paper: Property Prediction of Crystalline Solids from Composition and Crystal Structure                        #
# Authors: Bruno A. Calfa (bruno.calfa@wisc.edu), John R. Kitchin (jkitchin@andrew.cmu.edu)                      #
# Supplementary Material: ExhaustiveEnumerationSearch.R                                                          #
##################################################################################################################

# Clear all variables
rm(list = ls())

# Packages
library(R.utils)
library(foreach)

# Parallel support
if (.Platform$OS.type == "windows") { # Windows platform
  library(doParallel)
} else { # Unix-like platform
  library(doMC)
}

# Flags
do.multicore <- TRUE
# do.multicore.cores <- ifelse(do.parallel,2^seq.int(0,log(detectCores(),2)),1)
do.multicore.cores <- ifelse(do.multicore,parallel:::detectCores(),c(1))

# List of response variables
list_responses <- c('energy','formation_energy_per_atom','density',
                    'band_gap','efermi','total_magnetization')

# Load data
datafile <- "metal_oxides.txt"
data.df.full <- read.table(datafile,header=TRUE)
data.df.full <- data.df.full[, colSums(data.df.full != 0,na.rm=TRUE) > 0] # Remove columns with zeros

# List of predictors names
PredictorNames <- c()
for (column in names(data.df.full)) {
  if (!(column %in% list_responses) && (column != "pretty_formula") && (column != "task_id")) {
    PredictorNames <- c(PredictorNames, column)
  }
}

# Obtain mapping of unique (i,p) pairs
X <- data.df.full[,PredictorNames] # Matrix of predictor data

Xunique <- matrix(0,nrow=nrow(X),ncol=ncol(X))
colnames(Xunique) <- colnames(X)
Xunique[1,] <- 1 # First row elements are unique by convention
for (p in 1:ncol(Xunique)) {
  uniqueels <- X[1,p]
  for (i in 2:nrow(Xunique)) {
    if (!any(is.element(X[i,p],uniqueels))) {
      Xunique[i,p] <- 1
      uniqueels <- cbind(uniqueels,X[i,p])
    }
  }
}
IPunique <- matrix(1,nrow=nrow(X),ncol=ncol(X))
colnames(IPunique) <- colnames(X)
for (i in 2:nrow(IPunique)) {
  for (ii in 1:nrow(IPunique)) {
    for (p in 1:ncol(IPunique)) {
      if (Xunique[ii,p] == 1 && X[i,p] == X[ii,p]) {
        IPunique[i,p] <- ii
      }
    }
  }
}

# List of responses to match and their targets
list_responses.match <- c('band_gap','efermi','total_magnetization')
Y.match <- data.df.full[,list_responses.match]

for (ncores in do.multicore.cores) {
  if (.Platform$OS.type == "windows") { # Windows platform
    cl <- makeCluster(ncores,outfile="")
    expList <- list("X","Y","lambda","Xunique","IPunique","Y.match")
    clusterExport(cl,expList,environment())
    registerDoParallel(cl)
  } else { # Unix-like platform
    registerDoMC(ncores)
  }
  
  ee.pred <- list() # List of predicted responses from exaustive enumeration
  
  # For each data point row, loop over metals (with Xunique == 1) and form all combinations of metal oxides
  t <- system.time(ee.pred <- foreach (i = 1:nrow(Xunique),.combine=c) %dopar% {
    if (Xunique[i,"spacegroup.number"] == 1) {
      ee.pred.tmp <- list()
      ee.pred.ind <- 1
      source("regkernel.R")
      
      for (ii in 1:nrow(Xunique)) {
        Xunique.row <- Xunique[ii,]
        for (p in 1:ncol(Xunique)) {
          if (names(Xunique.row)[p] != "spacegroup.number" && names(Xunique.row)[p] != "O" && Xunique[ii,p] == 1 && X[ii,p] != 0) {
            for (iii in 1:nrow(Xunique)) {
              if (Xunique[iii,"O"] == 1) {
                z <- matrix(1,nrow=nrow(IPunique),ncol=ncol(IPunique))
                colnames(z) <- colnames(X)
                
                z[IPunique[,"spacegroup.number"] == IPunique[i,"spacegroup.number"],"spacegroup.number"] <- 0
                z[IPunique[,p] == IPunique[ii,p],p] <- 0
                z[IPunique[,"O"] == IPunique[iii,"O"],"O"] <- 0
                
                for (pp in 1:ncol(Xunique)) {
                  if (pp != p) {
                    z[X[,pp] == 0,pp] <- 0
                  }
                }
                
                Ytilde <- c()
                for (response in names(Y.match)) {
                  # Data for current response
                  Y <- as.matrix(Y.match[,response])
                  lambda <- as.matrix(readRDS(paste("npreg_categorical_",response,".rds",sep=""))$bw) # RDS files are generated after executing file MetalOxidesRegressionMPI.R
                  
                  Ytilde <- c(Ytilde,regkernel(Y,lambda,z,IPunique))
                }
                
                printf("i = %d, (ii,p) = (%d,%d), iii = %d: sg = %d %s %d %s %d\n",i,ii,p,iii,X[i,"spacegroup.number"],names(Xunique[ii,p]),X[ii,p],names(Xunique[iii,"O"]),X[iii,"O"])
                print(Ytilde)
                
                ee.pred.tmp[[ee.pred.ind]] <- list(Ytilde=Ytilde,i=i,ii=ii,p=p,iii=iii)
                ee.pred.ind <- ee.pred.ind + 1
              }
            }
          }
        }
      }
      
      ee.pred.tmp
    }
  })[3]

  if (.Platform$OS.type == "windows") { # Windows platform
    stopCluster(c1)
  }
}

## Check matchings
# Test 1
match_target <- 1.23
match_tol <- 1e-9
printf("|Ytilde(band_gap) - %g| <= %g\n",match_target,match_tol)
for (l in ee.pred) {
  if ((match_target - l$Ytilde[1])^2 <= match_tol) {
    printf("%s: Ytilde = %g, sg = %d, MO = %s %d O %d\n",list_responses.match[m],l$Ytilde[1],
           X[l$i,"spacegroup.number"],names(Xunique[l$ii,l$p]),X[l$ii,l$p],X[l$iii,"O"])
  }
}

# Test 2
match_targets <- c(0.5,1.5)
printf("min Ytilde(band_gap)\n%g <= Ytilde(efermi) <= %g",match_targets[1],match_targets[2])
l.match <- list()
m <- 1
for (l in ee.pred) {
  if ((match_targets[1] <= l$Ytilde[2]) && (l$Ytilde[2] <= match_targets[2])) {
    l.match[[m]] <- l
    m <- m + 1
  }
}
l.match.min <- c(Inf,0)
m <- 1
for (l in l.match) {
  if (l$Ytilde[2] <= l.match.min[1]) {
    l.match.min[1] <- l$Ytilde[2]
    l.match.min[2] <- m
  }
  m <- m + 1
}
l.match.res <- l.match[[l.match.min[2]]]
printf("Ytilde(band_gap) = %g, Ytilde(efermi) = %g, sg = %d, MO = %s %d O %d\n",l.match.res$Ytilde[1],l.match.res$Ytilde[2],
       X[l.match.res$i,"spacegroup.number"],names(Xunique[l.match.res$ii,l.match.res$p]),X[l.match.res$ii,l.match.res$p],X[l.match.res$iii,"O"])

# Test 3
match_targets <- c(2,1.5)
match_tol <- 1E-5
printf("Ytilde(band_gap) <= %g\n|Ytilde(efermi) - %g| <= %g",match_targets[1],match_targets[2],match_tol)
for (l in ee.pred) {
  if ((l$Ytilde[1] <= match_targets[1]) && (abs(l$Ytilde[2] - match_targets[2]) <= match_tol)) {
    printf("Ytilde(band_gap) = %g, Ytilde(efermi) = %g, sg = %d, MO = %s %d O %d\n",l$Ytilde[1],l$Ytilde[2],
           X[l$i,"spacegroup.number"],names(Xunique[l$ii,l$p]),X[l$ii,l$p],X[l$iii,"O"])
  }
}