##################################################################################################################
# Paper: Property Prediction of Crystalline Solids from Composition and Crystal Structure                        #
# Authors: Bruno A. Calfa (bruno.calfa@wisc.edu), John R. Kitchin (jkitchin@andrew.cmu.edu)                      #
# Supplementary Material: AllElasticityRegressionMPI.R                                                           #
##################################################################################################################

# Kernel types
ckertype <- "epanechnikov"
ukertype <- "liracine"

# List of response variables
list_responses <- c('energy','formation_energy_per_atom','density',
                    'G_Reuss','G_VRH','G_Voigt','K_Reuss','K_VRH','K_Voigt','poisson_ratio','universal_anisotropy')
list_responses.len <- length(list_responses)

# Load data
datafile <- "allelasticity.txt"
data.df.full <- read.table(datafile,header=TRUE)
data.df.full <- data.df.full[, colSums(data.df.full != 0,na.rm=TRUE) > 0] # Remove columns with zeros

# List of predictors names
PredictorNames <- c()
for (column in names(data.df.full)) {
  if (!(column %in% list_responses) && (column != "pretty_formula") && (column != "task_id")) {
    PredictorNames <- c(PredictorNames, column)
  }
}

# Timing variable
tnpyhat.fact <- numeric(list_responses.len)

# Initialize master and slaves
mpi.bcast.cmd(np.mpi.initialize(),caller.execute=TRUE)

# Broadcast options and parameters
mpi.bcast.cmd(options(np.tree=TRUE),caller.execute=TRUE)
mpi.bcast.Robj2slave(list_responses.len)

# Loop over all responses
for (i in 1:list_responses.len) {
  response <- list_responses[i]

	# Data set for current response and predictors
  data.df <- subset(data.df.full, select=-c(pretty_formula))
  data.df <- data.df[c(response,PredictorNames)]
  
  # Results files
  KR.outfile.bw <- paste("mp_allelasticity_KR_bw_",response,".txt",sep="")
  KR.outfile.resid <- paste("mp_allelasticity_KR_resid_",response,".txt",sep="")
  
  # Regression formula
  Formula <- formula(paste(response, "~", paste(PredictorNames, collapse=" + ")))
  
  # Factor predictors
  data.df.fact <- data.df
  for (column in names(data.df)) {
    if (!(column %in% list_responses) && (column != "pretty_formula")) {
      if (column != "spacegroup.number") {
        data.df.fact[column] <- factor(data.df.fact[column][,1])
      } else {
        data.df.fact[column] <- factor(data.df.fact[column][,1],levels=1:230)
      }
      contrasts(data.df.fact[column][,1]) <- contr.sum(length(levels(data.df.fact[column][,1]))) # Coding system
    }
  }
  
  # Broadcast data to all slaves (it will be known to the master node)
  mpi.bcast.Robj2slave(data.df.fact)
  mpi.bcast.Robj2slave(Formula)
  mpi.bcast.Robj2slave(ukertype)
  mpi.bcast.Robj2slave(ckertype)
  
  # Kernel regression with continuous predictors
  tnpyhat.fact[i] <- system.time(mpi.bcast.cmd(npyhat.fact <- npreg(Formula,data=data.df.fact,ukertype=ukertype,ckertype=ckertype,residuals=TRUE),caller.execute=TRUE))[3]
  print(summary(npyhat.fact))
  saveRDS(npyhat.fact,paste("npreg_categorical_",response,".rds",sep="")) # Serialize npreg object
  
  cat("Elapsed time =", tnpyhat.fact[i], "\n")
  
  # Write estimated coefficients and their confidence intervals to file
  if (file.exists(KR.outfile.bw)) {
    # Create backup copy if file already exists
    file.rename(KR.outfile.bw,paste(KR.outfile.bw,".bak",sep=""))
  }
  write(npyhat.fact$bw,KR.outfile.bw)
  if (file.exists(KR.outfile.resid)) {
    # Create backup copy if file already exists
    file.rename(KR.outfile.resid,paste(KR.outfile.resid,".bak",sep=""))
  }
  write.table(cbind(npyhat.fact$mean,npyhat.fact$resid),KR.outfile.resid,col.names=FALSE,row.names=FALSE,quote=FALSE,sep='\t')
	
	if (i == list_responses.len) break
}

# Print times
print(tnpyhat.fact)

# Clean up properly then quit()
mpi.close.Rslaves()
mpi.quit()
 
