library(stringr)

cmd <- commandArgs(trailingOnly = TRUE)

cmd[[4]] = str_trim(cmd[[4]])

set.seed(as.numeric(cmd[1]))

nrep = as.numeric(cmd[2])

p = as.numeric(cmd[3])

pDEG = as.numeric(cmd[4])

library(MOSim)

# Now we normalize per sample first and then per regulator
discretize_countsCOLROW <- function(df) {
  ## Transform into relative percentages per sample (sum of column is 1)
  df <- df/colSums(df)
  ## And now per gene
  df <- df/rowSums(df)
  # Threshold is the relative abundance if the genes are the same in the samples
  threshold <- 1/length(colnames(df))
  # If below the threshold, 0, if above 1
  for (i in 1:ncol(df)) {
    df[,i] <- ifelse(df[,i] < threshold, 0, 1)
  }
  return(df)
}

#Ask for the characteristics we want to simulate in each omic

omics_list <- c("RNA-seq",'ChIP-seq','DNase-seq', "miRNA-seq","Methyl-seq")
omics_options <- c(omicSim('DNase-seq', regulatorEffect = list('activator' = p,'repressor' = p,'NE' = 1-p-p )),
                   omicSim('miRNA-seq', regulatorEffect = list('activator' = 0,'repressor' = 2*p, 'NE' = 1-p-p)),
                   omicSim('ChIP-seq', regulatorEffect = list('activator' = p,'repressor' = p, 'NE' = 1-p-p)),
                   omicSim('Methyl-seq', regulatorEffect = list('activator' = p,'repressor' = p, 'NE'=1-p-p)))

ti = c(1:nrep) 
rnaseq_simulation <- mosim(omics = omics_list, omicsOptions = omics_options ,times = ti, depth = 30, numberReps = 1, numberGroups = 2, diffGenes=pDEG,TFtoGene = TRUE)
# Get the count tables
rnaseq_simulated <- omicResults(rnaseq_simulation, omics_list)
# get the settings used to generate each count table
all_settings <- omicSettings(rnaseq_simulation)
design_matrix <- experimentalDesign(rnaseq_simulation)

#Save the simulations
saveRDS(rnaseq_simulated, paste("~/MORE/data/my_simulated_bulkMOSim", 100*pDEG,100*p, cmd[1], nrep, "rep.rds", sep = '_'))
saveRDS(all_settings, paste("~/MORE/data/my_settings",100*pDEG,100*p, cmd[1], nrep, "rep.rds", sep = '_'))
saveRDS(rnaseq_simulation, paste("~/MORE/data/my_simulation_bulkMOSim",100*pDEG,100*p, cmd[1], nrep, "rep.rds", sep = '_'))
saveRDS(design_matrix, paste("~/MORE/data/edesign", 100*pDEG,100*p, cmd[1], nrep, "rep.rds",sep = '_'))

