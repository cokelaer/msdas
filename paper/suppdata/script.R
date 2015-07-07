# author T. Cokelaer@2014
# CNORode version:
# CellNOptR version:

print("Creating Fitness figure based on the median model")
library(CNORode)
source('../plotOptimResultsPan.R')


# READs PKN
pknmodel = readSIF("PKN-YEAST.sif")
# READ data
cnolist = CNOlist("MD-YEAST.csv")
# create set of parameters for the PKN
params = createLBodeContPars(pknmodel)

# load the decision vector, which correspond to one of the long simulation. 40
# simulations were run (boxplot of MSE provided in the paper). Here, we look at
# the median one. The decision vector was extracted manually and saved into
# data.R
source("data.R")
params$parValues = decision
params$verbose=F

# run the simulation for sanity check with a very smal number of iteration
# since we start we the parametesr that are already optimised
params = parEstimationLBodeSSm(cnolist, pknmodel, params, ndiverse=5, verbose=F,
                             dim_refset=5, maxeval=5)

# new_params and params should be identical. Here we print the MSE that should
# be around 0.06 or 0.010 depending on the model
print(paste("MSE = ", params$ssm_results$fbest, sep=""))

# saving parameters in a file
df = data.frame(names=params$parNames, values=params$parValues)
write.table(df, file="params.csv")


pdf("Fitness.pdf", width=18, height=8)
output = plotLBodeFitness(cnolist, pknmodel, params,
    plotParams=list(margin=0.1, cex=.8, rotation=90))
dev.off()

print("File Fitness.pdf created")


