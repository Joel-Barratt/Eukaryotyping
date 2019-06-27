require(stringr)
require(gtools)
require(parallel)

source("import_data_V2.r")

# only clean data
#data = cleandata
# only  first 10
#data = data[1:10,]
ids = data$ids
nids = length(ids)

source("euk_bayesian_fulldataset_V2.r")

source("euk_heuristic_fulldataset.r")

## average both methods

Bayesian_pairwisedistancematrix_norm_temp = (2-Bayesian_pairwisedistancematrix)
Bayesian_pairwisedistancematrix_norm_temp = Bayesian_pairwisedistancematrix_norm_temp - min(Bayesian_pairwisedistancematrix_norm_temp)
Bayesian_pairwisedistancematrix_norm = Bayesian_pairwisedistancematrix_norm_temp  / (max(Bayesian_pairwisedistancematrix_norm_temp ))

write.csv(Bayesian_pairwisedistancematrix_norm,"Bayesian_pairwisedistancematrix_norm.csv")


Heuristic_pairwisedistancematrix_norm_temp = (Heuristic_pairwisedistancematrix - min(Heuristic_pairwisedistancematrix,na.rm=TRUE)) 
Heuristic_pairwisedistancematrix_norm = Heuristic_pairwisedistancematrix_norm_temp  / (max(Heuristic_pairwisedistancematrix_norm_temp,na.rm=TRUE ))

D = sort(Heuristic_pairwisedistancematrix_norm)[rank(Bayesian_pairwisedistancematrix_norm,ties.method="random")]

ensemble_pairwisedistancematrix = matrix(NA,nids,nids)
ensemble_pairwisedistancematrix = list(Bayesian_pairwisedistancematrix_norm,Heuristic_pairwisedistancematrix_norm)

arraytemp = array(NA, dim=c(nids,nids,2))
arraytemp[,,1] = D
arraytemp[,,2] = Heuristic_pairwisedistancematrix_norm
ensemble_pairwisedistancematrix = apply(arraytemp,c(1,2),function (x) mean(x,na.rm=TRUE))
colnames(ensemble_pairwisedistancematrix) = colnames(Bayesian_pairwisedistancematrix_norm )
rownames(ensemble_pairwisedistancematrix) = rownames(Bayesian_pairwisedistancematrix_norm )

ensemble_pairwisedistancematrix = ensemble_pairwisedistancematrix/max(ensemble_pairwisedistancematrix)

write.csv(ensemble_pairwisedistancematrix,"pairwisedistancematrix_ensemble.csv")

