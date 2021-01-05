#### Calculate matrix using Barratt's Unsupervised Heuristic Mixture Model


alleles = list()
frequencies = list()

for (j in 1:nloci) {
	locicolumns = grepl(paste(locinames[j],"",sep=""),colnames(data))
	raw_alleles = c(as.matrix(data[,locicolumns]))
	raw_alleles[raw_alleles == "NA"] = NA
	raw_alleles[raw_alleles == 0] = NA
	alleles[[j]] = unique(raw_alleles[!is.na(raw_alleles)])
	frequencies[[j]] = sapply(alleles[[j]], function(x) sum(raw_alleles == x,na.rm=TRUE))
	frequencies[[j]] = frequencies[[j]] / sum(frequencies[[j]])
}


observeddatamatrix = list()
for (j in 1:nloci) {
	locus = locinames[j]
	locicolumns = grepl(paste(locus,"",sep=""),colnames(data))
	oldalleles = as.vector(data[,locicolumns])
	oldalleles [oldalleles == "NA"] = NA
	oldalleles [oldalleles == 0] = NA
	if (length(dim(oldalleles)[2]) == 0) {
		oldalleles = matrix(oldalleles,length(oldalleles),1)
	}
	observeddatamatrix[[j]] = oldalleles 
}

m <<- rep(1,nloci)

H_nu = sapply(1:nloci, function (j) -sum(frequencies[[j]] * logb(frequencies[[j]],2)))

sub_per_locus = function(isolate1,isolate2,j) {
		v1 = observeddatamatrix[[j]][isolate1,]
		v1 = v1[!is.na(v1)]
		p1 = frequencies[[j]][match(v1,alleles[[j]])]
		v2 = observeddatamatrix[[j]][isolate2,]
		v2 = v2[!is.na(v2)]
		p2 = frequencies[[j]][match(v2,alleles[[j]])]


		if (ploidy[j] > 1) {
			x = length(unique(v1)) + length(unique(v2))
			n = min(length(unique(v1)),length(unique(v2)))
			w = x * (n > 1) + 4 * (n == 1) * (x == 2) + (1+x) * (n== 1) *  (x > 2)
			jj = 2 * (m[j] == 1) + m[j] * (m[j] > 1)
			y = length(intersect(v1,v2))
			z = 3 * (((2*(n == 1) + 1 * (y == 1) * (x > 2)))==3) + 2 * (y *(n>1)*(jj>=y)+jj*(n > 1)* (y > jj) + jj * (n==1) * (x==2) * (y==1))
			delta_nu_raw = w * (y == 0) + 2 * jj * (y > 0) + sum(sapply(jj:(2*jj), function (ii) -ii*(z == ii)))
		
			shared_alleles = intersect(v1 , v2)
			if (length(shared_alleles) > 0) {
				temp_shared = sapply(1:nids, function (x) sum(shared_alleles %in% observeddatamatrix[[j]][x,]))
				notmissing = sapply(1:nids, function (x) sum(!is.na(observeddatamatrix[[j]][x,])))
			
				P_nu = (sum(temp_shared == length(shared_alleles)) / sum(notmissing != 0))^2
			} else {
				P_nu = 1
			}
			k = 1 * (y == 0) + P_nu * (y > 0)
			delta_nu = H_nu[j] * ( delta_nu_raw * (delta_nu_raw > 0) + P_nu * (delta_nu_raw == 0))*k
			delta = delta_nu
		} else {
			delta_ex_raw = 0
			x = length(unique(v1)) + length(unique(v2))
			y = length(intersect(v1,v2))
			delta_ex_raw = 2*x*(y == 0)
			shared_alleles = intersect(v1 , v2)
			if (length(shared_alleles) > 0) {
				temp_shared = sapply(1:nids, function (x) sum(shared_alleles %in% observeddatamatrix[[j]][x,]))
				notmissing = sapply(1:nids, function (x) sum(!is.na(observeddatamatrix[[j]][x,])))
				P_ex = (sum(temp_shared == length(shared_alleles)) / sum(notmissing != 0))^2
			} else {
				P_ex = 1
			}
			k = 1 * (y == 0) + P_ex * (y > 0)
			delta_ex = H_nu[j] * ( delta_ex_raw * (delta_ex_raw > 0) + P_ex * (delta_ex_raw == 0))*k
			delta = delta_ex
		}
		if (sum(!is.na(v1)) == 0 | sum(!is.na(v2)) == 0) { delta = NA }
		delta
}

pairwisedistance_heuristic = function(isolate1,isolate2){
	print(((isolate2-1)*nids+isolate1)/ (nids*nids))
	delta = sapply(1:nloci, function (x) sub_per_locus(isolate1,isolate2,x))
	c(delta,sum(delta))	
}

		       ####### MODIFY NUMBER OF CORES USED BELOW - mc.cores=##
		       
allpossiblepairs = expand.grid(1:nids,1:nids)
allpossiblepairs = unique(allpossiblepairs[allpossiblepairs[,1] <= allpossiblepairs[,2],])
pairwisedistancevector = do.call(cbind,mclapply(1:dim(allpossiblepairs)[1], function (x) pairwisedistance_heuristic(allpossiblepairs[x,1],allpossiblepairs[x,2]),mc.cores=12))

pairwisedistancematrix_components = list()
for (j in 1:(nloci+1)) { 
	pairwisedistancematrix_temp = matrix(NA,nids,nids)
	sapply(1:dim(allpossiblepairs)[1], function (x) pairwisedistancematrix_temp[allpossiblepairs[x,1],allpossiblepairs[x,2]] <<- pairwisedistancevector[j,x])
	sapply(1:dim(allpossiblepairs)[1], function (x) pairwisedistancematrix_temp[allpossiblepairs[x,2],allpossiblepairs[x,1]] <<- pairwisedistancevector[j,x])
	pairwisedistancematrix_components[[j]] = pairwisedistancematrix_temp
}


#### impute missing values
pairwisedistancematrix_components_imputed = pairwisedistancematrix_components
whichna = which(rowSums(is.na(pairwisedistancematrix_components[[nloci+1]])) == nids)

imputemissing = function(isolate1) {
	missingloci = which(sapply(1:nloci, function (j) sum(!is.na(observeddatamatrix[[j]][isolate1,]))) == 0)
	nonmissingloci = (1:nloci)[-missingloci]
	matchingsamples = which(rowSums(rbind(sapply(nonmissingloci, function (j) sapply(1:nids, function (x) (setequal(observeddatamatrix[[j]][x,],observeddatamatrix[[j]][isolate1,]))))))==length(nonmissingloci))
	matchingsamples = setdiff( matchingsamples , whichna)
	for (j in missingloci ) {
		if (length(matchingsamples) > 0) {
			sapply(1:nids, function (x) pairwisedistancematrix_components_imputed[[j]][isolate1,x] <<- mean(pairwisedistancematrix_components[[j]][x,matchingsamples],na.rm=TRUE))
			sapply(1:nids, function (x) pairwisedistancematrix_components_imputed[[j]][x,isolate1] <<- mean(pairwisedistancematrix_components[[j]][x,matchingsamples],na.rm=TRUE))
			pairwisedistancematrix_components_imputed[[j]][isolate1,isolate1] <<- mean(diag(pairwisedistancematrix_components_imputed[[j]])[matchingsamples],na.rm=TRUE)
		} else {
			pairwisedistancematrix_components_imputed[[j]][isolate1,] <<-mean(pairwisedistancematrix_components[[j]],na.rm=TRUE)
			pairwisedistancematrix_components_imputed[[j]][,isolate1] <<-mean(pairwisedistancematrix_components[[j]],na.rm=TRUE)
			pairwisedistancematrix_components_imputed[[j]][isolate1,isolate1] <<-mean(diag(pairwisedistancematrix_components_imputed[[j]]),na.rm=TRUE)
		}
	}
}
sapply(whichna, imputemissing)

temppairwisedistancematrix = matrix(0,nids,nids)
for (j in 1:(nloci)) { 
	temppairwisedistancematrix = temppairwisedistancematrix + pairwisedistancematrix_components_imputed[[j]]
}
whichna2 = which(rowSums(is.na(temppairwisedistancematrix )) != 0)

pairwisedistancematrix_components_imputed_secondpass = pairwisedistancematrix_components_imputed

imputemissing_secondpass = function(isolate1) {
	missingloci = which(sapply(1:nloci, function (j) sum(!is.na(observeddatamatrix[[j]][isolate1,]))) == 0)
	nonmissingloci = (1:nloci)[-missingloci]
	matchingsamples = which(rowSums(rbind(sapply(nonmissingloci, function (j) sapply(1:nids, function (x) (setequal(observeddatamatrix[[j]][x,],observeddatamatrix[[j]][isolate1,]))))))==length(nonmissingloci))
	matchingsamples = setdiff( matchingsamples , whichna)
	for (j in missingloci ) {
		if (length(matchingsamples) > 0) {
			sapply(1:nids, function (x) pairwisedistancematrix_components_imputed_secondpass[[j]][isolate1,x] <<- mean(pairwisedistancematrix_components_imputed[[j]][x,matchingsamples],na.rm=TRUE))
			sapply(1:nids, function (x) pairwisedistancematrix_components_imputed_secondpass[[j]][x,isolate1] <<- mean(pairwisedistancematrix_components_imputed[[j]][x,matchingsamples],na.rm=TRUE))
			pairwisedistancematrix_components_imputed_secondpass[[j]][isolate1,isolate1] <<- mean(diag(pairwisedistancematrix_components_imputed[[j]])[matchingsamples],na.rm=TRUE)
		} else {
			pairwisedistancematrix_components_imputed_secondpass[[j]][isolate1,] <<-mean(pairwisedistancematrix_components_imputed[[j]],na.rm=TRUE)
			pairwisedistancematrix_components_imputed_secondpass[[j]][,isolate1] <<-mean(pairwisedistancematrix_components_imputed[[j]],na.rm=TRUE)
			pairwisedistancematrix_components_imputed_secondpass[[j]][isolate1,isolate1] <<-mean(diag(pairwisedistancematrix_components_imputed[[j]]),na.rm=TRUE)
		}
	}
}

sapply(whichna2, imputemissing_secondpass)



# calculate final
finalpairwisedistancematrix = matrix(0,nids,nids)
for (j in 1:(nloci)) { 
	finalpairwisedistancematrix = finalpairwisedistancematrix + pairwisedistancematrix_components_imputed_secondpass[[j]]
}


#pairwisedistancematrix2 = sapply(1:nids, function (x) sapply(1:nids, function (y) pairwisedistance_heuristic(x,y)))
colnames(pairwisedistancematrix) = ids 
rownames(pairwisedistancematrix) = ids
#write.csv(finalpairwisedistancematrix,"pairwisedistancematrix_heuristic.csv")

Heuristic_pairwisedistancematrix = finalpairwisedistancematrix 


#normalized_finalpairwisedistancematrix <- finalpairwisedistancematrix/(max(finalpairwisedistancematrix))

#colnames(normalized_finalpairwisedistancematrix) <- ids
#rownames(normalized_finalpairwisedistancematrix) <- ids


#write.csv(normalized_finalpairwisedistancematrix,"Heuristic_pairwisedistancematrix_norm.csv")
print("Calculation of heuristic matrix complete")
			       
			       
			       
			       
