require(picante) #if not more...
#comdist: calculates inter-community mean pairwise distance
#comdist(comm, dis, abundance.weighted = FALSE)
#comm = community data matrix, dis = interspecific distance matrix

#comdistnt: calculates inter-community mean nearest taxon distance
#comdistnt(comm, dis, abundance.weighted = FALSE, exclude.conspecifics = FALSE)
#exclude.conspecifics: Should conspecific taxa in different communities be exclude from MNTD calculations? (default = FALSE)

commdata_array <- array(data = NA, dim=c(10,7,9))
mpd_beta <- rep(NA,length(DispV))
mntd_beta <- rep(NA,length(DispV))
for(i in 1:length(DispV)){ #go through all of the dispersal levels
commdata <- SIH_data[["Abund",i]] #extracts the abundance data at dispersal level i from SIHdata
commdata_array[,,i] <- t(commdata[400,,]) #creates a community data matrix from the community values at the last time step, rows = patches, columns = species - for all dispersal levels
commdata_matrix = commdata_array[,,i]
colnames(commdata_matrix) <- paste('species', 1:nspecies)
interspec_mat <- cophenetic(SIH_data[["phylo",i]])
colnames(interspec_mat) <- paste('species', 1:7)
rownames(interspec_mat) <- paste('species', 1:7)
mpd_beta[i] <- mean(comdist(commdata_matrix,interspec_mat,abundance.weighted=TRUE))

mntd_beta[i] <- mean(comdistnt(commdata_matrix,interspec_mat,abundance.weighted=TRUE,exclude.conspecifics=TRUE))
#gives this warning message: In FUN(newX[, i], ...) : no non-missing arguments to min; returning Inf
#this warning message apparently means that there's some vector of length 0?
}

#Hill Number variant
#AED = abundance-weighted evolutionary distinctiveness

require(pez)
#first, need to create a comparative.comm framework
#comparative.comm(phy, comm, traits = NULL, env = NULL, warn = TRUE,force.root = -1)
##need only supply comm and phy, nothing else is mandatory
##comm needs rownames and colnames, is a community data matrix
###it's telling me that the community data has no sites in common with the rest of the data...not really sure what to do about that ####:((####
##phy needs a phylogeny (in phylo format)

#Phylogenetic Shannon: - sum_i->nspecies ((AEDi/PD)ln(AEDi/PD))
#Hill's diversity number of order a: (sum_i((AEDi/PD)^a))^(1/(1-a))
##...I think? I substituted p_i for AEDi/PD as was done in the Shannon formula

##need to run the full model before running the code below
a = 2 #inverse of Gini-Simpson or something
commdata_array <- array(data = NA, dim=c(npatches,nspecies,length(DispV))) 
#phlgshannon_alpha <- rep(NA,npatches)
phlgshannon_alpha <- matrix(data = rep(0,npatches*length(DispV)), nrow = npatches, ncol = length(DispV))
shannonhillnum <- matrix(data = rep(0,npatches*length(DispV)), nrow = npatches, ncol = length(DispV))
generalhillnum <- matrix(data = rep(0,npatches*length(DispV)), nrow = npatches, ncol = length(DispV))
for(i in 1:length(DispV)){ #go through all of the dispersal levels
commdata <- SIH_data[["Abund",i]] #extracts the abundance data at dispersal level i from SIHdata
commdata_array[,,i] <- t(commdata[400,,]) #creates a community data matrix from the community values at the last time step, rows = patches, columns = species - for all dispersal levels
commdata_matrix = commdata_array[,,i]
rownames(commdata_matrix) <- paste('patch',1:10)
#colnames(commdata_matrix) <- paste('species', 1:nspecies)
colnames(commdata_matrix) <- 1:nspecies
#SIH_data[["phylo",i]]$tip.label <- paste('species', 1:nspecies)
SIH_data[["phylo",i]]$tip.label <- 1:nspecies
compcommdata <- comparative.comm(SIH_data[["phylo",i]],commdata_matrix)
for(j in 1:npatches){
	for(k in 1:nspecies){ #HAVE TO ADD ABUNDANCE IN b/c it's (AED_i*n_i)/PD
		pptn <- (.aed(compcommdata)[k,j])/pd(commdata_matrix,SIH_data[["phylo",i]])$PD[j]
		phlgshannon_alpha[j,i][is.na(phlgshannon_alpha[j,i])] <- 0 #checks if NA and replaces with 0
		shannonindex <- pptn*log(pptn)
		shannonindex[is.na(shannonindex)] <- 0
		phlgshannon_alpha[j,i] <- pptn*log(pptn) + phlgshannon_alpha[j,i]
		generalhillnum[j,i] <- (pptn^a)^(1/(1-a))
	}
shannonhillnum[j,i] <- exp(phlgshannon_alpha[j,i])
}
#.aed(comparativecommdata) #abundance-weighted evolutionary distinctiveness
	#if(i == 1){
		#listofcomms <- comparativecommdata
	#}
#listofcomms <- c(listofcomms,comparativecommdata)
}

##ignore - applies to ecoPD, not to pez
#evolutionary distinctiveness(ED) = evol.distinct(type="fair.proportion")
#evol.distinct(tree, type = c("equal.splits", "fair.proportion"), scale = FALSE, use.branch.lengths = TRUE)
##This function will return a vector of evolutionary distinctivenss for every species in the given tree.


#'raw' beta diversity
shan_alpha_hillnum <- matrix(data = rep(0,npatches*length(DispV)), nrow = npatches, ncol = length(DispV))
shannon_alpha <- matrix(data = rep(0,npatches*length(DispV)), nrow = npatches, ncol = length(DispV))
shannon_gamma <- rep(0, length = nspecies)
for(k in 1:npatches){
	for(m in 1:nspecies){
		relabund <- com_data[k,m]/sum(com_data[k,])
		shannon_alpha[k,i][is.na(shannon_alpha[k,i])] <- 0 #checks if NA and replaces with 0
		shannon_alpha[k,i] <- relabund*log(relabund) + shannon_alpha[k,i]
		gammaabund <- sum(com_data[,m])/sum(com_data)
		shannon_gamma[i][is.na(shannon_gamma[i])] <- 0
		shannon_gamma[i] <- gammaabund*log(gammaabund) + shannon_gamma[i]
	}
	shan_alpha_hillnum[k,i] <- exp(shannon_alpha[k,i])
	
}
shan_gamma_hillnum[i] <- exp(shannon_gamma[i])
