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
#Phylogenetic Shannon: - sum_i->nspecies ((AEDi/PD)ln(AEDi/PD))
#AED = abundance-weighted evolutionary distinctiveness
#evolutionary distinctiveness(ED) = evol.distinct(type="fair.proportion")
#evol.distinct(tree, type = c("equal.splits", "fair.proportion"), scale = FALSE, use.branch.lengths = TRUE)
##This function will return a vector of evolutionary distinctivenss for every species in the given tree.

beta_MPDabund =NA, beta_MNTDabund=NA

Mean_beta_MPDabund=mean(beta_MPDabund,na.rm=T), SD_beta_MPDabund=sd(beta_MPDabund,na.rm=T), 
Mean_beta_MNTDabund=mean(beta_MNTDabund,na.rm=T), SD_beta_MNTDabund=sd(beta_MNTDabund,na.rm=T)