#if rooted = F: get a trifurcating bottom, instead of a bifurcating one; the default is rooted = TRUE

#gives a 7x7 matrix, but i need a 7x2 matrix...
cophenetic(rcoal(nspecies)) 


phylo<-as.phylo(hclust(daisy(cbind(eOptimum,eff_values)),method ="ward"))

#Spectral Density
#first need to construct the modified graph laplacian
#degree matrix: the diagonal matrix where diagonal element i is the sum of the branch lengths from node i to all of the other nodes in the phylogeny
	##the rowsum or colsum of the distance matrix is each (i.i) element
#distance matrix: where element (i,j) is the negative of the branch-length b/w nodes i and j

#making the degree matrix
#...i don't think that this is feasible for me to do quickly because this wants all of the nodes - not just dealing with the tips..., wasting too much time
#branching.times() gives the distance from the root to the ith internal node, ordered from farthest from root to closest to root

#randomly generated phylogenies
nspecies <- 7
nreplicates <- 1000
fakecom <- matrix(Data = rep(1,nspecies), ncol = nspecies, nrow = 1)
colnames(fakecom)<-1:nspecies
mpdphylo <- rep(NA,nreplicates)
mntdphylo <- rep(NA,nreplicates)
phlgs <- array(data = NA, dim=c(nspecies, nspecies, nreplicates))
for(i in 1:nreplicates){
	placeholder = cophenetic(rcoal(nspecies))
	colnames(placeholder)<-1:nspecies
	rownames(placeholder)<-1:nspecies
	mpdphylo[i] = mpd(fakecom,placeholder, abundance.weighted = F)
	mntdphylo[i] = mntd(fakecom,placeholder, abundance.weighted = F)
	phlgs[,,i] <- placeholder
}
plot(mpdphylo,mntdphylo)

#randomly generated phylogenies, constrained within model parameters
nspecies <- 7
nreplicates <- 1000
eOptimum<-1-seq(0,eAMP, by=eAMP/(nspecies-1))
fakecom <- matrix(Data = rep(1,nspecies), ncol = nspecies, nrow = 1)
colnames(fakecom)<-1:nspecies
model_mpdphylo <- rep(NA,nreplicates)
model_mntdphylo <- rep(NA,nreplicates)
model_phlgs <- array(data = NA, dim=c(nspecies, nspecies, nreplicates))
for(i in 1:nreplicates){
	eff_values<-rnorm(nspecies,mean=0.2,sd=0.005)
	phylo<-as.phylo(hclust(daisy(cbind(eOptimum,eff_values)),method ="ward"))
	placeholder = cophenetic(phylo)
	colnames(placeholder)<-1:nspecies
	rownames(placeholder)<-1:nspecies
	model_mpdphylo[i] = mpd(fakecom,placeholder, abundance.weighted = F)
	model_mntdphylo[i] = mntd(fakecom,placeholder, abundance.weighted = F)
	model_phlgs[,,i] <- placeholder
}
plot(model_mpdphylo,model_mntdphylo)

#simulating traits from phylogenies
##use sim.char() need a trait variance-covariance matrix and a phylogeny
uncorr <- matrix(data = c(1,0,0,1), nrow = 2, ncol = 2)
phylo <- rcoal(nspecies)
sim.char(phylo,uncorr,nsim = 1, model = "BM", root = 1)

#simulating traits from restricted phylogenies, nope that's no good 
eAMP <- 1
nspecies <- 7
ntraits <- 2
nreplicates <- 100
model_traitvals <- array(data = NA, dim=c(nspecies, ntraits, nreplicates))
for(i in 1:nreplicates){
eOptimum<-1-seq(0,eAMP, by=eAMP/(nspecies-1))
eff_values<-rnorm(nspecies,mean=0.2,sd=0.005)
phylo<-as.phylo(hclust(daisy(cbind(eOptimum,eff_values)),method ="ward"))
uncorr <- matrix(data = c(1,0,0,1), nrow = 2, ncol = 2)	
model_traitvals[,,i] <- sim.char(phylo,uncorr,nsim = 1, model = "BM", root = 1)
}
