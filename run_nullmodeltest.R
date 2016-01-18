###null models###
#taxa.labels
#richness
#frequency
#sample.pool
#phylogeny.pool
#independentswap
#trialswap

require(geiger)
nspecies<-7 #the number of species
npatches<-10 #the number of patches
nreplicates<-100 #number of replicates
nfunctions<-1

DispV<-c(0.0001,0.0005,0.001,0.005,0.01,0.05,0.1,0.5,1) #the dispersal rates 

Data_storage<-data.frame(MPD_abund=NA,MNTD_abund=NA,sesMPDpp_abund_z = NA, sesMNTDpp_abund_z = NA, sesMPDtl_abund_z = NA, sesMNTDtl_abund_z = NA, sesMPDr_abund_z = NA, sesMNTDr_abund_z = NA, sesMPDf_abund_z = NA, sesMNTDf_abund_z = NA, sesMPDsp_abund_z = NA, sesMNTDsp_abund_z = NA, sesMPDis_abund_z = NA, sesMNTDis_abund_z = NA, sesMPDts_abund_z = NA, sesMNTDts_abund_z = NA,
Dispersal=rep(DispV,each=nreplicates),ReplicateNum=factor(1:nreplicates),Scale=rep(c("Local","Regional"),each=length(DispV)*nreplicates)) #building the data frame

MTraits<-t(matrix(1,nspecies,nfunctions))
fakecom <- matrix(1, ncol = nspecies, nrow = 1)
colnames(fakecom)<-1:nspecies
mpdphylo <- matrix(NA, nrow = length(DispV), ncol = nreplicates)
mntdphylo <- matrix(NA, nrow = length(DispV), ncol = nreplicates)
phlgs <- array(data = NA, dim=c(nspecies, nspecies, nreplicates, length(DispV)))

for(j in 1:nreplicates){
#runs the SIH model at all dispersal rates in DispV and saves the abundances and productivity in a list
set.seed(j)
eff_sd <- 0.005
eff_mean <- 0.2
SIH_data<-sapply(DispV,SIH_function,species=nspecies,patches=npatches,eff_vary=T,eff_sd=eff_sd, eff_mean = eff_mean)

for(i in 1:length(DispV)){
   #At the local scale...
    interspec_mat <- cophenetic(SIH_data[["phylo",i]])
     colnames(interspec_mat) <- 1:nspecies
     rownames(interspec_mat) <- 1:nspecies	
  com_data<-t(SIH_data[["Abund",i]][400,,])
  colnames(com_data)<-1:nspecies
  SIH_data[["phylo",i]]$tip.label <- 1:nspecies
  
  Data_storage$MPD_abund[Data_storage$Dispersal==DispV[i] & Data_storage$ReplicateNum==j & Data_storage$Scale == "Local"]<-mean(mpd(com_data,cophenetic(SIH_data[["phylo",i]]),abundance.weighted = T)) 
  Data_storage$MNTD_abund[Data_storage$Dispersal==DispV[i] & Data_storage$ReplicateNum==j & Data_storage$Scale == "Local"]<-mean(mntd(com_data,cophenetic(SIH_data[["phylo",i]]),abundance.weighted = T))
      
   Data_storage$sesMPDpp_abund_z[Data_storage$Dispersal==DispV[i] & Data_storage$ReplicateNum==j & Data_storage$Scale == "Local"]<-mean(ses.mpd(samp = com_data,dis = interspec_mat,null.model = "phylogeny.pool",abundance.weighted = T, runs = 999)$mpd.obs.z)

Data_storage$sesMNTDpp_abund_z[Data_storage$Dispersal==DispV[i] & Data_storage$ReplicateNum==j & Data_storage$Scale == "Local"]<-mean(ses.mntd(samp = com_data,dis = interspec_mat,null.model = "phylogeny.pool",abundance.weighted = T, runs = 999)$mntd.obs.z)

Data_storage$sesMPDtl_abund_z[Data_storage$Dispersal==DispV[i] & Data_storage$ReplicateNum==j & Data_storage$Scale == "Local"]<-mean(ses.mpd(samp = com_data,dis = interspec_mat,null.model = "taxa.labels",abundance.weighted = T, runs = 999)$mpd.obs.z)

Data_storage$sesMNTDtl_abund_z[Data_storage$Dispersal==DispV[i] & Data_storage$ReplicateNum==j & Data_storage$Scale == "Local"]<-mean(ses.mntd(samp = com_data,dis = interspec_mat,null.model = "taxa.labels",abundance.weighted = T, runs = 999)$mntd.obs.z)

Data_storage$sesMPDr_abund_z[Data_storage$Dispersal==DispV[i] & Data_storage$ReplicateNum==j & Data_storage$Scale == "Local"]<-mean(ses.mpd(samp = com_data,dis = interspec_mat,null.model = "richness",abundance.weighted = T, runs = 999)$mpd.obs.z)

Data_storage$sesMNTDr_abund_z[Data_storage$Dispersal==DispV[i] & Data_storage$ReplicateNum==j & Data_storage$Scale == "Local"]<-mean(ses.mntd(samp = com_data,dis = interspec_mat,null.model = "richness",abundance.weighted = T, runs = 999)$mntd.obs.z)

Data_storage$sesMPDf_abund_z[Data_storage$Dispersal==DispV[i] & Data_storage$ReplicateNum==j & Data_storage$Scale == "Local"]<-mean(ses.mpd(samp = com_data,dis = interspec_mat,null.model = "frequency",abundance.weighted = T, runs = 999)$mpd.obs.z)

Data_storage$sesMNTDf_abund_z[Data_storage$Dispersal==DispV[i] & Data_storage$ReplicateNum==j & Data_storage$Scale == "Local"]<-mean(ses.mntd(samp = com_data,dis = interspec_mat,null.model = "frequency",abundance.weighted = T, runs = 999)$mntd.obs.z)

Data_storage$sesMPDsp_abund_z[Data_storage$Dispersal==DispV[i] & Data_storage$ReplicateNum==j & Data_storage$Scale == "Local"]<-mean(ses.mpd(samp = com_data,dis = interspec_mat,null.model = "sample.pool",abundance.weighted = T, runs = 999)$mpd.obs.z)

Data_storage$sesMNTDsp_abund_z[Data_storage$Dispersal==DispV[i] & Data_storage$ReplicateNum==j & Data_storage$Scale == "Local"]<-mean(ses.mntd(samp = com_data,dis = interspec_mat,null.model = "sample.pool",abundance.weighted = T, runs = 999)$mntd.obs.z)

Data_storage$sesMPDis_abund_z[Data_storage$Dispersal==DispV[i] & Data_storage$ReplicateNum==j & Data_storage$Scale == "Local"]<-mean(ses.mpd(samp = com_data,dis = interspec_mat,null.model = "independentswap",abundance.weighted = T, runs = 999)$mpd.obs.z)

Data_storage$sesMNTDis_abund_z[Data_storage$Dispersal==DispV[i] & Data_storage$ReplicateNum==j & Data_storage$Scale == "Local"]<-mean(ses.mntd(samp = com_data,dis = interspec_mat,null.model = "independentswap",abundance.weighted = T, runs = 999)$mntd.obs.z)

Data_storage$sesMPDts_abund_z[Data_storage$Dispersal==DispV[i] & Data_storage$ReplicateNum==j & Data_storage$Scale == "Local"]<-mean(ses.mpd(samp = com_data,dis = interspec_mat,null.model = "trialswap",abundance.weighted = T, runs = 999)$mpd.obs.z)

Data_storage$sesMNTDts_abund_z[Data_storage$Dispersal==DispV[i] & Data_storage$ReplicateNum==j & Data_storage$Scale == "Local"]<-mean(ses.mntd(samp = com_data,dis = interspec_mat,null.model = "trialswap",abundance.weighted = T, runs = 999)$mntd.obs.z)

  #At the regional scale
  com_data<-matrix(colSums(t(SIH_data[["Abund",i]][400,,])),1,nspecies)
  colnames(com_data)<-1:nspecies
  
  Data_storage$MPD_abund[Data_storage$Dispersal==DispV[i] & Data_storage$ReplicateNum==j & Data_storage$Scale == "Regional"]<-mpd(com_data,cophenetic(SIH_data[["phylo",i]]),abundance.weighted = T)
  
  Data_storage$MNTD_abund[Data_storage$Dispersal==DispV[i] & Data_storage$ReplicateNum==j & Data_storage$Scale == "Regional"]<-mntd(com_data,cophenetic(SIH_data[["phylo",i]]),abundance.weighted = T)  
   
    #SES measures, taking the mean of the 3 seems more legit than just choosing the first? even though they all represent the exact same thing and were constructed from the same deterministic information
  reg_data<-matrix(rep(colSums(t(SIH_data[["Abund",i]][400,,])),each=3),3,nspecies)
  colnames(reg_data)<-1:nspecies
  
   Data_storage$sesMPDpp_abund_z[Data_storage$Dispersal==DispV[i] & Data_storage$ReplicateNum==j & Data_storage$Scale == "Regional"]<-mean(ses.mpd(samp = reg_data,dis = interspec_mat,null.model = "phylogeny.pool",abundance.weighted = T, runs = 999)$mpd.obs.z)

Data_storage$sesMNTDpp_abund_z[Data_storage$Dispersal==DispV[i] & Data_storage$ReplicateNum==j & Data_storage$Scale == "Regional"]<-mean(ses.mntd(samp = reg_data,dis = interspec_mat,null.model = "phylogeny.pool",abundance.weighted = T, runs = 999)$mntd.obs.z)

Data_storage$sesMPDtl_abund_z[Data_storage$Dispersal==DispV[i] & Data_storage$ReplicateNum==j & Data_storage$Scale == "Regional"]<-mean(ses.mpd(samp = reg_data,dis = interspec_mat,null.model = "taxa.labels",abundance.weighted = T, runs = 999)$mpd.obs.z)

Data_storage$sesMNTDtl_abund_z[Data_storage$Dispersal==DispV[i] & Data_storage$ReplicateNum==j & Data_storage$Scale == "Regional"]<-mean(ses.mntd(samp = reg_data,dis = interspec_mat,null.model = "taxa.labels",abundance.weighted = T, runs = 999)$mntd.obs.z)

Data_storage$sesMPDr_abund_z[Data_storage$Dispersal==DispV[i] & Data_storage$ReplicateNum==j & Data_storage$Scale == "Regional"]<-mean(ses.mpd(samp = reg_data,dis = interspec_mat,null.model = "richness",abundance.weighted = T, runs = 999)$mpd.obs.z)

Data_storage$sesMNTDr_abund_z[Data_storage$Dispersal==DispV[i] & Data_storage$ReplicateNum==j & Data_storage$Scale == "Regional"]<-mean(ses.mntd(samp = reg_data,dis = interspec_mat,null.model = "richness",abundance.weighted = T, runs = 999)$mntd.obs.z)

Data_storage$sesMPDf_abund_z[Data_storage$Dispersal==DispV[i] & Data_storage$ReplicateNum==j & Data_storage$Scale == "Regional"]<-mean(ses.mpd(samp = reg_data,dis = interspec_mat,null.model = "frequency",abundance.weighted = T, runs = 999)$mpd.obs.z)

Data_storage$sesMNTDf_abund_z[Data_storage$Dispersal==DispV[i] & Data_storage$ReplicateNum==j & Data_storage$Scale == "Regional"]<-mean(ses.mntd(samp = reg_data,dis = interspec_mat,null.model = "frequency",abundance.weighted = T, runs = 999)$mntd.obs.z)

Data_storage$sesMPDsp_abund_z[Data_storage$Dispersal==DispV[i] & Data_storage$ReplicateNum==j & Data_storage$Scale == "Regional"]<-mean(ses.mpd(samp = reg_data,dis = interspec_mat,null.model = "sample.pool",abundance.weighted = T, runs = 999)$mpd.obs.z)

Data_storage$sesMNTDsp_abund_z[Data_storage$Dispersal==DispV[i] & Data_storage$ReplicateNum==j & Data_storage$Scale == "Regional"]<-mean(ses.mntd(samp = reg_data,dis = interspec_mat,null.model = "sample.pool",abundance.weighted = T, runs = 999)$mntd.obs.z)

Data_storage$sesMPDis_abund_z[Data_storage$Dispersal==DispV[i] & Data_storage$ReplicateNum==j & Data_storage$Scale == "Regional"]<-mean(ses.mpd(samp = reg_data,dis = interspec_mat,null.model = "independentswap",abundance.weighted = T, runs = 999)$mpd.obs.z)

Data_storage$sesMNTDis_abund_z[Data_storage$Dispersal==DispV[i] & Data_storage$ReplicateNum==j & Data_storage$Scale == "Regional"]<-mean(ses.mntd(samp = reg_data,dis = interspec_mat,null.model = "independentswap",abundance.weighted = T, runs = 999)$mntd.obs.z)

Data_storage$sesMPDts_abund_z[Data_storage$Dispersal==DispV[i] & Data_storage$ReplicateNum==j & Data_storage$Scale == "Regional"]<-mean(ses.mpd(samp = reg_data,dis = interspec_mat,null.model = "trialswap",abundance.weighted = T, runs = 999)$mpd.obs.z)

Data_storage$sesMNTDts_abund_z[Data_storage$Dispersal==DispV[i] & Data_storage$ReplicateNum==j & Data_storage$Scale == "Regional"]<-mean(ses.mntd(samp = reg_data,dis = interspec_mat,null.model = "trialswap",abundance.weighted = T, runs = 999)$mntd.obs.z)

 

  }
  }
  


  

