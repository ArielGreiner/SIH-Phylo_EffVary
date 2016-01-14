nspecies<-7 #the number of species
npatches<-10 #the number of patches
nreplicates<-2 #number of replicates
nfunctions<-1

DispV<-c(0.0001,0.0005,0.001,0.005,0.01,0.05,0.1,0.5,1) #the dispersal rates 

Data_storage<-data.frame(SR=NA,Biomass=NA,Biomass_CV=NA,PD=NA,MPD_abund=NA,MPD_pa=NA,MNTD_abund=NA,beta_MPDabund =NA, beta_MNTDabund=NA, sesMPD_abund_z = NA, sesMNTD_abund_z = NA, sesMPD_abund_p = NA, sesMNTD_abund_p = NA,Dispersal=rep(DispV,each=nreplicates),ReplicateNum=factor(1:nreplicates),Scale=rep(c("Local","Regional"),each=length(DispV)*nreplicates)) #building the data frame

MTraits<-t(matrix(1,nspecies,nfunctions))

for(j in 1:nreplicates){
#runs the SIH model at all dispersal rates in DispV and saves the abundances and productivity in a list
set.seed(j)
eff_values<-rnorm(nspecies,mean=0.2,sd=0.005)
SIH_data<-sapply(DispV,SIH_function,species=nspecies,patches=npatches,eff_vary=T,eff_values=eff_values)


#MPD_abund = mean pairwise distance, abundance weighted (vs pa = presence absence)
#need to put the NAs in initially
#MNTD = mean nearest taxon index, measures tippiness of tree
#need to add 'ses' to some? of the phylogenetic functions to make them compare the observed with a null
	#right now none of them are doing that 
for(i in 1:length(DispV)){
 #BIOMASS
  LFunc_rate1<-apply(SIH_data[["Abund",i]],1,function(M){MTraits%*%M}) #calculates the amount of each function in each patch based on the abundances from the SIH model
  LFunc_rate<-array(t(LFunc_rate1),dim=c(nrow(SIH_data[["Abund",1]]),nfunctions,npatches)) #arranges the functions into an array that is easier to use
  Data_storage$Biomass[Data_storage$Dispersal==DispV[i] & Data_storage$ReplicateNum==j & Data_storage$Scale == "Local"]<-mean(apply(LFunc_rate,3,colMeans))
  Data_storage$Biomass[Data_storage$Dispersal==DispV[i] & Data_storage$ReplicateNum==j & Data_storage$Scale == "Regional"]<-colMeans(apply(LFunc_rate,2,rowSums)) #calculates and saves the mean regoinal rate for each function
  
	#calculate species richness at the local and regional scale	
  Data_storage$SR[Data_storage$Dispersal==DispV[i] & Data_storage$ReplicateNum==j & Data_storage$Scale == "Local"] <-mean(rowMeans(apply(SIH_data[["Abund",i]]>0,3,rowSums))) 
  Data_storage$SR[Data_storage$Dispersal==DispV[i] & Data_storage$ReplicateNum==j & Data_storage$Scale == "Regional"]<-mean(rowSums(apply(SIH_data[["Abund",i]],2,rowMeans)>0)) #regional SR of all species in each time step
  #At the local scale...	
  com_data<-t(SIH_data[["Abund",i]][400,,])
  colnames(com_data)<-1:nspecies
  Data_storage$PD[Data_storage$Dispersal==DispV[i] & Data_storage$ReplicateNum==j & Data_storage$Scale == "Local"]<-mean(pd(com_data,SIH_data[["phylo",i]])$PD) 
  Data_storage$MPD_abund[Data_storage$Dispersal==DispV[i] & Data_storage$ReplicateNum==j & Data_storage$Scale == "Local"]<-mean(mpd(com_data,cophenetic(SIH_data[["phylo",i]]),abundance.weighted = T)) 
  Data_storage$MPD_pa[Data_storage$Dispersal==DispV[i] & Data_storage$ReplicateNum==j & Data_storage$Scale == "Local"]<-mean(mpd(com_data,cophenetic(SIH_data[["phylo",i]]),abundance.weighted = F))
  Data_storage$MNTD_abund[Data_storage$Dispersal==DispV[i] & Data_storage$ReplicateNum==j & Data_storage$Scale == "Local"]<-mean(mntd(com_data,cophenetic(SIH_data[["phylo",i]]),abundance.weighted = T))
  
  #ses measures
  interspec_mat <- cophenetic(SIH_data[["phylo",i]])
     colnames(interspec_mat) <- 1:nspecies
     rownames(interspec_mat) <- 1:nspecies
Data_storage$sesMPD_abund_z[Data_storage$Dispersal==DispV[i] & Data_storage$ReplicateNum==j & Data_storage$Scale == "Local"]<-mean(ses.mpd(samp = com_data,dis = interspec_mat,null.model = "phylogeny.pool",abundance.weighted = T, runs = 999)$mpd.obs.z)

Data_storage$sesMPD_abund_p[Data_storage$Dispersal==DispV[i] & Data_storage$ReplicateNum==j & Data_storage$Scale == "Local"]<-mean(ses.mpd(samp = com_data,dis = interspec_mat,null.model = "phylogeny.pool",abundance.weighted = T, runs = 999)$mpd.obs.p)

Data_storage$sesMNTD_abund_z[Data_storage$Dispersal==DispV[i] & Data_storage$ReplicateNum==j & Data_storage$Scale == "Local"]<-mean(ses.mntd(samp = com_data,dis = interspec_mat,null.model = "phylogeny.pool",abundance.weighted = T, runs = 999)$mntd.obs.z)

Data_storage$sesMNTD_abund_p[Data_storage$Dispersal==DispV[i] & Data_storage$ReplicateNum==j & Data_storage$Scale == "Local"]<-mean(ses.mntd(samp = com_data,dis = interspec_mat,null.model = "phylogeny.pool",abundance.weighted = T, runs = 999)$mntd.obs.p)

  #beta measures
  #colnames(com_data) <- paste('species', 1:nspecies)
     #colnames(interspec_mat) <- paste('species', 1:nspecies)
     #rownames(interspec_mat) <- paste('species', 1:nspecies)
     #interspec_mat <- cophenetic(SIH_data[["phylo",i]])
     #colnames(interspec_mat) <- 1:nspecies
     #rownames(interspec_mat) <- 1:nspecies
  Data_storage$beta_MPDabund[Data_storage$Dispersal==DispV[i] & Data_storage$ReplicateNum==j & Data_storage$Scale == "Regional"] <- mean(comdist(com_data,interspec_mat,abundance.weighted=TRUE))
Data_storage$beta_MNTDabund[Data_storage$Dispersal==DispV[i] & Data_storage$ReplicateNum==j & Data_storage$Scale == "Regional"] <- mean(comdistnt(com_data,interspec_mat,abundance.weighted=TRUE,exclude.conspecifics=TRUE))
#^gives weird warning messages...

  #At the regional scale
  com_data<-matrix(colSums(t(SIH_data[["Abund",i]][400,,])),1,nspecies)
  colnames(com_data)<-1:nspecies
  Data_storage$PD[Data_storage$Dispersal==DispV[i] & Data_storage$ReplicateNum==j & Data_storage$Scale == "Regional"]<-pd(com_data,SIH_data[["phylo",i]])$PD
  Data_storage$MPD_abund[Data_storage$Dispersal==DispV[i] & Data_storage$ReplicateNum==j & Data_storage$Scale == "Regional"]<-mpd(com_data,cophenetic(SIH_data[["phylo",i]]),abundance.weighted = T)
  Data_storage$MPD_pa[Data_storage$Dispersal==DispV[i] & Data_storage$ReplicateNum==j & Data_storage$Scale == "Regional"]<-mpd(com_data,cophenetic(SIH_data[["phylo",i]]),abundance.weighted = F)
  Data_storage$MNTD_abund[Data_storage$Dispersal==DispV[i] & Data_storage$ReplicateNum==j & Data_storage$Scale == "Regional"]<-mntd(com_data,cophenetic(SIH_data[["phylo",i]]),abundance.weighted = T)  
  
  #SES measures, taking the mean of the 3 seems more legit than just choosing the first? even though they all represent the exact same thing and were constructed from the same deterministic information
  reg_data<-matrix(rep(colSums(t(SIH_data[["Abund",i]][400,,])),each=3),3,nspecies)
  colnames(reg_data)<-1:nspecies
  Data_storage$sesMPD_abund_z[Data_storage$Dispersal==DispV[i] & Data_storage$ReplicateNum==j & Data_storage$Scale == "Regional"]<-mean(ses.mpd(samp = reg_data,dis = cophenetic(SIH_data[["phylo",i]]),null.model = "phylogeny.pool",abundance.weighted = T, runs = 999)$mpd.obs.z)
  
    Data_storage$sesMPD_abund_p[Data_storage$Dispersal==DispV[i] & Data_storage$ReplicateNum==j & Data_storage$Scale == "Regional"]<-mean(ses.mpd(samp = reg_data,dis = cophenetic(SIH_data[["phylo",i]]),null.model = "phylogeny.pool",abundance.weighted = T, runs = 999)$mpd.obs.p)

Data_storage$sesMNTD_abund_z[Data_storage$Dispersal==DispV[i] & Data_storage$ReplicateNum==j & Data_storage$Scale == "Regional"]<-mean(ses.mntd(samp = reg_data,dis = cophenetic(SIH_data[["phylo",i]]),null.model = "phylogeny.pool",abundance.weighted = T, runs = 999)$mntd.obs.z)

Data_storage$sesMNTD_abund_p[Data_storage$Dispersal==DispV[i] & Data_storage$ReplicateNum==j & Data_storage$Scale == "Regional"]<-mean(ses.mntd(samp = reg_data,dis = cophenetic(SIH_data[["phylo",i]]),null.model = "phylogeny.pool",abundance.weighted = T, runs = 999)$mntd.obs.p)

  }
  }
require(dplyr) 
Data_storage_total<-summarise(group_by(Data_storage, Dispersal, Scale), Mean_SR=mean(SR,na.rm=T), SD_SR=sd(SR,na.rm=T), Mean_Biomass=mean(Biomass,na.rm=T), SD_Biomass=sd(Biomass,na.rm=T), Mean_PD=mean(PD,na.rm=T), SD_PD=sd(PD,na.rm=T), Mean_MPD_abund=mean(MPD_abund,na.rm=T), SD_MPD_abund=sd(MPD_abund,na.rm=T), Mean_MPD_pa=mean(MPD_pa,na.rm=T), SD_MPD_pa=sd(MPD_pa,na.rm=T), Mean_MNTD_abund=mean(MNTD_abund,na.rm=T), SD_MNTD_abund=sd(MNTD_abund,na.rm=T),Mean_beta_MPDabund=mean(beta_MPDabund,na.rm=T), SD_beta_MPDabund=sd(beta_MPDabund,na.rm=T), 
Mean_beta_MNTDabund=mean(beta_MNTDabund,na.rm=T), SD_beta_MNTDabund=sd(beta_MNTDabund,na.rm=T), Mean_sesMPD_abund_z = mean(sesMPD_abund_z,na.rm=T), SD_sesMPD_abund_z = sd(sesMPD_abund_z,na.rm=T),
Mean_sesMPD_abund_p = mean(sesMPD_abund_p,na.rm=T), SD_sesMPD_abund_p = sd(sesMPD_abund_p,na.rm=T), Mean_sesMNTD_abund_z = mean(sesMNTD_abund_z,na.rm=T), SD_sesMNTD_abund_z = sd(sesMNTD_abund_z,na.rm=T), Mean_sesMNTD_abund_p = mean(sesMNTD_abund_p,na.rm=T), SD_sesMNTD_abund_p = sd(sesMNTD_abund_p,na.rm=T))

###need to go to an older version of this to get the geom_errobar variant
#Plot species richness at different dispersal levels 
require(ggplot2) #need to define x and y within aes in ggplot
ggplot(Data_storage_total,aes(x=Dispersal,y=Mean_SR,color=Scale,group=Scale,fill=Scale,alpha=0.1))+
  geom_line(size=2)+ #plots data as lines
  geom_ribbon(aes(ymin=Mean_SR-SD_SR,ymax=Mean_SR+SD_SR),width=0.1)+
  scale_x_log10(breaks=DispV)+ #sets x axis to log10 scale
  theme_bw(base_size = 18)+ #gets rid of grey background
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) #removes grid lines
  
#Plot phylogenetic diversity at different dispersal levels 
ggplot(Data_storage_total,aes(x=Dispersal,y=Mean_PD,color=Scale,group=Scale,fill=Scale,alpha=0.1))+
  geom_line(size=2)+ #plots data as lines
  geom_ribbon(aes(ymin=Mean_PD-SD_PD,ymax=Mean_PD+SD_PD),width=0.1)+
  scale_x_log10(breaks=DispV)+ #sets x axis to log10 scale
  theme_bw(base_size = 18)+ #gets rid of grey background
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) #removes grid lines

#MPD_pa may not be robust to reps, need to check
#Plot mean pairwise distance at different dispersal levels  
ggplot(Data_storage_total,aes(x=Dispersal,y=Mean_MPD_pa,color=Scale,group=Scale,fill=Scale,alpha=0.1))+
  geom_line(size=2,na.rm=T)+ #plots data as lines
  geom_ribbon(aes(ymin=Mean_MPD_pa-SD_MPD_pa,ymax=Mean_MPD_pa+SD_MPD_pa),width=0.1)+
  scale_x_log10(breaks=DispV)+ #sets x axis to log10 scale
  facet_grid(.~Scale)+ #plots local and regional side-by-side 
  theme_bw(base_size = 18)+ #gets rid of grey background
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

#Plot abundance-weighted mean pairwise distance at different dispersal levels 
ggplot(Data_storage_total,aes(x=Dispersal,y=Mean_MPD_abund,color=Scale,group=Scale,fill=Scale,alpha=0.1))+
  geom_line(size=2,na.rm=T)+ #plots data as lines
  geom_ribbon(aes(ymin=Mean_MPD_abund-SD_MPD_abund,ymax=Mean_MPD_abund+SD_MPD_abund),width=0.1)+
  scale_x_log10(breaks=DispV)+ #sets x axis to log10 scale
  theme_bw(base_size = 18)+ #gets rid of grey background
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) #removes grid lines

#Plot abundance-weighted mean nearest taxon index at different dispersal levels 
ggplot(Data_storage_total,aes(x=Dispersal,y=Mean_MNTD_abund,color=Scale,group=Scale,fill=Scale,alpha=0.1))+
  geom_line(size=2,na.rm=T)+ #plots data as lines
  geom_ribbon(aes(ymin=Mean_MNTD_abund-SD_MNTD_abund,ymax=Mean_MNTD_abund+SD_MNTD_abund),width=0.1)+
  facet_grid(.~Scale)+
  scale_x_log10(breaks=DispV)+ #sets x axis to log10 scale
  theme_bw(base_size = 18)+ #gets rid of grey background
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) #removes grid lines

#Plot presence-absence mean pairwise distance at different species richnesses 
ggplot(Data_storage_total,aes(x=Mean_SR,y=Mean_MPD_pa,color=Scale,group=Scale,fill=Scale,alpha=0.1))+
  geom_line(size=2,na.rm=T)+ #plots data as lines
  geom_ribbon(aes(ymin=Mean_MPD_pa-SD_MPD_pa,ymax=Mean_MPD_pa+SD_MPD_pa),width=0.1)+
  #scale_x_log10(breaks=DispV)+ #sets x axis to log10 scale
  facet_grid(.~Scale)+
  theme_bw(base_size = 18)+ #gets rid of grey background
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  
#Plot Biomass
ggplot(Data_storage_total,aes(x=Dispersal,y=Mean_Biomass))+ #defines the variables that you are plotting
  geom_line(size=2)+ #plots data as lines
  facet_grid(Scale~.,scale="free")+ #separates out local and regional into two panels
  geom_errorbar(aes(ymin=Mean_Biomass-SD_Biomass,ymax=Mean_Biomass+SD_Biomass),width=0.1)+
  #facet_grid(ReplicateNume~.,scale='free')
  #geom_errorbar(data = Function_rates_sd)+
  scale_x_log10(breaks=DispV)+ #sets x axis to log10 scale
  theme_bw(base_size = 18)+ #gets rid of grey background
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) #removes grid lines
  
  #Could you produce some plots of PD vs SR (x axis). I curious about how far from linear this relationship deviates. 
  ggplot(Data_storage_total,aes(x=Mean_SR,y=Mean_PD,color=Scale,group=Scale))+
  geom_line(size=2)+ #plots data as lines
  geom_errorbar(aes(ymin=Mean_PD-SD_PD,ymax=Mean_PD+SD_PD),width=0.1)+
  scale_x_log10(breaks=DispV)+ #sets x axis to log10 scale
  theme_bw(base_size = 18)+ #gets rid of grey background
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) #removes grid lines

#example of how to plot data colour-coded and connected by dispersal
ggplot(Data_storage_total,aes(x=Mean_SR,y=Mean_Biomass,color=factor(Dispersal),group=Scale))+
geom_point()+ #plots data as points
geom_path()+
facet_grid(Scale~.,scale="free")+
theme_bw(base_size = 18)+ #gets rid of grey background
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

#plots data as points with error bars, colour-coded by scale and not by dispersal
ggplot(Data_storage_total,aes(x=Mean_SR,y=Mean_PD,color=Scale,group=Scale))+
geom_point()+ #plots data as points
geom_errorbar(aes(ymin=Mean_PD-SD_PD,ymax=Mean_PD+SD_PD),width=0.1)+
facet_grid(Scale~.,scale="free")+
theme_bw(base_size = 18)+ #gets rid of grey background
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  
ggplot(Data_storage_total,aes(x=Mean_SR,y=Mean_PD,color=factor(Dispersal),group=Scale))+
geom_point()+ #plots data as points
geom_path()+
geom_errorbar(aes(ymin=Mean_PD-SD_PD,ymax=Mean_PD+SD_PD),width=0.1)+
facet_grid(Scale~.,scale="free")+
theme_bw(base_size = 18)+ #gets rid of grey background
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 
 

##THESE SHOULD PROBABLY BE MOVED WITHIN THE FOR-LOOP
#plot phylo####
par(mfrow=c(3,3)) #local phylogenetic tree from whatever '3' is (see below) ##check
for(i in 1:length(DispV)){
com_data<-t(SIH_data[["Abund",i]][400,,])
plot(SIH_data[["phylo",i]],show.tip.label = F,main=paste("Dispersal = ",DispV[i]))
tiplabels(pch=22,bg=heat.colors(nspecies)[1:nspecies], cex=decostand(com_data[3,],method = 'range')*5)
}

par(mfrow=c(3,3)) #regional phylogenetic tree ##check
for(i in 1:length(DispV)){
  com_data<-matrix(colSums(t(SIH_data[["Abund",i]][400,,])),1,nspecies)
  plot(SIH_data[["phylo",i]],show.tip.label = F,main=paste("Dispersal = ",DispV[i]))
  tiplabels(pch=22,bg=heat.colors(nspecies)[1:nspecies], cex=com_data/max(com_data)*5)
}

for(i in 1:length(DispV)){
  com_data<-t(SIH_data[["Abund",i]][400,,])
  plot(SIH_data[["phylo",i]],show.tip.label = F,main=paste("Dispersal = ",DispV[i]))
  tiplabels(pch=22,bg=heat.colors(nspecies)[1:nspecies], cex=3*com_data[3,]>0)
}

###Plot 3x3 grid of possible species phylogenies
H<-seq(0,1,length=9)
par(mfrow=c(3,3)) #regional phylogenetic tree ##check
for(i in 1:9){
  e<-rnorm(nspecies,mean = 0.2,sd=0.005)
  phylo_e<-plot(hclust(daisy(cbind(H,e)),method ="single"))
    #mpd(com_data,cophenetic(phylo_e,abundance.weighted = F))
}
