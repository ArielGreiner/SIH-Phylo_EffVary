nspecies<-7 #the number of species
npatches<-10 #the number of patches

DispV<-c(0.0001,0.0005,0.001,0.005,0.01,0.05,0.1,0.5,1) #the dispersal rates 

#runs the SIH model at all dispersal rates in DispV and saves the abundances and productivity in a list
eff_values<-rnorm(nspecies,mean=0.2,sd=0.005)
SIH_data<-sapply(DispV,SIH_function,species=nspecies,patches=npatches,eff_vary=T,eff_values=eff_values)
#run until here first, this will take a few minutes

Data_storage<-data.frame(SR=NA,Biomass=NA,Biomass_CV=NA,PD=NA,MPD_abund=NA,MPD_pa=NA,MNTD_abund=NA,Dispersal=DispV,Scale=rep(c("Local","Regional"),each=length(DispV))) #building the data frame
#MPD_abund = mean pairwise distance, abundance weighted (vs pa = presence absence)
#need to put the NAs in initially
#MNTD = mean nearest taxon index, measures tippiness of tree
#need to add 'ses' to some? of the phylogenetic functions to make them compare the observed with a null
	#right now none of them are doing that 
for(i in 1:length(DispV)){
	#calculate species richness at the local and regional scale	
  Data_storage$SR[Data_storage$Dispersal==DispV[i] & Data_storage$Scale == "Local"] <-mean(rowMeans(apply(SIH_data[["Abund",i]]>0,3,rowSums))) 
  Data_storage$SR[Data_storage$Dispersal==DispV[i] & Data_storage$Scale == "Regional"]<-mean(rowSums(apply(SIH_data[["Abund",i]],2,rowMeans)>0)) #regional SR of all species in each time step
  #At the local scale...	
  com_data<-t(SIH_data[["Abund",i]][400,,])
  colnames(com_data)<-1:nspecies
  Data_storage$PD[Data_storage$Dispersal==DispV[i] & Data_storage$Scale == "Local"]<-mean(pd(com_data,SIH_data[["phylo",i]])$PD) 
  Data_storage$MPD_abund[Data_storage$Dispersal==DispV[i] & Data_storage$Scale == "Local"]<-mean(mpd(com_data,cophenetic(SIH_data[["phylo",i]]),abundance.weighted = T)) 
  Data_storage$MPD_pa[Data_storage$Dispersal==DispV[i] & Data_storage$Scale == "Local"]<-mean(mpd(com_data,cophenetic(SIH_data[["phylo",i]]),abundance.weighted = F))
  Data_storage$MNTD_abund[Data_storage$Dispersal==DispV[i] & Data_storage$Scale == "Local"]<-mean(mntd(com_data,cophenetic(SIH_data[["phylo",i]]),abundance.weighted = T))
  #At the regional scale
  com_data<-matrix(colSums(t(SIH_data[["Abund",i]][400,,])),1,nspecies)
  colnames(com_data)<-1:nspecies
  Data_storage$PD[Data_storage$Dispersal==DispV[i] & Data_storage$Scale == "Regional"]<-pd(com_data,SIH_data[["phylo",i]])$PD
  Data_storage$MPD_abund[Data_storage$Dispersal==DispV[i] & Data_storage$Scale == "Regional"]<-mpd(com_data,cophenetic(SIH_data[["phylo",i]]),abundance.weighted = T)
  Data_storage$MPD_pa[Data_storage$Dispersal==DispV[i] & Data_storage$Scale == "Regional"]<-mpd(com_data,cophenetic(SIH_data[["phylo",i]]),abundance.weighted = F)
  Data_storage$MNTD_abund[Data_storage$Dispersal==DispV[i] & Data_storage$Scale == "Regional"]<-mntd(com_data,cophenetic(SIH_data[["phylo",i]]),abundance.weighted = T)
  }

#Plot species richness at different dispersal levels 
require(ggplot2) #need to define x and y within aes in ggplot
ggplot(Data_storage,aes(x=Dispersal,y=SR,color=Scale,group=Scale))+
  geom_line(size=2)+ #plots data as lines
  scale_x_log10(breaks=DispV)+ #sets x axis to log10 scale
  theme_bw(base_size = 18)+ #gets rid of grey background
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) #removes grid lines

#Plot phylogenetic diversity at different dispersal levels 
ggplot(Data_storage,aes(x=Dispersal,y=PD,color=Scale,group=Scale))+
  geom_line(size=2)+ #plots data as lines
  scale_x_log10(breaks=DispV)+ #sets x axis to log10 scale
  theme_bw(base_size = 18)+ #gets rid of grey background
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) #removes grid lines

#MPD_pa may not be robust to reps, need to check
#Plot mean pairwise distance at different dispersal levels  
ggplot(Data_storage,aes(x=Dispersal,y=MPD_pa,color=Scale,group=Scale))+
  geom_line(size=2,na.rm=T)+ #plots data as lines
  scale_x_log10(breaks=DispV)+ #sets x axis to log10 scale
  facet_grid(.~Scale)+ #plots local and regional side-by-side 
  theme_bw(base_size = 18)+ #gets rid of grey background
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

#Plot abundance-weighted mean pairwise distance at different dispersal levels 
ggplot(Data_storage,aes(x=Dispersal,y=MPD_abund,color=Scale,group=Scale))+
  geom_line(size=2,na.rm=T)+ #plots data as lines
  scale_x_log10(breaks=DispV)+ #sets x axis to log10 scale
  theme_bw(base_size = 18)+ #gets rid of grey background
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) #removes grid lines

#Plot abundance-weighted mean nearest taxon index at different dispersal levels 
ggplot(Data_storage,aes(x=Dispersal,y=MNTD_abund,color=Scale,group=Scale))+
  geom_line(size=2,na.rm=T)+ #plots data as lines
  facet_grid(.~Scale)+
  scale_x_log10(breaks=DispV)+ #sets x axis to log10 scale
  theme_bw(base_size = 18)+ #gets rid of grey background
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) #removes grid lines

#Plot presence-absence mean pairwise distance at different species richnesses 
ggplot(Data_storage,aes(x=SR,y=MPD_pa,color=Scale,group=Scale))+
  geom_line(size=2,na.rm=T)+ #plots data as lines
  #scale_x_log10(breaks=DispV)+ #sets x axis to log10 scale
  facet_grid(.~Scale)+
  theme_bw(base_size = 18)+ #gets rid of grey background
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

#plot phylo####
par(mfrow=c(3,3)) #local phylogenetic tree from whatever '3' is ##check
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
