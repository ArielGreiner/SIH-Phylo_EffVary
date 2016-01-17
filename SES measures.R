 
sesMPD_abund_z = NA, sesMNTD_abund_z = NA, sesMPD_abund_p = NA, sesMNTD_abund_p = NA

#maybe add these two things within the for-loop
if(Data_storage$sesMPD_abund_p > 0.05){
	Data_storage$sesMPD_abund_p <- 0
}

if(Data_storage$sesMPD_abund_p <= 0.05){
	Data_storage$sesMPD_abund_p <- 1
}

#hopefully this plots sesMPD measures and colour codes the points by sesMPD_abund_p
ggplot(Data_storage_total,aes(x=Dispersal,y=Mean_sesMPD_abund_z,color=factor(Mean_sesMPD_abund_p),group=Scale))+
geom_point()+ #plots data as points
geom_line()+
facet_grid(Scale~.,scale="free")+
theme_bw(base_size = 18)+ #gets rid of grey background
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

#just checking things
ggplot(Data_storage_total,aes(x=Dispersal,y=Mean_sesMPD_abund_z,color=Scale,group=Scale,fill=Scale,alpha=0.1))+
  geom_line(size=2)+ #plots data as lines
  geom_ribbon(aes(ymin=Mean_sesMPD_abund_z-SD_sesMPD_abund_z,ymax=Mean_sesMPD_abund_z+SD_sesMPD_abund_z),width=0.1)+
  scale_x_log10(breaks=DispV)+ #sets x axis to log10 scale
  theme_bw(base_size = 18)+ #gets rid of grey background
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) #removes grid lines
  
 #plotting phylog evenness/clustering stuff
 ggplot(Data_storage_total,aes(x=Dispersal,y=Mean_sesMPD_abund_z,color=factor(sum_phylogeven),group=Scale))+
geom_point()+ #plots data as points
geom_line()+
facet_grid(Scale~.,scale="free")+
theme_bw(base_size = 18)+ #gets rid of grey background
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

#plotting phylog evenness/clustering stuff with error bars

ggplot(Data_storage_total,aes(x=Dispersal,y=Mean_sesMPD_abund_z,color=factor(sum_phylogeven),group=Scale))+
geom_point()+ #plots data as points
geom_errorbar(aes(ymin=Mean_sesMPD_abund_z-SD_sesMPD_abund_z,ymax=Mean_sesMPD_abund_z+SD_sesMPD_abund_z),width=0.1)+
geom_line()+
facet_grid(Scale~.,scale="free")+
theme_bw(base_size = 18)+ #gets rid of grey background
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
