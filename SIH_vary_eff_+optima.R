SIH_function<-function(dispersal=0.001,species=9,patches=30,eff_vary=F,eff_values=NA){ 
	#can define things here and then don't need to define them below 
  require(cluster)
  require(picante)
  
  #Constants####
  N<- matrix(10,ncol=species,nrow=patches) # Community x Species abundance matrix
  R<-rep(10*(species/10),patches) #Initial resources
  
  rInput<-150 #resource input
  rLoss<-10 #resource loss 
  mort<-0.2 #mortality
  Ext<- 0.1 #extinction Threshold
  eff<-0.2 #conversion efficiency
  
  if(eff_vary==T){eff<-eff_values} else{
   eff<-rep(eff,species) 
  }
  #always a vector, so that it doesn't mess up when eff_vary == F
  
  ePeriod<-40000 #period of env sinusoidal fluctuations
  eAMP<-1 #amplitude of envrionment sinusoidal fluctuations
  
  Tmax<-140000 #number of time steps in Simulation
  DT<- 0.08 # % size of discrete "time steps" - this is the Euler value
  
  #vectors####
  eOptimum<-1-seq(0,eAMP, by=eAMP/(species-1)) #species environmental optima
  #eOptimum<-runif(nspecies,min=0,max=1) #AG modif
  
  phylo<-as.phylo(hclust(daisy(cbind(eOptimum,eff_vary)),method ="ward")) #he had this set to ward.D
  plot(phylo,show.tip.label = F)
  tiplabels(pch=22,bg=heat.colors(species)[1:species]) #heat.colors = white -> red
  
  #dispersal conditions####
  dispersal_matrix<-matrix(1/(patches-1),patches,patches)
  diag(dispersal_matrix)<-0
  
  calc.immigration <- function(N,a,dispersal_matrix) dispersal_matrix%*%N*rep(a,each=patches)
  
  Prod<-matrix(NA,species*patches,40000)
  Abund<-Prod
  
  for(TS in 1:Tmax){
    Immigrants<-calc.immigration(N,dispersal,dispersal_matrix)
    envt.v<-0.5*eAMP*(sin((2*pi/ePeriod)*TS+1+(1:patches)*2*pi/patches)+1)
    consume <- 0.1*(1.5-abs(sapply(eOptimum,'-',envt.v)))
    Nt <- N*(1+DT*(rep(eff,each=patches)*R*consume - dispersal - mort)) + DT*Immigrants #abundance step
    Rt <- DT*rInput+R*(1-DT*(rLoss + rowSums(consume*N))) #resource step    
    N <- Nt * (Nt>Ext) # set to 0 if below extinction threshold
    R <- Rt
    
    if(TS>=100000){
      Prod[,(TS-100000)] <- c(t(eff*consume*R*N))
      Abund[,(TS-100000)] <- c(t(N))
    }
  } 
  
  Prod<-array(t(Prod),dim=c(40000,species,patches))
  Prod<-Prod[seq(1,40000,100),,]
  
  Abund<-array(t(Abund),dim=c(40000,species,patches))
  Abund<-Abund[seq(1,40000,100),,] #this only takes a subset, only every 100 time steps 
  return(list(Prod=Prod,Abund=Abund,phylo=phylo)) #will only return these things, nothing else
}


