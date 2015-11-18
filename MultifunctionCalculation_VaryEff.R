#source("SIH code_PLT.r") #the SIH model code
require(vegan)


nspecies<-9 #the number of species
npatches<-30 #the number of patches
nfunctions<-1 #the number of functions #CURRENTLY ONLY LOOKING AT PRODUCTION OF 1 FUNCTION
function_overlap<-0.5 #the amount of functional overlap
nreplicates<-5

DispV<-c(0.0001,0.0005,0.001,0.005,0.01,0.05,0.1,0.5,1) #the dispersal rates 

#runs the SIH model at all dispersal rates in DispV and saves the abundances and productivity in a list

Function_rates<-data.frame(Rate=NA,Function=factor(1:nfunctions),ReplicateNum=factor(1:nreplicates),Dispersal=rep(DispV,each=nfunctions*nreplicates),Scale=rep(c("Local","Regional"),each=nreplicates*nfunctions*length(DispV))) #storage dataframe for functional rates
Function_number<-data.frame(Number=NA,ReplicateNum=factor(1:nreplicates),Dispersal=rep(DispV,each=nreplicates),Scale=rep(c("Local","Regional"),each= nreplicates*length(DispV))) #storage dataframe for functional rates

for(j in 1:nreplicates){
eff_values<-rnorm(nspecies,mean=0.2,sd=0.005)
SIH_data<-sapply(DispV,SIH_function,species=nspecies,patches=npatches,eff_vary=T,eff_values=eff_values)	

#SIH_data<-sapply(DispV,SIH_function,species=nspecies,patches=npatches,eff_vary=T,eff_values=eff_values)
#SIH_data<-sapply(DispV,SIH_function,species=nspecies,patches=npatches)

#define traits#### - I have added this to make some but you can define your traits any way you like
MTraits<-t(matrix(runif(nspecies*nfunctions)*rbinom(nspecies*nfunctions,1,function_overlap),nspecies,nfunctions)) #no structure 
#MTraits[1,]<-c(1,1,1,0,0,0,0,0,0)
#MTraits[2,]<-c(0,0,0,1,1,1,0,0,0)
#MTraits[3,]<-c(0,0,0,0,0,0,1,1,1)
#print(MTraits) #just shows the traits

#Version that does not allow for replicates
#Function_rates<-data.frame(Rate=NA,Function=factor(1:nfunctions),Dispersal=rep(DispV,each=nfunctions),Scale=rep(c("Local","Regional"),Replicate Number=factor(1:nreplicates),each=nreplicates*nfunctions*length(DispV))) #storage dataframe for functional rates

#Version allowing for replicates (moved out of the loop)
#Function_rates<-data.frame(Rate=NA,Function=factor(1:nfunctions),ReplicateNum=factor(1:nreplicates),Dispersal=rep(DispV,each=nfunctions*nreplicates),Scale=rep(c("Local","Regional"),each=nreplicates*nfunctions*length(DispV))) #storage dataframe for functional rates
#Function_number<-data.frame(Number=NA,ReplicateNum=factor(1:nreplicates),Dispersal=rep(DispV,each=nreplicates),Scale=rep(c("Local","Regional"),each= nreplicates*length(DispV))) #storage dataframe for functional rates

for(i in 1:length(DispV)){ #loops through all dispersal rates
  LFunc_rate1<-apply(SIH_data[["Abund",i]],1,function(M){MTraits%*%M}) #calculates the amount of each function in each patch based on the abundances from the SIH model
  LFunc_rate<-array(t(LFunc_rate1),dim=c(nrow(SIH_data[["Abund",1]]),nfunctions,npatches)) #arranges the functions into an array that is easier to use
  Function_number[Function_number$Scale=="Local" & Function_number$Dispersal==DispV[i] & Function_number$ReplicateNum==j,"Number"]<-mean(apply(LFunc_rate>0,3,rowSums))
  Function_number[Function_number$Scale=="Regional" & Function_number$Dispersal==DispV[i] & Function_number$ReplicateNum==j,"Number"]<-mean(rowSums(apply(LFunc_rate,2,rowSums)>0))
  #Single function variant:
  Function_rates[Function_rates$Scale=="Local" & Function_rates$Dispersal==DispV[i] & Function_number$ReplicateNum==j,"Rate"]<-mean(apply(LFunc_rate,3,colMeans))
  #multifunction variant:
  #Function_rates[Function_rates$Scale=="Local" & Function_rates$Dispersal==DispV[i],"Rate"]<-rowMeans(apply(LFunc_rate,3,colMeans)) #calculates and saves the mean local rate for each function
  Function_rates[Function_rates$Scale=="Regional" & Function_rates$Dispersal==DispV[i] & Function_number$ReplicateNum==j,"Rate"]<-colMeans(apply(LFunc_rate,2,rowSums)) #calculates and saves the mean regoinal rate for each function
}
}

#Take the averages over all of the replicates
Function_rates_avg<-data.frame(Rate=NA,Function=factor(1:nfunctions),Dispersal=rep(DispV,each=nfunctions),Scale=rep(c("Local","Regional"),each=nfunctions*length(DispV))) #storage dataframe for functional rates
Function_number_avg<-data.frame(Number=NA,Dispersal=DispV,Scale=rep(c("Local","Regional"),each=length(DispV))) #storage dataframe for functional rates

for(i in 1:length(DispV)){ #loops through all dispersal rates
  #FuncNumbAvgLocal <- mean(Function_number$Number[Function_number$Scale=="Local" & Function_number$Dispersal==DispV[i]])
  Function_number_avg[Function_number_avg$Scale=="Local" & Function_number_avg$Dispersal==DispV[i],"Number"]<- mean(Function_number$Number[Function_number$Scale=="Local" & Function_number$Dispersal==DispV[i]])
  Function_number_avg[Function_number_avg$Scale=="Regional" & Function_number_avg$Dispersal==DispV[i],"Number"]<- mean(Function_number$Number[Function_number$Scale=="Regional" & Function_number$Dispersal==DispV[i]])
  #Single function variant:
  Function_rates_avg[Function_rates_avg$Scale=="Local" & Function_rates_avg$Dispersal==DispV[i],"Rate"]<- mean(Function_rates$Rate[Function_rates$Scale=="Local" & Function_rates$Dispersal==DispV[i]])
  #multifunction variant: (need to adjust still)
  #Function_rates[Function_rates$Scale=="Local" & Function_rates$Dispersal==DispV[i],"Rate"]<-rowMeans(apply(LFunc_rate,3,colMeans)) 
  Function_rates_avg[Function_rates_avg$Scale=="Regional" & Function_rates_avg$Dispersal==DispV[i],"Rate"]<- mean(Function_rates$Rate[Function_rates$Scale=="Regional" & Function_rates$Dispersal==DispV[i]]) 
}

#Take the standard deviation over all of the replicates
Function_rates_sd<-data.frame(Rate=NA,Function=factor(1:nfunctions),Dispersal=rep(DispV,each=nfunctions),Scale=rep(c("Local","Regional"),each=nfunctions*length(DispV))) #storage dataframe for functional rates
Function_number_sd<-data.frame(Number=NA,Dispersal=DispV,Scale=rep(c("Local","Regional"),each=length(DispV))) #storage dataframe for functional rates

for(i in 1:length(DispV)){ #loops through all dispersal rates
  Function_number_sd[Function_number_sd$Scale=="Local" & Function_number_sd$Dispersal==DispV[i],"Number"]<- sd(Function_number$Number[Function_number$Scale=="Local" & Function_number$Dispersal==DispV[i]])
  Function_number_sd[Function_number_sd$Scale=="Regional" & Function_number_sd$Dispersal==DispV[i],"Number"]<- sd(Function_number$Number[Function_number$Scale=="Regional" & Function_number$Dispersal==DispV[i]])
  #Single function variant:
  Function_rates_sd[Function_rates_sd$Scale=="Local" & Function_rates_sd$Dispersal==DispV[i],"Rate"]<- sd(Function_rates$Rate[Function_rates$Scale=="Local" & Function_rates$Dispersal==DispV[i]])
  #multifunction variant: (need to adjust still)
  #Function_rates[Function_rates$Scale=="Local" & Function_rates$Dispersal==DispV[i],"Rate"]<-rowMeans(apply(LFunc_rate,3,colMeans)) 
  Function_rates_sd[Function_rates_sd$Scale=="Regional" & Function_rates_sd$Dispersal==DispV[i],"Rate"]<- sd(Function_rates$Rate[Function_rates$Scale=="Regional" & Function_rates$Dispersal==DispV[i]]) 
}



require(ggplot2) #plotting package
ggplot(Function_rates_avg,aes(x=Dispersal,y=Rate,color=Function,group=Function,linetype=Function))+ #defines the variables that you are plotting
  geom_line(size=2)+ #plots data as lines
  facet_grid(Scale~.,scale="free")+ #separates out local and regional into two panels
  #facet_grid(ReplicateNume~.,scale='free')
  #geom_errorbar(data = Function_rates_sd)+
  scale_x_log10(breaks=DispV)+ #sets x axis to log10 scale
  theme_bw(base_size = 18)+ #gets rid of grey background
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) #removes grid lines

ggplot(Function_number_avg,aes(x=Dispersal,y=Number,group=Scale,color=Scale,linetype=Scale))+ #defines the variables that you are plotting
  geom_line(size=2)+ #plots data as lines
  scale_x_log10(breaks=DispV)+ #sets x axis to log10 scale
  theme_bw(base_size = 18)+ #gets rid of grey background
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) #removes grid lines
  