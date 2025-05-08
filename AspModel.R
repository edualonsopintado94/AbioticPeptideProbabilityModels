##Asp Model

##A quantitative investigation of the chemical space surrounding 
##amino acid alphabet formation
##Yi Lu, Stephen J. Freeland

##Constants
h<-0.87045479 #Frequency of previous 
i<-0.879428551 #Aspartic frequency + previous 
lambda<-5 #Average peptide length for the Poisson 
set.seed(1994)

##Variables
G<-c() #AA vector in each peptide
l<-0 #Initial length
Smax=150 #number of simulations
Ns=0 #Simulated peptide counter
Nse<-c() #Vector of successful simulations
NS=0 #Simulation counter
Lpep<-c() #length vector of each peptide
Asp<-0 #Asp counter in each peptide
AspC<-c() #vector of Asp number per peptide

##Simulation
repeat {
NS<-NS+1 #Simulation number is updated
Ns <- 0 #Simulated peptide number is reset to 0
  repeat {
  Ns<-Ns+1 #Simulated peptide number is updated
  lmax=1+rpois(1,lambda) #the peptide length is simulated
  Lpep[Ns]=lmax #length of each peptide is added to the length vector
  repeat {
    l=l+1 #if l<lmax an AA is added
    r<-runif(1,0,1) #random number between 0 and 1
    if (r>h & r<=i) AA="Asp" else AA="Other"  
    #determination of the amino acid based on the result of r
    if (r>h & r<=i) Asp<-Asp+1 
    #annotation of Asp number per peptide
    G[l]=AA #the new AA is added to the growing peptide
    if (l==lmax) break #when lmax is reached the loop stops
  }
  AspC[Ns]<-Asp #the amount of Asp per peptide is saved
  Asp<-0
  #Asp counter is set to zero for the next peptide
  if (AspC[Ns]/Lpep[Ns]==1 & Lpep[Ns]==2) print (G) #print successful peptide for visual checking
  if (AspC[Ns]/Lpep[Ns]==1 & Lpep[Ns]==2) print (NS) #print nº of simulation to check the progress
  G<-c() #the AAs vector is emptied into each peptide
  l<-0 #length is restored
  if (AspC[Ns]/Lpep[Ns]==1 & Lpep[Ns]==2) break 
  #If a Asp(6) is generated, the loop stops
  }
Nse[NS]<-Ns #The number of successful simulations is added to its vector
if(NS==Smax) break 
#When the stipulated number of successful simulations is reached, the loop stops
}

library(modeest)
mlv(1/Nse,method = "meanshift") 
#Mode of the probability of the appearance of Asp(2)
par(mar = c(5, 6, 4, 2)) 
hist(1/Nse,col="green", breaks=3000, xlim=c(0,mean(1/Nse)), main = "Asp(2)", 
     xlab = "Probability of the sequence\n(between 0 and 1)", 
     ylab = "Frequency of Observed Probabilities\n(n = 150)")
abline(v=mlv(1/Nse, method="meanshift"), col="black", lty=2)
#Histogram of the probability of Asp(2) appearing
