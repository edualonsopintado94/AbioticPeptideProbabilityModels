##Gly Model

##A quantitative investigation of the chemical space surrounding 
##amino acid alphabet formation
##Yi Lu, Stephen J. Freeland

##Constants
h<-0.209986001 #Frequency of previous 
i<-0.385871711 #Glycine frequency + previous 
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
Gly<-0 #Gly counter in each peptide
GlyC<-c() #vector of Gly number per peptide

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
    if (r>h & r<=i) AA="Gly" else AA="Other"  
    #determination of the amino acid based on the result of r
    if (r>h & r<=i) Gly<-Gly+1 
    #annotation of Gly number per peptide
    G[l]=AA #the new AA is added to the growing peptide
    if (l==lmax) break #when lmax is reached the loop stops
  }
  GlyC[Ns]<-Gly #the amount of Gly per peptide is saved
  Gly<-0
  #Gly counter is set to zero for the next peptide
  if (GlyC[Ns]/Lpep[Ns]==1 & Lpep[Ns]==4) print (G) #print successful peptide for visual checking
  if (GlyC[Ns]/Lpep[Ns]==1 & Lpep[Ns]==4) print (NS) #print nº of simulation to check the progress
  G<-c() #the AAs vector is emptied into each peptide
  l<-0 #length is restored
  if (GlyC[Ns]/Lpep[Ns]==1 & Lpep[Ns]==4) break
  #If a Gly(4) is generated, the loop stops
  }
Nse[NS]<-Ns #The number of successful simulations is added to its vector
if(NS==Smax) break 
#When the stipulated number of successful simulations is reached, the loop stops
}

library(modeest)
mlv(1/Nse,method = "meanshift") 
#Mode of the probability of the appearance of Gly(4)
par(mar = c(5, 6, 4, 2)) 
hist(1/Nse,col="red", breaks=1000, xlim=c(0,mean(1/Nse)), main = "Gly(4)", 
     xlab = "Probability of the sequence\n(between 0 and 1)", 
     ylab = "Frequency of Observed Probabilities\n(n = 150)")
abline(v=mlv(1/Nse, method="meanshift"), col="black", lty=2)
#Histogram of the probability of Gly(4) appearing
