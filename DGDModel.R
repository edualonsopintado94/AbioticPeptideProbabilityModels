##DGD Model

##A quantitative investigation of the chemical space surrounding 
##amino acid alphabet formation
##Yi Lu, Stephen J. Freeland

##Constants
h<-0.209986001 #Frequency of previous 
i<-0.385871711 #Glycine frequency + previous
j<-0.87045479 #Frequency of previous 
k<-0.879428551 #Aspartic frequency + previous 
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
Stop=0 #Stop signal

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
      if (r>h & r<=i) AA="Gly" else if (r>j & r<=k) AA="Asp" else AA="Other"  
      #determination of the amino acid based on the result of r
      G[l]=AA #the new AA is added to the growing peptide
      if (l==lmax) break #when lmax is reached the loop stops
    }
    if (lmax==3 & G[1]=="Asp" & G[2]=="Gly" & G[3]=="Asp") Stop=Stop+1
    #if DGD peptide is generated Stop signal is set to 1
    if (Stop==1) print (G) #print successful peptide for visual checking
    if (Stop==1) print (NS) #print nº of simulation to check the progress
    G<-c() #the AAs vector is emptied into each peptide
    l<-0 #length is restored
    if (Stop==1) break #when Stop signal is set to 1 the loop stops
  }
  Stop<-0 #Stop signal reset to 0
  Nse[NS]<-Ns #The number of successful simulations is added to its vector
  if(NS==Smax) break 
  #When the stipulated number of successful simulations is reached, the loop stops
}

library(modeest)
mlv(1/Nse,method = "meanshift") 
#Mode of the probability of the appearance of DGD
par(mar = c(5, 6, 4, 2)) 
hist(1/Nse,col="darkred", breaks=1500, xlim=c(0,mean(1/Nse)), main = "DGD", 
     xlab = "Probability of the sequence\n(between 0 and 1)", 
     ylab = "Frequency of Observed Probabilities\n(n = 150)")
abline(v=mlv(1/Nse, method="meanshift"), col="black", lty=2)
#Histogram of the probability of DGD appearing
