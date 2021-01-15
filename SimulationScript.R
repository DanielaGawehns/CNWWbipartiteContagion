
# Script for Simulation Study
# Based on CNWW Report ...... Link.... 

################################

# Total population size of each animal
#input: N (size of population), Place (number of locations), 
#probVec(vector of length Place with probabilities of an individual having an additional link to that place)

createFinalAMat<- function (N,Place,probVec) {
  
  #initalize empty matrix
  InitialA <- matrix(0,N,Place)
  
  #devide population equally across places:
  InitialA[1:(N/Place),1] <- 1
  
  InitialA[((N/Place)+1):(2*N/Place),2] <- 1
  
  InitialA[((2*N/Place)+1):N,3] <- 1

  #randomly assign individuals to additional places, depending on probVec input
  
  FirstA <- cbind(InitialA[1:(N/Place),1], matrix(rbinom(N,1,probVec [1]),(N/Place),2))
  
  SecondA <- cbind(matrix(rbinom(N,1,probVec [2]),(N/Place),1),InitialA[((N/Place)+1):(2*N/Place),2],matrix(rbinom(N,1,probVec [2]),(N/Place),1))
  
  ThirdA <- cbind(matrix(rbinom(N,1,probVec [3]),(N/Place),2), InitialA[((2*N/Place)+1):N,3])
  
  FinalA <- rbind(FirstA, SecondA, ThirdA)
  
  return(FinalA)
}

######################################

updateFunction <- function (Iterations, FinalA, N, contamVec, r, gamma,beta) {

  #Initialize a number of empty lists: 
  
  # Number of susceptible individuals 
  Suscep <- list()
  Suscep[[1]] <- matrix(1,N,1)

  # Number of infected individuals 
  Infect <- list()
  Infect[[1]] <- matrix(0,N,1)
  Infect[[1]][c(sample(1:N, 0.01*N, replace=FALSE)),] <- 1

  # Number of recovered individuals 
  Recover <- list()
  Recover[[1]] <- matrix(0,N,1)

  # Contamination level
  E <- list()
  E[[1]] <- matrix(contamVec,Place,1) 
  
  #initalize output:
  TotalSuscep <- numeric()
  TotalInfect <- numeric()
  TotalRecover <- numeric()
  
  #run for i Iterations
  for (i in 2:Iterations){
    
    Suscep[[i]] <- Suscep[[i-1]]*prod(1-FinalA%*%tanh(alpha*E[[i-1]]))
    
    Infect[[i]] <- (1-r)*Infect[[i-1]]+Suscep[[i-1]]*(1-prod(1-FinalA%*%tanh(alpha*E[[i-1]])))
    
    Recover[[i]] <- Recover[[i-1]] + r*Infect[[i-1]]
    
    E[[i]] <- E[[i-1]]*(1-gamma) + beta* t(FinalA)%*%Infect[[i]]
    
    TotalSuscep [i] <- mean(Suscep[[i]] )
    TotalInfect [i] <- mean(Infect[[i]] )
    TotalRecover[i] <- mean(Recover[[i]])
  }
  
  return(list(TotalSuscep,TotalInfect,TotalRecover))
}

###############################
###############################

##Parameters to be set for simulation: 

# Recovery rate parameter
r <- 0.8

# Environment virus decay 
gamma <- 0.1

# Environment infection parameter
alpha <- ((r*gamma)/(N+1))

#shedding rate per location
beta<- c(0.1,0.4,0.8)

# Number of iterations
Iterations <- 20

#size of population
N <- 300

#number of locations
Place <- 3

#vector with probabilities of having additional link to a given location
# Choose the proportion of second degree 
# some individuals contact several places; 
#probVec[1]: probability of having a second degree at location1; other indices respectively
#sums to 1
#probVec <- c( 0.4, 0.4, 0.2)
probVec <- c(0,0,0)

#contamination levels
contamVec<- c(5,10,15)

#####################################

#run Simulation
FinalA <- createInitialAMat (N,Place,probVec)
simulationOutput<- updateFunction (Iterations, FinalA, N, contamVec, r, gamma)

###################################

#Plot Simulation

plot(seq(1,Iterations,1), simulationOutput[[1]], type="l", col=c("blue"))

lines(simulationOutput [[2]], col=c("red"))

lines(simulationOutput [[3]], col=c("green"))


###################

#300 indiv. 3 environment
#change % of people in 3 environment

#First condition:
#1/3, 1/3, 1/3 - everyone 1 connection

#different steps/percentages also in different environment



####################
#make ggplot graph

install.packages("ggplot2")
library("ggplot2")

data<- data.frame(Iterations = rep(seq(1,Iterations,1),3), 
                  Vals = c(simulationOutput [[1]],simulationOutput [[2]],simulationOutput [[3]] ),
                  Type = rep(c("TotalSusceptible", "TotalInfected", "TotalRecovered"), each=20) )

ggplot(data=data, aes(x=Iterations, y=Vals, group=Type)) +
  geom_line(aes(linetype=Type, color = Type))+
  geom_point()


