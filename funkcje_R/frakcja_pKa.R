library(openxlsx)

DS <- read.xlsx("G:/Model_ACN_1026/1_dane/database_Stan_pelna.xlsx")
PKA <- read.csv("G:/Model_ACN_1026/isocratic_database_ACD_data2.csv")

# z matlaba
fr <- function(pKas, pHo){
  n = length(pKas)
  cpKas = cumsum(pKas)
  x <- c()
  x[1] = 1
  for(i in 1:n){
    x[i+1] <- 10^(i*pHo-cpKas[i])
  }
  return(x/sum(x))
}





fr <- function(pKas, pHo){
  n = length(pKas)
  cpKas = cumsum(sort(pKas))
  x <- c()
  x[n+1] = 1
  for(i in 1:n){
    x[i] <- 10^(i*pHo-cpKas[i])
  }
  return(x/sum(x))
}



x <- matrix(rep(NA,5*140),5,140)
for(i in 1:140){
  x[,i] <- fr(PKA$pKa_ACD[which(PKA$Compound==1)],i*0.1)
}

plot(seq(0.1,14,by=0.1),x[1,],type="l",ylab="fractions",xlab="pH",ylim=c(0,1),main=paste(unique(DS$Names[which(DS$ID==1)])))
lines(seq(0.1,14,by=0.1),x[2,],col=2)
lines(seq(0.1,14,by=0.1),x[3,],col=3)
lines(seq(0.1,14,by=0.1),x[4,],col=4)
lines(seq(0.1,14,by=0.1),x[5,],col=5)

x <- matrix(rep(NA,3*140),3,140)
for(i in 1:140){
  x[,i] <- fr(PKA$pKa_ACD[which(PKA$Compound==2)],i*0.1)
}
plot(seq(0.1,14,by=0.1),x[1,],type="l",ylab="fractions",xlab="pH",ylim=c(0,1),main=paste(unique(DS$Names[which(DS$ID==2)])))
lines(seq(0.1,14,by=0.1),x[2,],col=2)
lines(seq(0.1,14,by=0.1),x[3,],col=3)



DS <- read.xlsx("G:/Model_ACN_1026/1_dane/database_Stan_pelna.xlsx")
DS <- read.xlsx("E:/Model_ACN_1026/1_dane/database_Stan_pelna.xlsx")
DS_PKA <- DS[!duplicated(DS$ID),-c(3,4)]


library(tidyr)
library(dplyr)


PKA <- DS_PKA %>% 

  gather(PKA_ozn, PKA, c(pKa1_ACD, pKa2_ACD, pKa3_ACD, pKa4_ACD, pKa5_ACD, pKa6_ACD, pKa7_ACD, pKa8_ACD, pKa9_ACD),
         -ID,-Names,-logP_ACD,-MW_ACD,-PSA_ACD) %>%
  gather(PKA_ozn_type, PKA_type, c(pKa1_ACD_type, pKa2_ACD_type, pKa3_ACD_type, pKa4_ACD_type, pKa5_ACD_type, pKa6_ACD_type,
                              pKa7_ACD_type, pKa8_ACD_type, pKa9_ACD_type),-ID,-Names,-logP_ACD,-MW_ACD,-PSA_ACD) %>%
  select(ID, Names, PKA_ozn, PKA,PKA_ozn_type, PKA_type)

PKA %>% sort(ID)
arrange(PKA,ID)


PKA <- DS_PKA %>% 
  gather(PKA_ozn, PKA, c(pKa1_ACD, pKa2_ACD, pKa3_ACD, pKa4_ACD, pKa5_ACD, pKa6_ACD, pKa7_ACD, pKa8_ACD, pKa9_ACD))

PKA <- PKA %>%
  gather(PKA_ozn_type, PKA_type, c(pKa1_ACD_type, pKa2_ACD_type, pKa3_ACD_type, pKa4_ACD_type, pKa5_ACD_type, pKa6_ACD_type,
                                   pKa7_ACD_type, pKa8_ACD_type, pKa9_ACD_type),-PKA_ozn,-PKA) %>%
  select(ID, Names, PKA_ozn, PKA,PKA_ozn_type, PKA_type)
