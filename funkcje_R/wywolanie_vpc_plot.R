setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()

load("E:/Model_1024/AA/Fit_sample.Rsave")
DS <- read.csv(file="E:/Model_1024/AA/database_Stan.csv", header=TRUE, sep=',',dec = ",")


sim <- as.data.frame(t(fit_sample$logkPred))
sim <- gather(sim, key = "REP", V1:V4000)
sim$concentration <- rep(DS$concentration,4000)
colnames(sim)[2] <- 'logk'

obs <- DS

fi <- seq(0,0.99,by=0.099)


sim <- as.data.frame(t(fit_sample$logkPred))
sim <- gather(sim, key = "REP", V1:V4000)
sim$concentration <- rep(DS$concentration,4000)
colnames(sim)[2] <- 'logk'

obs <- DS

fi <- seq(0,0.99,by=0.099)

