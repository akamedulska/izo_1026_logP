vpc_plot <- function(concentration_Obs,logk_Obs,logk_Pred,fi){
  
  library(ggplot2)
  library(lattice)
  library(Cairo)
  # zapewnia pliki bitmap PNG, JPEG, TIFF, pliki PDF o wysokiej rozdzielczoœci
  library(scales)
  # poprawia estetykê grafiki i zapewnia metody automatycznego okreœlania przerw, etykiet dla osi i legend.
  library(gridExtra)
  library(dplyr)
  library(tidyr)
  
  # data preparation
  sim <- as.data.frame(t(logk_Pred))
  sim <- gather(sim, key = "REP")
  sim$concentration <- rep(concentration_Obs,4000)
  colnames(sim)[2] <- 'logk'
  
  DS1 <- as.data.frame(cbind(concentration_Obs,logk_Obs))
  colnames(DS1) <- c("concentration","logk")
  
  # observed
  DS1 <- mutate(DS1,
                bin = cut(DS1$concentration,   
                          breaks = fi, 
                          labels=1:(length(fi)-1)))
  
  by_obs <- group_by(DS1, bin)
  
  sum_obs <- summarise(by_obs, 
                       mediana = median(logk), 
                       q_0.05 = quantile(logk, probs=0.05),
                       q_0.95 = quantile(logk, probs=0.95))
  
  # simulated
  sim <- mutate(sim, bin = cut(sim$concentration, 
                               breaks = fi, 
                               labels=1:(length(fi)-1)))
  
  by_sim <- group_by(sim, bin, REP)
  
  sum_sim <- summarise(by_sim, 
                       mediana = median(logk), 
                       q_0.05 = quantile(logk, probs=0.05),
                       q_0.95 = quantile(logk, probs=0.95))
  
  CI <- summarise(sum_sim, 
                  CI_med_low = quantile(mediana, probs=0.05),
                  CI_med_up = quantile(mediana, probs=0.95),
                  CI_q_0.05_low = quantile(q_0.05, probs=0.05),
                  CI_q_0.05_up = quantile(q_0.05, probs=0.95),
                  CI_q_0.95_low = quantile(q_0.95, probs=0.05),
                  CI_q_0.95_up = quantile(q_0.95, probs=0.95))
  
  #œrodki binów
  middle <-  (fi[1:(length(fi)-1)] + fi[2:(length(fi))]) / 2 #zmienna niezale¿na, czas dla œrodków binów
  
  #zmienna niezale¿na - doklejam do df srodki binów
  sum_obs$x <- middle
  
  ###########################################################################
  # Figures -----------------------------------------------------------------
  ###########################################################################
  
  plot <- ggplot() +  
    #wygl¹d wykresu
    ggtitle(NULL) +
    theme_classic() +
    theme(axis.text=element_text(size=12),
          axis.title=element_text(size=14)) +
    xlab(expression(varphi)) + ylab("Log k") +
    
    #w³aœciwy wykres  
    geom_rect(aes(xmin = fi[1:(length(fi)-1)], 
                  xmax = fi[2:length(fi)],   
                  ymin = CI$CI_med_low, 
                  ymax = CI$CI_med_up), 
              fill = "turquoise", 
              alpha=0.5) +
    geom_rect(aes(xmin = fi[1:(length(fi)-1)], 
                  xmax = fi[2:(length(fi))],   
                  ymin = CI$CI_q_0.05_low[1:(length(fi)-1)], 
                  ymax = CI$CI_q_0.05_up[1:(length(fi)-1)]), 
              fill = "blue", 
              alpha=0.2) +
    geom_rect(aes(xmin = fi[1:(length(fi)-1)], 
                  xmax = fi[2:(length(fi))],   
                  ymin = CI$CI_q_0.95_low[1:(length(fi)-1)], 
                  ymax = CI$CI_q_0.95_up[1:(length(fi)-1)]), 
              fill = "blue", 
              alpha=0.2) +
    geom_line(aes(x=sum_obs$x, y=sum_obs$mediana), size=0.6) +
    geom_line(aes(x=sum_obs$x, y=sum_obs$q_0.95), size=0.6, linetype = 2) +
    geom_line(aes(x=sum_obs$x, y=sum_obs$q_0.05), size=0.6, linetype =2) +
    
    geom_point(aes(x=sum_obs$x, y=sum_obs$mediana), shape=21, size=1, fill="black") +
    geom_point(aes(x=sum_obs$x, y=sum_obs$q_0.95), shape=21, size=1, fill="black") +
    geom_point(aes(x=sum_obs$x, y=sum_obs$q_0.05), shape=21, size=1, fill="black") 
  
  #geom_point(aes(x=obs$TIME, y=obs$DV), size=3, colour="grey42", alpha = 0.5) +
  
  #ggsave(filename = "VPC.png", plot=plot, width=10, height=5, type = "cairo-png")
  #line types 0 = blank, 1 = solid, 2 = dashed, 3 = dotted, 4 = dotdash, 5 = longdash, 6 = twodash
  
return(plot)
  
}

