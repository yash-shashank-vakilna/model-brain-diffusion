#setup data
setwd("~/brian/final_plots/select/csv_for_selected")
csv_table = read.csv("~/brian/final_plots/select/csv_for_selected/102513_gmm_clus.csv",header = FALSE)   #eg of clustered data
#csv_table=csv_table_full[sample(nrow(csv_table_full),replace=F,size=5e4),]
csv_table[csv_table$V3==5,]$V3=4
csv_df = data.frame(lp=csv_table[,1],log.lt=csv_table[,2],clus=factor(csv_table[,3]))
rm("csv_table")

#summarise calculate FA, exponentiate log.lt and calculate means
library(dplyr)
calc_FA <- function(lp,lt){
  lm = lm=(lp+2*lt)/3
  fa = sqrt(1.5*((lp =lp-lm)^2+2*(lt-lm)^2)
            /((lp)^2+(lt)^2))
}
csv_df=mutate(csv_df,lt = exp(log.lt), fa = calc_FA(lp,lt))
clus_select = group_by(csv_df,clus)
summarise(clus_select,mn_lp = mean(lp), sd_lp = sd(lp), mean_lt = mean(lt), sd_lt = sd(lt))

#plot histogram
library(ggplot2)

plot = ggplot(csv_df,aes(x = lp, y = ((log.lt))))
                +geom_point(aes(col=clus),alpha  = 1/3) #, shape = 20

ticks = c(0.02,0.05,0.1,0.2,0.5,1,2)
exp_formatter <- function(x){    #to convert to log scale
 {
    exp_l = exp(x)
    lab = sprintf('%g',exp_l)
  }
}
log.plot = plot + scale_y_continuous(label=exp_formatter,breaks=log(c(0.02,0.05,0.1,0.2,0.5,1,2)),
                                     name = expression(paste("Transverse diffusion coefficient ",lambda["\u27C2"]," [" ,"\u03BC","m"^{2},"/ms]")),
                                     limits = log(c(0.01,2.3)))   #x-label
log.plot = log.plot + scale_x_continuous(breaks = seq(0,3.3,0.5),
                                         name = expression(paste("Longitudinal diffusion coefficient ",lambda["\u2225"]," [" ,"\u03BC","m"^{2},"/ms]"))
                                         )    #y-klabel
log.plot+ scale_color_manual(values=c("red","yellow", "green", "blue")) 


lp = c(2.98,2.07,2.86,2.77)
lt = c(1.34,0.314,0.395,0.138)
