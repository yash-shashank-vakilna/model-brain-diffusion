
csv_table = read.csv("~/brian/pgms/r_for_gmm/output/102513_smt.csv",
#                     header = TRUE)   #eg of clustered data
csv_table = read.csv("~/brian/dhcp_data/csv_smt_output/sub-CC00069XX12_ses-26300_smt.csv",
                     header = T)
csv_table=csv_table[sample(nrow(csv_table),replace=F,size=10e4),]
csv_df = data.frame(lp=csv_table[,1],log.lt=log(csv_table[,2]))

#summarise calculate FA, exponentiate log.lt and calculate means
library(dplyr)
calc_FA <- function(lp,lt){
  lm = lm=(lp+2*lt)/3
  fa = sqrt(1.5*((lp =lp-lm)^2+2*(lt-lm)^2)
            /((lp)^2+(lt)^2))
}
csv_df=mutate(csv_df, fa = calc_FA(lp,exp(log.lt)))


#diff_model = Mclust(csv_df, modelNames = "VVV", G =1:5)
#log-GMM
library(mclust)
diff_model = Mclust(csv_df, modelNames = "VVV", G =10)
to.plot = data.frame(csv_df,class = factor(diff_model$classification),index = csv_table[,3])
library(ggplot2)
ggplot(to.plot, aes(lp,((log.lt)), col = class))+geom_point(size = 0.001)
write.csv(to.plot,file="102513_rclus.csv",col.names = FALSE)

#log-kmeans
library(fpc)
km.fit.log = kmeans(csv_df, 10, nstart = 20, iter.max = 40)
to.plot = data.frame(csv_df, cluster = km.fit.log$cluster)
to.save = data.frame(index = csv_table[3], class = km.fit.log$cluster)
library(ggplot2)
ggplot(to.plot,aes(lp,log.lt,col = factor(cluster)))+geom_point(size = 0.001)
write.table(to.save,file="~/brian/dhcp_data/csv_smt_output/sub-CC00069XX12_ses-26300_rclus.csv",col.names = FALSE, row.names = FALSE)
