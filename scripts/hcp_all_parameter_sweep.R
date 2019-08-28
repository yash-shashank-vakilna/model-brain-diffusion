
wd = "~/brian/data/hcp_output/rclus_output/"
smt_file="102513_smt.csv"
setwd(wd)
csv_table = read.csv(paste0("../smt_output/csv_file/",smt_file),header = T)
csv_table=csv_table[sample(nrow(csv_table),replace=F,size=20e4),]
csv_df = data.frame(lp=csv_table[,1],log.lt=log(csv_table[,2]))

library(ggplot2)
source("~/brian/scripts/clus_func.R")
for(if_fa in c(1,0)){
  for(algo_type in c('kmeans','gmm')){
    for (no_col in 10:20){
      #print(no_col)
      #algo_type = 'kmeans' #either 'gmm' or 'kmeans'
      to.plot = cluster_data(csv_df,algo_type,no_col,if_fa)
      # plotting
      plot_fn = gsub("_smt.csv",paste0("_",no_col,"_",algo_type,"_FA_",if_fa,
                                       ".png"),smt_file)
      csv_fn = gsub("_smt.csv",paste0("_",no_col,"_",algo_type,"_FA_",if_fa,
                                      "_rclus.csv"),smt_file)
      plot=ggplot(to.plot,aes(lp,log.lt,col = factor(cluster)))+geom_point(size = 0.001)
      sprintf("%s saved",plot_fn)
      to.save = data.frame(index = csv_table[3], class = as.numeric(to.plot$cluster))
      ggsave(plot_fn, plot,width = 12, height = 7.2)
      write.table(to.save,file=paste0(wd,csv_fn),col.names = FALSE, row.names = FALSE, quote = F)
      print(no_col)
    }}}

