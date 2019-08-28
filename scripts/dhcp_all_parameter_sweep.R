
setwd("~/brian/dhcp_data/rclus_output/")
smt_file="sub-CC00069XX12_ses-26300_smt.csv"
name = paste0("~/brian/dhcp_data/smt_output/csv_file/",smt_file)
csv_table = read.csv(name,header = T)
csv_table=csv_table[sample(nrow(csv_table),replace=F,size=10e4),]
csv_df = data.frame(lp=csv_table[,1],log.lt=log(csv_table[,2]))

library(ggplot2)
for (no_col in 5:10){
  #print(no_col)
  algo_type = 'gmm'
  to.plot = cluster_data(csv_df,algo_type,no_col,0)
  # plotting
  plot_fn = gsub("_smt.csv",paste0("_",no_col,"_",algo_type,".png"),smt_file)
  csv_fn = gsub("_smt.csv",paste0("_",no_col,"_",algo_type,"_rclus.csv"),smt_file)
  plot=ggplot(to.plot,aes(lp,log.lt,col = factor(cluster)))+geom_point(size = 0.001)
  sprintf("%s saved",plot_fn)
  to.save = data.frame(index = csv_table[3], class = as.numeric(to.plot$cluster))
  ggsave(paste0("~/brian/dhcp_data/rclus_output/",plot_fn), plot,width = 12, height = 7.2)
  write.table(to.save,file=paste0("~/brian/dhcp_data/rclus_output/",csv_fn),col.names = FALSE, row.names = FALSE, quote = F)
  no_col
}
