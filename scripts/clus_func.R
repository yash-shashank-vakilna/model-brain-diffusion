cluster_data <- function (csv_df, algo = 'gmm', nc = 4, if_fa = 1)
{
  if(if_fa == 1){
    require(dplyr)
    calc_FA <- function(lp,lt){
      lm = lm=(lp+2*lt)/3
      fa = sqrt(1.5*((lp =lp-lm)^2+2*(lt-lm)^2)
                /((lp)^2+(lt)^2))
    }
    csv_df=mutate(csv_df, fa = calc_FA(lp,exp(log.lt)))
  }
  #diff_model = Mclust(csv_df, modelNames = "VVV", G =1:5)
  #log-GMM
  if(algo == 'gmm'){
    require(mclust)
    diff_model = Mclust(csv_df, modelNames = "VVV", nc)
    to.plot = data.frame(csv_df,cluster = factor(diff_model$classification))
  }
  #log-kmeans
  else if(algo == 'kmeans'){
    require(fpc)
    km.fit.log = kmeans(csv_df, nc, nstart = 20, iter.max = 40)
    to.plot = data.frame(csv_df, cluster = factor(km.fit.log$cluster))
  }
    return(to.plot)
}
