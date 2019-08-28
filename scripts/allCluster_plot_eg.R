#load data  
setwd('/home/yashvakilna/brian/pgms/print2r/class_plots/')
diff = read.csv("102513_r.csv")
diff_sub=diff[sample(nrow(diff),replace=F,size=40e4),]
diff_sub_log = data.frame(lp=diff_sub$lp, "log.lt"=(log(diff_sub$lt)))
diff_log = data.frame(lp=diff$lp, "log.lt"=(log(diff$lt)))

require(factoextra)
#perform EM
library(mclust)
#plot(diff_sub_log, pch=".")
diff_model = Mclust(diff_sub, modelNames = "VVV", G =5)
to.plot = data.frame(diff_sub,class = factor(diff_model$classification))
library(ggplot2)
ggplot(to.plot, aes(lp,log(lt), col = class))+geom_point()

diff_model = Mclust(diff_sub_log, modelNames = "VVV", G =5)
data.plot = data.frame(diff_sub_log, class = factor(as.vector(diff_model$classification)))
ggplot(data.plot, aes(lp,log.lt, col = class))+geom_point()
library(ggplot2)
diff_model_log = Mclust(diff_sub_log, modelNames = "VVV", G =1:5)
fviz_cluster(diff_model_log,diff_sub_log, ellipse =TRUE, geom = "point")

library(fpc)
km.fit = kmeans(diff_sub, 5, nstart = 40)
fviz_cluster(km.fit,diff_sub_log, ellipse =TRUE, geom = "point")
km.fit.log = kmeans(diff_sub_log, 10, nstart = 20)
fviz_cluster(km.fit.log,diff_sub_log, ellipse =TRUE, geom = "point")

library(dbscan)
kNNdistplot(as.matrix(diff_sub), k=5)
abline(0.02,0, col = "red")
db.fit=dbscan(diff_sub, 0.02)
fviz_cluster(db.fit,diff_sub_log, ellipse =TRUE, geom = "point")

#just density plot of smt data
library(ggplot2)
p1 = ggplot(diff,aes(x = lp, y = (log(lt))))
plot(p1)
p2 = ggplot(diff_sub,aes(x = lp, y = (log(lt))))+ stat_bin2d( bins = 512)

m = ggplot(diff_sub, aes(x = lp, y = (log(lt/b0)/b1))) +
  geom_point() 
  
m + geom_density_2d()+ stat_density_2d(geom = "raster", aes(fill = ..density..), contour = FALSE, n = 512)
