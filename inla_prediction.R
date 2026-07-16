library(INLA)
library(sf)
library(dplyr)
# read township point
setwd("./dataset")
dsn <- "Stack_Shanghai.gdb"
lay <- "HSSES_final_2020_point_albers"
system.time(df_Shanghai_point <- st_read(dsn = dsn, layer = lay))
df_Shanghai_point <- st_zm(df_Shanghai_point,drop = T,what = "ZM")
# get lon and lat
coordinates <- st_coordinates(df_Shanghai_point) 
meta_data <- cbind(meta_data,coordinates)

set.seed(1002)
coords = coordinates
st_crs(df_Shanghai_point)
#plot(df_Shanghai_point)
Shanghai_bdry_albers <- st_as_sf(Shanghai_bdry_albers)
Shanghai_bdry_albers <- st_simplify(Shanghai_bdry_albers, dTolerance = 50)
bdry <- inla.sp2segment(Shanghai_bdry_albers)
plot(bdry)
bdry$loc <- inla.mesh.map(bdry$loc)


# establishing mesh
dis <- diff(range(st_coordinates(df_Shanghai_point)[,1])) 
max.edge <- c(dis / 30, dis / 10 )
offset <- dis / 8
mesh0 <- inla.mesh.2d(loc = coords, boundary = bdry, max.edge=max.edge,offset = offset,cutoff = dis / 100)
plot(mesh0)
points(coords, col = 2, pch = 16, cex = 0.1)

A<-inla.spde.make.A(mesh=mesh0,loc=as.matrix(coords));dim(A)
spde <- inla.spde2.matern(mesh0, alpha=2)
iset <- inla.spde.make.index(name = "spatial.field", spde$n.spde)
stk <- inla.stack(data=list(y=meta_data$meta_y), #the response
                  
                  A=list(A,1),  #the A matrix; the 1 is included to make the list(covariates)
                  
                  effects=list(c(list(Intercept=1), #the Intercept
                                 iset),  #the spatial index
                               #the covariates
                               list(RF = meta_data$RF,
                                    SVM = meta_data$SVM,
                                    ET = meta_data$ET
                                    )
                  ), 
                  
                  #this is a quick name so you can call upon easily
                  tag='dat')



###model selection with WAIC (other criteria can be used)
for(i in 1:3){
  
  f1 <- as.formula(paste0("y ~ -1 + Intercept + f(spatial.field, model=spde) + ", paste0(c("RF", "SVM", "ET")[1:i], collapse = " + ")))
  
  model1<-inla(f1, data=inla.stack.data(stk,spde=spde),family= 'gaussian', 
               control.predictor=list(A=inla.stack.A(stk),compute=TRUE),
               control.compute = list(dic = TRUE, cpo=TRUE, waic = TRUE)) #verbose=TRUE,
  
  model_selection <- if(i==1){rbind(c(model = paste(c("RF", "SVM", "ET")[1:i]),waic = model1$waic$waic))}else{rbind(model_selection,c(model = paste(c("RF", "SVM", "ET")[1:i],collapse = " + "),waic = model1$waic$waic))
  }
  
}


model_selection

formula0<-y ~ -1 + Intercept + RF + SVM + ET + f(spatial.field, model=spde) 
inla.model0<-inla(formula0, #the formula
                  data=inla.stack.data(stk,spde=spde),  #the data stack
                  family= 'gaussian',   #which family the data comes from
                  control.predictor=list(A=inla.stack.A(stk),compute=TRUE),  #compute gives you the marginals of the linear predictor
                  control.compute = list(dic = TRUE, waic = TRUE, config = TRUE), #model diagnostics and config = TRUE gives you the GMRF
                  verbose = FALSE) #can include verbose=TRUE to see the log of the model runs

inla.model0$summary.fix
inla.model0$summary.hyperpar
##observe the plots for fixed parameters
par(mfrow=c(2,3))
plot(inla.model0$marginals.fixed[[1]], ty = "l", xlab = expression(beta[0]), ylab = "Density") 
plot(inla.model0$marginals.fixed[[2]], ty = "l", xlab = expression(beta[RF]), ylab = "Density") 
plot(inla.model0$marginals.fixed[[3]], ty = "l", xlab = expression(beta[SVM]), ylab = "Density") 
plot(inla.model0$marginals.fixed[[4]], ty = "l", xlab = expression(beta[ET]), ylab = "Density") 

inla.model0.res<-inla.spde2.result(inla.model0, 'spatial.field', spde, do.transf=TRUE)
#observe the plots for hyper parameters
par(mfrow=c(1,3))
plot(inla.model0.res$marginals.var[[1]], ty = "l", xlab = expression(sigma[randomfield]^2), ylab = "Density") 
plot(inla.model0.res$marginals.kap[[1]], type = "l", xlab = expression(kappa), ylab = "Density")
plot(inla.model0.res$marginals.range[[1]], type = "l", xlab = "range nominal", ylab = "Density")

#looking at the spatial field and what it looks like
gproj <- inla.mesh.projector(mesh0)
g.mean <- inla.mesh.project(gproj, inla.model0$summary.random$spatial.field$mean)
g.sd <- inla.mesh.project(gproj, inla.model0$summary.random$spatial.field$sd)
library(gridExtra)
grid.arrange(levelplot(g.mean, scales=list(draw=F), xlab='', ylab='', main='mean',col.regions = heat.colors(16)),
             levelplot(g.sd, scal=list(draw=F), xla='', yla='', main='sd' ,col.regions = heat.colors(16)), nrow=1)

xmin <- min(gproj$x)  # x min
xmax <- max(gproj$x)  # x max
ymin <- min(gproj$y)  # y min
ymax <- max(gproj$y)  # y max
ncols <- length(gproj$x)  # columns
nrows <- length(gproj$y)  # lines
r <- rast(xmin = xmin, xmax = xmax,
          ymin = ymin, ymax = ymax,
          ncols = ncols, nrows = nrows)

# fill raster
values(r) <- as.vector(g.sd[, ncols:1])
crs(r) <- crs(building)
plot(r)
print(r)
#writeRaster(r,"sptialrandomfieldSD.tif")
# Prediction
#remake the A matrix for prediction
Aprediction <- inla.spde.make.A(mesh = mesh0, loc = coords);
dim(Aprediction)
stk.pred <- inla.stack(data=list(y=NA), 
                       A=list(Aprediction,1), 
                       effects=list(c(list(Intercept=1)
                                      ,iset),
                                    list(RF = meta_data$RF,
                                         SVM = meta_data$SVM,
                                         ET = meta_data$ET
                                         )
                       ), 
                       tag='pred')

#join the prediction stack with the one for the full data
stk.full <- inla.stack(stk, stk.pred)
p.res.pred<-inla(formula0, data=inla.stack.data(stk.full,spde=spde),
                 family= 'gaussian', quantiles = NULL,
                 control.predictor=list(link = 1, A=inla.stack.A(stk.full),compute=FALSE),  #compute gives you the marginals of the linear predictor
                 control.compute = list(config = TRUE), #model diagnostics and config = TRUE gives you the GMRF
                 control.inla(strategy = 'simplified.laplace', huge = TRUE),  #this is to make it run faster
                 verbose = FALSE) 
# get index
index.pred <- inla.stack.index(stk.full, tag = "pred")$data

# get non-spatial parts results
linear.predictions <- p.res.pred$summary.linear.predictor[index.pred, ]
meta_data <- cbind(meta_data, linear.predictions$mean)




library(Metrics)
library(stats)

rf_pred <- meta_data$RF
svm_pred <- meta_data$SVM
et_pred <- meta_data$ET
meta_pred <- meta_data$meta_model
inla_pred <- meta_data$`linear.predictions$mean`
true_values <- meta_data$meta_y

# RF
rf_r2 <- cor(true_values, rf_pred)^2
rf_rmse <- rmse(true_values, rf_pred)
rf_mae <- mae(true_values, rf_pred)

# SVM
svm_r2 <- cor(true_values, svm_pred)^2
svm_rmse <- rmse(true_values, svm_pred)
svm_mae <- mae(true_values, svm_pred)

# ET
et_r2 <- cor(true_values, et_pred)^2
et_rmse <- rmse(true_values, et_pred)
et_mae <- mae(true_values, et_pred)

# OLS meta model
meta_r2 <- cor(true_values, meta_pred)^2
meta_rmse <- rmse(true_values, meta_pred)
meta_mae <- mae(true_values, meta_pred)

#inla meta model
inla_r2 <- cor(true_values, inla_pred)^2
inla_rmse <- rmse(true_values, inla_pred)
inla_mae <- mae(true_values, inla_pred)

# to a table
results <- data.frame(
  Model = c("Random Forest", "SVM", "ET", "meta_model", "inla_model"),
  R2 = c(rf_r2, svm_r2, et_r2, meta_r2, inla_r2),
  RMSE = c(rf_rmse, svm_rmse, et_rmse, meta_rmse,inla_rmse),
  MAE = c(rf_mae, svm_mae, et_rmse, meta_mae, inla_mae)
  
)

print(results)

# test set
## using results from Model1
model = inla.model0
## recall:: formula<-y ~ -1 +Intercept + f(spatial.field, model=spde) + Access + Elevation + EVI + LST_day

# Covariates for prediction points

RF<- accuracy.basemodel$RF
SVM <-  accuracy.basemodel$SVM
ET <- accuracy.basemodel$ET


#create the spatial structure
sfield_nodes <- model$summary.random$spatial.field['mean']
Aprediction <- inla.spde.make.A(mesh = mesh0, loc = reference.coordinates);
field <- (Aprediction %*% as.data.frame(sfield_nodes)[, 1])

#make empty matrix to fill predictions
pred <- matrix(NA, nrow = dim(Aprediction)[1], ncol = 1)

## Calculate Predicted values using regression formula
pred <- model$summary.fixed['Intercept', 'mean'] + 
  model$summary.fixed['RF', 'mean'] * RF +
  model$summary.fixed['SVM', 'mean'] * SVM +
  model$summary.fixed['ET', 'mean'] * ET +
  field 

# as df
inla.results <- as.data.frame(as.matrix(exp(pred)))
colnames(inla.results) <- "inla_exp"
accuracy.inlamodel <- cbind(accuracy.basemodel,inla.results)


inla_sum <- aggregate(inla_exp ~ county_id, data = accuracy.inlamodel, sum)
colnames(inla_sum) <- c("county_id", "inla_sum")
accuracy.inlamodel <- merge(accuracy.inlamodel, inla_sum, by = "county_id", all.x = TRUE)
accuracy.inlamodel$inla_pop <- (accuracy.inlamodel$inla_exp/accuracy.inlamodel$inla_sum)*accuracy.inlamodel$county_pop

inla.plot <- accuracy.inlamodel[, c("ID","RF", "RF_pop","SVM", "SVM_pop", "ET", "ET_pop", 
                                    "meta_model", "meta_pop", "inla_exp","inla_pop")]
inla.plot <- inla.plot[order(inla.plot$ID), ]
inla.plot <-cbind(inla.plot,reference.sf)
inla_plot <- rasterize(st_as_sf(inla.plot), building, field = "inla_exp")
plot(inla_plot)
inla_pop_plot <- rasterize(st_as_sf(inla.plot), building, field = "inla_pop")
plot(inla_pop_plot, breaks = c(0, 10, 25, 50, 100, 250, 500))
#writeRaster(inla_pop_plot, "inla_RDCV.tif")

inla_p_town <- aggregate(inla_pop ~ town_id, data = accuracy.inlamodel, sum)
town_first <- merge(town_first, inla_p_town, by = "town_id")

# get township level predictions
rf_pred <- town_first$RF_pop
svm_pred <- town_first$SVM_pop
et_pred <- town_first$ET_pop
meta_pred <- town_first$meta_pop
inla_pred <- town_first$inla_pop
true_values <- town_first$town_pop

# rf
rf_r2 <- cor(true_values, rf_pred)^2
rf_rmse <- rmse(true_values, rf_pred)
rf_mae <- mae(true_values, rf_pred)

# svm
svm_r2 <- cor(true_values, svm_pred)^2
svm_rmse <- rmse(true_values, svm_pred)
svm_mae <- mae(true_values, svm_pred)

# et
et_r2 <- cor(true_values, et_pred)^2
et_rmse <- rmse(true_values, et_pred)
et_mae <- mae(true_values, et_pred)

# linear meta
meta_r2 <- cor(true_values, meta_pred)^2
meta_rmse <- rmse(true_values, meta_pred)
meta_mae <- mae(true_values, meta_pred)

# inla meta
inla_r2 <- cor(true_values, inla_pred)^2
inla_rmse <- rmse(true_values, inla_pred)
inla_mae <- mae(true_values, inla_pred)

# to a table
results <- data.frame(
  Model = c("Random Forest", "SVM", "ET", "meta_model", "inla"),
  R2 = c(rf_r2, svm_r2, et_r2, meta_r2, inla_r2),
  RMSE = c(rf_rmse, svm_rmse, et_rmse,  meta_rmse, inla_rmse),
  MAE = c(rf_mae, svm_mae, et_mae,  meta_mae,inla_mae)
  
)
print(results)

#write.csv(results, file = "randomCV.csv")

#write.csv(town_first, file = "randomCVscatterplot.csv")


