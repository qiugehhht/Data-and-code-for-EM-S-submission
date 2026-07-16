library(caret)
library(randomForest)
library(xgboost)
library(extraTrees)
library(glmnet)
library(e1071)
#model_list <- modelLookup()
#unique(model_list$model)
X_train <- as.matrix(df[, -16])
y_train <- df$y_data

set.seed(888)
# base models' prediction as Meta Model's inputs
cv_folds <- createFolds(y_train, k = 10)
meta_features <- data.frame()
meta_y <- c()
sample_indices <- c()  # indexing

for (fold in names(cv_folds)) { 
  train_idx <- setdiff(1:nrow(X_train), cv_folds[[fold]])
  vali_idx <- cv_folds[[fold]]
  
  # training base models
  rf_fit <- randomForest(X_train[train_idx, ], y_train[train_idx], mtry = 5, ntree = 100)

  et_fit <- extraTrees(X_train[train_idx, ], y_train[train_idx], 
                       mtry = 4, 
                       numRandomCuts = 1, 
                       ntree = 100)
  svm_fit <- svm(X_train[train_idx, ], y_train[train_idx],
                 cost = 1,
                 gamma = 0.2)
  
  # base model prediction
  rf_pred <- predict(rf_fit, X_train[vali_idx, ])
  et_pred <- predict(et_fit, X_train[vali_idx, ])
  svm_pred <- predict(svm_fit, X_train[vali_idx, ])
  
  # collecting base model prediction
  fold_meta <- data.frame(RF = rf_pred, SVM = svm_pred, ET=et_pred)
  meta_features <- rbind(meta_features, fold_meta)
  meta_y <- c(meta_y, y_train[vali_idx])
  sample_indices <- c(sample_indices, vali_idx)  # indexing
}

# order by indexing
sorted_indices <- order(sample_indices)
meta_features <- meta_features[sorted_indices, ]
meta_y <- meta_y[sorted_indices]


meta_data <- cbind(meta_features,meta_y)
colnames(meta_data)[3] <- "ET"

#  Meta Model (OLS) training
meta_fit <- lm(meta_y ~ RF+SVM+ET, data = meta_data)
summary(meta_fit)
meta_data$meta_model <- predict(meta_fit, meta_data)

library(Metrics)  # RMSE and MAE
library(stats)    # R2

# get base model predictions
rf_pred <- meta_data$RF
svm_pred <- meta_data$SVM
et_pred <- meta_data$ET
meta_pred <- meta_data$meta_model
true_values <- meta_data$meta_y

# get Random Forest model acc matrix
rf_r2 <- cor(true_values, rf_pred)^2
rf_rmse <- rmse(true_values, rf_pred)
rf_mae <- mae(true_values, rf_pred)

# get svm model acc matrix
svm_r2 <- cor(true_values, svm_pred)^2
svm_rmse <- rmse(true_values, svm_pred)
svm_mae <- mae(true_values, svm_pred)

# get extreme trees acc matrix
et_r2 <- cor(true_values, et_pred)^2
et_rmse <- rmse(true_values, et_pred)
et_mae <- mae(true_values, et_pred)

# get meta model (OLS) acc matrix
meta_r2 <- cor(true_values, meta_pred)^2
meta_rmse <- rmse(true_values, meta_pred)
meta_mae <- mae(true_values, meta_pred)


# to a table
results <- data.frame(
  Model = c("Random Forest", "SVM","ET", "meta_model"),
  R2 = c(rf_r2, svm_r2, et_r2, meta_r2),
  RMSE = c(rf_rmse, svm_rmse, et_rmse, meta_rmse),
  MAE = c(rf_mae, svm_mae, et_mae,meta_mae)
  
)
print(results)

# get acc in the test set 
town_id <- "./town_id.tif"
town_id <- rast(town_id)
town_id <- crop(town_id, ext(building))
town_id <- resample(town_id, building, method = "near")
town_id <- mask(town_id, building)
plot(town_id)
town_pop <- "./town_pop.tif"
town_pop <- rast(town_pop)
town_pop <- crop(town_pop, ext(building))
town_pop <- resample(town_pop, building, method = "near")
town_pop <- mask(town_pop, building)
plot(town_pop)
county_id <- "./county_id.tif"
county_id <- rast(county_id)
county_id <- crop(county_id, ext(building))
county_id <- resample(county_id, building, method = "near")
county_id <- mask(county_id, building)
plot(county_id)
county_pop <- "./county_pop.tif"
county_pop <- rast(county_pop)
county_pop <- crop(county_pop, ext(building))
county_pop <- resample(county_pop, building, method = "near")
county_pop <- mask(county_pop, building)
plot(county_pop)
# get other benchmark data for comparison
WorldPop <- "./WorldPop_Shanghai_albers.tif"
WorldPop <- rast(WorldPop)
WorldPop <- mask(WorldPop, Shanghai_bdry_albers)
print(WorldPop)
plot(WorldPop)

GPW <- "./GPW_albers.tif"
GPW <- rast(GPW)
GPW <- mask(GPW, Shanghai_bdry_albers)
print(GPW)
plot(GPW)
GHP <- "./GHP_albers.tif"
GHP <- rast(GHP)
GHP <- mask(GHP, Shanghai_bdry_albers)
print(GHP)
plot(GHP)
LandScan <- "./LandScan_albers.tif"
LandScan <- rast(LandScan)
LandScan <- mask(LandScan, Shanghai_bdry_albers)
print(LandScan)
plot(LandScan)



st_crs(building)
in.predict <- which(!is.na(values(building)))
reference.coordinates <- xyFromCell(building, which(!is.na(values(building))))

#make these into points and extract covariates for prediction grid
reference.sf <- st_as_sf(as.data.frame(reference.coordinates), 
                         coords = c("x", "y"), 
                         crs = crs(building))  #  building  CRS

pred.covs <- terra::extract(raster_stack, reference.sf, df=T)
ref_stack <- c(town_id, town_pop, county_id, county_pop)
names(ref_stack) <- c("town_id", "town_pop", "county_id", "county_pop")
pop.ref <- terra::extract(ref_stack, reference.sf, df=T)
sum(is.na(pred.covs))
# NA fill
for (i in 1:ncol(pred.covs)) {
  if (any(is.na(pred.covs[, i]))) {  
    mean_value <- mean(pred.covs[, i], na.rm = TRUE)  # mean to NA
    pred.covs[is.na(pred.covs[, i]), i] <- mean_value  
  }
}
pred.covs <- as.matrix(pred.covs[, -1])
pred.covs <- pred.covs[, colnames(X_train), drop = FALSE]
cat("features:", colnames(pred.covs), "\n")
cat("model features:", colnames(X_train) , "\n")

# base models for testing
rf_pred.test <- predict(rf_fit, pred.covs)

svm_pred.test <- predict(svm_fit, pred.covs)

et_pred.test <- predict(et_fit, pred.covs)


# base model accuracy data frame
accuracy.basemodel <- cbind(rf_pred.test,svm_pred.test,et_pred.test, pop.ref)
colnames(accuracy.basemodel) <- c("RF", "SVM", "ET",   "ID","town_id","town_pop","county_id","county_pop")

meta_pred.test <- predict(meta_fit, accuracy.basemodel)

accuracy.basemodel <- cbind(accuracy.basemodel,meta_pred.test)

accuracy.basemodel$RF_exp <- exp(accuracy.basemodel$RF)
accuracy.basemodel$SVM_exp <- exp(accuracy.basemodel$SVM)
accuracy.basemodel$ET_exp <- exp(accuracy.basemodel$ET)
accuracy.basemodel$meta_model <- exp(accuracy.basemodel$meta_pred.test)

RF_sum <- aggregate(RF_exp ~ county_id, data = accuracy.basemodel, sum)
colnames(RF_sum) <- c("county_id", "RF_sum")
SVM_sum <- aggregate(SVM_exp ~ county_id, data = accuracy.basemodel, sum)
colnames(SVM_sum) <- c("county_id", "SVM_sum")
ET_sum <- aggregate(ET_exp ~ county_id, data = accuracy.basemodel, sum)
colnames(ET_sum) <- c("county_id", "ET_sum")
meta_sum <- aggregate(meta_model ~ county_id, data = accuracy.basemodel, sum)
colnames(meta_sum) <- c("county_id", "meta_sum")
accuracy.basemodel <- merge(accuracy.basemodel, RF_sum, by = "county_id", all.x = TRUE)
accuracy.basemodel <- merge(accuracy.basemodel, SVM_sum, by = "county_id", all.x = TRUE)
accuracy.basemodel <- merge(accuracy.basemodel, ET_sum, by = "county_id", all.x = TRUE)
accuracy.basemodel <- merge(accuracy.basemodel, meta_sum, by = "county_id", all.x = TRUE)

accuracy.basemodel$RF_pop <- (accuracy.basemodel$RF_exp/accuracy.basemodel$RF_sum)*accuracy.basemodel$county_pop
accuracy.basemodel$SVM_pop <- (accuracy.basemodel$SVM_exp/accuracy.basemodel$SVM_sum)*accuracy.basemodel$county_pop
accuracy.basemodel$ET_pop <- (accuracy.basemodel$ET_exp/accuracy.basemodel$ET_sum)*accuracy.basemodel$county_pop
accuracy.basemodel$meta_pop <- (accuracy.basemodel$meta_model/accuracy.basemodel$meta_sum)*accuracy.basemodel$county_pop

base.plot <- accuracy.basemodel[, c("ID","RF_exp", "RF_pop","SVM_exp", "SVM_pop", "ET_exp", "ET_pop",  "meta_model", "meta_pop")]
base.plot <- base.plot[order(base.plot$ID), ]
base.plot <-cbind(base.plot,reference.sf)
rf_plot <- rasterize(st_as_sf(base.plot), building, field = "RF_pop")
plot(rf_plot)
#writeRaster(rf_plot, "rf_plot_RDCV.tif")
et_plot <- rasterize(st_as_sf(base.plot), building, field = "ET_pop")
plot(et_plot)
#writeRaster(et_plot, "et_plot_RDCV.tif")
svm_plot <- rasterize(st_as_sf(base.plot), building, field = "SVM_pop")
plot(svm_plot)
#writeRaster(svm_plot, "svm_plot_RDCV.tif")
meta_plot <- rasterize(st_as_sf(base.plot), building, field = "meta_pop")
plot(meta_plot)
#writeRaster(meta_plot, "meta_plot_RDCV.tif")

RF_p_town <- aggregate(RF_pop ~ town_id, data = accuracy.basemodel, sum)
SVM_p_town <- aggregate(SVM_pop ~ town_id, data = accuracy.basemodel, sum)
ET_p_town <- aggregate(ET_pop ~ town_id, data = accuracy.basemodel, sum)
meta_p_town <- aggregate(meta_pop ~ town_id, data = accuracy.basemodel, sum)

town_first <- accuracy.basemodel[!duplicated(accuracy.basemodel$town_id), ]
town_first <- na.omit(town_first[, c("town_id", "town_pop")])

town_first <- merge(town_first, RF_p_town, by = "town_id")
town_first <- merge(town_first, SVM_p_town, by = "town_id")
town_first <- merge(town_first, ET_p_town, by = "town_id")
town_first <- merge(town_first, meta_p_town, by = "town_id")

WorldPop_acc <- terra::extract(WorldPop, df_Shanghai_albers, fun = sum, na.rm = TRUE)
colnames(WorldPop_acc) <- c("ID", "WorldPop_acc")
PopSE_acc <- terra::extract(PopSE, df_Shanghai_albers, fun = sum, na.rm = TRUE)
colnames(PopSE_acc) <- c("ID", "PopSE_acc")

GPW_acc <- terra::extract(GPW, df_Shanghai_albers, fun = sum, na.rm = TRUE)
colnames(GPW_acc) <- c("ID", "GPW_acc")
GHP_acc <- terra::extract(GHP, df_Shanghai_albers, fun = sum, na.rm = TRUE)
colnames(GHP_acc) <- c("ID", "GHP_acc")
LandScan_acc <- terra::extract(LandScan, df_Shanghai_albers, fun = sum, na.rm = TRUE)
colnames(LandScan_acc) <- c("ID", "LandScan_acc")

other_product_acc <- cbind(df_Shanghai_albers,WorldPop_acc,PopSE_acc, GPW_acc, GHP_acc, LandScan_acc)

WorldPop_pred <- other_product_acc$WorldPop_acc
PopSE_pred <- other_product_acc$PopSE_acc
GPW_pred <- other_product_acc$GPW_acc
GHP_pred <- other_product_acc$GHP_acc
LandScan_pred <- other_product_acc$LandScan_acc

true_values <- other_product_acc$age_total


WorldPop_r2 <- summary(lm(true_values ~ WorldPop_pred))$r.squared
WorldPop_rmse <- rmse(true_values, WorldPop_pred)
WorldPop_mae <- mae(true_values, WorldPop_pred)

PopSE_r2 <- summary(lm(true_values ~ PopSE_pred))$r.squared
PopSE_rmse <- rmse(true_values, PopSE_pred)
PopSE_mae <- mae(true_values, PopSE_pred)

GPW_r2 <- summary(lm(true_values ~ GPW_pred))$r.squared
GPW_rmse <- rmse(true_values, GPW_pred)
GPW_mae <- mae(true_values, GPW_pred)

GHP_r2 <- summary(lm(true_values ~ GHP_pred))$r.squared
GHP_rmse <- rmse(true_values, GHP_pred)
GHP_mae <- mae(true_values, GHP_pred)

LS_r2 <- summary(lm(true_values ~ LandScan_pred))$r.squared
LS_rmse <- rmse(true_values, LandScan_pred)
LS_mae <- mae(true_values, LandScan_pred)

# MAPE and fit line with out intercept
MAPE = function(m, o){
  mean(abs((m - o) / o)) * 100
}
WorldPop_r2_no_intercept <- summary(lm(true_values ~ WorldPop_pred-1))$r.squared
WorldPop_MAPE <- MAPE(true_values, WorldPop_pred)
print(WorldPop_r2_no_intercept)
print(WorldPop_MAPE)

GPW_r2_no_intercept <- summary(lm(true_values ~ GPW_pred-1))$r.squared
GPW_MAPE <- MAPE(true_values, GPW_pred)
print(GPW_r2_no_intercept)
print(GPW_MAPE)

GHP_r2_no_intercept <- summary(lm(true_values ~ GHP_pred-1))$r.squared
GHP_MAPE <- MAPE(true_values, GHP_pred)
print(GHP_r2_no_intercept)
print(GHP_MAPE)

LS__no_intercept <- summary(lm(true_values ~ LandScan_pred-1))$r.squared
LS_MAPE <- MAPE(true_values, LandScan_pred)
print(LS__no_intercept)
print(LS_MAPE)

results <- data.frame(
  Model = c("WorldPop", "PopSE", "GPW", "GHP", "LS"),
  R2 = c(WorldPop_r2, PopSE_r2, GPW_r2, GHP_r2, LS_r2),
  RMSE = c(WorldPop_rmse, PopSE_rmse, GPW_rmse, GHP_rmse, LS_rmse),
  MAE = c(WorldPop_mae, PopSE_mae, GPW_mae, GHP_mae, LS_mae)
  
)
print(results)
# get township level results
rf_pred <- town_first$RF_pop
svm_pred <- town_first$SVM_pop
et_pred <- town_first$ET_pop
meta_pred <- town_first$meta_pop
true_values <- town_first$town_pop

# testing acc matrix RF
rf_r2 <- summary(lm(true_values ~ rf_pred))$r.squared
rf_rmse <- rmse(true_values, rf_pred)
rf_mae <- mae(true_values, rf_pred)

# testing acc matrix SVM
svm_r2 <- summary(lm(true_values ~ svm_pred))$r.squared
svm_rmse <- rmse(true_values, svm_pred)
svm_mae <- mae(true_values, svm_pred)

# testing acc matrix ET
et_r2 <- summary(lm(true_values ~ et_pred))$r.squared
et_rmse <- rmse(true_values, et_pred)
et_mae <- mae(true_values, et_pred)

# testing acc matrix OLS meta
meta_r2 <- summary(lm(true_values ~ meta_pred))$r.squared
meta_rmse <- rmse(true_values, meta_pred)
meta_mae <- mae(true_values, meta_pred)

# to table
results <- data.frame(
  Model = c("Random Forest", "SVM","ET", "meta_model"),
  R2 = c(rf_r2, svm_r2, et_r2, meta_r2),
  RMSE = c(rf_rmse, svm_rmse, et_rmse, meta_rmse),
  MAE = c(rf_mae, svm_mae,  et_mae,meta_mae)
  
)
print(results)



