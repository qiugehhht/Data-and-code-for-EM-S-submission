library(sf)
library(dplyr)
library(moments)

# 读取township-level数据
setwd("H:/INLA/Stack/Stack_Shanghai")
dsn <- "Stack_Shanghai.gdb"
lay <- "HSSES_final_2020_albers"
system.time(df_Shanghai_albers <- st_read(dsn = dsn, layer = lay))
df_Shanghai_albers <- st_zm(df_Shanghai_albers,drop = T,what = "ZM")
# 检查当前的CRS（坐标参考系统）
st_crs(df_Shanghai_albers)

# 投影变换为 Albers Equal Area Conic
#albers_proj <- "+proj=aea +lat_1=25 +lat_2=47 +lat_0=40 +lon_0=110 +datum=WGS84 +units=m +no_defs"
#df_Shanghai_albers <- st_transform(df_Shanghai, crs = albers_proj)

# 检查 age_0 到 age_100 的列类型
df_Shanghai_albers <- df_Shanghai_albers %>%
  mutate(across(starts_with("age_"), as.numeric)) %>%
  mutate(age_total = rowSums(across(starts_with("age_"))))

# check population data quality
#df_Shanghai_albers$sex_total <- df_Shanghai_albers$male + df_Shanghai_albers$female
#inconsistent_indices <- which(df_Shanghai_albers$age_total != df_Shanghai_albers$sex_total)
#print(inconsistent_indices)
#print(result)
#county <- aggregate(age_total ~ district, data = df_Shanghai_albers, sum)


df_Shanghai_albers$area_km2 <- as.numeric(st_area(df_Shanghai_albers)) / 1e6
df_Shanghai_albers$pop_den <- df_Shanghai_albers$age_total/df_Shanghai_albers$area_km2
skewness(df_Shanghai_albers$pop_den, na.rm = TRUE)
kurtosis(df_Shanghai_albers$pop_den, na.rm = TRUE)
df_Shanghai_albers$lnpop <- log(df_Shanghai_albers$pop_den)
skewness(df_Shanghai_albers$lnpop, na.rm = TRUE)
kurtosis(df_Shanghai_albers$lnpop, na.rm = TRUE)
#st_write(df_Shanghai_albers, "H:/INLA/Stack/Stack_Shanghai/df_Shanghai_albers.shp",
#         driver = "ESRI Shapefile", layer_options = "ENCODING=UTF-8")

library(ggplot2)
ggplot(df_Shanghai_albers) +
  geom_sf(aes(fill = pop_den)) +
  scale_fill_viridis_c()

# 读取渔网
#setwd("H:/INLA/Stack/Stack_Shanghai")
#dsn <- "Stack_Shanghai.gdb"
#lay <- "fishnet40m_point"
#system.time(fishnet40m_point <- st_read(dsn = dsn, layer = lay))
#fishnet40m_point <- st_zm(fishnet40m_point,drop = T,what = "ZM")
#st_crs(fishnet40m_point)

# 读取边界
setwd("H:/INLA/Stack/Stack_Shanghai")
dsn <- "Stack_Shanghai.gdb"
lay <- "Shanghai_bdry_albers"
system.time(Shanghai_bdry_albers <- st_read(dsn = dsn, layer = lay))
Shanghai_bdry_albers <- st_zm(Shanghai_bdry_albers,drop = T,what = "ZM")
st_crs(Shanghai_bdry_albers)

library(terra)
# 读取raster
building <- "H:/INLA/building/building_agg40m.tif"
building <- rast(building)
print(building)
Shanghai_bdry_albers <- vect(Shanghai_bdry_albers)
building <- mask(building, Shanghai_bdry_albers)
print(building)
plot(building)
ntl_rgb <- "H:/INLA/raster_clip/RGB_merge_shanghai_proj.tif"
ntl_rgb <- rast(ntl_rgb)
ntl_rgb <- crop(ntl_rgb, ext(Shanghai_bdry_albers))
ntl_rgb <- mask(ntl_rgb, Shanghai_bdry_albers)
print(ntl_rgb)
plot(ntl_rgb)

# 提取 R、G、B 波段
red_band <- ntl_rgb[[1]]  # 第一波段（Red）
green_band <- ntl_rgb[[2]]  # 第二波段（Green）
blue_band <- ntl_rgb[[3]]  # 第三波段（Blue）
# 辐射定标
red_band <- red_band*0.0000102744+0.0000099253
green_band <- green_band*0.0000041779+0.0000060840
blue_band <- blue_band*0.0000070119+0.0000136754
plot(blue_band)

red_band <- resample(red_band, building, method = "bilinear")
green_band <- resample(green_band, building, method = "bilinear")
blue_band <- resample(blue_band, building, method = "bilinear")


red_band <- mask(red_band, building)  # 使用边界掩膜
green_band <- mask(green_band, building)
blue_band <- mask(blue_band, building)

plot(blue_band)
# 计算 Brightness
brightness <- 0.2989 * red_band + 0.5870 * green_band + 0.1141 * blue_band
plot(brightness)
# 计算 NDIGR (Green-Red)/(Green+Red)
ndigr <- (green_band - red_band) / (green_band + red_band)

# 计算 NDIBG (Blue-Green)/(Blue+Green)
ndibg <- (blue_band - green_band) / (blue_band + green_band)

ndirb <- (red_band - blue_band) / (red_band + blue_band)

plot(ndirb)
# is zero denominator
den_igr <- green_band + red_band
den_ibg <- blue_band + green_band
den_irb <- red_band + blue_band
global(den_igr == 0, "sum", na.rm = TRUE)
global(den_ibg == 0, "sum", na.rm = TRUE)
global(den_irb == 0, "sum", na.rm = TRUE)

# 计算 SONDI (NDIGR + NDIBG)
#sondi <- ndigr + ndibg
TIS <- "H:/INLA/raster_clip/TIS_shanghai_proj.tif"
TIS <- rast(TIS)
TIS <- mask(TIS, Shanghai_bdry_albers)
print(TIS)
plot(TIS)

# 提取3 波段
TIS_band1 <- TIS[[1]]  # 第一波段
TIS_band2 <- TIS[[2]]  # 第二波段
TIS_band3 <- TIS[[3]]  # 第三波段

# 辐射定标
TIS_band1 <- TIS_band1*0.003947+0.167126
TIS_band2 <- TIS_band2*0.003946+0.124622
TIS_band3 <- TIS_band3*0.005329+0.222530

TIS_band1 <- resample(TIS_band1, building, method = "bilinear")
TIS_band2 <- resample(TIS_band2, building, method = "bilinear")
TIS_band3 <- resample(TIS_band3, building, method = "bilinear")

TIS_band1 <- mask(TIS_band1, building)  # 使用边界掩膜
TIS_band2 <- mask(TIS_band2, building)
TIS_band3 <- mask(TIS_band3, building)
plot(TIS_band1)

library(ppcor)
# land use for partial correlation index
LULC_level1 <- "H:/INLA/GAIA/EULUC_level1.tif"
LULC_level2 <- "H:/INLA/GAIA/EULUC_level2.tif"
LULC_level1 <- rast(LULC_level1)
print(LULC_level1)
LULC_level2 <- rast(LULC_level2)
print(LULC_level2)

LULC_level1 <- crop(LULC_level1, ext(building))

LULC_level1 <- resample(LULC_level1, building, method = "near")
LULC_level1 <- mask(LULC_level1, building)
plot(LULC_level1)

LULC_level2 <- crop(LULC_level2, ext(building))

LULC_level2 <- resample(LULC_level2, building, method = "near")
LULC_level2 <- mask(LULC_level2, building)

# 创建一个新的栅格，将 LULC_level2 为 NA 且 building 不为 NA 的区域标记为 "Others"
others1 <- is.na(LULC_level1) & !is.na(building)
others2 <- is.na(LULC_level2) & !is.na(building)
# 将 Others 分类赋值为一个新的类别（如 99 表示 Others）
LULC_level1[others1] <- 99
plot(LULC_level1)
LULC_level2[others2] <- 999
plot(LULC_level2)

extracted1 <- terra::extract(LULC_level1, df_Shanghai_albers, cells = FALSE, df = TRUE)
extracted1 <- extracted1[complete.cases(extracted1), ]

# 统计每个区域内的栅格数量和比例
result1 <- extracted1 %>%
  group_by(ID, EULUC_level1) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(ID) %>%
  mutate(proportion = count / sum(count))

library(tidyr)
# 转换为宽格式
result_wide1 <- result1 %>%
  pivot_wider(
    names_from = EULUC_level1,  # 基于土地利用类型生成新列
    values_from = c(count, proportion),  # 获取计数和占比
    names_glue = "lv1_{EULUC_level1}_{.value}"  # 定义新列名称
  )
result_wide1 <- result_wide1 %>%
  mutate_all(~replace(., is.na(.), 0))

extracted2 <- terra::extract(LULC_level2, df_Shanghai_albers, cells = FALSE, df = TRUE)
extracted2 <- extracted2[complete.cases(extracted2), ]

# 统计每个区域内的栅格数量和比例
result2 <- extracted2 %>%
  group_by(ID, EULUC_level2) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(ID) %>%
  mutate(proportion = count / sum(count))

# 转换为宽格式
result_wide2 <- result2 %>%
  pivot_wider(
    names_from = EULUC_level2,  # 基于土地利用类型生成新列
    values_from = c(count, proportion),  # 获取计数和占比
    names_glue = "lv2_{EULUC_level2}_{.value}"  # 定义新列名称
  )

result_wide2 <- result_wide2 %>%
  mutate_all(~replace(., is.na(.), 0))
pcortest <- cbind(df_Shanghai_albers[,132:136],result_wide1[,2:13],result_wide2[,2:25])
pcortest <- as.data.frame(pcortest)

# 根据偏相关系数计算lv1
#pcortotal1 <- pcortest[,2:8]
#pcortotal1 <- pcor(pcortotal1)
#print(pcortotal1)

# 只计算lv1就好，lv2太过复杂且不合理
#pcortotal2 <- cbind(pcortest$age_total,pcortest[,15:26])
#pcortotal2 <- pcor(pcortotal2)
#print(pcortotal2)

pcorpercent1 <- cbind(pcortest$lnpop,pcortest[,11:15])
pcorpercent1 <- pcor(pcorpercent1)
print(pcorpercent1)
pcorpercent1 <- cbind(pcortest$lnpop,pcortest[,12:16])
pcorpercent1 <- pcor(pcorpercent1)
print(pcorpercent1)
#pcorpercent2 <- cbind(pcortest$age_total,pcortest[,27:38])
#pcorpercent2 <- pcor(pcorpercent2)
#print(pcorpercent2)

weights <- data.frame(
  land_use_type = c(1, 2, 3, 4, 5, 99),  # 土地利用类型
  weight = c(0.8053113, 0.71897312, 0.4098747, 0.026517986, 0.342286718, -0.8053113)  # 对应权重
)

#weights <- data.frame(
#  land_use_type = c(1, 2, 3, 4, 5, 99),  # 土地利用类型
#  weight = c(0.51269794, 0.61167080, -0.25635526, -0.1683228, -0.10234842, -0.57369563)  # 对应权重
#)

#weights <- data.frame(
#  land_use_type = c(1, 2, 3, 4, 5, 99),  # 土地利用类型
#  weight = c(0.4658723, 0.43421035, -0.2291584, -0.07952953, -0.0007614426, -0.46587232)  # 对应权重
#)

# 检查栅格唯一值
unique_values <- unique(LULC_level1[])
print(unique_values)  # 检查栅格中的土地利用类型

# 创建分类矩阵
weight_matrix <- as.matrix(weights)

# 替换土地利用类型为权重
LULC_lv1_weighted <- classify(LULC_level1, weight_matrix)
plot(LULC_lv1_weighted)
brightness_weighted <- LULC_lv1_weighted*brightness
plot(brightness_weighted)
blue_band_weighted <- LULC_lv1_weighted*blue_band
green_band_weighted <- LULC_lv1_weighted*green_band
red_band_weighted <- LULC_lv1_weighted*red_band
# 读取NDVI
ndvi_SDG <- "H:/INLA/raster_clip/MII_shanghai_proj.tif"
ndvi_SDG <- rast(ndvi_SDG)
print(ndvi_SDG)
MII_red <- ndvi_SDG[[5]]
MII_NIR <- ndvi_SDG[[7]] 
ndvi_SDG <- (MII_red - MII_NIR)/(MII_red + MII_NIR)
#template <- rast(ncol=ncol(ndvi_landsat_albers), nrow=nrow(ndvi_landsat_albers), resolution=c(40, 40),
#                 ext=ext(ndvi_landsat_albers), crs=crs(ndvi_landsat_albers))

# 重采样：将NDVI图像重采样到40米分辨率并与NTL数据对齐
ndvi_SDG <- resample(ndvi_SDG, building, method = "bilinear")
ndvi_SDG <- mask(ndvi_SDG, building)

plot(ndvi_SDG)
print(ndvi_SDG)
vanui <- (1-ndvi_SDG)*brightness
vanui_weighted <- (1-ndvi_SDG)*brightness*LULC_lv1_weighted
building_weighted <- building*LULC_lv1_weighted

#ndigr_weighted <- (green_band_weighted - red_band_weighted) / (green_band_weighted + red_band_weighted)
#ndibg_weighted <- (blue_band_weighted - green_band_weighted) / (blue_band_weighted + green_band_weighted)
#ndirb_weighted <- (red_band_weighted - blue_band_weighted) / (red_band_weighted + blue_band_weighted)



raster_stack <- c(building,building_weighted, red_band, green_band, blue_band, brightness,ndigr, ndibg,ndirb, vanui,brightness_weighted,
                  blue_band_weighted, green_band_weighted, red_band_weighted, vanui_weighted)
names(raster_stack) <- c("building","building_weighted","red_band", "green_band", "blue_band", "brightness",
                         "ndigr", "ndibg","ndirb", "vanui", "brightness_weighted",
                         "blue_band_weighted", "green_band_weighted", "red_band_weighted","vanui_weighted")
extracted_values <- terra::extract(raster_stack, df_Shanghai_albers, fun = mean, na.rm = TRUE)

# 将提取的值添加到 sf 对象中
df_Shanghai_albers <- cbind(df_Shanghai_albers, extracted_values)
ggplot(df_Shanghai_albers) +
  geom_sf(aes(fill = vanui)) +
  scale_fill_viridis_c()

library(randomForest)


x <- as.data.frame(cbind(df_Shanghai_albers$building,df_Shanghai_albers$building_weighted, df_Shanghai_albers$red_band, df_Shanghai_albers$green_band, df_Shanghai_albers$blue_band, 
                         df_Shanghai_albers$ndigr, df_Shanghai_albers$ndibg,df_Shanghai_albers$ndirb, 
                         df_Shanghai_albers$vanui,df_Shanghai_albers$brightness_weighted,df_Shanghai_albers$brightness,
                         df_Shanghai_albers$blue_band_weighted, df_Shanghai_albers$green_band_weighted,df_Shanghai_albers$red_band_weighted,
                         df_Shanghai_albers$vanui_weighted))
colnames(x) <- c("building","building_weighted","red_band", "green_band", "blue_band", 
                 "ndigr", "ndibg", "ndirb", "vanui","brightness_weighted","brightness",
                 "blue_band_weighted", "green_band_weighted", "red_band_weighted","vanui_weighted")
x_data = x
y_data = df_Shanghai_albers$lnpop
#x_data

##  Now we are going tune our randomForest demographic regression: 
init_popfit = tuneRF(
  x=x_data, y=y_data, plot=TRUE, mtryStart=length(x_data)/3, 
  ntreeTry=length(y_data)/20, improve=0.0001, stepFactor=1.20, trace=TRUE,
  doBest=TRUE, nodesize=length(y_data)/1000, na.action=na.omit, 
  importance=TRUE, proximity=TRUE, sampsize=length(y_data), replace=TRUE
) 
init_popfit
varImpPlot(init_popfit, color = par("fg"))
##  Get list of covariates that have an importance score greater than 0:
importance_scores <- 
  randomForest::importance(init_popfit)[order(randomForest::importance(init_popfit)[,1], decreasing=TRUE),]
pos_importance <- rownames(importance_scores)[importance_scores[,1] > 0]
pos_importance

popfit <- randomForest(x_data,y_data,importance = T,mtry = init_popfit$mtry, ntree = init_popfit$ntree)
popfit
varImpPlot(popfit, color = par("fg"))

library(Boruta)
library(ranger)
#traindata <- mid_build_Urban_pt_adj_cls_jicheng_
df <- cbind(x_data,y_data)
df <- as.data.frame(df)
traindata <- df[, !(colnames(df) %in% c("geometry","cityclass"))]
str(traindata)
names(traindata) <- gsub("_", "", names(traindata))
summary(traindata)
traindata[traindata == ""] <- NA
set.seed(2)
boruta.train <- Boruta(ydata~., data = traindata, doTrace = 2)
print(boruta.train)
plot(boruta.train, ylab = '', xlab = "Importance scores",las =2, horizontal = T)
legend("bottomright", inset=.02,
       c("Important","Tentative","Unimportant","Shadow features"), fill=c("green", "yellow", "red", "blue"), horiz=F, cex=0.95)
#lz<-lapply(1:ncol(boruta.train$ImpHistory),function(i)boruta.train$ImpHistory[is.finite(boruta.train$ImpHistory[,i]),i])
#names(lz) <- colnames(boruta.train$ImpHistory)
#Labels <- sort(sapply(lz,median))
#axis(side = 1,las=2,labels = names(Labels),at = 1:ncol(boruta.train$ImpHistory), cex.axis = 0.7)
final.boruta <- TentativeRoughFix(boruta.train)
print(final.boruta)
getSelectedAttributes(final.boruta, withTentative = T)
boruta.df <- attStats(final.boruta)
class(boruta.df)

library(corrplot)

# 示例数据：假设 x_data 已存在，包含数值型变量
# 如果需要筛选数值型变量，请先转换或过滤
# x_data <- x_data[, sapply(x_data, is.numeric)]

# 计算相关系数矩阵
cor_matrix <- cor(x_data, use = "pairwise.complete.obs")

# 绘制相关性热图并添加相关性数值
corrplot(cor_matrix, method = "color", type = "upper",
         tl.col = "black", tl.srt = 45, 
         col = colorRampPalette(c("blue", "white", "red"))(200),
         addCoef.col = "black", # 在图上添加相关性数值，颜色为黑色
         number.cex = 0.7)      # 调整相关性数值的字体大小

library(car)

lmdata <- cbind(x_data,y_data)
lmdata <- as.data.frame(lmdata)
model <- lm(y_data ~ brightness+ndigr+ndibg, data = lmdata)
summary(model)
alias(model)
vif(model)

#library(h2o)
#h2o.init()
summary(df)
str(df)
library(e1071)
mean_values <- apply(df, 2, mean, na.rm = TRUE)
sd_values <- apply(df, 2, sd, na.rm = TRUE)
median_values <- apply(df, 2, median, na.rm = TRUE)
min_values <- apply(df, 2, min, na.rm = TRUE)
max_values <- apply(df, 2, max, na.rm = TRUE)

# Skewness 和 kurtosis
skewness_values <- apply(
  df, 2,
  function(x) e1071::skewness(x, na.rm = TRUE, type = 2)
)

# e1071::kurtosis() 返回 excess kurtosis，因此加 3 转为普通 kurtosis
kurtosis_values <- apply(
  df, 2,
  function(x) e1071::kurtosis(x, na.rm = TRUE, type = 2) + 3
)

# 汇总成一个数据框
descriptive_stats <- data.frame(
  Mean = mean_values,
  `Standard Deviation` = sd_values,
  Median = median_values,
  Minimum = min_values,
  Maximum = max_values,
  Skewness = skewness_values,
  Kurtosis = kurtosis_values
)

# 查看结果
print(descriptive_stats)


