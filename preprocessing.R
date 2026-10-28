library(sf)
library(dplyr)
library(moments)

# read township-level census
setwd("H:/INLA/Stack/Stack_Shanghai")
dsn <- "Stack_Shanghai.gdb"
lay <- "HSSES_final_2020_albers"
system.time(df_Shanghai_albers <- st_read(dsn = dsn, layer = lay))
df_Shanghai_albers <- st_zm(df_Shanghai_albers,drop = T,what = "ZM")
# check CRS
st_crs(df_Shanghai_albers)

# aggregated by age
df_Shanghai_albers <- df_Shanghai_albers %>%
  mutate(across(starts_with("age_"), as.numeric)) %>%
  mutate(age_total = rowSums(across(starts_with("age_"))))

df_Shanghai_albers$area_km2 <- as.numeric(st_area(df_Shanghai_albers)) / 1e6
df_Shanghai_albers$pop_den <- df_Shanghai_albers$age_total/df_Shanghai_albers$area_km2
df_Shanghai_albers$lnpop <- log(df_Shanghai_albers$pop_den)
skewness(df_Shanghai_albers$lnpop, na.rm = TRUE)
kurtosis(df_Shanghai_albers$lnpop, na.rm = TRUE)

# plot township level census total
library(ggplot2)
ggplot(df_Shanghai_albers) +
  geom_sf(aes(fill = pop_den)) +
  scale_fill_viridis_c()


# Read city-level boundary for clipping ALAN data
setwd("H:/INLA/Stack/Stack_Shanghai")
dsn <- "Stack_Shanghai.gdb"
lay <- "Shanghai_bdry_albers"
system.time(Shanghai_bdry_albers <- st_read(dsn = dsn, layer = lay))
Shanghai_bdry_albers <- st_zm(Shanghai_bdry_albers,drop = T,what = "ZM")
st_crs(Shanghai_bdry_albers)

library(terra)
# read raster
building <- "H:/INLA/building/building_agg40m.tif" # this is rasterized building footprint data
building <- rast(building)
print(building)
Shanghai_bdry_albers <- vect(Shanghai_bdry_albers)
building <- mask(building, Shanghai_bdry_albers) # clip building footprint data by city boundary
print(building)
plot(building)
# read multi-spectral ALAN data and plot
ntl_rgb <- "H:/INLA/raster_clip/RGB_merge_shanghai_proj.tif" # this is SDGSAT-1 nighttime light data, which are subject to third-party licensing restrictions and available at: https://www.sdgsat.ac.cn.
ntl_rgb <- rast(ntl_rgb)
ntl_rgb <- crop(ntl_rgb, ext(Shanghai_bdry_albers))
ntl_rgb <- mask(ntl_rgb, Shanghai_bdry_albers)
print(ntl_rgb)
plot(ntl_rgb)

# get R, G, B bands
red_band <- ntl_rgb[[1]]  # Red
green_band <- ntl_rgb[[2]]  # Green
blue_band <- ntl_rgb[[3]]  # Blue
# radiometric calibration
red_band <- red_band*0.0000102744+0.0000099253
green_band <- green_band*0.0000041779+0.0000060840
blue_band <- blue_band*0.0000070119+0.0000136754
plot(blue_band)
# RGB to building footprint areas
red_band <- resample(red_band, building, method = "bilinear")
green_band <- resample(green_band, building, method = "bilinear")
blue_band <- resample(blue_band, building, method = "bilinear")
red_band <- mask(red_band, building)
green_band <- mask(green_band, building)
blue_band <- mask(blue_band, building)

plot(blue_band)
# get greyscale Brightness
brightness <- 0.2989 * red_band + 0.5870 * green_band + 0.1141 * blue_band
plot(brightness)
# get NDIGR (Green-Red)/(Green+Red)
ndigr <- (green_band - red_band) / (green_band + red_band)

# get NDIBG (Blue-Green)/(Blue+Green)
ndibg <- (blue_band - green_band) / (blue_band + green_band)

# get NDIRB
ndirb <- (red_band - blue_band) / (red_band + blue_band)
plot(ndirb)

# find if zero denominator exist
den_igr <- green_band + red_band
den_ibg <- blue_band + green_band
den_irb <- red_band + blue_band
global(den_igr == 0, "sum", na.rm = TRUE)
global(den_ibg == 0, "sum", na.rm = TRUE)
global(den_irb == 0, "sum", na.rm = TRUE)

library(ppcor)
# land use for partial correlation coefficient
LULC_level1 <- "H:/INLA/GAIA/EULUC_level1.tif"
LULC_level1 <- rast(LULC_level1)
print(LULC_level1)

LULC_level1 <- crop(LULC_level1, ext(building))

LULC_level1 <- resample(LULC_level1, building, method = "near")
LULC_level1 <- mask(LULC_level1, building)
plot(LULC_level1)


# make NA values in LULC_level2 as "Others"
others1 <- is.na(LULC_level1) & !is.na(building)

# assign "Others" as a type
LULC_level1[others1] <- 99
plot(LULC_level1)


extracted1 <- terra::extract(LULC_level1, df_Shanghai_albers, cells = FALSE, df = TRUE)
extracted1 <- extracted1[complete.cases(extracted1), ]

# calculate Land use percentages in each township
result1 <- extracted1 %>%
  group_by(ID, EULUC_level1) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(ID) %>%
  mutate(proportion = count / sum(count))

library(tidyr)
result_wide1 <- result1 %>%
  pivot_wider(
    names_from = EULUC_level1,
    values_from = c(count, proportion),
    names_glue = "lv1_{EULUC_level1}_{.value}" 
  )
result_wide1 <- result_wide1 %>%
  mutate_all(~replace(., is.na(.), 0))


pcortest <- cbind(df_Shanghai_albers[,132:136],result_wide1[,2:13])
pcortest <- as.data.frame(pcortest)
pcorpercent1 <- cbind(pcortest$lnpop,pcortest[,11:15])
pcorpercent1 <- pcor(pcorpercent1)
print(pcorpercent1)

weights <- data.frame(
  land_use_type = c(1, 2, 3, 4, 5, 99),  # Land use types
  weight = c(0.8053113, 0.71897312, 0.4098747, 0.026517986, 0.342286718, -0.8053113)  # corresponding weights from pcorpercent1
)


# check unique
unique_values <- unique(LULC_level1[])
print(unique_values)

# create matrix for weighting ALAN
weight_matrix <- as.matrix(weights)
LULC_lv1_weighted <- classify(LULC_level1, weight_matrix)
plot(LULC_lv1_weighted)
# weighting ALAN by partial correlation coefficient
brightness_weighted <- LULC_lv1_weighted*brightness
plot(brightness_weighted)
blue_band_weighted <- LULC_lv1_weighted*blue_band
green_band_weighted <- LULC_lv1_weighted*green_band
red_band_weighted <- LULC_lv1_weighted*red_band

# get NDVI data from SDGSAT-1 daytime data
ndvi_SDG <- "H:/INLA/raster_clip/MII_shanghai_proj.tif" # this is SDGSAT-1 daytime data, which are subject to third-party licensing restrictions and available at: https://www.sdgsat.ac.cn.
ndvi_SDG <- rast(ndvi_SDG)
print(ndvi_SDG)
MII_red <- ndvi_SDG[[5]] # red band
MII_NIR <- ndvi_SDG[[7]] # near-infrared band
ndvi_SDG <- (MII_red - MII_NIR)/(MII_red + MII_NIR)

# resample to 40m resolution to align with ALAN data
ndvi_SDG <- resample(ndvi_SDG, building, method = "bilinear")
ndvi_SDG <- mask(ndvi_SDG, building)

# get VANUI and land-use weighted VANUI
vanui <- (1-ndvi_SDG)*brightness
vanui_weighted <- (1-ndvi_SDG)*brightness*LULC_lv1_weighted
# get land-use weighted building height
building_weighted <- building*LULC_lv1_weighted
# get all variable in a stack for later processing
raster_stack <- c(building,building_weighted, red_band, green_band, blue_band, brightness,ndigr, ndibg,ndirb, vanui,brightness_weighted,
                  blue_band_weighted, green_band_weighted, red_band_weighted, vanui_weighted)
names(raster_stack) <- c("building","building_weighted","red_band", "green_band", "blue_band", "brightness",
                         "ndigr", "ndibg","ndirb", "vanui", "brightness_weighted",
                         "blue_band_weighted", "green_band_weighted", "red_band_weighted","vanui_weighted")
# extract variable values to the township level
extracted_values <- terra::extract(raster_stack, df_Shanghai_albers, fun = mean, na.rm = TRUE)

# to sf 
df_Shanghai_albers <- cbind(df_Shanghai_albers, extracted_values)
ggplot(df_Shanghai_albers) +
  geom_sf(aes(fill = vanui)) +
  scale_fill_viridis_c()


# get x and y data
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

# run Boruta variable selection
library(Boruta)
library(ranger)
df <- cbind(x_data,y_data)
df <- as.data.frame(df)
traindata <- df[, !(colnames(df) %in% c("geometry","cityclass"))]
str(traindata)
names(traindata) <- gsub("_", "", names(traindata))
summary(traindata)
traindata[traindata == ""] <- NA
# set a random seed for reproducibility
set.seed(2)
boruta.train <- Boruta(ydata~., data = traindata, doTrace = 2)
print(boruta.train)
plot(boruta.train, ylab = '', xlab = "Importance scores",las =2, horizontal = T)
legend("bottomright", inset=.02,
       c("Important","Tentative","Unimportant","Shadow features"), fill=c("green", "yellow", "red", "blue"), horiz=F, cex=0.95)
final.boruta <- TentativeRoughFix(boruta.train)
print(final.boruta)
getSelectedAttributes(final.boruta, withTentative = F)
boruta.df <- attStats(final.boruta)
class(boruta.df)

# plot a heatmap
library(corrplot)
cor_matrix <- cor(x_data, use = "pairwise.complete.obs")
corrplot(cor_matrix, method = "color", type = "upper",
         tl.col = "black", tl.srt = 45, 
         col = colorRampPalette(c("blue", "white", "red"))(200),
         addCoef.col = "black", 
         number.cex = 0.7) 

library(car)

# get collinearity
lmdata <- cbind(x_data,y_data)
lmdata <- as.data.frame(lmdata)
model <- lm(y_data ~ brightness+ndigr+ndibg, data = lmdata)
summary(model)
alias(model)
vif(model)

summary(df)
str(df)
library(e1071)
mean_values <- apply(df, 2, mean, na.rm = TRUE)
sd_values <- apply(df, 2, sd, na.rm = TRUE)
median_values <- apply(df, 2, median, na.rm = TRUE)
min_values <- apply(df, 2, min, na.rm = TRUE)
max_values <- apply(df, 2, max, na.rm = TRUE)

# Skewness and kurtosis
skewness_values <- apply(
  df, 2,
  function(x) e1071::skewness(x, na.rm = TRUE, type = 2)
)

# e1071::kurtosis() return excess kurtosis, plus 3 to normal kurtosis
kurtosis_values <- apply(
  df, 2,
  function(x) e1071::kurtosis(x, na.rm = TRUE, type = 2) + 3
)

# get descriptive statistics
descriptive_stats <- data.frame(
  Mean = mean_values,
  `Standard Deviation` = sd_values,
  Median = median_values,
  Minimum = min_values,
  Maximum = max_values,
  Skewness = skewness_values,
  Kurtosis = kurtosis_values
)

print(descriptive_stats)


