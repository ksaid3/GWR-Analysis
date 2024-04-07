easypackages::packages ("sf", "GWmodel", "sp", "spdep","spatialreg" ,"spatialsample", "spgwr", "tmap", "mapview", "mapview", "lmtest",  "dplyr", "ggplot2", "caret", "tidymodels","spatialsample","mgwrsar","blockCV", "mice","RColorBrewer")

cali <- st_read("C:/Users/Kalee/Desktop/SSML/Californiahealth/Calicountydata.gpkg")
cali



columns_to_convert <- c("Hispanic__", "White____", "African_Am", "Native_Ame", "Asian_Amer", "Other____")
# Convert columns to NUM
cali[columns_to_convert] <- lapply(cali[columns_to_convert], function(x) as.numeric(as.character(x)))
cali[is.na(cali)] <- 0

cali <- cali %>%
  mutate(Non_White = rowSums(across(c("Hispanic__", "African_Am", "Native_Ame", "Asian_Amer", "Other____")), na.rm = TRUE))

cali
#Select only county of california
selected_counties <- subset(cali, County %in% c("Los Angeles"))

selected_counties
#selected variables used for analysis
selected_vars <- c( "Asthma", "Ozone", "PM2_5", "Diesel_PM", "DrinkingWa", "Pesticides", "Traffic", "Poverty", "Education", "Children__", "Elderly___","Non_White")


selected_counties <- subset(selected_counties, select = selected_vars)
selected_counties


#add longitude and latitude
selected_counties$centroid <- st_centroid(selected_counties$geom)
selected_counties$latitude <- st_coordinates(selected_counties$centroid,crs = 4326)[, "Y"]
selected_counties$longitude <- st_coordinates(selected_counties$centroid,crs = 4326)[, "X"]
selected_counties <- selected_counties %>% select(-centroid)

county_data<-selected_counties



#Finding missing data
non_spatial_cols <- setdiff(names(county_data), c("geom"))
missingdata <- rep(FALSE, nrow(county_data))


for (i in 1:nrow(county_data)) {
  # Check if any of the non-spatial columns in this row have -999
  missingdata[i] <- any(unlist(county_data[i, non_spatial_cols]) == -999)
}

# subset the data that has missing values
country_data <- county_data[missingdata, ]
country_data


#judging by values, missing data is not random and occurs more in specific columns

# change -999 to na value instead
county_data <- county_data %>%
  mutate(across(all_of(non_spatial_cols), ~na_if(., -999)))


# apply Multiple Imputation

geom_data <- county_data$geom
county_data_no_geom <- st_drop_geometry(county_data)
imputed_data <- mice(county_data_no_geom, m = 5, method = 'pmm')

completed_data <- complete(imputed_data, 1)
county_data_imputed <- st_sf(completed_data, geometry = geom_data)
summary(county_data_imputed)

county_data<-county_data_imputed

#Density plot for Asthma
Dens_Asthma <- ggplot(data = county_data) +
  geom_density(alpha=0.8, colour="black", fill="lightblue", aes(x = Asthma)) +
  theme_classic()
Dens_Asthma

#Asthma and air pollutants visualization
tm_shape(county_data) +
  tm_polygons(c("Asthma", "PM2_5")) + 
  tm_layout(legend.position = c("left", "bottom"), 
            legend.text.size = 0.4, legend.title.size = 0.5)


tm_shape(county_data) +
  tm_polygons(c("Ozone","Diesel_PM")) + 
  tm_layout(legend.position = c("left", "bottom"), 
            legend.text.size = 0.4, legend.title.size = 0.5)

#####running GWR

#Non-spatial CV

# new equation based on previous regression and its significant variables
formula <- Asthma ~ Ozone+Diesel_PM +Pesticides + Poverty + Education + Children__+DrinkingWa+Elderly___+Traffic+PM2_5+ Non_White


# Generate spatial coordinates
coords <- cbind(county_data$latitude, county_data$longitude)

#drop geometry and lat/lon
geometry_data <- st_geometry(county_data)
non_geometry_data <- st_drop_geometry(county_data)
df <- subset(non_geometry_data, select = -c(longitude, latitude))
df
#prepare format for cv
Y <- county_data$Asthma
X <- as.matrix(df)
# Add a column of ones to X for the intercept
X <- cbind(1, X)

# Define a range of bandwidths to test
bandwidths <- seq(from = 10, to =20, by = 1)
# vector to store CV scores
cv_scores <- numeric(length(bandwidths))
#create distance matrix
dMat <- gw.dist(dp.locat = coords)

# Loop over bandwidths to compute CV scores(sse)
for (i in seq_along(bandwidths)) {
  bw <- bandwidths[i]
  cv_scores[i] <- gwr.cv(bw = bw,X=X, Y=Y, kernel= "gaussian", dp.locat =coords, adaptive= TRUE,verbose=TRUE,dMat = dMat)
}

# Identify the optimal bandwidth
optimal_bw_gauss <- bandwidths[which.min(cv_scores)]
optimal_bw_gauss


#CV for biqsuare kernel
for (i in seq_along(bandwidths)) {
  bw <- bandwidths[i]
  cv_scores[i] <- gwr.cv(bw = bw,X=X, Y=Y, kernel= "bisquare", dp.locat =coords, adaptive= TRUE,verbose=TRUE,dMat = dMat)
}
optimal_bw_bi <- bandwidths[which.min(cv_scores)]
optimal_bw_bi


#CV for exponential kernel
for (i in seq_along(bandwidths)) {
  bw <- bandwidths[i]
  cv_scores[i] <- gwr.cv(bw = bw,X=X, Y=Y, kernel= "exponential", dp.locat =coords, adaptive= TRUE,verbose=TRUE,dMat = dMat)
}
optimal_bw_ex <- bandwidths[which.min(cv_scores)]
optimal_bw_ex


county_data_sp<-as_Spatial(county_data)

#Run GWR
gwr <- gwr.basic(formula=formula, #the equation
                 adaptive = TRUE,
                 kernel="gaussian", #indicate the Kernel again
                 bw =20, #give the optimal bandwidth we found in the last stage
                 data=county_data_sp) 

#print the model result
gwr


####spatial cv


#create spatial cross validation folds
sb1 <- cv_spatial(x = county_data,
                  k = 10, # number of folds
                  size = 35000, # size of the blocks in metres
                  selection = "random") # random blocks-to-fold



###Looping through bandwith search grid and the spatial cv

set.seed(999)

bandwidths <- seq(from = 10, to = 110, by = 10)

# list to store results for each bandwidth
overall_results <- list()
mean_results <- data.frame(bandwidth = numeric(), mean_rmse = numeric())

# loop for bandwidths
for (bw in bandwidths) {
  # results for each bandwidth
  results <- data.frame()
  
  # loop for folds
  for (fold in 1:length(sb1$folds_list)) {
    train_indices <- sb1$folds_list[[fold]][[1]]
    test_indices <- sb1$folds_list[[fold]][[2]]
    
    # created training and testing indices applied on data 
    train_data <- as.data.frame(county_data_sp[train_indices, ])
    test_data <- as.data.frame(county_data_sp[test_indices, ])
    
    #coordinates extracted
    coords_train <- cbind(train_data$latitude, train_data$longitude)
    coords_test <- cbind(test_data$latitude, test_data$longitude)
    
    # fit GWR model on training data using each bandwidth
    mgwrsar_model <- MGWRSAR(formula = formula, data = train_data,
                             coords = coords_train, kernels = c('gauss'), 
                             H = bw, Model = 'GWR', control = list(adaptive = TRUE))
    
    # predict on test data
    predictions <- predict_mgwrsar(model = mgwrsar_model, newdata = test_data,
                                   newdata_coords = coords_test, method_pred = 'kernel_model')
    
    
    # Calculate RMSE
    rmse <- sqrt(mean((test_data$Asthma - predictions)^2))
    
    # record results for each fold
    results <- rbind(results, data.frame(fold = fold, rmse = rmse))
  }
  # Calculate mean RMSE for each bandwidth
  mean_rmse <- mean(results$rmse)
  
  # Store mean RMSE in the mean_results data frame
  mean_results <- rbind(mean_results, data.frame(bandwidth = bw, mean_rmse = mean_rmse))
  
  # Store the fold-wise results for the current bandwidth in the list
  overall_results[[as.character(bw)]] <- results
  
}

mean_results


#Running final model with bandwith that was found with spatial cv
gwrmodel <- gwr.basic(formula=formula, #the equation
                      adaptive = TRUE,
                      kernel="gaussian", #indicate the Kernel again
                      bw = 80, #give the optimal bandwidth we found in the last stage
                      data=county_data_sp) 
gwrmodel


#### spatial cv result visualization 
gwr_sf = st_as_sf(gwrmodel$SDF)

#produce the maps of the most significant variables


agwrM1 <- tm_shape(gwr_sf) +
  tm_fill("Diesel_PM") + 
  tm_borders() +
  tm_layout(legend.position = c(0.03, 0.01), legend.width = 0.5)

# Formatting agwrM2
agwrM2 <- tm_shape(gwr_sf) +
  tm_fill("Ozone") + 
  tm_borders() +
  tm_layout(legend.position = c(0.03, 0.01), legend.width = 0.5)

# Formatting agwrM3
agwrM3 <- tm_shape(gwr_sf) +
  tm_fill("Traffic") + 
  tm_borders() +
  tm_layout(legend.position = c(0.03, 0.01), legend.width = 0.5)

# Formatting agwrM4
agwrM4 <- tm_shape(gwr_sf) +
  tm_fill("Pesticides") + 
  tm_borders() +
  tm_layout(legend.position = c(0.03, 0.01), legend.width = 0.5)

# Formatting agwrM5
agwrM5 <- tm_shape(gwr_sf) +
  tm_fill("DrinkingWa") + 
  tm_borders() +
  tm_layout(legend.position = c(0.03, 0.01), legend.width = 0.5)

# Formatting agwrM6
agwrM6 <- tm_shape(gwr_sf) +
  tm_fill("Elderly___") + 
  tm_borders() +
  tm_layout(legend.position = c(0.03, 0.01), legend.width = 0.5)

# Formatting agwrM7
agwrM7 <- tm_shape(gwr_sf) +
  tm_fill("Children__") + 
  tm_borders() +
  tm_layout(legend.position = c(0.03, 0.01), legend.width = 0.5)

agwrM8 <- tm_shape(gwr_sf) +
  tm_fill("Education") + 
  tm_borders() +
  tm_layout(legend.position = c(0.03, 0.01), legend.width = 0.5)

agwrM9 <- tm_shape(gwr_sf) +
  tm_fill("Poverty") + 
  tm_borders() +
  tm_layout(legend.position = c(0.03, 0.01), legend.width = 0.5)

agwrM10 <- tm_shape(gwr_sf) +
  tm_fill("PM2_5") + 
  tm_borders() +
  tm_layout(legend.position = c(0.03, 0.01), legend.width = 0.5)

agwrM11 <- tm_shape(gwr_sf) +
  tm_fill("Non_White") + 
  tm_borders() +
  tm_layout(legend.position = c(0.03, 0.01), legend.width = 0.5)
# Formatting agwrRrsq
agwrRrsq <- tm_shape(gwr_sf) +
  tm_fill("Local_R2") + 
  tm_borders() +
  tm_layout(legend.position = c(0.03, 0.01), legend.width = 0.5)



tmap_arrange(agwrM1, agwrM2,agwrM3)
tmap_arrange(agwrM4,agwrM5,agwrM6)
tmap_arrange(agwrM7,agwrM8,agwrM9)
tmap_arrange(agwrM10,agwrM11,agwrRrsq)

agwrRrsq

#Visualizing significant polygons

visualize_variable <- function(variable_name) {
  # Data is always gwr_sf
  data <- gwr_sf
  
  # Create a new variable for t-values
  t_value_var <- paste0(variable_name, "_TV")
  
  # Select the t-value variable and drop geometry
  variable_viz <- data %>% 
    dplyr::select(all_of(t_value_var)) %>% 
    st_drop_geometry()
  
  # identify significant values
  variable_sig <- variable_viz < -1.96 | variable_viz > 1.96
  
  
  tmap_mode("view")
  
  
  tm_shape(data) +
    tm_fill(variable_name, midpoint = 0) + 
    tm_style("col_blind") +
    tm_basemap("OpenStreetMap") +
    tm_layout(legend.position = c("right", "top")) +
    # Add the layer for significant values
    tm_shape(data[variable_sig, ]) + 
    tm_borders()
}

#Statistically significant areas for each variable

visualize_variable("Diesel_PM")
visualize_variable("Ozone")
visualize_variable("Pesticides")
visualize_variable("Traffic")
visualize_variable("DrinkingWa")
visualize_variable("Elderly___")
visualize_variable("Children__")
visualize_variable("PM2_5")
visualize_variable("Education")
visualize_variable("Poverty")







####Testing spatial autocorrelation

#creating adjacency matrix, queen: Considers units as neighbors if they share a border or a point 
cali_nbq <- poly2nb(county_data, queen=TRUE) #Queen’s Contiguity neighborhood, if we put queen= FALSE, it will use Rook’s case contiguity
#Print the summary of the matrix
summary(cali_nbq)

#convert to list
cali_nbq_w<- nb2listw(cali_nbq, style="W", zero.policy = TRUE)


#Due to irregular distribution of polygons, monte carlo morans I instead of regular is use
mc_global <- moran.mc(county_data$Asthma, cali_nbq_w , 999, alternative="greater") 
mc_global
#The p value is significant:reject the null-hypothesis that there is no significant autocorrelations in the variable


#Look at local variations
moran_local<- localmoran(county_data$Asthma,cali_nbq_w )
summary(moran_local)


#Visualization of local moran I and significance
county_data$moran_local <- moran_local[,1] 
# extract p-values
county_data$moran_local_p <- moran_local[,5] 

#the local Moran's I with t-map+ map significancy of clusters
map_local <- tm_shape(county_data) + 
  tm_polygons(col= "moran_local", title= "Local Moran’s I", midpoint=0,
              palette = "RdYlBu", breaks= c(-10, -5, 0, 5, 10, 20)) 
map_local_p <- tm_shape(county_data) + 
  tm_polygons(col= "moran_local_p", title= "p-values",
              breaks= c(0, 0.01, 0.05, 1), palette = "Reds") 

tmap_arrange(map_local, map_local_p)


###Initial linear model

#OLS

#run the model
linearMod <- lm (formula, data = county_data) 

#get summary
summary(linearMod)

#test morans I
lm_ols <- moran.mc(linearMod$residuals, cali_nbq_w, 999, zero.policy= TRUE, alternative="greater")
lm_ols ###there is spatial autocorrelation in the errors


#test to see autocorrelation
test <- lm.LMtests(linearMod, listw= cali_nbq_w, test="all") 
tres <- t(sapply(test, function(x) c(x$statistic, x$parameter, x$p.value))) #create a table format from the results of lm.LMtests

colnames(tres) <- c("Statistic", "df", "p-value") #rename the columns

#print the table
tres 

#both RSerr and RSlag (robust versions) have low p-values, indicating strong evidence of spatial autocorrelation
##all the p values are below <0.05, indicating the errors of the model and independent variables have spatial autocorrelation.


#Therefore spatial model needed

#lag model
spa_lagmodel = lagsarlm (formula, data = county_data, listw= cali_nbq_w, zero.policy = TRUE)
summary(spa_lagmodel, Nagelkerke=T)
lag_moran <-moran.mc(spa_lagmodel$residuals, cali_nbq_w, 999, alternative="greater")
#not statistically significant, spatial autocorrelation 
#lag has accounted for spatial dependance
lag_moran 


#SEM model
spa_errmodel <- errorsarlm(formula, data = county_data, listw = cali_nbq_w, zero.policy = TRUE)
summary(spa_errmodel, Nagelkerke=T)


sem_moran <-moran.mc(spa_errmodel$residuals, cali_nbq_w, 999, alternative="greater")
#not statistically significant, spatial autocorrelation 
#sem has accounted for spatial dependance
sem_moran 




#visualizations

tmap_mode("plot")

map<-tm_shape(county_data) +
  tm_polygons(c("Asthma"), palette="GnBu", border.col="black", border.alpha=0.2) +
  tm_borders(alpha = 0) +
  tm_scale_bar(position = c("left", "BOTTOM"),width=0.2) +
  tm_compass(type = "arrow", position = c("center", "top"), size = 2) +
  tm_layout(
    legend.position = c("left", "center"),  
    inner.margins = c(0, 0.35, 0.30, 0.25),
    legend.text.size = 1.25,
    legend.title.size = 1.75,
    frame = FALSE,
    legend.frame = TRUE, 
    legend.frame.lwd = 0.3,
    legend.bg.color = "white", 
    legend.bg.alpha = 0.75
  )

tmap_save(map, "wider_map.png", width = 11, height = 8.5) 


gwr1 <- tm_shape(gwr_sf) +
  tm_polygons(c("DrinkingWa"), palette="GnBu", border.col="black", border.alpha=0.2,title="Drinking Water") + 
  tm_borders(alpha = 0) + # Set the alpha for borders
  tm_scale_bar(position = c("left", "bottom"), width=0.2) + 
  tm_compass(type = "arrow", position = c("right", "center"), size = 2) + 
  tm_layout(main.title = "Coefficient map: Drinking Water Contaminant Index",
            main.title.position = "center",
            main.title.size = 1,
    legend.position = c("left", "center"), 
    inner.margins = c(0, 0.35, 0.30, 0.25), 
    legend.text.size = 1, 
    legend.title.size = 1.5, 
    frame = FALSE, 
    legend.frame = TRUE, 
    legend.frame.lwd = 0.3, 
    legend.bg.color = "white", 
    legend.bg.alpha = 0.75 
  )

gwr1




gwr2 <- tm_shape(gwr_sf) +
  tm_polygons(c("Diesel_PM"), palette="GnBu", border.col="black", border.alpha=0.2,title="Diesel PM") +
  tm_borders(alpha = 0) + # Set the alpha for borders
  tm_scale_bar(position = c("left", "bottom"), width=0.2) + 
  tm_compass(type = "arrow", position = c("right", "center"), size = 2) + 
  tm_layout(main.title = "Coefficient map: Diesel PM emissions(kg/day)",
            main.title.position = "center",
            main.title.size = 1,
            legend.position = c("left", "center"), 
            inner.margins = c(0, 0.35, 0.30, 0.25), 
            legend.text.size = 1, 
            legend.title.size = 1.5, 
            frame = FALSE, 
            legend.frame = TRUE, 
            legend.frame.lwd = 0.3,
            legend.bg.color = "white", 
            legend.bg.alpha = 0.75 
  )

gwr2

gwr3 <- tm_shape(gwr_sf) +
  tm_polygons(c("Education"), palette="GnBu", border.col="black", border.alpha=0.2,title="Educational Attainment") + 
  tm_borders(alpha = 0) + 
  tm_scale_bar(position = c("left", "bottom"), width=0.2) + # Add a scale bar
  tm_compass(type = "arrow", position = c("right", "center"), size = 2) + # Add a compass
  tm_layout(main.title = "Coefficient map:Educational attainment",
            main.title.position = "center",
            main.title.size = 1,
    legend.position = c("left", "center"), 
    inner.margins = c(0, 0.35, 0.30, 0.25), 
    legend.text.size = 1, 
    legend.title.size = 1.5, 
    frame = FALSE, 
    legend.frame = TRUE, 
    legend.frame.lwd = 0.3, 
    legend.bg.color = "white", 
    legend.bg.alpha = 0.75 
  )
gwr3



gwr4 <- tm_shape(gwr_sf) +
  tm_polygons(c("Non_White"), palette="GnBu", border.col="black", border.alpha=0.2,title="Non white population") + 
  tm_borders(alpha = 0) + 
  tm_scale_bar(position = c("left", "bottom"), width=0.2) +
  tm_compass(type = "arrow", position = c("right", "center"), size = 2) + 
  tm_layout(main.title = "Coefficient map: Non white population",
            main.title.position = "center",
            main.title.size = 1,
    legend.position = c("left", "center"), 
    inner.margins = c(0, 0.35, 0.30, 0.25), 
    legend.text.size = 1, 
    legend.title.size = 1.5, 
    frame = FALSE, 
    legend.frame = TRUE, 
    legend.frame.lwd = 0.3, 
    legend.bg.color = "white", 
    legend.bg.alpha = 0.75 
  )
gwr4


breaks <- seq(0, 1, by = 0.1) 

# Set midpoint for diverging scale
midpoint <- 0

agwrsquare <- tm_shape(gwr_sf) +
  tm_polygons("Local_R2", palette="Spectral", border.col="black", border.alpha=0.2, breaks=breaks, midpoint=midpoint) +
  tm_shape(gwr_sf) + 
  tm_borders(alpha = 0) +
  tm_scale_bar(position = c("left", "bottom"), width=0.2) +
  tm_compass(type = "arrow", position = c("right", "top"), size = 2) +
  tm_layout(
    legend.position = c("left", "center"),
    inner.margins = c(0, 0.35, 0.30, 0.25),
    legend.text.size = 0.4,
    legend.title.size = 0.8,
    frame = FALSE,
    legend.frame = TRUE,
    legend.frame.lwd = 0.3,
    legend.bg.color = "white",
    legend.bg.alpha = 0.75
  )

# Display the map
agwrsquare

two_maps<-tmap_arrange(gwr1,gwr2)
two_maps
tmap_save(two_maps, "two_map.png", width = 13, height = 8.5) # width and height in inches


two_maps2<-tmap_arrange(gwr3,gwr4)
two_maps2
tmap_save(two_maps2, "two_map2.png", width = 13, height = 8.5) # width and height in inches
