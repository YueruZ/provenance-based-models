---
#Manuscript title: "Predicting the global fundamental climate niche of lodgepole pine based on provenance trials"
#Authors: Yueru Zhao and Tongli Wang
---
################################################## 
## code for each model mentioned in the manuscript (Table1 and Table2)
##################################################
## Loading packages
library(Matrix);library(raster);library(lme4)
## dat2: provenance data used for model training
## ht20: 20 year tree height
## mat_s: mean annual temperature for site
## mat_p: mean annual temperature for provenance
## mat_d: transfer distance for mean annual temerature
mat_d = mat_s - mat_p
## lahm_s: log-transformed annual heat-moisture index for site
## lahm_p: log-transformed annual heat-moisture index for peovenance
## lahm_d: transfer distance for log-transformed annual heat-moisture index
lahm_d = lahm_s - lahm_p

## simple regression models
##linear URF
urf <- lm(ht20 ~ mat_s + I(mat_s^2) + mat_p + I(mat_p^2) + mat_s*mat_p, data=dat2);

##linear mixed URF
urf_mix <- lmer(ht20 ~ mat_s + I(mat_s^2) + mat_p + I(mat_p^2) + mat_s*mat_p + (mat_s|site) +
                  (1|prov), REML = FALSE,data=dat2)

##linear UTF
utf <- lm(ht20 ~ mat_p + I(mat_p^2) + mat_d + I(mat_d^2) + mat_p*mat_d, data=dat2)

##linear mixe UTF
utf_mix <- lmer(ht20 ~ mat_d + mat_p + I(mat_d^2) + I(mat_p^2) + mat_d*mat_p + (mat_d|site) +
                  (1|prov), REML = FALSE, data = dat2)


## multiple regression model models
## linear URF, it is also the final URF
urf_2 <- lm(ht20~mat_s + I(mat_s^2) + mat_p + I(mat_p^2) + lahm_s + I(lahm_s^2) + lahm_p + 
             I(lahm_p^2) + mat_s*mat_p + I(mat_s^2)*mat_p + I(mat_s^2)*I(mat_p^2), data = dat2)

##linear mixed URF
urf2_mix <- lmer( ht20 ~ mat_s + mat_p + I(mat_s^2) + I(mat_p^2) + lahm_s + I(lahm_s^2) + lahm_p +
                    I(lahm_p^2) + mat_s*mat_p + lahm_s*lahm_p + (lahm_p + mat_s|site) + (1|prov), REML = FALSE,
                  data=dat2)

##linear UTF
utf2 <- lm( ht20 ~ mat_d + mat_p + I(mat_d^2) + I(mat_p^2) + mat_d*mat_p + lahm_d + I(lahm_d^2) +
              lahm_p + I(lahm_p^2) + lahm_d*lahm_p, data=dat2)

##linear mixed UTF
utf2_mix <- lmer( ht20 ~ mat_d + mat_p + I(mat_d^2) + I(mat_p^2) + mat_d*mat_p + lahm_d +
                    I(lahm_d^2) + lahm_p + I(lahm_p^2) + lahm_d*lahm_p + (lahm_p + mat_d|site) + (1|prov), REML =
                    FALSE, data=dat2)


##Figure 1 in manuscript used urf and utf mentioned above and was generated in Excel.


##################################################
##code for random sampling from full data and generating Figure 2 in manuscript. Here using urf and urf_mix with different provenance sizes for demonstration.
##################################################
prov_result <- data.frame(matrix(nrow = 50, ncol = 38))
colnames(prov_result) <- c("pduu5","pdl5","pduu10", "pdl10","pduu15", "pdl15", "pduu20", "pdl20","pduu25", "pdl25", "pduu30", "pdll30","pduu35", "pdll35", "pduu40", "pdll40","pduu45", "pdll45", "pduu50", "pdll50", "pduu55", "pdll55", "pduu60", "pdll60", "pduu65", "pdll65", "pduu70", "pdll70","pduu75", "pdll75", "pduu80", "pdll80","pduu85", "pdll85", "pduu90", "pdll90","pduu95", "pdll95","pduu100", "pdll100")
prov_result

i <- 0
repeat 
{ 
  sample <- data.frame(rdSample(prov.list,n=42));sample
  names(sample)[c(1)] <- c('prov')
  sp2 <- dtJoin(sample,dat2, by = 'prov', type = "left");sp2
  
  ##training both models with sampled data
  urf <- lm(ht20~mat_s+I(mat_s^2)+mat_p+I(mat_p^2)+mat_s*mat_p,data=sp2)
  urf_mix <- lmer(ht20 ~ mat_s + I(mat_s^2) + mat_p + I(mat_p^2) + mat_s*mat_p + (1|site) + (1|prov), REML = FALSE,data=sp2)
  
  obs <- c(dat2$ht20);
  ##testing both models with 
  purf <- predict(urf_mat, dat2, re.form=NA); 
  purf_mix <- predict(urf_mix, dat2, re.form=NA); 
  #writing prediction errors with different sample size into dataframe prov_result
  prov_result[i,1] <- sqrt(mean((obs - purf)^2))
  prov_result[i,2] <- sqrt(mean((obs - purf_mix)^2))
  
  ## incrementing the iteration variable 
  i = i + 1
  
  ## checking the stop condition 
  if (i == 51) 
  { 
    ## using break statement 
    ## to terminate the loop       
    break
  } 
} 
prov_result

##################################################
## code for generating maps used in the manuscript (Figure 3 and 4)
## Raster stack of Western North America or world climate data (here using world data for demonstration)
## Function Predict_map() for predicting spatial distribution of fundamental niche and productivity using climate data in raster format
## period: climate period. Can be "current" or "future".
## type: type of seed source. Can be "local" or "optimal". 
##################################################
Predict_map <- function(period, type) {
  if(period=='current'){
    stk_w <- stack(mat_s_current,mat_p_current,lahm_s_current,lahm_p_current);stk # stack containing raster layers of each climate variable
    names(stk_w) <- c('mat_s','mat_p', 'lahm_s', 'lahm_p');stk_w
  }
  if(period=='future'){
    stk_w <- stack(mat_s_future,mat_p_future,mat_p_future,lahm_p_future);stk
    names(stk_w) <- c('mat_s','mat_p', 'lahm_s', 'lahm_p');stk_w
  }
  if(type=='local'){
    p_w <- predict(stk_w,urf_f);p_w[p_w < 3] <- 0
    plot(p_w)
  }
  if(type=='optimal'){
    mat_p= (0.206168 - 0.078614*mat_s + 0.032213*mat_s^2)/(0.09582 + 0.002164*(mat_s^2));
    lahm_p = 3.1621
    y = -20.921887 + 2.025553*mat_s - 0.294167*mat_s^2 + 0.206168*mat_p - 0.047910*mat_p^2 + 
      16.004807*lahm_s  -3.039169*lahm_s^2 + 2.812786*lahm_p - 0.444765*lahm_p^2 - 
      0.078614*mat_s*mat_p + 0.032213*(mat_s^2)*mat_p - 0.001082*(mat_s^2)*(mat_p^2)
    y[y<3]=0;plot(y);
  }
}

p <- Predict_map("current","local"); writeRaster(p, 'fURF_w.tif')
##.tif files of prediction were saved. Maps used in manuscript were generated in ArcGIS pro.
## difference maps were generated by using future maps (either with local or optimal seed source) minus current ditribution with local seed source.


##################################################
## code for spatial block cross-validation R^2 amd RMSE in the manuscript (Table 1 and 2)
##################################################
#list of linear model
regression_lm <- list(
  c("ht20~mat_s+I(mat_s^2)+mat_p+I(mat_p^2)+mat_s*mat_p"),
  c("ht20 ~ mat_p + I(mat_p^2) + mat_d + I(mat_d^2) + mat_p*mat_d"),
  c("ht20~mat_s+I(mat_s^2)+mat_p+I(mat_p^2)+lahm_s+I(lahm_s^2)+lahm_p+I(lahm_p^2)+mat_s*mat_p+lahm_p*lahm_s+I(mat_s^2)*mat_p+I(mat_s^2)*I(mat_p^2)"), 
  c("ht20~mat_s+I(mat_s^2)+mat_p+I(mat_p^2)+lahm_s+I(lahm_s^2)+lahm_p+I(lahm_p^2)+mat_s*mat_p+I(mat_s^2)*mat_p+I(mat_s^2)*I(mat_p^2)"),
  c("ht20~mat_d+I(mat_d^2)+mat_p+I(mat_p^2)+lahm_d+I(lahm_d^2)+lahm_p+I(lahm_p^2)+mat_d*mat_p+lahm_p*lahm_d+I(mat_d^2)*mat_p+I(mat_d^2)*I(mat_p^2)"))

#list of linear mixed model
regression_lme <- list(
  c("ht20 ~ mat_s + mat_p + mat_s*mat_p + I(mat_s^2) + I(mat_p^2) + (mat_s|site) + (1|prov)"),
  c("ht20 ~ mat_d + mat_p + mat_d*mat_p + I(mat_d^2) + I(mat_p^2) + (mat_d|site) + (1|prov)"),
  c("ht20 ~ mat_s + mat_p + I(mat_s^2) + I(mat_p^2) + lahm_s + I(lahm_s^2) + lahm_p + I(lahm_p^2) + mat_s*mat_p + lahm_s*lahm_p + (lahm_p + mat_s|site) + (1|prov)"), 
  c("ht20 ~ mat_d + mat_p + I(mat_d^2) + I(mat_p^2) +  mat_d*mat_p + lahm_d + I(lahm_d^2) + lahm_p + I(lahm_p^2) + lahm_d*lahm_p + (lahm_p + mat_d|site) + (1|prov)"),
  c("ht20 ~ mat_s + mat_p + I(mat_s^2) + I(mat_p^2) + lahm_s + I(lahm_s^2) + lahm_p + I(lahm_p^2) + mat_s*mat_p + lahm_s*lahm_p + (lahm_p + mat_s|site) + I(mat_s^2)*mat_p + I(mat_s^2)*I(mat_p^2) + (1|prov)"))


#plotting test site data points on British Columbia map
mp <- raster("BC mapping/climatedata/Normal_1961_1990SY/mat.tif")
dat2$lon_s <- -dat2$lon_s
pt <- st_as_sf(dat2, coords = c("lon_s", "lat_s"), crs = crs(mp))
plot(mp)
points(dat2$lon_s, dat2$lat_s, pch = ".")

metrics <- data.frame(matrix(nrow = 5, ncol = 10));
metrics$model <- c("URF_sgl","UTF_sgl","URF_mul","fURF","UTF_mul"); metrics_rmse <- metrics
metrics_x <- data.frame(matrix(nrow = 5, ncol = 10)); metrics_x$model <- c("URFm_sgl","UTFm_sgl","URFm_mul","UTFm_mul","URFm_f"); metrics_x_rmse <- metrics_x
for (t in 1:10){
  sb <- spatialBlock(speciesData = pt,
                     species = "lpp",
                     rasterLayer = mp,
                     theRange = 300000, # size of the blocks
                     k = 5,
                     selection = "random",
                     iteration = 100, # find evenly dispersed folds
                     biomod2Format = TRUE,
                     xOffset = 0, # shift the blocks horizontally
                     yOffset = 0)
  sb$plots + geom_sf(data = pt, alpha = 0.5)
  dat2$fold <- sb$foldID;folds <- dat2$fold
  
  cv <- data.frame(matrix(nrow = 5, ncol = 3)); cvx <- cv
  
  for (i in 1:length(regression_lm)){
    for (k in c(1:5)){
      trainSet <- which(folds != k) # training set indices
      testSet <- which(folds == k) # testing set indices
      #build model
      lm <- lm(as.formula(regression_lm[[i]]), dat2[trainSet,])
      lme <- lmer(as.formula(regression_lme[[i]]), REML = FALSE, dat2[trainSet,])
      # Make predictions and compute the R2, RMSE and MAE
      predictions <- predict(lm, dat2[testSet,],re.form=NA)
      cv[k,1]<- R2(predictions, dat2[testSet,]$ht20)
      cv[k,2]<- RMSE(predictions, dat2[testSet,]$ht20)
      cv[k,3]<- MAE(predictions, dat2[testSet,]$ht20)
      predictions_x <- predict(lme, dat2[testSet,],re.form=NA)
      cvx[k,1]<- R2(predictions_x, dat2[testSet,]$ht20)
      cvx[k,2]<- RMSE(predictions_x, dat2[testSet,]$ht20)
      cvx[k,3]<- MAE(predictions_x, dat2[testSet,]$ht20)
    }
    metrics[i,t] <- mean(cv$X1); metrics_rmse[i,t] <- mean(cv$X2)
    metrics_x[i,t] <- mean(cvx$X1); metrics_x_rmse[i,t] <- mean(cvx$X2)
  }
}
data.frame(metrics$model,rowMeans(metrics[1:5,1:10]));
data.frame(metrics_x$model,rowMeans(metrics_x[1:5,1:10]))
data.frame(metrics_rmse$model,rowMeans(metrics_rmse[1:5,1:10]));
data.frame(metrics_x_rmse$model,rowMeans(metrics_x_rmse[1:5,1:10]))

##################################################
## Citations
##Robert J. Hijmans (2021). raster: Geographic Data Analysis
## and Modeling. R package version 3.4-13.
##  https://CRAN.R-project.org/package=raster

##Douglas Bates and Martin Maechler (2021). Matrix: Sparse
##  and Dense Matrix Classes and Methods. R package version
##  1.3-4. https://CRAN.R-project.org/package=Matrix

##Douglas Bates, Martin Maechler, Ben Bolker, Steve Walker
##  (2015). Fitting Linear Mixed-Effects Models Using lme4.
##  Journal of Statistical Software, 67(1), 1-48.
##  doi:10.18637/jss.v067.i01.
##################################################