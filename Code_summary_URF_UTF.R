---
  title: "final models"
output:
  html_document: default
pdf_document: default
---
  
# code for each model mentioned in the manuscript (Table1 and Table2)
  ## Loading packages
library(Matrix);library(sp);library(raster);library(lme4);library(nlme);library(r2glmm);library(MuMIn)
library(raster)

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


  ##multiple regression model models
  ##linear URF
urf2 <- lm( ht20 ~ mat_s + mat_p + I(mat_s^2) + I(mat_p^2) + lahm_s + I(lahm_s^2) + lahm_p +
              I(lahm_p^2) + mat_s*mat_p + lahm_s*lahm_p, data=dat2)

  ##linear mixed URF
urf2_mix <- lmer( ht20 ~ mat_s + mat_p + I(mat_s^2) + I(mat_p^2) + lahm_s + I(lahm_s^2) + lahm_p +
                    I(lahm_p^2) + mat_s*mat_p + lahm_s*lahm_p + (lahm_p + mat_s|site) + (1|prov), REML = FALSE,
                  data=dat2)

  ##final URF
urf_f <- lm ( ht20 ~ mat_s + I(mat_s^2) + mat_p + I(mat_p^2) + lahm_s + I(lahm_s^2) + lahm_p +
                I(lahm_p^2) + mat_s*mat_p + I(mat_s^2)*mat_p + I(mat_s^2)*I(mat_p^2),data=dat2)

  ##linear UTF
utf2 <- lm( ht20 ~ mat_d + mat_p + I(mat_d^2) + I(mat_p^2) + mat_d*mat_p + lahm_d + I(lahm_d^2) +
              lahm_p + I(lahm_p^2) + lahm_d*lahm_p, data=dat2)

  ##linear mixed UTF
utf2_mix <- lmer( ht20 ~ mat_d + mat_p + I(mat_d^2) + I(mat_p^2) + mat_d*mat_p + lahm_d +
                    I(lahm_d^2) + lahm_p + I(lahm_p^2) + lahm_d*lahm_p + (lahm_p + mat_d|site) + (1|prov), REML =
                    FALSE, data=dat2)


##Figure 1 in manuscript used urf and utf mentioned above and was generated in Excel.

#code for random sampling from full data and generating Figure 2 in manuscript. Here using urf and urf_mix with different provenance sizes for demonstration.
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

#code for generating maps used in the manuscript (Figure 3 and 4)
## Raster stack of Western North America or world climate data (here using world data for demonstration)
## Function Predict_map() for predicting spatial distribution of fundamental niche and productivity
## period: climate period. Can be "current" or "future".
## type: type of seed source. Can be "local" or "optimal". 
Predict_map <- function(period, type) {
  if(period=='current'){
    stk_w <- stack(mat_s_current,mat_p_current,lahm_s_current,lahm_p_current);stk
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
    y = -56.20 + 2.140*mat_s - 0.2832*mat_s^2 + 0.223*mat_p - 0.0481*mat_p^2 + 39.70*lahm_s + -7.068*lahm_s^2 + 2.698*lahm_p - 0.3887*lahm_p^2 - 0.09199*mat_s*mat_p + 0.03632*(mat_s^2)*mat_p - 0.001548*(mat_s^2)*(mat_p^2)
    lahm_p = 3.4705;
    mat_p= (0.2230 - 0.09199*mat_s + 0.03632*mat_s^2)/(0.0962 + 0.003096*(mat_s^2))
    y[y<3]=0;plot(y);
  }
}

p <- Predict_map("current","local"); writeRaster(p, 'urf3_w.tif')
##.tif files of prediction were saved. Maps used in manuscript were generated in ArcGIS pro.
## difference maps were generated by using future maps (either with local or optimal seed source) minus current ditribution with local seed source.


#code for cross-validation R^2 amd RMSE in the manuscript (Table 1 and 2)
## using simple regression URF here for demonstration

cv <- data.frame(matrix(nrow = 10, ncol = 3))
## Split the data into training and test set
dat2$id <- c(1:4583)
##rdSample(site.list,n=15)
##set.seed(123)
for (i in c(1:10)){
  training.samples <- dat2$id %>%
    createDataPartition(p = 0.75, list = FALSE)
  train.data  <- dat2[training.samples, ]
  test.data <- dat2[-training.samples, ]
  ## training model
  urf <- lm(ht20 ~ mat_s + I(mat_s^2) + mat_p + I(mat_p^2) + mat_s*mat_p, data=train.data)
  ## Make predictions and compute the R2, RMSE and MAE
  predictions <- predict(urf, test.data)
  cv[i,1]<- R2(predictions, test.data$y)
  cv[i,2]<- RMSE(predictions, test.data$y)
  cv[i,3]<- MAE(predictions, test.data$y)
  i=i+1
  mean(cv$X1);mean(cv$X2)
}
