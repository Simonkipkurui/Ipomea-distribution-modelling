library(biomod2)
library(raster)
library(terra)
library(usdm)
library(tidyverse)
library(sp)
library(ggplot2)
library(corrplot)
library(dplyr)
library(ggcorrplot)
library(data.table)

##loading the predictor variable==============================================================================
######Bio climatic variables
Bioclimatic<-list.files(path ="D:/Mr kerich/Bioclimatic_Variables" ,pattern=".tif$" ,full.names = TRUE)
Bioclimatic_variable<-stack(Bioclimatic)
plot(Bioclimatic_variable)
crs(Bioclimatic_variable)
######Soil variables  
Soil_bioclimatic<-stack(list.files(path="D:/Mr kerich/GEE_Exports/Gee_masked",pattern =".tif$",full.names = TRUE ))
plot(Soil_bioclimatic)
soil<-stack(list.files(path = "D:/Mr kerich/Soil/soil_projected",pattern = "\\.tif$",full.names = TRUE))
soil<-extend(soil, Bioclimatic_variable)
short_names<-c( "SBIO1", "SBIO10",
 "SBIO11", "SBIO2" ,                 
 "SBIO3","SBIO4" ,            
 "SBIO5","SBIO6",    
 "SBIO7","SBIO8", 
 "SBIO9")
names(Soil_bioclimatic)
names(Soil_bioclimatic)<-short_names
soil<-resample(soil,Bioclimatic_variable)
soil<-projectRaster(soil,Bioclimatic_variable,method="bilinear")
Soil_bioclimatic<-resample(Soil_bioclimatic,Bioclimatic_variable)
######Topography variables
topographic_variables<-stack(list.files(path ="D:/Mr kerich/Dem" ,pattern ="\\.tif$" ,full.names = TRUE))

topographic_variables<-resample(topographic_variables,Bioclimatic_variable)
res(topographic_variables)
crs(topographic_variables)
##Ipomea occurrence data =====================================================================================
ipomea_occurrence<-read.csv("D:/Mr kerich/Ipomea Occurence data/Ipomea_plant.csv")

predictor<-stack(topographic_variables,soil,Soil_bioclimatic,Bioclimatic_variable)
future_stacked_soil_bio<-dropLayer(Soil_bioclimatic,c("SBIO1","SBIO10","SBIO11","SBIO2","SBIO5","SBIO6","SBIO8" ,"SBIO9"))
plot(future_stacked_soil_bio)
names(Soil_bioclimatic)

##multicollnearity analysis====================================================================================
predictor_multi<-terra::extract(predictor,ipomea_occurrence[,c("longitude","latitude")])
collinear_data<-as.data.frame(predictor_multi)
collinear_data<-scale(collinear_data)
########filling missing value with regression anlaysis==============================================

filled_data <- as.data.frame(predictor_multi) 
colnames(filled_data)

for (col_name in colnames(filled_data)) {
  if (any(is.na(filled_data[[col_name]]))) {
    
    train_data <- filled_data[!is.na(filled_data[[col_name]]), ]
    test_data <- filled_data[is.na(filled_data[[col_name]]), ]
    
    # Only use complete predictors (drop columns with NA in train)
    train_data <- train_data[, colSums(is.na(train_data)) == 0]
    
    # Proceed only if enough predictors remain
    if (ncol(train_data) > 1) {
      lm_formula <- as.formula(paste(col_name, "~ ."))
      lm_model <- lm(lm_formula, data = train_data)
      test_data_sub <- test_data[, colnames(train_data)[-1], drop = FALSE]
      complete_cases <- complete.cases(test_data_sub)
      if (any(complete_cases)) {
        predicted <- predict(lm_model, newdata = test_data_sub[complete_cases, ])
        filled_data[is.na(filled_data[[col_name]]), ][complete_cases, col_name] <- predicted
      }
    }
    filled_data[[col_name]][is.na(filled_data[[col_name]])] <- mean(filled_data[[col_name]], na.rm = TRUE)
  }
}

sum(is.na(filled_data)) 
str(filled_data)   
collinear_data <- scale(filled_data)
##plotting Pearson correlation
ggplot(cor_melt, aes(Var1, Var2, fill = Freq)) +
  geom_tile(color = "white") + 
  geom_text(aes(label = round(Freq, 2)), size = 2.5, color = "black") +  # Add correlation values
  scale_fill_gradient2(
    low = "#00796B",     
    mid = "white",     
    high = "#F57C00",    
    midpoint = 0,
    limit = c(-1, 1),
    name = "Pearson\nCorrelation"
  )+#00796B

  theme_minimal() +
  
  theme(
    axis.text.x = element_text(angle = 90,vjust = 0.5, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10), 
    axis.title = element_blank(), 
    plot.title = element_text(size = 12,face = "bold.italic"),
    panel.grid = element_blank(),
    plot.margin = unit(c(1, 1, 1,1), "cm"),
  ) + 
  labs(title = "Pearson Correlation Matrix between predictor variable")
ggsave("correlation_matrix.png", width = 12, height = 10, dpi = 300)
#########Variance inflation factor
result<-vifcor(collinear_data, th=0.80)###Variance Inflation factor using vifcor function (usdm) package
print(result)
######adding presences data to Ipomea occurrence
ipomea_occurrence$presence<-1
######removing highly correlated variable
predictor_modelling<-dropLayer(predictor,c("BIO11","BIO10","BIO17",
                               "BIO16","SBIO1"
                               ,"SBIO10",
                               "BIO9", "BIO6", "BIO5", "SBIO5", 
                               "SBIO11", "BIO13","BIO14", 
                               "SBIO8", "SBIO6",
                               "SBIO9", "BIO19",
                               "SBIO2",) )

names(predictor_modelling)
##Preparing  data for biomod2 , pseudo absence data=============================================================
Biomod_data<- BIOMOD_FormatingData(
  resp.name = "Ipomea",
  resp.var =ipomea_occurrence$presence,
  expl.var = predictor_modelling,  
  resp.xy = ipomea_occurrence[,c("longitude","latitude")],
  PA.nb.rep = 7,  #number of times pseudo absence data is to drawn
  PA.nb.absence =45, #number of pseudo absence data to be generated
  PA.strategy = "random",## randomly selected 
  na.rm = TRUE,
  filter.raster = TRUE,
)
Biomod_data
summary(Biomod_data)
######biomod2  modelling option
opt<- bm_ModelingOptions(data.type = "binary",
                         models = c("RF","MAXNET","XGBOOST"),
                         strategy = "bigboss",
                         bm.format=Biomod_data)

##Random forest,gradient boost machine and Maxent modelling===================================================
Ipomea_model<-BIOMOD_Modeling(bm.format =Biomod_data,
                           models = c("RF","XGBOOST","MAXNET"),
                           CV.strategy = "kfold",
                           CV.nb.rep = 2,
                           CV.k=7,
                           bm.options =opt,
                           metric.eval = c("TSS","ROC"),
                           var.import = 3,
                           do.progress = TRUE,
                           CV.do.full.models = FALSE,
                           nb.cpu=1,
                           seed.val = 24)
                           
                           

str(Ipomea_model)

bm_PlotEvalMean(bm.out = Ipomea_model, dataset ="validation",main="Model evaluation by algorithm on the test data")
Ipomea_model
##variable importance=======================================================================================
ipomea_variable_importance<-get_variables_importance(Ipomea_model)
class(ipomea_variable_importance)
str(ipomea_variable_importance)
ipomea_variable_importance$algo <- factor(ipomea_variable_importance$algo)
unique(ipomea_variable_importance$algo)
                            
######plotting 2d plots of variable importance for the model
ggplot(ipomea_variable_importance, aes(x = expl.var, y = var.imp, fill = algo)) +
  geom_bar(stat = "identity",position="dodge" ) +
  theme_minimal()+
  labs(title = "Variable Importance by Algorithm",x = "Explanatory Variable",y = "Relative Importance") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = c("RF"="blue", "XGBOOST" = "red", "MAXNET" = "green"))
######Single modell projection
projection <- BIOMOD_Projection(
  bm.mod = Ipomea_model,
  new.env= predictor_modelling,  ######can feed new data
  proj.name = "current_distribution_ipomea",
  models.chosen=c("all" ),
  metric.binary = c("TSS"),
  metric.filter = c("TSS"),
  metric.threshold=c(0.57),
  build.clamping.mask=TRUE,######uncertainty 
  do.stack=TRUE,
  nb.cpu=1,
  seed.val=5678
)
plot(projection)
####Ensemble all the models================================================================================
ensembled <- BIOMOD_EnsembleModeling(
  bm.mod = Ipomea_model,
  models.chosen = "all",
  em.by = 'PA+run',##Ensemble by number of run and pseudo absence 
  em.algo ='EMmedian' ,##ensemble by median probabilities of the models
  metric.select = c("TSS","ROC"), 
  metric.select.thresh= c(0.65,0.70),#excluding model <than 
  metric.eval = c("TSS","ROC"),####model performance using ROC and TSS metric from biomode2
  metric.select.dataset = "validation", #ensembles by the test data during training and testing
  var.import = 0,
  do.progress = TRUE,
  nb.cpu = 1,
  seed.val = 156
)
bm_PlotEvalMean(bm.out = ensembled, dataset = 'validation')#####plotting model performance 
############Current and future prediction distribution by Ensemble model #####
ensembel_bar <-BIOMOD_EnsembleForecasting(
  bm.em=ensembled,
  proj.name = "New_current habitat_Kajiado_ready",
  new.env = predictor_modelling,
  new.env.xy = NULL,
  models.chosen = "all",
  metric.binary = c("TSS"),
  
  build.clamp = TRUE,
  nb.cpu = 1,
  na.rm = TRUE,
  do.stack=TRUE,
  )
ensembel_bar2 <-BIOMOD_EnsembleForecasting(
  bm.em=ensembled,
  proj.name = "New_current habitat_Kajiado_ready",
  new.env = future_stacked_ssp,
  new.env.xy = NULL,
  models.chosen = "all",
  metric.binary = c("TSS"),
  
  build.clamp = TRUE,
  nb.cpu = 1,
  na.rm = TRUE,
  do.stack=TRUE,
)
#######fUTURE PREDICTED ===================================================================================
###############future prediction in 2050s data################
#########scenario_data  SSP 126 data #########
bioclimatic_future<-stack(list.files(path = "D:/Mr kerich/Bioclimatic_Variables/Future_bioclimatic/bioclimatic_126",pattern =".tif$",full.names = TRUE))
future_stacked_ssp<-stack(bioclimatic_future,topographic_variables,future_stacked_soil_bio,soil)
#########scenario_data  SSP 585 data #########
bioclimatic_future_585<-stack(list.files(path="D:/Mr kerich/Bioclimatic_Variables/Future_bioclimatic/bioclimatic_585",pattern =".tif$",full.names = TRUE))
future_stacked_ssp_585<-stack(bioclimatic_future_585,topographic_variables,future_stacked_soil_bio,soil)
################CHANGE DISTRIBUTION ANALYSIS
CurrentProj <- get_predictions(ensembel_bar,
                               metric.binary = "TSS",
                               model.as.col = TRUE)
FutureProj <- get_predictions(myBiomodProjectionFuture,
                              metric.binary = "TSS",
                              model.as.col = TRUE)
# Compute differences