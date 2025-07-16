# 🌿 Species Distribution Modelling of *Ipomoea* Species in Kenya's🌍
# 🌍 Species Distribution Modelling of Aloe Species in Kenya's ASAL Regions  

---

## 📘 Abstract

This project applies advanced species distribution modelling (SDM) techniques to assess the historical, current, and future potential distribution of **Ipomea plant** in **Kajiado**. Ipomea  plants are invasive threathening other speceis 

We leverage open-source geospatial data (bioclimatic, edaphic, and topographic) and occurrence records from global biodiversity databases. Machine learning models such as **Random Forest (RF)** , **MaxEnt** and **XGboost**are employed to understand the environmental determinants of Ipomea habitat suitability, generate predictive distribution maps, and forecast changes under future climate scenarios in 2040s under climate scenarion SSP126 and SSP 585.

---

## 🎯 Project Objectives

1. **Identify and quantify the environmental variables** influencing ipomea distributions .
2. **Model the spatial distribution of ipomea plant** between 1981 and 2011 using machine learning algorithms.
3. **Forecast the potential future distribution** of ipomea species by integrating projected climate data to inform conservation and adaptation strategies.

---

## 🔬 Methodological Workflow

1. **Data Acquisition**  
   - Bioclimatic variables from [CHELSA v2.1](https://chelsa-climate.org/)
   - Soil data from [SoilGrids](https://soilgrids.org/)
   - Elevation and terrain derivatives from SRTM DEM
   - Species occurrence data from [GBIF](https://www.gbif.org/) and [iNaturalist](https://www.inaturalist.org/)

2. **Data Pre-processing**  
   - Occurrence data cleaning: removing duplicates, spatial thinning
   - Environmental raster alignment and resampling
   - Variable selection using correlation and Variance Inflation Factor (VIF)

3. **Modeling**  
   - Algorithms used: MaxEnt, Random Forest (via `biomod2`)
   - Data split into training (70%) and testing (30%)
   - Evaluation metrics: ROC and  TSS

4. **Prediction and Mapping**  
   - Habitat suitability mapping across timeframes
   - Future projections using downscaled climate scenarios
   - Map visualization with `ggplot2` and `tmap`

5. **Interpretation and Reporting**  
   - Variable importance plots
   - Response curves
   - Interpretation of temporal and spatial patterns

---

📊 Outputs
✔️ High-resolution raster maps of predicted Aloe suitability

✔️ Model performance plots (AUC curves, TSS scores)

✔️ Variable response plots

✔️ Comparative maps showing change between historical, current, and future distributions


🌐 Data Sources
Dataset	Source	Description
CHELSA Bioclim	https://chelsa-climate.org	Climate variables BIO1–BIO19
SoilGrids	https://soilgrids.org	Soil pH, organic carbon, nitrogen, etc.
GBIF	https://gbif.org	Aloe species presence records
SRTM	NASA	Digital elevation data (DEM)

📍 Study Area
Region: Narok County, Kenya

Coordinates: approx. 1.08° S, 35.86° E

Ecozone: Arid and Semi-Arid Lands (ASALs)

Significance: Region of high  biodiversity and vulnerable to climate shifts

💡 Motivation & Impact
Understanding invasive ipomea  respond to environmental changes is critical for:

Biodiversity conservation

Climate-resilient land management

Ethnobotanical value preservation in pastoral communities

This project contributes to climate adaptation planning by delivering actionable geospatial insights for policymakers and conservationists.

🙋 Author
Simon Kipkurui
GIS & Remote Sensing | Earth Observation & AI for Climate | SDM & Environmental Intelligence

🔗 LinkedIn

✉️ simonkipkurui759@gmail.com



