## Display Heat Maps
rm(list = ls())
cat("\014")
graphics.off

#unique to each computer, change for new 
setwd("~/Desktop/Files/School/USC/Radiomics Lab/Code") 

## package install
#install.packages("readxl")
#install.packages("phytools")
#install.packages("hrbrthemes")
#install.packages("viridis")
#install.packages("phyloseq")
#install.packages("plotly")
#install.packages(patchwork)

## libraries
library(ggplot2)
library(hrbrthemes)
library(plotly)
library(readr)
library(readxl)
library(tidyr)
library(tibble)
library(hrbrthemes)
library(dplyr)
library(data.table)
library(patchwork)

sortAndPlotHM <- function(fileName,graphTitle)
{
  dataTableForm <- read_csv(fileName,show_col_types = FALSE)
  framedData = data.frame(dataTableForm)
  setorder(framedData,featureNum)
  FeatureLabel = framedData[,1]
  ROI = framedData[,2]
  FeatureValue = framedData[,3]
  FeatureNumber = framedData[,4]
  FeatureGroup = framedData[,5]
  framedData <- framedData %>%
    mutate(tooltiptext = paste0("Category: ", ROI, "\n", "Image Feature: ", FeatureLabel, "\n", "Value: ",round(FeatureValue,2)))
  caps <- quantile(FeatureValue, probs=c(.05, .95), na.rm = T)
  #subSetLower <- framedData[framedData$dataColumn < quantile(framedData$dataColumn, 0.05)]
  #subSetUpper <- framedData[framedData$dataColumn > quantile(framedData$dataColumn,0.95)]

  FeatureValue[FeatureValue < caps[1]] <- caps[1] # if cutting out outliers
  FeatureValue[FeatureValue > caps[2]] <- caps[2] # if cutting out outliers
  statPlotHM = ggplot(framedData, aes(x = ROI, y = FeatureNumber, fill= FeatureValue),show.legend = FALSE) + 
    geom_tile() +
    scale_y_continuous(breaks = c(1, 20, 50, 206,302, 350, 398, 413))+
    annotate(geom = "text", x = 6, y = c(10, 35, 128, 254, 326, 374,405,438),
                       label = c("Intensity","Histogram","GLCM","GLRLM","GLSZM","GLDZM","NGTDM","NGLDM"))+
              
    coord_cartesian(xlim = c(0,7), expand = FALSE, clip = "off") +
    ggtitle(graphTitle) +
    xlab("ROI") +
    ylab("Image Feature")+
    scale_fill_gradient2(low="red", mid = "white", high="blue")+ 
    guides(size = "none",colour = guide_colourbar(title.position = "right")) +
    theme(legend.key.height = unit(2.5, "cm"),
          legend.title = element_text(size = 12, angle = 90),
          legend.title.align = 0.5,
          legend.direction = "vertical"
    )
  
  ggplotly(statPlotHM, hoverinfo="tooltiptext")
  #
}


sortAndPlotArtPhase <- function(fileName,graphTitle)
{
  dataTableForm <- read_csv(fileName,show_col_types = FALSE)
  framedData = data.frame(dataTableForm)
  colnames(framedData)[4]="featureNum"
  setorder(framedData,classColumn)
  for (kernelNum in 1:nrow(framedData))
    {actualKernel = framedData[kernelNum,2]
      if (actualKernel == 1 |actualKernel == 2 | actualKernel ==3 | actualKernel == 4
          | actualKernel == 5 | actualKernel == 6 | actualKernel == 7 |
          actualKernel == 8 | actualKernel == 9)
        {actualKernel = paste0("0",actualKernel)}
      framedData[kernelNum,2]=paste0("FC_",actualKernel)
      }
  FeatureLabel = framedData[,1]
  numRows = nrow(framedData)
  Kernel = framedData[,2]
  FeatureValue = framedData[,3]
  FeatureNumber = framedData[,4]
  FeatureGroup = framedData[,5]
  
  framedData <- framedData %>%
    mutate(tooltiptext = paste0("Category: ", Kernel, "\n", "Image Feature: ", FeatureNumber, "\n", "Value: ",round(FeatureValue,2)))
  statPlotheatMap = ggplot(framedData, aes(x = Kernel, y = FeatureNumber, fill= FeatureValue)) + 
    geom_tile() +
    scale_y_continuous(breaks = c(1, 20, 50, 206,302, 350, 398, 413))+
    annotate(geom = "text", x = 58, y = c(10, 35, 128, 254, 326, 374,405,438),
             label = c("Intensity","Histogram","GLCM","GLRLM","GLSZM","GLDZM","NGTDM","NGLDM"))+
    coord_cartesian(xlim = c(1,65), expand = FALSE, clip = "off") +
    ggtitle(graphTitle) +
    xlab("Image Kernel") +
    ylab("Image Feature")+
    scale_fill_gradient2(low="red", mid = "white", high="blue")+ 
    guides(size = "none",colour = guide_colourbar(title.position = "right")) +
    theme(legend.key.height = unit(2.5, "cm"),
          legend.title = element_text(size = 12, angle = 90),
          legend.title.align = 0.5,
          axis.text.x=element_text(angle=-90),
          legend.direction = "vertical"
    )

  ggplotly(statPlotheatMap, hoverinfo="tooltiptext")
  #
}


sortAndPlotVenPhase <- function(fileName,graphTitle)
{
  dataTableForm <- read_csv(fileName,show_col_types = FALSE)
  framedData = data.frame(dataTableForm)
  colnames(framedData)[4]="featureNum"
  setorder(framedData,classColumn)
  for (kernelNum in 1:nrow(framedData))
  {actualKernel = framedData[kernelNum,2]
  if (actualKernel == 1 |actualKernel == 2 | actualKernel ==3 | actualKernel == 4
      | actualKernel == 5 | actualKernel == 6 | actualKernel == 7 |
      actualKernel == 8 | actualKernel == 9)
  {actualKernel = paste0("0",actualKernel)}
  framedData[kernelNum,2]=paste0("FC_",actualKernel)
  }
  FeatureLabel = framedData[,1]
  numRows = nrow(framedData)
  Kernel = framedData[,2]
  FeatureValue = framedData[,3]
  FeatureNumber = framedData[,4]
  FeatureGroup = framedData[,5]
  
  framedData <- framedData %>%
    mutate(tooltiptext = paste0("Category: ", Kernel, "\n", "Image Feature: ", FeatureNumber, "\n", "Value: ",round(FeatureValue,2)))
  statPlotheatMap = ggplot(framedData, aes(x = Kernel, y = FeatureNumber, fill= FeatureValue)) + 
    geom_tile() +
    scale_y_continuous(breaks = c(1, 20, 50, 206,302, 350, 398, 413))+
    annotate(geom = "text", x = 65, y = c(10, 35, 128, 254, 326, 374,405,438),
             label = c("Intensity","Histogram","GLCM","GLRLM","GLSZM","GLDZM","NGTDM","NGLDM"))+
    coord_cartesian(xlim = c(1,70), expand = FALSE, clip = "off") +
    ggtitle(graphTitle) +
    xlab("Image Kernel") +
    ylab("Image Feature")+
    scale_fill_gradient2(low="red", mid = "white", high="blue")+ 
    guides(size = "none",colour = guide_colourbar(title.position = "right")) +
    theme(legend.key.height = unit(2.5, "cm"),
          legend.title = element_text(size = 12, angle = 90),
          legend.title.align = 0.5,
          axis.text.x=element_text(angle=-90),
          legend.direction = "vertical"
    )
  
  ggplotly(statPlotheatMap, hoverinfo="tooltiptext")
  #
}

sortAndPlotHM("venSD.csv","Venous Standard Deviation")
sortAndPlotHM("artSD.csv","Arterial Standard Deviation")

sortAndPlotArtPhase("complexCystZ.csv","Complex Cyst Z Scores (Arterial)")
sortAndPlotArtPhase("fatZ.csv","Fat ROI Z Scores (Arterial)")
sortAndPlotArtPhase("muscleZ.csv","Muscle ROI Z Scores (Arterial)")
sortAndPlotArtPhase("liverZ.csv","Liver ROI Z Scores (Arterial)")
sortAndPlotArtPhase("spleenZ.csv","Spleen ROI Z Scores (Arterial)")

sortAndPlotVenPhase("complexCystZven.csv","Conplex Cyst Z Scores (Venous)")
sortAndPlotVenPhase("fatZven.csv","Fat Z Scores (Venous)")
sortAndPlotVenPhase("muscleZven.csv","Muscle ROI Z Scores (Venous)")
sortAndPlotVenPhase("liverZven.csv","Liver ROI Z Scores (Venous)")
sortAndPlotVenPhase("spleenZven.csv","Spleen ROI Z Scores (Venous)")