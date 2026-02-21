# Glacial-Retreat-Remote-Sensing-Application
Detection of the glacial retreat for the Fox Glacier in New Zealand, using satellite remote sensing (Landsat 4-8, ERA5-Land) with Google Earth Engine and statistical analysis in R.

This Repository contains the complete analysis workflow for investigating the Glacial Retreat, with an Focus on area-loss, temperature dependency and redistribution of snowmass trough an hypsometry-analysis using satellite based remote sensing data. 

This study is Part of the Course: MNF-Geogr-333: Remote Sensing Applications WiSe25/26, 
Author: Justin Lingg-Laham, Liev Paustian, Lisa Wohlgemuth

# Project Overview
The aim of this study was to investigate the evolution of Fox Glacier in New Zealand by analysing time series of glacier area, snowline altitude, and summer air temperatures for the period 1989–2025 (January to March). In addition, correlation analyses were conducted to examine the relationship between glacier area and snowline variability, while combined analyses explored the interaction between glacier area and snowfall–precipitation patterns.

## Methods
The Workflow consists of 2 Main components

### Part 1: Google Earth Engine 

1.1 Define Region of Interest: Fox Glacier, New Zealand
1.2 Glacial Analysis - Get Data about Glacial condition
1.3 Climate Data Extraction
1.4 Hypsometry

### Part 2: R

2.1 Load Data from GEE
2.2 Preprocess Statistics (trends, Correlation)
2.3 Plots: Glacial conditions
2.4 Climate Analysis (Preciptation, Temperature)
2.5 Plots: Climate TimeSeries
2.6 Plots: PReciptation + Temperature
2.7 Plots: Anomaly TimeSeries
2.8 Hypsometrie


