---
title: "README.md"
author: "Anikó B. Tóth"
date: "15/11/2018"
output: html_document
---

# Megafauna
Data and code for manuscript "End-Pleistocene extinction caused a fundamental shift in mammal survivor community structure" by Toth et al. 

### Data 

All data needed to reproduce this project can be found in the Data folder. 
Data is provided as .RData objects.
The data provided is pre-cleaned, but code is provided for data formatting
Calibrated radiocarbon dates and climate estimates are provided.

### Code

The main script is in Analysis_Script.R
This file includes a source call to the helper functions file (Helper_functions.R).
To reproduce the plots in the manuscript, use Plots_Script.R. 
The Data objects needed for plotting are produced by Analysis_Script.R, so it needs to be run first.
