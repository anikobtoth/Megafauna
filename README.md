---
title: "README.md"
author: "Anikó B. Tóth"
date: "15/11/2018"
output: html_document
---

# Megafauna
Data and code for manuscript "Reorganization of surviving mammal communities after end-Pleistocene megafaunal extinction" by Toth et al. 

### Data 

All data needed to reproduce this project can be found in the Data folder. 
Data is provided as .RData objects.
The data provided is pre-cleaned, but code is provided for data formatting.
Calibrated radiocarbon dates and climate estimates are provided in the pre-cleaned dataset.

### Code

The main script is in Analysis_Script.R
This file includes a source call to the helper functions file (Helper_functions.R).
Hypervolumes are calculated in a separate script, which is sourced in Analysis_Script.R.
  The source call to Hypervolume_script.R is commented out by because it takes a long time to run. Uncomment to run. 
To reproduce the plots in the manuscript, use Plots_Script.R. 
The Data objects needed for plotting are produced by Analysis_Script.R and Hypervolume_script.R, so they need to be run first.
It is a good idea to run the code in order, as temporary object names are recycled. 
