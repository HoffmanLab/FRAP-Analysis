clear
close all
clc

% get general info about the data sets
Get_FRAP_Info

% for each data set, check if need to generate polygons & parameters -> do it
Init_Poly_Gen

% run through FA_FRAP & FRAP_FIT for each -> take out parts that generate
% initial polygons
% FRAP_Analysis

% do additional FRAP analysis
FRAP_PostAnalysis