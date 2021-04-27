Analysis files for Schandry, et al., 2021: Plant-derived benzoxazinoids act as antibiotics and shape bacterial communities
================

# About

This repository contains the code used for the analysis presented in Schandry, et al., 2021.

# Contents

The script Figures_only.Rmd contains the most important steps of the analysis used for Figure creation.

The files Reads_in.Rmd contains read preprocessing.
Growth_trees.Rmd contains the script to create Figure 1.
Community_Analysis.Rmd contains code for all other analysis included in the manuscript.

The folders _files_publication_ and _rds_ contain rds files that are created at various stages of the scripts. This should enable running each of the .Rmd files independently of each other.

The folder _functions_ contains custom functions that were used for the analysis, which are sourced in the respective scripts. These functions are not documented in detail, but there is some minimal information in roxygen format included. 

A renv environment snapshot is included and it should be possible to restore the analysis environment using renv::restore(). 