# FMEMs-quantitative-genetics
This is the PhD project at KCL about applying functional mixed-effect models in quantitative genetics.
## Research Summary
Outline of the project research and expected outcomes:	In evolutionary biology, function-valued traits, such as growth trajectories, are phenotypes of living organisms whose value can be described by a function of some continuous index. These traits can be assessed for continuous genetic variation using a quantitative genetic approach where one of the primary aims is to decompose the genetic and environmental variations. While existing literature has explored mixed effects models for such traits, this research uniquely focuses on using modern functional data analysis methodologies in quantitative genetics. The primary objective is to extend functional mixed-effect models tailored for quantitative genetics experiments. A key secondary objective is to conduct an in-depth theoretical analysis, specifically addressing time and phase variation within these models. This involves advancing our understanding of the dynamic changes in traits over time and capturing variations that may occur due to genetic factors or environmental influences.
## Contents
The pdf file `FMEMs in Quantitative Genetics` contains an overview of this research project, including the mathematical formulation of the functional mixed-effect model,  model representations and fitting procedures. This document will be updated periodically to reflect the progress made in the research project.
### R Code
- `TRFUN25PUP4.DAT`: containing 873 individuals and 6860 measurements (body mass), with an average of approximately 8 measurements per individual. Sampling points are not taken at fixed times as they vary in number and location, and the range of days measured is 1-25 days.
- `convert_to_basisfunctions.R`: converting an eigenvector to an eigenfunction using interpolation.
- `fit_genetic_fmm.R`:  fitting a functional mixed-effect model to genetic data using the `lme4` package in R. 
- `TC_alignment_fda.R`: smoothing data using penalised smoothing spline and align curves by continuous registration from the `fda` package.
- `TC_smoothing_FPCA.R`: smoothing, aligning data and running functional principal component analysis using the `fdasrvf` package.
- `TC_computingZ.R`: using principal components as basis functions to fit random-effect model. Here only focus on computing the random-effect design matrices. We test our code on a subset consisting of 3 subjects from the original dataset and examine the structure of $Z^G$ and $Z^E$.
- `Compare_genetic_fmm.R`: (up to Marth 7th) This R script fits functional mixed-effect models with genetic and environmental random effects, using the principal components obtained from running FPCA on the aligned curves. Here we consider two ways to fit the fixed effect and compare the corresponding results. We also construct the genetic, environmental, and phenotypic covariance functions from the fitted mixed-effect models and visualise them using 3D surface plots.

### R package help documents
- `fdasrvf`: perform alignment, PCA, and modelling of multidimensional and unidimensional functions using the square-root velocity framework.
- `pedigreemm`: fit pedigree-based mixed-effects models.
