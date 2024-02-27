# FMEMs-quantitative-genetics
This is the PhD project about applying functional mixed-effect models in quantitative genetics.
## Research outline
In evolutionary biology, function-valued traits, such as growth trajectories, are phenotypes of living organisms whose value can be described by a function of some continuous index. These traits can be assessed for continuous genetic variation using a quantitative genetic approach where one of the primary aims is to decompose the genetic and environmental variations. While existing literature has explored mixed effects models for such traits, this research uniquely focuses on using modern functional data analysis methodologies in quantitative genetics. The primary objective is to extend functional mixed-effect models tailored for quantitative genetics experiments. A key secondary objective is to conduct an in-depth theoretical analysis, specifically addressing time and phase variation within these models. This involves advancing our understanding of the dynamic changes in traits over time and capturing variations that may occur due to genetic factors or environmental influences.
## Contents
### R codes
- `TRFUN25PUP4.DAT`: contain 873 individuals and 6860 measurements (body mass), with an average of approximately 8 measurements per individual. Sampling points are not taken at fixed times as they vary in number and location, and the range of days measured is 1-25 days.
- `TC_alignment_fda.R`: smooth data using penalised smoothing spline and align curves by continuous registration from the `fda` package.
- `TC_smoothing_FPCA.R`: smooth, align data and run functional principal component analysis using the `fdasrvf` package.
- `convert_to_basisfunctions.R`: convert an eigenvector to an eigenfunction using interpolation.
- `TC_computingZ.R`: use principal components as basis functions and compute the design matrix $Z$ for the mixed-effects model.
### R package help documents
- `fdasrvf`: perform alignment, PCA, and modelling of multidimensional and unidimensional functions using the square-root velocity framework.
- `pedigreemm`: fit pedigree-based mixed-effects models.
