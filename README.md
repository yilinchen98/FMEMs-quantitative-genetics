# FMEMs-quantitative-genetics
This is the PhD project about applying functional mixed-effect models in quantitative genetics.
## Research outline
In evolutionary biology, function-valued traits, such as growth trajectories, are phenotypes of living organisms whose value can be described by a function of some continuous index. These traits can be assessed for continuous genetic variation using a quantitative genetic approach where one of the primary aims is to decompose the genetic and environmental variations. While existing literature has explored mixed effects models for such traits, this research uniquely focuses on using modern functional data analysis methodologies in quantitative genetics. The primary objective is to extend functional mixed-effect models tailored for quantitative genetics experiments. A key secondary objective is to conduct an in-depth theoretical analysis, specifically addressing time and phase variation within these models. This involves advancing our understanding of the dynamic changes in traits over time and capturing variations that may occur due to genetic factors or environmental influences.
## Contents
### R codes
- `TC_smoothing_FPCA.R`: smoothes data and runs functional principal component analysis using the `fdasrvf` package.
- `convert_to_basisfunctions.R`: converts an eigenvector to an eigenfunction using interpolation.
- `TC_computingZ.R`: uses principal components as basis functions and computes the design matrix $Z$ for the mixed-effects model.
### R package help documents
- `fdasrvf`: performs alignment, PCA, and modeling of multidimensional and unidimensional functions using the square-root velocity framework.
- `pedigreemm`: fits pedigree-based mixed-effects models.
