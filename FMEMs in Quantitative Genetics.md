# FMEMs in Quantitative Genetics: Research Project Outline
## Genetic Framework
### Quantiative Genetics
* Quantitative genetics studys the quantitative characters that are characterised by the presence of continuous type of variation in a population. 
* Function-valued (FV) traits, such as growth trajectories are phenotypes of living organisms whose value can be described by a function of some continuous index. This functional nature suggests a potential of using modern functional data analysis in quantitative genetics. 
* A key concept here is that the variation between individuals can be decomposed into genetic composition and environmental effects.
$$ Var(y) = V_G + V_E$$
In particular, the genetic variation is estimated by comparing phenotypical variation in related individuals where the genetic relationship is modeled by the additive relationship matrix $A$. It is a positive definite matrix with each term $A_{ij}$ represents how closely related between individuals $i$ and $j$. An example of $A$:
$$
\begin{pmatrix}
1 & 0.5 & 0.125 & 0 & \dots \\
0 & 1 & 0.25 & 0 &\dots \\
0 & 0.5 & 1 & 0.25 & \dots \\
\vdots
\end{pmatrix}
$$
* In a functional setting, the change in mean phenotype in one generration of selection, $ \Delta \bar{y}(t)$, can be modeled by 
 $$ \Delta \bar{y}(t) = \int_\mathbb{R} G(s,t)\beta(s) \: ds$$ where $G(s,t)$ is the additive genetic covariance function and $\beta(s)$ is the selection gradient function. The selection gradient in the direction of the first eigenfunction of the additive genetic covariance would result in the maximum response to selection. Therefore, our main goal is to estimate the additive genetic covariance function.
### Data Description
In this project, we will be working on the $\textit{Triboleum Castaneum}$ dataset, $\textrm{TRFUN25PUP4}$. All $\textit{Tribolium}$ individuals were derived from a stock population of the cSM++ strain, with approximately 300 individuals initially divided into 17 populations of about 20 individuals each and allowed to breed freely. Offspring were collected from these populations as they pupated, their sex was determined and recorded, and they were isolated into separate vials and used to create a stock of isolated, virgin adults. 

A half-sib/full-sib breeding design was used to facilitate quantitative genetic analysis. From the stock of isolated virgin adults, one randomly chosen male was mated with five randomly chosen females, none of which were his siblings. This was repeated 30 times, using a total of 30 males and 150 females, thereby producing 30 half-sib families and 150 full-sib families. The sire, dam and subjects' ids in the the dataset will be used to build the pedigree and extract the additive genetic relationship matrix $A$. 

The trait reported in this dataset is the body mass measured for each lava at different ages during the larval period and mass at pupation was included as the final mass measure for each growth curve. It contains 873 individuals and 6860 measurements, with an average of approximately 8 measurements per individual. Sampling points are not taken at fixed times as they vary in number and location, the range of days measured is 1-25 days. 

## Functional Mixed-Effect Model
We are interested in fitting a functional mixed-effect model
$$
\underbrace{Y_{ij}}_{\text{trait}} = \underbrace{u(t_{ij})}_{\text{fixed effect: population mean}} + \underbrace{\sum_{k=1}^K \alpha_{ik}\phi(t_{ij})}_{\text{genetic random effect}} + \underbrace{\sum_{k=1}^K \gamma_{ik}\phi(t_{ij})}_{\text{environmental random effect}} + \underbrace{\epsilon_{ij}}_{\text{measurement error}}
$$ 
In matrix form: $$\mathbf{Y} = \bm{X\beta} + \bm{Zu} + \bm{\epsilon}$$ where $\bm{X}$ is the predictor matrix; $\bm{\beta}$ is the fixed-effect coefficients; $\bm{Z} = [\bm{Z^G}, \bm{Z^E}]$ is the design matrix for both genetic and environmental random effects; $\bm{u} = [\bm{\alpha}, \bm{\gamma}] ^\mathit{T}$ is the random-effect vector; $\bm{\epsilon}$ is the meansurement error vector. 

The random terms are assumed to have the following distributions: 

\begin{align*}
\bm{\alpha} & \thicksim N(\bm{0}, \bm{A} \otimes \bm{C^G}) \\
\bm{\gamma} & \thicksim N(\bm{0}, \bm{I} \otimes \bm{C^E}) \\
\bm{\epsilon} & \thicksim N(\bm{0}, \sigma^2_{res} \otimes \bm{I})
\end{align*}
$$

The covariance is 
$$
\begin{equation*}
Cov(\mathbf{Y}_{ij}, \mathbf{Y}_{ij'}) = \sum_{k = 1}^K \sum_{l=1}^K \phi_k(t_{ij})\phi_l(t_{ij'})Cov(\alpha_{ik}, \alpha_{il}) + \sum_{k = 1}^K \sum_{l=1}^K \phi_k(t_{ij})\phi_l(t_{ij'})Cov(\gamma_{ik}, \gamma_{il}) + Cov(\epsilon_{ij},\epsilon_{ij'})
\end{equation*}
$$

### Model Reparameterisation
### Curve Registration
## Model-Fitting Procedure
### 1. Data Smoothing
### 2. Curve Alignment
### 3. Functional Principal Component Analysis
### 4. Pedigree and Genetic Relationship Matrix
### 5. Mixed-Effect Model Formula
#### Fit Fixed Effect
#### Fit Random Effect
### 6. `fit_genetic_fmm()`
### 7. Covariance Functions
## What to do next?