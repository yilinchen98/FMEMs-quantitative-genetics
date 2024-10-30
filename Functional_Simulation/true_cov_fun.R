## set up true covariance
N <- 873 # total number of individuals
n <- 7 # number of measurements per individual
nbasis <- 5 # number of basis
ngroups <- 50

### use b-spline basis
basisObj <- create.bspline.basis(c(0,1), nbasis = nbasis, norder = 4)

### true covariance matrix
C_true <- matrix(c(750, 10 ,130, 80, 250,
                  10, 800, 30, 15, 40,
                  130, 30, 700, 50, 130,
                  80, 15, 50, 420, 50,
                  250, 40, 130, 50, 330), nrow = 5, byrow = T)

### model reform
time_rang <- seq(0,1,length = 14)
timefine <- seq(0,1, length =100)
basis <- eval.basis(time_rang, basisObj)

C_fun_true <- basis %*% C_true %*% t(basis) # true genetic covariance function

# Create a 3D surface plot using plotly and remove the color bar
plot_ly(x = ~time_rang, y = ~time_rang, z = ~C_fun_true, type = "surface", showscale = FALSE) %>%
  layout(scene = list(
    xaxis = list(title = "Time"),
    yaxis = list(title = "Time"),
    zaxis = list(title = "Covariance Function", range = c(0,800)))
    )

plot_ly(x = time_rang,
        y = time_rang,
        z = ~C_fun_true, type = "heatmap", colorscale = "Viridis",colorbar = list(title = "Covariance")) %>%
  layout(title = "",
         xaxis = list(title = "Time", range = time_rang),
        yaxis = list(title = "Time", range = time_rang))

