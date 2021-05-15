library('ggplot2')
library('rgl')

m_estimate = 1.773521e+68

# #########################################################################################
# Stuff already in Jupyter Notebook
# #########################################################################################

dados <- read.csv(file = 'data_Antonio.R',sep=' ', header=FALSE)
names(dados) = c('freq','intensity')

Nu.data <- dados[1]
G.data <- dados[2]

sigma <- 0.05
N <- 101

nu.min <- Nu.data$freq[1]
nu.max <- tail(Nu.data$freq,n=1)
Delta.nu <- nu.max - nu.min

dnu = Nu.data$freq[2]-Nu.data$freq[1]
g.min = dnu/2
g.max = sqrt(log(2)/pi)/sigma
Delta.g = g.max - g.min

nu_i = Nu.data$freq
d_i = G.data$intensity
S_1 = sum((d_i-1)^2)

# #############################################################################

f_i <- function(nu,g){exp(-log(2)*(nu_i-nu)^2/g^2) * sqrt(log(2)/pi) / g}

# #############################################################################

log.Likelihood <- function(nu,g){
  fi = f_i(nu,g)
  
  (-0.5/sigma^2) * ( S_1 + sum( fi^2 - 2*(d_i-1)*fi ) )
}

unnormalized.posterior <- function(nu,g){
  
  prior = 1 / ( Delta.nu * log(g.max/g.min) * g)
  L.likelihood = log.Likelihood(nu,g) 
  
  prior * exp(L.likelihood) / (sigma * sqrt(2*pi))^N
}

vec_unnormalized.posterior <- function(Nu,G){Vectorize(unnormalized.posterior,vectorize.args = c('nu','g'))(Nu,G)}

log10.Gp.posterior <- function(nu,g){ log( vec_unnormalized.posterior(nu,g) / m_estimate , base=10) } 

# #########################################################################################
# #########################################################################################

# #########################################################################################
# 3D plotting
# #########################################################################################
#x = seq(1417,1419,length.out = 100)
#y = seq(1.4,2.8,length.out = 100)

x = seq(1414,1422,length.out = 200)
y = seq(1.2,3,length.out = 200)
z = outer(x,y,log10.Gp.posterior)
z[!is.finite(z)] = NA
open3d()

persp3d(x = x,
        y = y,
        z = z,
        col = hsv(0.5, .35, seq(.45,.90,length.out = 12)),
        xlab = "v0 (MHz)",
        ylab = "g",
        zlab = "log10 [ p(v0,g|D,M,I) ]"
)