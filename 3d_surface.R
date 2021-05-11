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

Gp.posterior <- function(nu,g){ log(vec_unnormalized.posterior(nu,g) / m_estimate ,base=10) } 

# #########################################################################################
# #########################################################################################

# #########################################################################################
# 3D plotting
# #########################################################################################
#x = seq(1417,1419,length.out = 100)
#y = seq(1.4,2.8,length.out = 100)

x = seq(1400,1440,length.out = 100)
y = seq(0.2,9.39,length.out = 100)

z = outer(x,y,Gp.posterior)

open3d()

persp3d(x = x,
        y = y,
        z = z,
        col = gray.colors(100, start =1, end = 0),
        xlab = "v0 (MHz)",
        ylab = "g",
        zlab = "p(v0,g|D,M,I)"
)