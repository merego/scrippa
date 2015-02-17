library("ggplot2")



nx <- 21
ny <- 21
nt <- 50
nit <- 50

c <- 1
dx <- 2.0/(nx-1)
dy <- 2.0/(ny-1)
x <- seq(0, 2.0, length.out = nx)
y <- seq(0, 2.0, length.out = ny)
XY <- expand.grid(x=x,y=y)


rho <- 1
nu <- .05
sigma <- 0.25
#dt <- .001
#dt <- .025
dt <- sigma*dx*dy/nu
c <- 1.0

u <- matrix(0,nx,ny)
v <- matrix(0,nx,ny)
p <- matrix(0,nx,ny)
b <- matrix(0,nx,ny)


u <- matrix(1,nx,ny)      #numpy function ones()
un <- matrix(1,nx,ny)      #numpy function ones()
u [ c((0.5/dx):(1/dx+1)),  c((0.5/dx):(1/dx+1)) ] <- 2  #setting u = 2 between 0.5 and 1 as per our I.C.s
XY$z <- matrix(u)
qplot(x, y, data = XY, fill = z, geom = "raster")

nt <- 50
for (n in seq(1,nt)) { 
  un <- u
  u[2:(nx-1),2:(ny-1)] <- un[2:(nx-1),2:(ny-1)] + 
                          nu * (dt/dx)^2 * ( un[3:nx,2:(ny-1)] - 2* un[2:(nx-1),2:(ny-1)] + un[1:(nx-2),2:(ny-1)] ) +  
                          nu * (dt/dy)^2 * ( un[2:(nx-1),3:ny] - 2* un[2:(nx-1),2:(ny-1)] + un[2:(nx-1),1:(ny-2)] )  
  u[1,] <- 1
  u[nx,] <- 1
  u[,1] <- 1
  u[,ny] <- 1
}
XY$z <- matrix(u,byrow = TRUE)
qplot(x, y, data = XY, fill = z, geom = "raster")

ggplot() +
  geom_raster(data=XY,aes(x=x,y=y,fill=z),hjust = 0, vjust = 0) +
  stat_contour(data=XY,aes(x=x,y=y,z=z))
  

