
setwd("~/Desktop/stat222/")
library(dplyr)
library(stringr)


###########################################
#### Read in data

## Downloaded from https://bpd-transparency-initiative-berkeleypd.hub.arcgis.com/datasets
cases <- read.csv("cases.csv")
cases <- cases[cases$Statute == "487 (A)", ]
cases$Year <- as.numeric(substring(cases$Occurred_Datetime, 1, 4))

cases <- (cases %>% filter(Year >= 2020, Incident_Type == "Theft Felony (over $950)")
          %>% select(Case_Number, Occurred_Datetime, Block_Address, ZIP_Code, Year))


## Get time stamps as numeric
myorigin <- as.POSIXct("2020/01/01 00:00:00")
pdates <- as.POSIXct(cases$Occurred_Datetime)
times <- as.numeric(pdates) - as.numeric(myorigin)

tstamps <- format(pdates, format = "%H:%M:%S")
hourCounts <- table(as.numeric(str_extract(tstamps, "^\\d{2}")))
barplot(hourCounts)

dates <- lubridate::date(pdates)
countByDay <- aggregate(cases$Case_Number, by = list(dates), length)
plot(countByDay$Group.1, ksmooth(countByDay$Group.1, countByDay$x, b = 10, x.points = countByDay$Group.1)$y, 
     t = "l", xlab = "", ylab = "Number of thefts", main = "Daily thefts over time")

### Geographic location
# match addresses from crime report datato lon/lat coords in another police dataset
locs <- read.csv("heatmap_data.csv") # links address to geocoded coords
locs <- locs %>% select(Block_Location, BLKADDR)
locs$BLKADDR <- gsub("BLOCK ", "", locs$BLKADDR)
locs$BLKADDR <- gsub("/ ", "&", locs$BLKADDR) 

checkMatch <- list(old = c(), new = c())
case.coords <- list(x = c(), y = c())
KnownLocs <- unique(locs$BLKADDR)
for (addy in cases$Block_Address) {
  lonlat <- "0"
  if (addy %in% KnownLocs) {
    lonlat <- locs[locs$BLKADDR == addy, "Block_Location"][1]
  } else { # if unknown, use the location of nearest known address
    stnum <- str_extract(addy, "^\\d+")
    stname <- gsub("^\\d+\\s", "", addy)
    lastword <- str_extract(stname,"\\w+$")
    
    ## some of the entries are messy and omit st/ave/etc, only use street name
    if (lastword %in% c("AVE", "ST", "RD", "LN", "WAY", "BLVD", "DR", "HWY", "CIR", "PL", "CT")) {
      stname <- gsub(lastword, "", stname)
      stname <- trimws(stname, which = "right")
    }
    
    neighbors <- KnownLocs[grepl(stname, KnownLocs)] # look for matches on the same street
    
    # if you have a st number and neighbors, choose the nearest building number
    if ((!is.na(stnum)) & (length(neighbors) > 0)) { 
      stnum <- as.numeric(stnum)
      neighborNums <- as.numeric(str_extract(neighbors, "^\\d+"))
      # if a neighbor is missing building number give it heavy negative weight
      neighborNums <- ifelse(is.na(neighborNums), -999, neighborNums)
      nearestNeighbor <- which.min(abs(stnum - neighborNums))
      addy2 <- neighbors[nearestNeighbor]
      lonlat <- locs[locs$BLKADDR == addy2, "Block_Location"][1]
      checkMatch$old <- c(checkMatch$old, addy) # track where buildings are matched to
      checkMatch$new <- c(checkMatch$new, addy2)
    } 
  }
  lon <- str_match(lonlat, "\\d{2}\\.\\d+")[,1]
  lat <- str_match(lonlat, "-\\d{3}\\.\\d+")[,1]
  case.coords$x <- c(case.coords$x, lon)
  case.coords$y <- c(case.coords$y, lat)
}

length(case.coords$y) == nrow(cases)
mean(!is.na(case.coords$y))

# percent exact match
1 - (length(checkMatch[[1]]) + sum(is.na(case.coords$y))) / nrow(cases)
# percent matched
length(checkMatch[[1]]) / nrow(cases)
# percent missing
sum(is.na(case.coords$y)) / nrow(cases)

missingLoc <- is.na(case.coords$y)

library(sf)
coordDf <- (data.frame(case.coords)[!missingLoc,] 
            # %>% filter(y < -122.31)
            %>% st_as_sf(coords = c("y", "x"), crs = 4326)
            %>% st_transform(crs = 32610)) # project lon/lat coordinates

x <- st_coordinates(coordDf)[,1]
y <- st_coordinates(coordDf)[,2]
times <- times[!missingLoc]

coordDf2 <- data.frame(case.coords)[!missingLoc,] %>% st_as_sf(coords = c("y", "x"), crs = 4326)
xlon <- st_coordinates(coordDf2)[,1]
ylat <- st_coordinates(coordDf2)[,2]

distMat <- st_distance(coordDf, coordDf)
maxD <- as.numeric(max(distMat))

### Convex boundary
pts <-list(x = x, y = y)
hpts <- chull(pts)
hpts <- c(hpts, hpts[1]) # makes closed hull
conv.bdry <- st_coordinates(coordDf[hpts,])


#####################################
###### Analysis

library(spatstat)
library(splancs)
library(spatial)

#### Map / kernel smooth
n <- length(x)
stddist = sqrt(((n-1) / n) * (var(x) + var(y))) ## standard distance 
ds = sqrt((x - mean(x))^2 + (y - mean(y))^2) ## distances to mean 
dm = median(ds) 
bdw = .9 * min(stddist, sqrt(1/log(2)) * dm)*(n^(-.2))
bdw <- sqrt(bw.nrd(x)^2 + bw.nrd(y)^2)  

berkeley <- st_read("Land Boundary/")
berkeley <- st_transform(berkeley, crs = 32610)
bowin <- as.owin(berkeley)

bdry <- st_coordinates(st_as_sfc(st_bbox(berkeley)))[,c("X", "Y")]
b1 <- st_coordinates(coordDf)


#### KDE
z = kernel2d(st_coordinates(coordDf), conv.bdry, bdw)
image(z, col=gray((64:20)/64), xlab="x", ylab="y", 
      main = "Catalytic Converter Thefts, 2020-23")
plot(bowin, add = T)
points(b1, lwd = 0.5, cex = 0.5, pch = 4)
lines(st_coordinates(coordDf[hpts,]), lty = 2, col = "red")

## K function: identifies whether clustering exists at different scales
s = seq(0, maxD, length = 100)
k4 = khat(b1, bdry, s, checkpoly = TRUE)
plot(s, k4 / 1e5, xlab="distance", ylab = "K(h) * 1e-5", pch = "*", 
     ylim = range(k4, pi * s^2)/1e5, main = "K(h)")
lines(s, k4 / 1e5)
lines(s, pi * s^2 / 1e5, lty=2)
legend("topleft", lty = c(1, 2), legend = c("Khat(r)", "Kpois(r)"))


## L function
L4 = sqrt(k4/pi) - s
plot(c(0, maxD), range(L4), type="n",xlab = "lag, h",ylab = "L(h) - h")
points(s, L4, pch = "*")
lines(s, L4)
abline(h = 0, lty = 2)


## F function (empty-space function): 
## The cumulative distribution function (cdf), F, 
## of the distance from a fixed location to the nearest point of X. 
## Lower F indicates clustering. 
## If F(0.2) = 0.4, for instance, then 
## 40% of locations are within distance 0.2 of a point of the process. 
b2 = as.ppp(coordDf[!duplicated(coordDf$geometry),]) 
par(mfrow=c(1,1)) 
f4 = Fest(b2, correction = "best", r = seq(0, maxD, length = 4000)) 
plot(f4, main = "F(h)") 
max(f4$theo)

#### G function: nearest neighbor
g4 <- Gest(b2, r = seq(0, 300, length = 1000)) 
plot(g4$r, g4$km, xlab = "h", ylab = "G(h)", type = "l", lty = 1, main = "G(h)") 
lines(g4$r, g4$theo, lty=2) 
legend("topleft", lty = c(1,2), legend = c("data","Poisson")) 

#### J-function: 
## J(r) = (1-G(r))/(1-F(r)). 
## J = 1 corresponds to a stationary Poisson process. 
## J < 1 indicates clustering. J > 1 indicates inhibition. 
j4 = Jest(b2, r = seq(0, 400, length = 1000), correction = "best") 
plot(j4) 


#### Run MISD
times2 <- sort(times) / 3600 # converted to hours
tdiffs <- c(dist(times2, method = "manhattan")) # in units of hours

hist(tdiffs[tdiffs < 24 * 3], breaks = 100,
     xlab = "Hours between reported thefts", 
     main = "Time between events")

source("misd.R")
grange <- seq(0, 24 * 7)
tmisd <- misd(times2, grange, 
              # tot_time = 365 * 24,
              tot_time = max(times2),
              num_iter = 1000)

tmisd$mu
gsm <- ksmooth(grange[-1], tmisd$g / (tmisd$delta_t * sum(tmisd$g)),
               x.points = grange[-1], b = bw.nrd(grange))
bestFitExp <- function(rate) {
  g <- tmisd$g / (tmisd$delta_t * sum(tmisd$g))
  return(sqrt(mean((g - dexp(grange[-1], rate = 1/rate))^2)))
}
myrate <- optimize(f = bestFitExp, interval = c(1, 50))

plot(gsm, t = "l", xlab = "Hours between events", xaxt = "n",
     ylab = "g(t)", main = "MISD Estimated g(t)")
lines(grange[-1], dexp(grange[-1], rate = 1/myrate$minimum), t = "l", col = "blue")
axis(1, at = seq(0, 24 * 7, by = 24))
legend("topright", lty = c(1, 1), col = c("black", "blue"),
        legend = c("MISD", "Exp(33)"))


### Maximum Likelihood
# z <- list(t = times2, lon = xlon[order(times)],
#           lat = ylat[order(times)], n = length(times))
# xtrans <- (x - min(x)) / max(x - min(x))
# ytrans <- (y - min(y)) / max(y - min(y))
xtrans <- xlon - min(xlon)
ytrans <- ylat - min(ylat)
z <- list(t = times2,  n = length(times),
          lon = xtrans[order(times)],
          lat = ytrans[order(times)])
T <- max(times2)
X1 <- max(xtrans)
Y1 <- max(ytrans)

##### This is for fitting a Hawkes model with no magnitudes. 
##### lambda(t,x,y) = mu rho(x,y) + K SUM gt(t-t_i) gxy(x-xi,y-yi),  
##### with rho(x,y) = 1/(X1Y1), 
##### gt(t) = beta e^(-beta t),
##### g(x,y) = alpha/pi exp(-alpha r^2), with x^2+y^2=r^2,
loglhawk = function(theta, draw = 0){
  K = theta[1]; alpha = theta[2];# beta = theta[4] 
  beta <- 1/33
  mu = tmisd$mu
  
  if(min(mu,K,alpha,beta) < 0.000000001) return(99999) 
  if(K > .99999) return(99999)
  
  sumlog = log(mu/X1/Y1) 
  intlam = mu * T + K * z$n
  const = K * alpha/pi * beta
  for(j in 2:(z$n)){
    gij = 0
    for(i in 1:(j-1)){
      r2 = (z$lon[j]-z$lon[i])^2 + (z$lat[j] - z$lat[i])^2
      gij = gij + exp(-beta*(z$t[j]-z$t[i]) - alpha * r2)
    }
    lamj = mu / X1 / Y1 + const * gij
    if(lamj < 0){
      cat("lambda ",j," is less than 0.")
      return(99999)
    }
    sumlog = sumlog + log(lamj)
  }
  loglik = sumlog - intlam
  # cat("loglike is ", loglik, ". sumlog = ", sumlog,". integral = ", intlam,".\n")
  return(-1.0*loglik)
}



start <- Sys.time()
theta1 = c(.5,.5) 
b1 = optim(theta1,loglhawk)
end1 <- Sys.time()
difftime(end1, start, units = "mins")
b2 = optim(b1$par,loglhawk,hessian=T)
difftime(Sys.time(), end1, units = "mins")

theta2 = b2$par
sqrt(diag(solve(b2$hess)))

theta2 <- c(theta2, "beta" = 0.03)
round(theta2, 3)

#####################################
##################### Stoyan Method
xtrans <- xlon - min(xlon)
ytrans <- ylat - min(ylat)
z <- list(t = times2,  n = length(times),
          lon = xtrans[order(times)],
          lat = ytrans[order(times)])
T <- max(times2) + 1e-3
X1 <- max(xtrans)
Y1 <- max(ytrans)

## bin the points
wbin = list()
for(i in 1:1000) wbin[[i]] = c(0)
for(m in 1:z$n) {
  gridindex = 10 * 10 * floor(z$t[m]*10/T) 
  gridindex <- gridindex + 10 * floor(z$lon[m]*10 / X1) + ceiling(z$lat[m] * 10 / Y1)
  wbin[[gridindex]] = c(wbin[[gridindex]], m)
}
for(i in 1:1000) wbin[[i]] = wbin[[i]][-1]


sumsqstoyan = function(theta){
  mu = theta[1]; K = theta[2]; alpha = theta[3]; 
  # beta = theta[4]
  beta = 1/33
  
  if(min(mu,K,alpha,beta) < 0.000000001) return(99999) 
  if(K > .99999) return(99999)

  const = K*alpha/pi*beta
  b = T * X1 * Y1/10/10/10
  mysum = rep(b, 1000)
  for(i in 1:1000){ ## i is the bin index. 
    if(length(wbin[[i]]) > .5){
      mysum[i] = 0
      for(j in wbin[[i]]){ ## j is the index of a point in bin i. 
        gkj = 0
        if(j>1) for(k in 1:(j-1)){ ## k are indices of previous points. 
          r2 = (z$lon[j]-z$lon[k])^2+(z$lat[j]-z$lat[k])^2
          gkj = gkj + exp(-beta*(z$t[j]-z$t[k])-alpha*r2)
        }
        lambdaj = mu/X1/Y1 + const*gkj
        if(lambdaj < 0){
          return(99999)
        }
        mysum[i] = mysum[i] + 1/lambdaj
      }
    }
  }
  sum((mysum - b)^2)
}

theta1 = c(.08,.75, 2.5)/2
start <- Sys.time()
b1 = optim(theta1,sumsqstoyan,control=list(maxit=100))
end1 <- Sys.time()
difftime(end1, start, units = "mins")
b2 = optim(b1$par,sumsqstoyan,hessian=T,control=list(maxit=100))
difftime(Sys.time(), end1, units = "mins")

theta = b1$par
theta <- c(theta, "beta" = 0.03)
round(theta, 3)
sqrt(diag(solve(b2$hess))) ## for SEs 


### g(t)
plot(gsm, t = "l", xlab = "Hours between events", xaxt = "n",
     ylab = "g(t)", main = "Estimated g(t)")
axis(1, at = seq(0, 24 * 7, by = 24))
lines(grange, dexp(grange, 0.03), col = "blue")
# lines(grange, dexp(grange, theta2[4]), col="blue") 
legend("topright", col = c("black", "blue"), 
       legend = c("MISD", "Exp(0.03)"), lty = 1)


# #### g(x,y)
# library(MASS)
# g4 = expxy(100000, 1, theta = list(alpha = theta2[2]))
# g5 = kde2d(g4[,1], g4[,2],lims=c(-1, 1, -1, 1) / 2)
# image(g5, main="estimated")
# contour(g5, add=T)

### Estimated lam
mu = theta[1]; K = theta[2]; alpha = theta[3]; beta = theta[4]
lambda = rep(mu/X1/Y1, z$n)
const = K*alpha/pi*beta
for(j in 2:(z$n)){
  gij = 0
  for(i in 1:(j-1)){
    r2 = (z$lon[j]-z$lon[i])^2+(z$lat[j]-z$lat[i])^2
    gij = gij + exp(-beta*(z$t[j]-z$t[i])-alpha*r2)
  }
  lambda[j] = mu / X1 / Y1 + const*gij
}

## Rate vs time
plot(sort(dates[!missingLoc]), ksmooth(times2, lambda, b = 1000)$y, t = "l", xaxt = "n",
     xlab = "", ylab = expression(lambda), main = expression(lambda))
axis(1, at = seq(min(dates), max(dates), length = 25),
     labels = as.character( seq(min(dates), max(dates), length = 25)))

## Rate over space
library(wesanderson)
library(ggplot2)
casesF <- cases[!missingLoc, ]
pal <- wes_palette("Zissou1", 10, type = "continuous")
pal <- adjustcolor(pal, alpha.f = 0.5)

(ggplot() 
  + geom_sf(data = berkeley)
  + geom_point(aes(x = X, y = Y, col = lambda, size = lambda),
               data = data.frame(st_coordinates(coordDf)) %>% mutate(lambda = lambda) 
  )
  + scale_fill_gradientn(colours = pal, aesthetics = "col")
  + xlab("") + ylab("") + ggtitle("") 
)


#############################
########### Super thinning
f = function(t,x,y,z){
  ## compute lambda(t,x,y) given data, z. 
  const = K * alpha/pi * beta
  gij = 0
  j = 0
  if(t > z$t[1]) j = max(c(1:z$n)[z$t < t])
  if(j > 0) {
    for (i in 1:j){
    r2 = (x-z$lon[i])^2+(y-z$lat[i])^2
    gij = gij + exp(-beta*(t-z$t[i])-alpha*r2)
    }
  }
  return(mu / X1 / Y1 + const * gij)
}

supthin = function(z,lambda, f, b = mean(lambda)){
  ## z = data, lambda = conditional intensity at pts, 
  # f = function to compute lambda, 
  ## and b = resulting rate.
  ## First thin, then superpose
  keepz = list()
  for(i in 1:z$n){
    if(runif(1) < b/lambda[i]){
      keepz$t = c(keepz$t,z$t[i])
      keepz$lon = c(keepz$lon,z$lon[i])
      keepz$lat = c(keepz$lat,z$lat[i])
    }
  }
  candn = rpois(1, b * X1 * Y1 * T) # simulate pois process w rate b
  candt = sort(runif(candn) * T)
  candx = runif(candn) * X1
  candy = runif(candn) * Y1
  for(i in 1:candn) { # thin the candidate process
    v = f(candt[i], candx[i], candy[i], z)
    if(v < b){
      if(runif(1) < (b - v) / b){
        keepz$t = c(keepz$t,candt[i])
        keepz$lon = c(keepz$lon,candx[i])
        keepz$lat = c(keepz$lat,candy[i])
      }}
  }
  keepz$lon = keepz$lon[order(keepz$t)]
  keepz$lat = keepz$lat[order(keepz$t)]
  keepz$t = sort(keepz$t)
  keepz$n = length(keepz$t)
  return(keepz)
}

s = supthin(z, lambda, f, b = 20)


plot(s$lon,s$lat,pch=3,cex=.5, xlab="lon", ylab="lat",
     main = "Superthinned Points")

thinsf <- (data.frame(x = s$lon + min(xlon), y = s$lat + min(ylat)) 
            %>% st_as_sf(coords = c("x", "y"), crs = 4326)
            %>% st_transform(crs = 32610))

distMat <- st_distance(thinsf, thinsf)
maxD <- as.numeric(max(distMat))
drange <- seq(0, maxD, length = 100)

## K function
kthin = khat(st_coordinates(thinsf), bdry, drange, checkpoly = TRUE)
plot(drange, kthin/1e5, xlab="distance", ylab = "K(h) * 1e-5", pch = "*", 
     ylim = range(kthin, pi * drange^2)/1e5, main = "K(h)")
lines(drange, kthin / 1e5)
lines(drange, pi * drange^2 / 1e5, lty=2)
legend("topleft", lty = c(1, 2), legend = c("Khat(r)", "Kpois(r)"))


## F function (empty-space function)
b2 = as.ppp(thinsf) 
par(mfrow=c(1,1)) 
f4 = Fest(b2, correction = "best", r = seq(0, maxD, length = 4000)) 
plot(f4, main = "F(h)") 
max(f4$theo)

#### G function: nearest neighbor
g4 <- Gest(b2, r = seq(1, 300, length = 1000)) 
plot(g4$r, g4$km, xlab = "h", ylab = "G(h)", type = "l", lty = 1, main = "G(h)") 
lines(g4$r, g4$theo, lty=2) 
legend("topleft", lty = c(1,2), legend = c("data","Poisson")) 
