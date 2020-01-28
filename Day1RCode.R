# Solar Radiation
# Day 1

library(ggplot2) # plotting package
library(tidyr)

#####
# Lecture 1 - Rad Basics
#####

# 1. calculate frequency with given wavelength
# f = c/lambda

lambda <- 0.55 * 10^(-6) # enter wavelength in um
c <- 3*10^8 # approximation to speed of light
f <- c/lambda # calculate frequency
print(f)

# calculate energy of a photon
# e = hf

h <- 6.63 * 10^(-34) # Planck's constant
e <- h*f # calculate energy per photon
print(e)

# how does the frequency, wavelength, and energy per photon change when
# we move to the infrared (e.g. 4 um) waveband? How about the
# ultraviolet (e.g. 0.1 um)?

# 2. calculate total energy reaching the surface of 
#  Freiburg over a year

avgIrradiance <- 130 # W/m2
sPerYear <- 60*60*24*365 # number of seconds per year
joules <- avgIrradiance * sPerYear # joules at Freiburg surface

GJ <- joules / 1000000000 # convert to gigajoules
print(GJ)

kWh <- joules / (3600*1000) # convert to kWh
print(kWh)

# The flux density of solar radiation at the edge of Earth's atmosphere
# normal to the equator is ~1361 W/m2. How many GJ is this over a year?
# How many kWh?

# 3. Relate the sun's surface temperature (5778 K) to its wavelength of
#  maximum emission

b <- 2.8978e-3 # Constant
TSun <- 5778 # Kelvin
lambdaMax <- b/TSun
lambdaMax <- lambdaMax * 1000000 # convert to um
print(lambdaMax)

# Seen from space Earth's radiative (effective) temperature is roughly
#  255 K (-18 deg.C). What is Earth's wavelength of max emission?


######
# Lecture 2 -> Solar Position
######

# 1. Calculating solar declination

# function for calculating the solar declination
getSolarDec <- function(yd) {
  gamma <- (2*pi)/365 * (yd-1)
  delta <- 0.006918 - 0.399912*cos(gamma) + 0.070257*sin(gamma) - 
    0.006758*cos(2*gamma) + 0.000907*sin(2*gamma) - 
    0.002697*cos(3*gamma) + 0.00148*sin(3*gamma)
  return(delta)
}

# convert radians to degrees for coherence and vice-verse
radToDeg <- function(r) {return(r * (180/pi))}
degToRad <- function(deg) {return(deg*(pi/180))}

# calculate the solar declination on, e.g., your birthday
yd <- 55 # adjust to your birthday

delta <- getSolarDec(yd)
print(radToDeg(delta))

# how does solar declination vary through the year?
# -> perform the above steps for each day of the year
# -> plot the results

ydVector <- seq(1,365,1) # list of integer days
deltaVector <- radToDeg(getSolarDec(ydVector)) # corresponding declinations

# put in a plottable form and plot
plotFrame <- as.data.frame(cbind(ydVector,deltaVector)) 
ggplot(plotFrame, aes(x=ydVector, y=deltaVector)) + geom_line(col='blue') + theme_linedraw() +
  labs(x="YD", y="Solar Declination (degrees)") + geom_hline(yintercept=0) +
  geom_hline(yintercept=23.5) + geom_hline(yintercept=-23.5) +
  annotate("text", x=50, y=20, label="Max. ~23.5 deg") +
  annotate("text", x=50, y=-20, label="Min. ~ -23.5 deg")

# Note where the extreme values occur -> do not line up with calendar year!
# What does a comparison to the simplified calculation look like?


# 2. It's 14:30 on March 23 2021 in Freiburg. Where's the sun?
# Functions transcribed from NOAA Low Accuracy Equations

# functions
getGamma <- function(yd) {return((2*pi)/365 * (yd-1))} # returned in radians!

getDelta <- function(gamma) {
  delta <- 0.006918 - 0.399912*cos(gamma) + 0.070257*sin(gamma) - 
    0.006758*cos(2*gamma) + 0.000907*sin(2*gamma) - 
    0.002697*cos(3*gamma) + 0.00148*sin(3*gamma)
  return(delta) # returned in radians!
}

getEqnofTime <- function(gamma) {
  dTlat <- 229.18 * (0.000075 + 0.001868*cos(gamma) -
                        0.032077*sin(gamma) - 0.014615*cos(2*gamma)-
                        0.040849*sin(2*gamma))
  return(dTlat) # returned in MINUTES!
}

getTOffset <- function(dTlat, lambda, epsilon) {
  To <- dTlat + 4*lambda - 60*epsilon
  return(To) # returned in MINUTES!
}

getTlat <- function(lst, To) {
  H <- as.numeric(unlist(strsplit(lst,"\\:"))[1])
  M <- as.numeric(unlist(strsplit(lst,"\\:"))[2])
  S <- as.numeric(unlist(strsplit(lst,"\\:"))[3])
  Tsol <- H*60 + M + S/60 + To # returned in MINUTES!
}

getHourAngle <- function(Tsol) {
  return((Tsol/4)-180) # returned in DEGREES!
}

getZenith <- function(phi, delta, h) {
  phi <- degToRad(phi)
  h <- degToRad(h)
  cosZ <- sin(phi)*sin(delta)+cos(phi)*cos(delta)*cos(h)
  zenith <- radToDeg(acos(cosZ))
  return(zenith) # returned in DEGREES!
}

getAzimuth <- function(phi, Z, delta) {
  # I never got to this - give it a go!
  }

# these values are looked up / determined from Google
yd <- 82 # year day
lambda <- 7.8 # longitude
phi <- 48 # latitude
epsilon <- 1
lst <- "12:30:00"

# these values are calculated as intermediate steps
gamma <- getGamma(yd) # radians
delta <- getDelta(gamma) # radians
dTlat <- getEqnofTime(gamma) # minutes
To <- getTOffset(dTlat, lambda, epsilon) # minutes
Tlat <- getTlat(lst, To) # minutes
h <- getHourAngle(Tlat) # degrees

# determine solar angles
zenithAngle <- getZenith(phi, delta, h) # in deg.
azimuthAngle <- getAzimuth(phi, Z, delta)

######
# Lecture 3
# Extraterrestrial Solar Rad
######

# 1. How much radiation is emitted by the Sun?
sunTemp <- 5778 # Kelvin
sigma <- 5.67e-8 # S-B Constant
totalEmit <- sigma * sunTemp ** 4 # calculate integral emittance
print(totalEmit)

# The star 63 Ophiuchi has a measured surface temperature of approximately 
# 34,000 K -> what is the radiation flux density at the surface of this star?

# 2. How much of the solar radiation reaches the tropopause?
# ensure totalEmit is still calculated for temp of 5778 K

r1 <- 694250 # radius of the Sun in km
r2 <- 149600000 # ~mean radius of Earth's orbit
E2 <- totalEmit * (r1/r2)**2
print(E2)

# Mercury has a mean distance of 58,000,000 km from the Sun. What is its solar constant?

# 3. Determine the variation in insolation caused by Earth's eccentric orbit

# functions
getGamma <- function(yd) {return((2*pi)/365 * (yd-1))} # returned in radians!
orbitRatio <- function(gamma, solarConst) {
  # equates the mean to real orbit using the inverse square law
  ratio <- 1.00011 + 0.034221*cos(gamma) + 0.001280*sin(gamma) +
    0.000719*cos(2*gamma) + 0.000077*sin(2*gamma)
  insol <- solarConst * ratio
  return(insol)
}

# inputs
solarConst <- 1361.012
days <- seq(1,365,1) # simple list of days
gammaDays <- getGamma(days) # calculate fractional year
solarRad <- orbitRatio(gammaDays, solarConst) # apply the inverse square law to each fractional day

# for plotting
insolation <- as.data.frame(cbind(days,solarRad))

# determine perihelion and aphelion as days w/ max and min insolation
perihelion <- insolation$days[insolation$solarRad == max(insolation$solarRad)]
aphelion <- insolation$days[insolation$solarRad == min(insolation$solarRad)]
  
# plot the results
ggplot(insolation, aes(x=days, y=solarRad)) + geom_line(col="red") + theme_linedraw() +
  labs(x="YD", y="Solar constant (W/m2)") + geom_vline(xintercept=perihelion) + geom_vline(xintercept=aphelion) +
  geom_hline(yintercept=mean(insolation$solarRad), col="green") +
  annotate("text", x=15, y=insolation$solarRad[perihelion]-7, label="Perihelion") +
  annotate("text", x=175, y=insolation$solarRad[aphelion]+7, label="Aphelion") +
  annotate("text", x=205, y=mean(insolation$solarRad)+5, label="Mean ~ 1361", col="green")

# 4. Examine how radiation is distributed between the poles and equator
# Assume: mean insolation at solar noon on an equinox

# functions
degToRad <- function(deg) {return(deg*(pi/180))}

cosineLaw <- function(z, solarConst) { # implement cosine law of illumination
  return(solarConst*cos(degToRad(z)))
}

latitudes <- seq(-90,90,1) # list of latitudes
solarConst <- 1361 # in W/m2
insolationDistribution <- cosineLaw(latitudes, solarConst) # apply cosine law to each latitude

# plot
distribution <- as.data.frame(cbind(latitudes,insolationDistribution))
ggplot(distribution, aes(x=latitudes, y=insolationDistribution)) + geom_line() + theme_linedraw() +
  labs(x="Latitude (deg)", y="Insolation (W/m2)")

# How would this distribution change at the perihelion and aphelion?

# Get the solar declination at perihelion and aphelion
perihelionDeclination <- radToDeg(getDelta(getGamma(perihelion)))
aphelionDeclination <- radToDeg(getDelta(getGamma(aphelion)))

# get the insolation at tropopause for perihelion and aphelion
perihelionKParallel <- max(insolation$solarRad)
aphelionKParallel <- min(insolation$solarRad)

# calculate cosine law
perihelionDistribution <- cosineLaw(latitudes-perihelionDeclination, perihelionKParallel)
aphelionDistribution <- cosineLaw(latitudes-aphelionDeclination, aphelionKParallel)

perihelionDistribution[perihelionDistribution < 0] <- NaN # handle negative values
aphelionDistribution[aphelionDistribution < 0] <- NaN

# organising
equinoxDistribution <- insolationDistribution
distributionTwo <- as.data.frame(cbind(latitudes,equinoxDistribution,perihelionDistribution,aphelionDistribution))
distributionTwo <- distributionTwo %>% gather("OrbitPeriod", "FluxDensity", -latitudes)

# plot
ggplot(distributionTwo, aes(x=latitudes, y=FluxDensity, col=OrbitPeriod)) + geom_line() + theme_linedraw() +
  labs(x="Latitude", y="Insolation (W/m2)") +
  annotate("text", x=0, y=550, label="Perihelion ~Jan 5", col="blue") +
  annotate("text", x=0, y=500, label="Equinox ~March 20 and September 22", col="green") +
  annotate("text", x=0, y=450, label="Aphelion ~July 4", col="red")