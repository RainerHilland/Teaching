######
# Procrastination in the form of solar modelling
######

library(tidyr)
library(ggplot2)

#####
# For each day we
# 1. calculate fractional year and solar declination
# 2. calculate solar constant due to orbital effects
# 3. For each time step: calculate solar zenith (elevation) and azimuth
# 4. Calculate irradiance on a flat surface or slope
# 5. keep track of our data for plotting/analysis later

#####
# because we aren't concerned about specific times we assume we're always on local solar time
# and just count minutes in to the day

##### 
# The big list of assumptions:
# 1. atmospheric transmission = 1, weather = 0
# 2. perfectly horizontal surface with
# 3. no horizon obstructions or anything like that
# 4. no atmosphere (no diffuse solar rad)

#####
# TODOs
# 2. implement the slope functions


#####
# User parameters
#####

# Site and model specifications
latitude <- 48 # site latitude in degrees. 48 ~= Freiburg (positive N, negative S)
timestep <- 30 # model timestep in minutes
Io <- 1361 # solar constant
startday <- 1 # YD to begin model
ndays <- 365 # 1 year

#####


#####
# FUNCTIONS
#####

# simple functions
radToDeg <- function(r) {return(r * (180/pi))}
degToRad <- function(deg) {return(deg*(pi/180))}


# sun location
fractionalYear <- function(yd) {return((2*pi)/365 * (yd-1))} # gamma, returned in RADIANS

solarDeclination <- function(gamma) { # delta / solar declination
  delta <- 0.006918 - 0.399912*cos(gamma) + 0.070257*sin(gamma) - 
    0.006758*cos(2*gamma) + 0.000907*sin(2*gamma) - 
    0.002697*cos(3*gamma) + 0.00148*sin(3*gamma)
  return(delta) # returned in radians!
}

EqnofTime <- function(gamma) { # equation of time IN MINUTES
  dTlat <- 229.18 * (0.000075 + 0.001868*cos(gamma) -
                       0.032077*sin(gamma) - 0.014615*cos(2*gamma)-
                       0.040849*sin(2*gamma))
  return(dTlat) # returned in MINUTES!
}

hourAngle <- function(minutes) { # simplified hour angle calculation assuming we're always in local solar time
  return(degToRad((1*(minutes-720))/4)) # returned IN RADIANS
}

solarElevation <- function(latitudeRadians, declination, hourAngle) { # determine solar elevation
  # assumes all input variables are already in radians
  elevation <- asin(sin(declination)*sin(latitudeRadians)+cos(declination)*cos(latitudeRadians)*cos(hourAngle))
  return(radToDeg(elevation)) # returned IN DEGREES
}

solarZenith <- function(elevation) { # returns solar zenith in degrees, requires elevation be passed in degrees
  return(90-elevation)
}

solarAzimuth <- function(latitudeRadians, declination, hourAngle, elevationRadians) { # returns solar azimuth
  # assumes all input variables are already in radians
  
  oldw <- getOption("warn") # suppress warnings here because something's not working 100% right
  options(warn = -1)
  
  azimuth <- acos( (sin(declination)*cos(latitudeRadians)-cos(declination)*sin(latitudeRadians)*cos(hourAngle))/cos(elevationRadians) )

  
  if (hourAngle > 0) { # past noon we do 360-azimuth
    az <- 360-radToDeg(azimuth)
  } else {
    az <- radToDeg(azimuth) # returned IN DEGREES
  }
  
  if (is.nan(az)) { # handling errors I don't understand
    az <- 180
  }
  
  options(warn=oldw) # turn warnings back on
  
  return(az)
}

slopeGeometry <- function(slopeAngleRad, solarElevRad, solarAzRad, slopeAzRad) {# implementing slopes
  # all input variables need to be in radians, lads
  relativeAngle <- acos( cos(slopeAngleRad)*sin(solarElevRad) + sin(slopeAngleRad)*cos(solarElevRad)*cos(solarAzRad-slopeAzRad) )

  return(radToDeg(relativeAngle))
}

# radiation stuff
solarConst <- function(gamma, solarConst) { # determine variation in solar constant due to orbit
  # equates the mean to real orbit using the inverse square law
  ratio <- 1.00011 + 0.034221*cos(gamma) + 0.001280*sin(gamma) +
    0.000719*cos(2*gamma) + 0.000077*sin(2*gamma)
  insol <- solarConst * ratio
  return(insol)
}

cosineLaw <- function(z, KdownTOA) { # implement cosine law of illumination
  # takes z as DEGREES, returns W/m2 at horizontal surface!
  return(KdownTOA*cos(degToRad(z)))
}

cosineLawSlope <- function(zRad, kdownTOA, slopeGeomRad) {
  
  return((kdownTOA/cos(zRad))*cos(slopeGeomRad)) # returns W/m2
  
}
#####

#####
# MODEL
#####

# data storage
# i, day, minutes, gamma, declination, kTOA, hourangle, solarElev, solarAz, solarZenith, KHoriz, JoulesHoriz
bigFrame <- matrix(nrow=((1440/timestep)+1)*ndays, ncol=12)

# day, gamma, declination, kTOA, maxSolarElev, meanSolarElev, meanKhoriz, cumEhoriz
smallFrame <- matrix(nrow=ndays, ncol=8)


latitudeRadians <- degToRad(latitude) # for simplicity
i <- 0 # tracks iteration

# model start
for (day in startday:ndays) { # day loop (for leap year need to modify gamma too)
  
  gamma <- fractionalYear(day) # calculate gamma (fractional year)
  declination <- solarDeclination(gamma) # calculate solar declination for the day
  kDownTOA <- solarConst(gamma, Io) # determine Kdown at Tropopause due to orbital variations
  
  minutes <- 0 # start the day
  
  while (minutes <= 1440) {
    
    i <- i + 1
    
    # solar position calculations
    hourAngleRads <- hourAngle(minutes) # get hour angle
    solarElevDeg <- solarElevation(latitudeRadians, declination, hourAngleRads) # elevation
    
    if (solarElevDeg > 0) {
      
      solarAzDeg <- solarAzimuth(latitudeRadians, declination, hourAngleRads, degToRad(solarElevDeg)) # get azimuth
      solarZenDeg <- solarZenith(solarElevDeg) # get Zenith
      
      Khoriz <- cosineLaw(solarZenDeg, kDownTOA) # determine incident K at a horizontal surface (W/m2)
      Jhoriz <- Khoriz * 60 * timestep
      
      
      
    } else {
      solarAzDeg <- NaN
      solarZenDeg <- NaN
      Khoriz <- 0
      Jhoriz <- 0
    }
    
    bigFrame[i,] <- c(i, day, minutes, gamma, radToDeg(declination), kDownTOA, radToDeg(hourAngleRads), solarElevDeg, solarAzDeg, solarZenDeg, Khoriz, Jhoriz)
    
    minutes <- minutes + timestep
  }
  
  # day, gamma, declination, kTOA, maxSolarElev, meanSolarElev, meanKsurf, cumJoules
  dayData <- bigFrame[which(bigFrame[,2]==day),]

  smallFrame[day,] <- c(day, gamma, radToDeg(declination), kDownTOA, max(dayData[,8]), mean(dayData[,8]), mean(dayData[,11]), sum(dayData[,12]))
  
} # model end
#####

# cleaning the data
bigFrame <- as.data.frame(bigFrame)
colnames(bigFrame) <- c("i","day","minute","fractionalYear_rad", "solarDeclination_deg", "KdownTOA_Wm2", "hourAngle_deg", "solarElevation_deg",
                        "solarAzimuth_deg", "solarZenith_deg", "Ksurface_Wm2")

smallFrame <- as.data.frame(smallFrame)
colnames(smallFrame) <- c("day", "fractionalYear_rad", "solarDeclination_deg", "KdownTOA_Wm2", "maxSolarElevation_deg", "meanSolarElevation_deg", "meanKsurface_Wm2", "cumulativeEnergy_J")

# plot
ggplot(smallFrame, aes(x=day, y=cumulativeEnergy_J)) + geom_line() + theme_linedraw()

