# This code is example of analysing a dataset composed of several point pattern
# https://www.google.com/url?sa=t&source=web&rct=j&opi=89978449&url=https://cran.r-project.org/web/packages/spatstat/vignettes/replicated.pdf&ved=2ahUKEwiQl7zoz72TAxV3U6QEHfvnCy0QFnoECBgQAQ&usg=AOvVaw1UpcanKOmgMRX7Z7IAAqtA

#############
#load packages
#############
library(here)
library(sf)
library(spatstat)
library(spatstat.model)
library(dplyr)
library(ggplot2)
library(ggspatial)
library(fields)

#############
#load data
#############
demo <- demohyper

list_PPP <- readRDS(paste(
  here(),"/via3_data_exploration/Data/processed/list_PPP_tuyau_2025.rds",
  sep=""))
PPPs <- hyperframe(Pattern=list_PPP)

########################
# Start with example
########################
plot(demo)
H <- hyperframe(Bugs=waterstriders)
plot(H, quote(plot(Kest(Bugs))), marsize=1)
with(H, npoints(Bugs)) # number of points
with(H, distmap(Bugs)) 
K <- with(H, Kest(Bugs)) # K-functions
plot(K)
with(H, nndist(Bugs)) # nearest neighbour distances
with(H, min(nndist(Bugs))) # minimum interpoint distance in each pattern
# For a new window (on Windows)
windows()
plot(H, quote(plot(density(Bugs), main="")), nrows=2)

rhos <- with(demohyper, rhohat(Points, Image))# Covariate effects
windows()
plot(rhos)

# Fitting models of spatial trend
mppm(Bugs~x+y,H,Poisson())
mppm(Bugs~x+y,H)

mppm(Points ~ Image, data=demohyper)
mppm(Points ~ offset(log(Image)), data=demohyper)
mppm(Points ~ log(Image), data=demohyper)

#Interaction depending on sampling design
fit <- mppm(Points ~ 1, simba, Strauss(0.07), iformula = ~Interaction*group)
fit
coef(fit)

# Fitting a mixed-effects model
mppm(Bugs ~ 1, H, random=~1|id)
# Studying the fitted model
subfits(fit) # Fits for each pattern
with(H, ppm(Bugs))# Fitting separately to each pattern for comparison
#check residual
fit <- mppm(Bugs ~ x, H)
res <- residuals(fit)
windows();plot(res)
smor <- with(hyperframe(res=res), Smooth(res, sigma=4))
windows();plot(smor) #smoothed residual field
totres <- sapply(res, integral.msr) # Sums of residuals

#In designed experiments we can plot these total residuals against the design covariates:
fit <- mppm(Points~Image, data=demohyper)
resids <- residuals(fit, type="Pearson")
totres <- sapply(resids, integral.msr)
areas <- with(demohyper, area.owin(as.owin(Points)))
df <- as.data.frame(demohyper[, "Group"])
df$resids <- totres/areas
plot(resids~Group, df)

# Four-panel diagnostic plots
fit <- mppm(P ~ 1, hyperframe(P=waterstriders))
sub <- hyperframe(Model=subfits(fit))
plot(sub, quote(diagnose.ppm(Model)))

# Residuals of the parameter estimates
H <- hyperframe(P = waterstriders)
fitall <- mppm(P ~ 1, H)
together <- subfits(fitall)
separate <- with(H, ppm(P))
Fits <- hyperframe(Together=together, Separate=separate)
dr <- with(Fits, unlist(coef(Separate)) - unlist(coef(Together)))
dr # objective is to minimize dr values
exp(dr)

# Goodness-of-fit tests
#quadrat count test for a fitted Poisso, point process only
H <- hyperframe(X=waterstriders)
# Poisson with constant intensity for all patterns
fit1 <- mppm(X~1, H)
quadrat.test(fit1, nx=2)
# uniform Poisson with different intensity for each pattern
fit2 <- mppm(X ~ id, H)
quadrat.test(fit2, nx=2)
anova(fit1, fit2, test="LR")


########################
# Test on real data
########################
#windows()
#plot(PPPs[c(1:5),])
windows()
plot(PPPs, quote(plot(Kest(Pattern))), marsize=1)
with(PPPs, npoints(Pattern))
with(PPPs, distmap(Pattern))
with(PPPs, Kest(Pattern))
with(PPPs, nndist(Pattern))
windows()
plot(PPPs, quote(plot(density(Pattern), main="")), nrows=8)
