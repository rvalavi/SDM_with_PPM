## This file constructs all the data frames needed for modeling.
## Note that the result of this code is saved at the end as allData.RData, so there is
##    no need to execute this script if you just want to replicate the modeling in modeling.R

# the data provided in this archive are described in Online Appendix C of Fithian et al. (2014).
# The presence-only species data are sourced from Atlas of Living Australia and Atlas of NSW Wildlife,
# Office of Environment and Heritage (OEH), both publicly available. The presence-absence data were downloaded
# from the Flora Survey Module of the Atlas of NSW Wildlife, Office of Environment and Heritage (OEH),
# and we thank them for permission to archive the data here. Any further use of these data should
# cite Fithian et al. (2014) and acknowledge the data sources.

library(maptools)
library(rgdal)
library(dismo)

wd <- "~/Desktop/biasCorrection"
grid.dir <- "data/grids"

species <- sort(c("eucadalr","eucadeld","eucadive","eucaniph","eucaovat","eucapauc","eucaaggl","eucablax","eucacype","eucafast","eucafrax","eucaobli","eucaobst","eucaparv","eucapilu","eucapipe","eucarobu","eucasieb","eucasten","eucatric","eucaaggr","eucacine","eucadean","eucamolu","eucaorea","eucapunc","eucaquad","eucaross","eucagreg","eucalueh","eucasqua","angobake","coryexim","corymacu","eucaparr", "eucadalh"))

## presence-absence and presence-only data
load(file.path(wd, "moddat.RData"))

ibra <- readShapeSpatial(file.path(wd,"ibraone"))


## construct stack from grid files
if(file.exists("nsw.stack.RData")) {
    load("nsw.stack.RData")
} else {
    vars <- list.files(grid.dir) #this contains data -  I'm sending via Cloudstor
    vars <- vars[vars!="info"] #the info folder is an ARCGIS format thing..
    nsw.stack <- stack(paste(grid.dir, vars, sep="/"))
    nsw.stack$iswater <- 1*is.na(nsw.stack$elev)
    nsw.stack$landNA <- calc(nsw.stack$iswater, function(x) ifelse(x==1,1,NA))
    nsw.stack$fieldNA <- calc(nsw.stack$d2r, function(x) ifelse(is.na(x),1,NA))
    nsw.stack$d2coast <- distance(nsw.stack$landNA)
    save(nsw.stack,file="nsw.stack.RData")
}

# check for different CRS in rasters
for(i in vars){
    r <- raster(paste(grid.dir, i, sep="/"))
    print(paste(i, ":", crs(r)))
    # plot(r)
}

writeRaster(nsw_stack, filename = paste0(names(nsw_stack), ".tif") , bylayer = TRUE, overwrite = TRUE)

vars <- c("bc04", "rsea", "bc33", "bc12", "rjja", "bc02", "bc05", "bc14", "bc21", "bc32", "mvbf", "rugg", "subs", "twmd", "twmx")
nsw.stack <- stack(paste(grid.dir, vars, sep="/"))
crs(nsw.stack) <- "+proj=longlat +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +no_defs"
plot(nsw.stack)

## sample random background points
if(file.exists("rpoints.RData")) {
    load("rpoints.RData")
} else {
    ## using d2r ensures they are only in region of interest - many
    ## of the grids extend west of there..
    set.seed(1)
    rpoints<- randomPoints(nsw.stack[["d2r"]],40000)
    save(rpoints,file="rpoints.RData")
}


## extract data for random background points
if(file.exists("samp.RData")) {
    load("samp.RData")
} else {
    samp <- extract(nsw.stack, rpoints)
    save(samp,file="samp.RData")
}

## get grids of varying fineness, mostly used for plotting in modeling.R
grid.pts <- expand.grid(x=seq(147.58,153.65,by=.01),y=seq(-37.53,-28.14,by=.01))
huge.samp <- extract(nsw.stack,grid.pts)

grid.pts.med <- expand.grid(x=seq(147.58,153.65,by=.05),y=seq(-37.53,-28.14,by=.05))
med.samp <- extract(nsw.stack,grid.pts.med)

grid.pts.coarse <- expand.grid(x=seq(147.58,153.65,by=.1),y=seq(-37.53,-28.14,by=.1))
coarse.samp <- extract(nsw.stack,grid.pts.coarse)

rast.fine <- rasterFromXYZ(cbind(grid.pts,1))
rast.med <- rasterFromXYZ(cbind(grid.pts.med,1))
rast.coarse <- rasterFromXYZ(cbind(grid.pts.coarse,1))

## Get data in format needed for model. Convert things to data frames, add extra derived variables,
PA <- survey_moddat
PA$isPO <- 0
PO.list <- lapply(as.list(species), function(sp) eval(parse(text=paste(sp,"_moddat",sep=""))))
names(PO.list) <- species
background <- as.data.frame(cbind(rpoints,samp))
huge.BG <- as.data.frame(cbind(grid.pts,huge.samp))
med.BG <- as.data.frame(cbind(grid.pts.med,med.samp))

PA$x <- PA$long
PA$y <- PA$lat
names(nsw.stack)
## Get data in format needed for model.
## Note that this is more complex than needed if you were starting from scratch – here we have some variables
##    already sampled (in the files that came as moddat.Rdata) and need to add others that are
##    now in the stack. However we’ve supplied the code since it may be helpful to people wanting
##    examples of use of the package raster.
PA.extr1 <- extract(nsw.stack,cbind(PA$x,PA$y),layer=12,nl=17)
PA.extr2 <- extract(nsw.stack,cbind(PA$x,PA$y),layer=48,nl=1)
PA$d2r <- PA.extr1[,"d2r"]
PA$d2t <- PA.extr1[,"d2t"]
PA$d2coast <- PA.extr2[,"d2coast"]
PA$dist2Sydne <- PA.extr1[,"dist2Sydne"]
PA$dist2Canbe <- PA.extr1[,"dist2Canbe"]
PA$dist2Armid <- PA.extr1[,"dist2Armid"]
PA$dist2Newca <- PA.extr1[,"dist2Newca"]

## coast.trans computes principal components of all the background points' latitude and longitude.
##   Because of the geography of our study region, the first principle component points in a direction that runs
##   along the coastline, from Southwest to Northeast.
coast.trans <- eigen(cov(background[,c("x","y")]))$vectors
construct.derived.vars <- function(dat) {
    z.dat <- as.matrix(dat[,c("x","y")]) %*% coast.trans
    dat$alongCoast <- z.dat[,1]
    dat$ld2coast <- log(1000+dat$d2coast)
    dat$ld2r <- log1p(dat$d2r) # log-transform d2r because it is highly skewed
    dat$ld2Sydne <- log(1000+dat$dist2Sydne)
    dat$ld2Canbe <- log(1000+dat$dist2Canbe)
    dat$ld2Newca <- log(1000+dat$dist2Newca)
    dat$ld2Armid <- log(1000+dat$dist2Armid)
    dat
}

PA <- construct.derived.vars(PA)
PA$isPO <- 0
background <- construct.derived.vars(background)
background$isPO <- 1
med.BG <- construct.derived.vars(med.BG)
med.BG$isPO <- 1
huge.BG <- construct.derived.vars(huge.BG)
huge.BG$isPO <- 1

## Here we are constructing the variable survey.bg: essentially, it is a variable that
##   measures how many presence-absence sites are nearby to a given site s.
##   We use it as a covariate when we model sampling effort for the presence-only data.
pa.site.rast <- rasterize(PA[,c("x","y")],rast.coarse,field=rep(1,nrow(PA)),
                          fun=function(zz,na.rm=TRUE) log1p(sum(zz,na.rm=na.rm)))

pa.site.rast$denom <- rasterize(background[,c("x","y")],rast.coarse,field=rep(1,nrow(background)),
                                fun=function(zz,na.rm=TRUE) log1p(sum(zz,na.rm=na.rm)))
pa.site.rast$value <- calc(pa.site.rast$layer, function(x) ifelse(is.na(x),0,x)) -
    calc(pa.site.rast$denom,function(x) ifelse(is.na(x),0,x))

PA$survey.bg <- extract(pa.site.rast,PA[,c("x","y")],layer=3,nl=1)
background$survey.bg <- extract(pa.site.rast,background[,c("x","y")],layer=3,nl=1)
med.BG$survey.bg <- extract(pa.site.rast,med.BG[,c("x","y")],layer=3,nl=1)
huge.BG$survey.bg <- extract(pa.site.rast,huge.BG[,c("x","y")],layer=3,nl=1)

## Add the extra variables and derived variables for the presence-only data.
for(sp in names(PO.list)) {
    PO.list[[sp]]$x <- PO.list[[sp]]$long
    PO.list[[sp]]$y <- PO.list[[sp]]$lat
    xy <- PO.list[[sp]][,c("x","y")]
    PO.list[[sp]]$d2coast <- extract(nsw.stack,xy,layer=match("d2coast",names(nsw.stack)),nl=1)
    PO.list[[sp]]$dist2Sydne <- extract(nsw.stack,xy,layer=match("dist2Sydne",names(nsw.stack)),nl=1)
    PO.list[[sp]]$dist2Armid <- extract(nsw.stack,xy,layer=match("dist2Armid",names(nsw.stack)),nl=1)
    PO.list[[sp]]$dist2Canbe <- extract(nsw.stack,xy,layer=match("dist2Canbe",names(nsw.stack)),nl=1)
    PO.list[[sp]]$dist2Newca <- extract(nsw.stack,xy,layer=match("dist2Newca",names(nsw.stack)),nl=1)
    PO.list[[sp]] <- construct.derived.vars(PO.list[[sp]])
    PO.list[[sp]]$isPO <- 1
    PO.list[[sp]]$survey.bg <- extract(pa.site.rast,PO.list[[sp]][,c("x","y")],layer=3,nl=1)
}

a <- area(nsw.stack$d2r,na.rm=TRUE) #don’t worry about the warning – our data ARE correct for this.
region.size <- sum(as.data.frame(a),na.rm=TRUE)

towns <- readShapeSpatial(file.path(wd, "/towns/ecologist_towns"))
rast <- aggregate(raster(paste(wd,"/towns1k/bio5clip", sep="")),fact=10)

## Save the entire workspace, for use in modeling.R
save.image("allData.RData")

