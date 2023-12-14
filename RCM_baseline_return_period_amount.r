library(chron)
library(RColorBrewer)
library(lattice)
library(ncdf4)
library(Lmoments)
library(distillery)
#library(extRemes)
library(evd)
library(lmomRFA)
library(raster)
library(rgdal)
#library(oceanmap)

################
#This script takes the highest annual values for each model (out of the 1951-1980 historical period) and
#fits the samples to a generalized extreme value distribution (GEVD) using the method of L-moments method
#Regional Frequency Analysis is used
#A total of 150,000 random values are generated with the GEVD. The output is a lon,lat,#ofmodels file
#of the 99th percentile of the 150,000 random values


#loop through files#
years_period <- "2070-2090"
scenario = 'ssp585'
out_percentile <- 0.99

all_files <- Sys.glob(paste("/Users/ddusseau/Documents/precip/CMIP6_BASD/annual_max_periods/pr*",scenario,"*",years_period,".nc",sep=""))


modelnames <- array(NA, length(all_files))
temp <- nc_open(all_files[1])
lat <- ncvar_get(temp,"lat")
lon <- ncvar_get(temp,"lon")
percent_99 <- array(NA, c(dim(lon),dim(lat),length(all_files)))


#add region masks file
mask <- raster("/Users/ddusseau/Documents/precip/CMIP6_BASD/final_pr_regions_01192021_ddedit_02172021_globalgrid.tif") #

regions <- raster::freq(mask)
regions <- regions[-nrow(regions),]

mask_nc <- nc_open("/Users/ddusseau/Documents/precip/CMIP6_BASD/final_pr_regions_01192021_ddedit_02172021_globalgrid.nc") #
mask <- as.matrix(ncvar_get(mask_nc, "Band1"))

k=1
for (file in all_files) {
  model <- nc_open(file)
  #get model names
  model_name <- unlist(strsplit(file,"_"))[7] #change indexing if necessary
  print(model_name)
  modelnames[k] = model_name

  pr <- ncvar_get(model, "annual_max")


  for (p in 1:length(regions[,"value"])){
    region_data <- matrix(,nrow=dim(pr)[3],ncol=regions[p,"count"]) ##create matrix to hold data series for each pixel in region
    index <- 1
    #loop through each gridpoint skip if any pixel isn't in the region
    for (i in 1:dim(lat)) {
      for (j in 1:dim(lon)) {
        if (is.na(mask[j,i])){
          next
        }
        
        if (mask[j,i] == regions[p,"value"]){ ##check if pixel is in region 
#            if (any(is.na(pr[j,i,]))){
#             next
#           } else {
          region_data[,index] = pr[j,i,] #append data series for pixel
          index = index+1
#           }
        }
      }
    }
    
    if (all(is.na(region_data))){
      next
    }
    momen <- regsamlmu(region_data,nmom=5,sort.data = TRUE, lcv = TRUE) #calculate the lmoments for each site
    params <- regfit(momen,"gev") #generates the GEV regional distribution
    quantfunc <- regqfunc(params) #generates a function to compute the quantiles
    index_mm <- params$index #indexes for the site-specific scale factors. These are multipliers for the regional frequency distribution
    rm(region_data)

    floodindex<-1 #index to loop through floodindex values in fitting output
    for (i in 1:dim(lat)) {
      for (j in 1:dim(lon)) {
        if (is.na(mask[j,i])){
          next
        }
        if (mask[j,i] == regions[p,"value"]){
          percentile_value <- quantfunc(out_percentile)*index_mm[floodindex] #set which percentile you want, 0.98 is 1-in-50  ##returns value at percentile, multiplied by scale factor
          percent_99[j,i,k] <- percentile_value
          floodindex=floodindex+1
        }
      }
    }
  }
  k=k+1
  rm(pr)
}

lat <- ncvar_get(model,"lat")
lon <- ncvar_get(model,"lon")

dimModel <- ncdim_def(name='model',units='name of model',longname='model',vals=array(1:(k-1)))

dimLon <- ncdim_def(name='lon',units='degrees_east',longname='longitude',vals=lon)
dimLat <- ncdim_def(name='lat',units='degrees_north',longname='latitude',vals=lat)

dimnchar <- ncdim_def('nchar', "",1:max(nchar(modelnames)),create_dimvar=FALSE)
dimModels <- ncdim_def(name='nmodels',units='',vals=array(1:(k-1)),create_dimvar=FALSE)

var_model <- ncvar_def("models","",list(dimnchar,dimModels),prec="char")

var_pr99 <- ncvar_def(name='pr99',units='mm/day',dim=list(dimLon,dimLat,dimModel),longname='99th percentile of bootstrapped values')

outputfile <- paste('/Users/ddusseau/Documents/precip/CMIP6_BASD/amounts_mm/RFA_pr_',scenario,"_",years_period,'_99percentile_Lmoments.nc',sep="") #make sure to change output file name
ncnew <- nc_create(outputfile, list(var_pr99,var_model))

ncvar_put(ncnew,var_model,modelnames)
ncvar_put(ncnew,var_pr99,percent_99)

nc_close(ncnew)






