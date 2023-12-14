#library(chron)
#library(RColorBrewer)
#library(lattice)
library(ncdf4)
#library(Lmoments)
#library(distillery)
#library(extRemes)
#library(evd)
library(lmomRFA)
library(raster)
#library(rgdal)

#library(oceanmap)

#library(foreach)
#library(doSNOW)
#library(doMC)

#cores <- 6
#registerDoMC(cores)
#cl<-makeCluster(cores, outfile="")
#registerDoSNOW(cl)


#library(doFuture)
#registerDoFuture()
#plan(multiprocess)

################
#This script takes the highest annual values for each model (out of the 1951-1980 historical period) and
#fits the samples to a generalized extreme value distribution (GEVD) using the method of L-moments method
#Regional Frequency Analysis is used
#A total of 150,000 random values are generated with the GEVD. The output is a lon,lat,#ofmodels file
#of the 99th percentile of the 150,000 random values

#loop through files#

years_period <- "2040-2060"
hist_period <-"2000-2020"
percentile_value_file <- "99"
scenario = 'ssp585'


pattern <- paste("/Users/ddusseau/Documents/precip/CMIP6_BASD/annual_max_periods/pr*",scenario,"*",years_period,".nc",sep="")
future_files <- Sys.glob(pattern)

hist_boot <- nc_open(paste("/Users/ddusseau/Documents/precip/CMIP6_BASD/amounts_mm/RFA_pr_",scenario,"_",hist_period,'_99percentile_Lmoments.nc',sep=""))
hist99 <- ncvar_get(hist_boot, paste("pr",percentile_value_file,sep=""))

models <- ncvar_get(hist_boot, "models")

modelnames <- array(NA, length(future_files))
temp <- nc_open(future_files[1])
lat <- ncvar_get(temp,"lat")
lon <- ncvar_get(temp,"lon")
percentile_output <- array(NA, c(dim(lon),dim(lat),length(future_files)))

#add region masks file
mask <- raster("/Users/ddusseau/Documents/precip/CMIP6_BASD/final_pr_regions_01192021_ddedit_02172021_globalgrid.tif") #

regions <- raster::freq(mask)
regions <- regions[-nrow(regions),]

mask_nc <- nc_open("/Users/ddusseau/Documents/precip/CMIP6_BASD/final_pr_regions_01192021_ddedit_02172021_globalgrid.nc") #
mask <- as.matrix(ncvar_get(mask_nc, "Band1"))



#foreach(file=future_files,.packages=c('ncdf4','lmomRFA')) %dopar% {
for (file in future_files) {
  #print(file)
  model <- nc_open(file)
  pr <- ncvar_get(model, "annual_max")
  
  #find the index for the corresponding model in the historical file
  for (p in 1:length(models)) {
    #print(models[p])
    model_name <- unlist(strsplit(file,"_"))[7] #change indexing if necessary
    #print(model_name)
    if (models[p]==model_name){ #change substring index if necessary
      hist_model = p
      modelnames[p] = model_name
      print(model_name)
      print(Sys.time())
      break
    }
  }

  for (z in 1:length(regions[,1])){
    #print(z/length(regions[,1]))
    region_data = matrix(,nrow=dim(pr)[3],ncol=regions[z,"count"])
    index <- 1
    #loop through each gridpoint skip if any pixel isn't in the region
    for (i in 1:dim(lat)) {
      for (j in 1:dim(lon)) {
        if (is.na(mask[j,i])){
          next
        }
        if (mask[j,i] == regions[z,1]){
          histperct <- hist99[j,i,hist_model]
          if (is.na(histperct)){
            next
          }
          region_data[,index] = pr[j,i,]
          index = index+1
        }
      }
    }
    
    if (all(is.na(region_data))){
      next
    }
  
    momen <- regsamlmu(region_data,nmom=5,sort.data = TRUE, lcv = TRUE)
    params <- regfit(momen,"gev")
    quantfunc <- regqfunc(params)
    index_mm <- params$index
    rm(region_data)
    
    floodindex <- 1 #index to loop through floodindex values in fitting output
    for (i in 1:dim(lat)) {
      for (j in 1:dim(lon)) {
        if (is.na(mask[j,i])){
          next
        }
        if (mask[j,i] == regions[z,1]){ #checks if the region mask value for that cell equals the region mask value that the loop is currently going through
          #generate random numbers from the GEV distribution
          histperct <- hist99[j,i,hist_model]
          if (is.na(histperct)){
            next
          }
          
          quantiles <- seq(0.000001,0.999999,0.000001)
          percentile_values <- quantfunc(quantiles)*index_mm[floodindex]
          
          #finds the index where the 99th historical value is closest to a value in the quantile function
          percentile_index <- which(abs(percentile_values-histperct)==min(abs(percentile_values-histperct)))
          percentile_output[j,i,p] <- quantiles[percentile_index]
          floodindex=floodindex+1
        }
      }
    }
  }
  rm(pr)
  
}

dimModel <- ncdim_def(name='model',units='name of model',longname='model',vals=array(1:length(models)))

dimLon <- ncdim_def(name='lon',units='degrees_east',longname='longitude',vals=lon)
dimLat <- ncdim_def(name='lat',units='degrees_north',longname='latitude',vals=lat)

dimnchar <- ncdim_def('nchar', "",1:max(nchar(modelnames)),create_dimvar=FALSE)
dimModels <- ncdim_def(name='nmodels',units='',vals=array(1:length(models)),create_dimvar=FALSE)

var_model <- ncvar_def("models","",list(dimnchar,dimModels),prec="char")

var_pr99 <- ncvar_def(name=paste('pr',percentile_value_file,sep=''),units='percentile',dim=list(dimLon,dimLat,dimModel),longname=paste('Future (',years_period,') Percentile of Historical ',percentile_value_file,'th Percentile',sep=''))
outputfile <- paste('/Users/ddusseau/Documents/precip/CMIP6_BASD/future_RP/RFA_pr_',scenario,'_',years_period,'_vs_',hist_period,'_',percentile_value_file,'percentile_Lmoments.nc',sep='')

ncnew <- nc_create(outputfile, var_pr99)
ncvar_put(ncnew,var_pr99,percentile_output)
  
nc_close(ncnew)





