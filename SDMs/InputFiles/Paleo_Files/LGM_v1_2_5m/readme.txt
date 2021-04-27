Paleoclim.org Bioclmate Variables
________________________________________________________________________________________________________________________________
Layers:CHELSA Current Bioclims, Version 1.2, 5/29/2018


By downloading the data you agree to cite the following peer reviewed article:

*Karger, D.N., Conrad, O., Böhner, J., Kawohl, T., Kreft, H., Soria-Auza, R.W., Zimmermann, N.E., Linder, H.P. & Kessler, M. (2017) Climatologies at high resolution for the earth’s land surface areas. Scientific Data 4, 170122.



Distributed, downsclaed and renamed by PaleoClim.org
             ___.....___
       ,..-.=--.-.       "".
     .{_..        `        ,. .
   .'     \      /        | ,'.\`.
  /        :   ;'          `____> \
 :         `. (           /       :
 |           `>\_         \      r|
             /   \         `._   \
 |          |      `          ;   |
  :          \     /          `   ;
   \          \.  '            ` /
     `.        | /             .'
        `      `/          . '
           `---'.....---''

    +-+-+-+-+-+-+-+-+-+-+-+-+-+
    |P|a|l|e|o|c|l|i|m|.|o|r|g|
    +-+-+-+-+-+-+-+-+-+-+-+-+-+



CHELSA README

CHELSA – Climatologies at high resolution for the Earth land surface areas. Version 1.2


CHELSA (http://chelsa-climate.org/) is a high resolution (30 arc sec, ~1 km) climate data set for the earth land surface areas. It includes monthly and annual mean temperature and precipitation patterns for the time period 1979-2013. CHELSA_v1 is based on a quasi-mechanistical statistical downscaling of the ERA interim global circulation model
(http://www.ecmwf.int/en/research/climate-reanalysis/era-interim) with a GPCC (https://www.dwd.de/EN/ourservices/gpcc/gpcc.html) and GHCN (https://www.ncdc.noaa.gov/ghcnm/) bias correction.

CHELSA Version 1.2 is licensed under a Creative Commons Attribution 4.0 International License.

Specifications:
High resolution (30 arcsec, ~1 km)
Precipitation & Temperature
Climatologies for the years 1979 – 2013
Incorporation of topoclimate (e.g. orographic rainfall & wind fields). 

All products of CHELSA are in a geographic coordinate system referenced to the WGS 84 horizontal datum, with the horizontal coordinates expressed in decimal degrees. The CHELSA layer extents (minimum and maximum latitude and longitude) are a result of the coordinate system inherited from the 1-arc-second GMTED2010 data which itself inherited the grid extent from the 1-arc-second SRTM data.

Note that because of the pixel center referencing of the input GMTED2010 data the full extent of each CHELSA grid as defined by the outside  edges of the pixels differs from an integer value of latitude or  longitude by 0.000138888888 degree (or 1/2 arc-second). Users of products based on the legacy GTOPO30 product should note that the coordinate referencing of CHELSA (and GMTED2010) and GTOPO30 are not the same. In GTOPO30, the integer lines of latitude and longitude fall directly on the edges of a 30-arc-second pixel. Thus, when overlaying CHELSA with products based on GTOPO30 a slight shift of 1/2 arc-second will be observed between the edges of corresponding 30-arc-second pixels.

The dataset is in GEOtiff format. GEOtiff can be viewed using standard GIS software such as:

Naming convention:

CHELSA_<variable>_<month>_<Version>_land.tif

variable:	prec 	= precipitation [mm/month]
		temp	= monthly mean of daily mean temperature [°C*10]
		tmax	= monthly mean of daily maximum temperature [°C*10]
		tmin	= monthly mean of daily minimum temperature [°C*10]
		bio1	= Annual Mean Temperature [°C*10]
		bio2	= Mean Diurnal Range [°C]
		bio3	= Isothermality
		bio4	= Temperature Seasonality [standard deviation]
		bio5	= Max Temperature of Warmest Month [°C*10]
		bio6	= Min Temperature of Coldest Month [°C*10]
		bio7	= Temperature Annual Range [°C*10]
		bio8	= Mean Temperature of Wettest Quarter [°C*10]
		bio9	= Mean Temperature of Driest Quarter [°C*10]
		bio10	= Mean Temperature of Warmest Quarter [°C*10]
		bio11	= Mean Temperature of Coldest Quarter [°C*10]
		bio12	= Annual Precipitation [mm/year]
		bio13	= Precipitation of Wettest Month [mm/month]
		bio14	= Precipitation of Driest Month [mm/month]
		bio15	= Precipitation Seasonality [coefficient of variation]
		bio16	= Precipitation of Wettest Quarter [mm/quarter]
		bio17	= Precipitation of Driest Quarter [mm/quarter]
		bio18	= Precipitation of Warmest Quarter [mm/quarter]
		bio19 	= Precipitation of Coldest Quarter [mm/quarter]
