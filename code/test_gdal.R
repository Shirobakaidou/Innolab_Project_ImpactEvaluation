library(sf)
library(terra)
library(raster)
#library(rgdal)
library(gdalUtilities)
#library(vapour)
library(geodata)
library(wdpar)
library(tidyr)

# Ref:
# About gdal virtual file system:
# https://gis.stackexchange.com/questions/266838/gdal-warp-to-in-memory-raster
# https://gis.stackexchange.com/questions/397451/why-does-vsimem-have-a-path
# Read vrt files stored online using gdal virtual file system:
# https://gis.stackexchange.com/questions/361086/reading-vrt-from-soilgrids-in-r

country = "Madagascar"
wdir = file.path(getwd(), 'data')

# Vector: Reproject GADM ####
gadm <- gadm(country = country, resolution = 1, level = 0, path = wdir) %>%
  st_as_sf()

# Find UTM zone of the country centroid
centroid = st_coordinates(st_centroid(gadm))
lonlat2UTM = function(lonlat) {
  utm = (floor((lonlat[1] + 180) / 6) %% 60) + 1
  if (lonlat[2] > 0) {
    utm + 32600
  } else{
    utm + 32700
  }
}
utm_code = lonlat2UTM(centroid)

# Reproject GADM
gadm_prj = gadm %>% st_transform(crs = utm_code)

# Export projected GADM for use in GDAL
#st_write(gadm_prj, 
#         dsn = file.path(getwd(), "country_prj.gpkg"), delete_dsn = TRUE) 


# Vector: Groups (funded PA, non-funded PA, buffer, control candidates) ####

# Get the PA polygons/points of the specified country; 
# They're downloaded to the working directory.
wdpa = wdpa_fetch(country, wait = TRUE, download_dir = wdir)
# If the PA file already exists, it can be loaded in this way
wdpa = wdpa_read(paste0(wdir, '/WDPA_Dec2022_MDG-shapefile.zip'))


# PAs are projected, and column "geometry_type" is added
wdpa_prj = wdpa_clean(wdpa, geometry_precision = 1000) %>%
  # Remove the PAs that are only proposed, or have geometry type "point"
  filter(STATUS != "Proposed") %>%
  filter(GEOMETRY_TYPE != "POINT") %>%
  # Project PA polygons to the previously determined UTM zone
  st_transform(crs = utm_code) 


# Make Buffers of 10km around all protected areas
buffer <- st_buffer(wdpa_prj, dist = buffer_m) %>% 
  # Assign an ID "3" to the buffer group
  mutate(group=3)


# Separate funded and non-funded protected areas
wdpaID_funded = paid
wdpa_funded = wdpa_prj %>% filter(WDPAID %in% wdpaID_funded) %>%
  mutate(group=1) # Assign an ID "1" to the funded PA group
wdpa_nofund = wdpa_prj %>% filter(!WDPAID %in% wdpaID_funded) %>% 
  mutate(group=2) # Assign an ID "2" to the non-funded PA group


# Merge the dataframes of funded PAs, non-funded PAs and buffers
wdpa_groups = rbind(wdpa_funded, wdpa_nofund, buffer)


# Co-var: Soilgrids
files = list.files(path = file.path(getwd(), "GEE"), 
                   full.names = TRUE) %>% head(., -1)

fInfo = gdal_utils(util = "info",
                   source = files[1])

gdalwarp(srcfile = files[1],
         dstfile = file.path(getwd(), "GDAL_R", "warpSum.tiff"), # Output file in GDAL virtual file system
         t_srs = paste0("EPSG:", utm_code),
         tr = c("1000", "1000"),
         r = "sum",
         ot = "UInt16")
         #cutline = file.path(getwd(), "country_prj.gpkg"),
         #crop_to_cutline = TRUE
    

# Rasters --> Dataframe
#gdal_utils(util = "info", "/vsimem/warpSum.vrt")
r.warpSum = rast(file.path(getwd(), "GDAL_R", "warpSum.tiff"))
#as.data.frame(r, na.rm = TRUE)[1,]

df = as.data.frame(r.warpSum, na.rm=TRUE, xy = TRUE, centroids = TRUE)

# DF to sf object
gdf = st_as_sf(df, 
               coords = c("x", "y"), 
               crs = crs(gadm_prj))

gdf_aoi = gdf %>% st_join(gadm_prj) %>% drop_na(COUNTRY)
