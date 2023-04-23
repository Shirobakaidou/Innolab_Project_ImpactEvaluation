library(sf)
library(terra)
library(raster)
library(gdalUtilities)
library(geodata)
library(wdpar)
library(tidyr)
library(dplyr)

country = "Madagascar"
wdir = file.path(getwd(), 'data')

# Specify buffer width in meter
buffer_m = 10000

# Specify the grid cell size in meter
gridSize = 1000

# Specify a list of WDPA IDs of funded protected areas (treated areas)
paid = c(5037, 2303, 5024, 303698, 303702, 26072, 5028, 303695, 20272, 20287, 20273)


#--------------------------------------------------------------
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
st_write(test, 
         dsn = file.path(wdir, "group_prj.gpkg"), delete_dsn = TRUE) 

# Make bounding box of projected country polygon
bbox = st_bbox(gadm_prj) %>% st_as_sfc() %>% st_as_sf() 

# Make a Grid to the extent of the bounding box
grid = st_make_grid(bbox, cellsize = c(gridSize,gridSize))

# Crop Grid to the extent of country boundary by
# subsetting to the grid cells that intersect with the country
grid.sub = grid %>% 
  st_intersects(gadm_prj, .) %>% 
  unlist()
# Filter the grid to the subset
grid = grid[sort(grid.sub)] %>%
  st_as_sf() %>%
  mutate(gridID = seq(1:nrow(.))) # Add id for grid cells
rm(grid.sub)
# Vector: Groups (funded PA, non-funded PA, buffer, control candidates) ####

# Get the PA polygons/points of the specified country; 
# They're downloaded to the working directory.
wdpa = wdpa_fetch(country, wait = TRUE, download_dir = wdir)
# If the PA file already exists, it can be loaded in this way
#wdpa = wdpa_read(paste0(wdir, '/WDPA_Dec2022_MDG-shapefile.zip'))


# PAs are projected, and column "geometry_type" is added
wdpa_prj = wdpa_clean(wdpa, geometry_precision = 1000) %>%
  # Remove the PAs that are only proposed, or have geometry type "point"
  filter(STATUS != "Proposed") %>%
  filter(GEOMETRY_TYPE != "POINT") %>%
  # Project PA polygons to the previously determined UTM zone
  st_transform(crs = utm_code) 
rm(wdpa)


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

rm(buffer, wdpa_funded, wdpa_nofund)


# Initialize an empty raster to the spatial extent of the country
r.ini = raster()
extent(r.ini) = extent(bbox)
# Specify the raster resolution as same as the pre-defined 'gridSize'
res(r.ini) = gridSize

# Assign the raster pixels with "Group" values, 
# Take the minial value if a pixel is covered by overlapped polygons, so that PA Group ID has higher priority than Buffer ID.
# Assign value "0" to the background pixels (control candidates group)
r.group = rasterize(wdpa_groups, r.ini, field="group", fun="min", background=0) %>%
  mask(., gadm_prj)
# Rename Layer
names(r.group) = "group"

# Rasterize wdpaid
r.wdpaid = rasterize(wdpa_prj, r.ini, field="WDPAID", fun="first", background=0) %>%
  mask(., gadm_prj)
names(r.wdpaid) = "wdpaid"


#-------------------------------------------------------------
# GDAL WARP ####
files = list.files(path = file.path(getwd(), "data", "GEE", "fc"), 
                   full.names = TRUE)

fInfo = gdal_utils(util = "info",
                   source = files[1])


gdalwarp(srcfile = files[8],
         dstfile = file.path(getwd(), "data", "GDAL", "fc", "fcSum_1km_8.tiff"), # Output file in GDAL virtual file system
         t_srs = paste0("EPSG:", utm_code),
         tr = c("1000", "1000"),
         r = "sum",
         ot = "UInt16")
         #cutline = file.path(getwd(), "country_prj.gpkg"),
         #crop_to_cutline = TRUE

# Travel Time
files_travel = list.files(path = file.path(getwd(), "data", "GEE", "travelTime"), 
                          full.names = TRUE)

gdalwarp(srcfile = files_travel[1],
         dstfile = file.path(getwd(), "data", "GDAL", "travelTime", "travelMean_1km.tiff"), # Output file in GDAL virtual file system
         t_srs = paste0("EPSG:", utm_code),
         tr = c("1000", "1000"),
         r = "average",
         ot = "UInt16")

# Clay Content
files_clay = list.files(path = file.path(getwd(), "data", "GEE", "clay"),
                        full.names = TRUE)

gdalwarp(srcfile = files_clay[1],
         dstfile = file.path(getwd(), "data", "GDAL", "clay", "clayMean_1km.tiff"), # Output file in GDAL virtual file system
         t_srs = paste0("EPSG:", utm_code),
         tr = c("1000", "1000"),
         r = "average",
         ot = "UInt16")

# DEM and TRI
files_dem = list.files(path = file.path(getwd(), "data", "GEE", "dem"),
                        full.names = TRUE)

gdalwarp(srcfile = files_dem[1],
         dstfile = file.path(getwd(), "data", "GDAL", "dem", "demMean_1km.tiff"), 
         t_srs = paste0("EPSG:", utm_code),
         tr = c("1000", "1000"),
         r = "average",
         ot = "UInt16")

gdaldem(mode = "TRI", input_dem = files_dem[1],
        output_map = file.path(getwd(), "data", "GEE", "dem", "tri.tiff"))

gdalwarp(srcfile = files_dem[2],
         dstfile = file.path(getwd(), "data", "GDAL", "dem", "triMean_1km.tiff"), 
         t_srs = paste0("EPSG:", utm_code),
         tr = c("1000", "1000"),
         r = "average",
         ot = "UInt16")

#------------------------------------------------------------
#------------------------------------------------------------

# Rasters --> Dataframe
files_fc = list.files(file.path(getwd(), "data", "GDAL", "fc"),
                      full.names = TRUE)
files_travel = list.files(file.path(getwd(), "data", "GDAL", "travelTime"),
                      full.names = TRUE)
files_clay = list.files(file.path(getwd(), "data", "GDAL", "clay"),
                          full.names = TRUE)
files_dem = list.files(file.path(getwd(), "data", "GDAL", "dem"),
                          full.names = TRUE)
files_tri = list.files(file.path(getwd(), "data", "GDAL", "tri"),
                       full.names = TRUE)

rToDF = function(file){
  r = rast(file)
  df = as.data.frame(r, na.rm=TRUE, xy=TRUE, centroids=TRUE)
  }

df_fc = do.call("rbind", lapply(files_fc, FUN = rToDF)) %>%
  # Transform warped rasters into dataframes, merge dfs by row
  
  # Turn df into projected geodataframe
  st_as_sf(., coords = c("x", "y"),
           crs = crs(gadm_prj)) %>%
  
  # Spatial join gdf with grid gdf, match unique gridID and observation
  st_join(grid) %>% 
  
  # Remove observations not covered by grid
  drop_na(gridID) %>%
  
  # Drop geometry, turn gdf into df, drop unneeded column
  mutate(geometry = NULL) %>%
  as.data.frame() %>%
  select(-c(centroids)) %>%
  
  # Summarise observations with duplicated gridID
  group_by(gridID) %>% 
  summarise(across(everything(), list(sum)))
  # Alternative: aggregate(.~gridID, gdf_fc, sum)


df_travel = do.call("rbind", lapply(files_travel, FUN = rToDF)) %>%
  st_as_sf(., coords = c("x", "y"),
           crs = crs(gadm_prj)) %>%
  st_join(grid) %>% 
  drop_na(gridID) %>%
  mutate(geometry = NULL) %>%
  as.data.frame() %>%
  select(-c(centroids)) %>%
  group_by(gridID) %>% 
  summarise(across(everything(), list(mean)))

df_clay = do.call("rbind", lapply(files_clay, FUN = rToDF)) %>%
  st_as_sf(., coords = c("x", "y"),
           crs = crs(gadm_prj)) %>%
  st_join(grid) %>% 
  drop_na(gridID) %>%
  mutate(geometry = NULL) %>%
  as.data.frame() %>%
  select(-c(centroids)) %>%
  group_by(gridID) %>% 
  summarise(across(everything(), list(mean)))

df_dem = do.call("rbind", lapply(files_dem, FUN = rToDF)) %>%
  st_as_sf(., coords = c("x", "y"),
           crs = crs(gadm_prj)) %>%
  st_join(grid) %>% 
  drop_na(gridID) %>%
  mutate(geometry = NULL) %>%
  as.data.frame() %>%
  select(-c(centroids)) %>%
  group_by(gridID) %>% 
  summarise(across(everything(), list(mean)))

df_tri = do.call("rbind", lapply(files_tri, FUN = rToDF)) %>%
  st_as_sf(., coords = c("x", "y"),
           crs = crs(gadm_prj)) %>%
  st_join(grid) %>% 
  drop_na(gridID) %>%
  mutate(geometry = NULL) %>%
  as.data.frame() %>%
  select(-c(centroids)) %>%
  group_by(gridID) %>% 
  summarise(across(everything(), mean))

df_group = rToDF(r.group) %>%
  st_as_sf(., coords = c("x", "y"),
           crs = crs(gadm_prj)) %>%
  st_join(grid) %>%
  drop_na(gridID) %>%
  mutate(geometry = NULL) %>%
  as.data.frame() %>%
  select(-c(centroids)) %>%
  # duplicated gridID hardly happens; only for just in case.
  group_by(gridID) %>% 
  summarise(across(everything(), max))

df_wdpaID = rToDF(r.wdpaid) %>%
  st_as_sf(., coords = c("x", "y"),
           crs = crs(gadm_prj)) %>%
  st_join(grid) %>%
  drop_na(gridID) %>%
  mutate(geometry = NULL) %>%
  as.data.frame() %>%
  select(-c(centroids)) %>%
  group_by(gridID) %>% 
  summarise(across(everything(), max))

df_grid = df_group %>% inner_join(df_wdpaID)
rm(df_group, df_wdpaID)

mf = grid %>%
  inner_join(df_grid) %>%
  inner_join(df_fc) %>%
  inner_join(df_dem) %>%
  inner_join(df_tri) %>%
  inner_join(df_travel) %>%
  inner_join(df_clay) #%>%
  inner_join(grid) %>%
  st_as_sf(., coords)

# Export projected GADM for use in GDAL
st_write(mf, 
         dsn = file.path(wdir, "mf_rBase_1km.gpkg"), 
         delete_dsn = TRUE) 
  