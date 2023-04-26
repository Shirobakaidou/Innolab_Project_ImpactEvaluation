library(sf)
library(terra)
library(raster)
library(gdalUtilities)
library(geodata)
library(wdpar)
library(tidyr)
library(dplyr)

country = "Cameroon"
path.input = file.path(getwd(), 'data', 'input', 'Cameroon')
path.output = file.path(getwd(), 'data', 'output', 'Cameroon')

# Specify buffer width in meter
buffer_m = 10000

# Specify the grid cell size in meter
gridSize = 1000

# Specify a list of WDPA IDs of funded protected areas (treated areas)
#paid = c(5037, 2303, 5024, 303698, 303702, 26072, 5028, 303695, 20272, 20287, 20273)
# IDs of funded PAs in Cameroon
paid = c(20058, 555547996, 555547994, 20112)

#--------------------------------------------------------------
# Vector: Reproject GADM ####
gadm <- gadm(country = country, resolution = 1, level = 0, path = path.input) %>%
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
st_write(gadm_prj, 
         dsn = file.path(path.input, "gadm_prj.gpkg"), delete_dsn = TRUE) 

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
wdpa = wdpa_fetch(country, wait = TRUE, download_dir = path.input)
# If the PA file already exists, it can be loaded in this way
#wdpa = wdpa_read(paste0(path.input, '/WDPA_Dec2022_MDG-shapefile.zip'))


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

# Warp Forest Cover
f.fc = Sys.glob(file.path(path.input, 'fc*'))

fun.warp.fc = function(file){
  newname = paste0("warped-", tail(strsplit(file, '/')[[1]], 1))
  gdalwarp(srcfile = file,
         dstfile = file.path(path.output, newname), # Output file in GDAL virtual file system
         t_srs = paste0("EPSG:", utm_code),
         tr = c(as.character(gridSize), as.character(gridSize)),
         r = "sum",
         ot = "UInt16")}

lapply(f.fc, fun.warp.fc)

#----------------------------------
# Terrain Ruggedness Index (TRI)
f.dem = Sys.glob(file.path(path.input, 'srtm*'))

fun.tri = function(file){
  gdaldem(mode = "TRI", 
        input_dem = file,
        output_map = file.path(path.input, "tri.tif"))
}

lapply(f.dem, fun.tri)


#--------------------------------
# Warp Other Variables: travel time, elevation, tri, clay content
f.others = Sys.glob(file.path(path.input, '[!f]*.tif'))

fun.warp.other = function(file){
  newname = paste0("warped-", tail(strsplit(file, '/')[[1]], 1))
  gdalwarp(srcfile = file,
         dstfile = file.path(path.output, newname), # Output file in GDAL virtual file system
         t_srs = paste0("EPSG:", utm_code),
         tr = c(as.character(gridSize), as.character(gridSize)),
         r = "average",
         ot = "UInt16")
}

lapply(f.others, fun.warp.other)


#------------------------------------------------------------
# Rasters --> Dataframe

rToDF = function(file){
  r = rast(file)
  df = as.data.frame(r, na.rm=TRUE, xy=TRUE, centroids=TRUE)
}

# Forest Cover
f.fc = Sys.glob(file.path(path.output, '*fc*'))
df.fc = do.call("rbind", lapply(f.fc, FUN = rToDF)) %>%
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


# Travel Time
f.travel = Sys.glob(file.path(path.output, '*travel*'))
df.travel = do.call("rbind", lapply(f.travel, FUN = rToDF)) %>%
  st_as_sf(., coords = c("x", "y"),
           crs = crs(gadm_prj)) %>%
  st_join(grid) %>% 
  drop_na(gridID) %>%
  mutate(geometry = NULL) %>%
  as.data.frame() %>%
  select(-c(centroids)) %>%
  group_by(gridID) %>% 
  summarise(across(everything(), list(mean)))


# Clay Content
f.clay = Sys.glob(file.path(path.output, '*clay*'))
df.clay = do.call("rbind", lapply(f.clay, FUN = rToDF)) %>%
  st_as_sf(., coords = c("x", "y"),
           crs = crs(gadm_prj)) %>%
  st_join(grid) %>% 
  drop_na(gridID) %>%
  mutate(geometry = NULL) %>%
  as.data.frame() %>%
  select(-c(centroids)) %>%
  group_by(gridID) %>% 
  summarise(across(everything(), list(mean)))


# DEM
f.dem = Sys.glob(file.path(path.output, '*srtm*'))
df.dem = do.call("rbind", lapply(f.dem, FUN = rToDF)) %>%
  st_as_sf(., coords = c("x", "y"),
           crs = crs(gadm_prj)) %>%
  st_join(grid) %>% 
  drop_na(gridID) %>%
  mutate(geometry = NULL) %>%
  as.data.frame() %>%
  select(-c(centroids)) %>%
  group_by(gridID) %>% 
  summarise(across(everything(), list(mean)))


# TRI
f.tri = Sys.glob(file.path(path.output, '*tri*'))
df.tri = do.call("rbind", lapply(f.tri, FUN = rToDF)) %>%
  st_as_sf(., coords = c("x", "y"),
           crs = crs(gadm_prj)) %>%
  st_join(grid) %>% 
  drop_na(gridID) %>%
  mutate(geometry = NULL) %>%
  as.data.frame() %>%
  select(-c(centroids)) %>%
  group_by(gridID) %>% 
  summarise(across(everything(), mean))


# Groups and WDPA-ID
df.group = rToDF(r.group) %>%
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

df.wdpaid = rToDF(r.wdpaid) %>%
  st_as_sf(., coords = c("x", "y"),
           crs = crs(gadm_prj)) %>%
  st_join(grid) %>%
  drop_na(gridID) %>%
  mutate(geometry = NULL) %>%
  as.data.frame() %>%
  select(-c(centroids)) %>%
  group_by(gridID) %>% 
  summarise(across(everything(), max))

df.grid = df.group %>% inner_join(df.wdpaid)
rm(df.group, df.wdpaid)

mf = grid %>%
  inner_join(df.grid) %>%
  inner_join(df.fc) %>%
  inner_join(df.dem) %>%
  inner_join(df.tri) %>%
  inner_join(df.travel) %>%
  inner_join(df.clay) 

# Export projected GADM for use in GDAL
st_write(mf, 
         dsn = file.path(path.output, "mf_rBase_1km.gpkg"), 
         delete_dsn = TRUE) 
  