library(sf)
library(terra)
library(raster)
library(gdalUtilities)
library(geodata)
library(wdpar)
library(tidyr)
library(dplyr)

country = "MDG"
path.input = file.path(getwd(), 'data', 'input', "Madagascar")
path.output = file.path(getwd(), 'data', 'test')

# Specify buffer width in meter
buffer_m = 10000
# Specify the grid cell size in meter
gridSize = 1000
# Specify a list of WDPA IDs of funded protected areas (treated areas)
paid = c(5037, 2303, 5024, 303698, 303702, 26072, 5028, 303695, 20272, 20287, 20273)


# If the country polygon file already exists in working directory, it can be loaded into work space by
gadm = readRDS(file.path(path.input, paste0('gadm41_',country,'_0_pk.rds'))) %>% st_as_sf()
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
# Make bounding box of projected country polygon
bbox_prj = st_bbox(gadm_prj) %>% st_as_sfc() %>% st_as_sf()

###=====================WDPA===============================###
# If the PA file already exists, it can be loaded in this way
wdpa = wdpa_read(Sys.glob(file.path(path.input, "WDPA*")))

# PAs are projected, and column "geometry_type" is added
wdpa_prj = wdpa %>%
  
  wdpa_clean(retain_status = NULL,
             erase_overlaps = FALSE,
             exclude_unesco = FALSE,
             verbose = TRUE) %>% 
  # Drop geometry POINT
  filter(GEOMETRY_TYPE != "POINT") %>%
  # Remove the PAs that are only proposed, or have geometry type "point"
  filter(STATUS != "Proposed") %>%
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
# CAREFUL : the order of the arguments does matter. 
## During rasterization, in case a cell of the raster is on both funded analysed and non-funded, 
# we want the cell to take the WDPAID of the funded analysed.
wdpa_groups = rbind(wdpa_funded, wdpa_nofund, buffer)
# Subset to polygons that intersect with country boundary
wdpa.sub = wdpa_groups %>% 
  st_intersects(gadm_prj, .) %>% 
  unlist()
# Filter the PA+buffer to the subset
wdpa_groups = wdpa_groups[sort(wdpa.sub),] %>%
  st_as_sf()

#------------------------------------------------------------
# Initialize an empty raster to the spatial extent of the country
r.ini = raster()
# Assign CRS to raster
crs(r.ini) = CRS(paste0('+init=EPSG:',utm_code))
# Define Extent
extent(r.ini) = extent(bbox_prj)
# Specify the raster resolution as same as the pre-defined 'gridSize'
res(r.ini) = gridSize

# Assign the raster pixels with "Group" values, 
# Take the minimal value if a pixel is covered by overlapped polygons, so that PA Group ID has higher priority than Buffer ID.
# Assign value "0" to the background pixels (control candidates group)
# fun = "min" can lead to bad group assignment. This issue is developed and tackled below
r.group = rasterize(wdpa_groups, r.ini, field="group", fun="min", background=0) %>%
  mask(., gadm_prj)
# Rename Layer
names(r.group) = "group"

# Rasterize wdpaid
## CAREFUL : as stated above, the wdpa_groups raster is ordered so that the first layer is the one of funded, analyzed PA. Thus one needs to have fun = "first"
r.wdpaid = rasterize(wdpa_groups, r.ini, field="WDPAID", fun="first", background=0) %>%
  mask(., gadm_prj)
names(r.wdpaid) = "wdpaid"

# Export Rasters of "group" and "wdpaid"
writeRaster(r.group, file.path(path.output,'r_group_test.tif'))
writeRaster(r.wdpaid, file.path(path.output, 'r_wdpaid_test.tif'))

###====================================================###

# GDAL: Rasterize Bounding Box
gdal_rasterize(src_datasource = file.path(path.output, "bboxMDG_prj.gpkg"), # Input
               dst_filename = file.path(path.output, "r_bbox1.tif"), # Output
               of = 'GTiff', # Format
               ot = "Byte", # Output Type
               co = c("COMPRESS=DEFLATE"),
               a_nodata = NA,
               tr = c(as.character(gridSize), as.character(gridSize)), # Pixel Resolution
               #at = TRUE, # all touched
               burn = 1)

# GDAL: Rasterize Attribute "Group"
gdal_rasterize(src_datasource = file.path(path.output, "wdpa_groups_MDG.gpkg"), # Input
               dst_filename = file.path(path.output, "r_group.tif"), # Output
               of = 'GTiff', # Format
               ot = "Byte", # Output Type
               co = c("COMPRESS=DEFLATE"),
               a = "group", # The attribute to burn value
               a_nodata = 0, # Assign specific no data value: 0 for control group
               te = c(311443.29363404075, 7168029.484370736, 1095443.2936340407, 8678029.484370736),
               tr = c(as.character(gridSize), as.character(gridSize)) # Pixel Resolution
)

# GDAL: Rasterize Attribute "WDPAID"
gdal_rasterize(src_datasource = file.path(path.output, "wdpa_groups_MDG.gpkg"), # Input
               dst_filename = file.path(path.output, "r_wdpaid.tif"), # Output
               of = 'GTiff', # Format
               ot = "Byte", # Output Type
               co = c("COMPRESS=DEFLATE"),
               a = "WDPAID", # The attribute to burn value
               a_nodata = 0, # Assign specific no data value: 0 for control group
               te = c(311443.29363404075, 7168029.484370736, 1095443.2936340407, 8678029.484370736),
               tr = c(as.character(gridSize), as.character(gridSize)) # Pixel Resolution
)


# GDAL WARP ####
#jsonlite::fromJSON(gdalinfo(file.path(path.output, "r_bbox.tif"),
#         json = TRUE))

#------------------------------------------------------------
# GDAL: Build VRT
gdalbuildvrt(gdalfile = Sys.glob(file.path(path.input, 'fc*.tif')),
             output.vrt = file.path(path.output, "VRT_fc.vrt"))
gdalbuildvrt(gdalfile = Sys.glob(file.path(path.input, 'srtm*.tif')),
             output.vrt = file.path(path.output, "VRT_srtm.vrt"))
gdalbuildvrt(gdalfile = Sys.glob(file.path(path.input, 'travel*.tif')),
             output.vrt = file.path(path.output, "VRT_traveltime.vrt"))
gdalbuildvrt(gdalfile = Sys.glob(file.path(path.input, 'clay*.tif')),
             output.vrt = file.path(path.output, "VRT_clay.vrt"))


# Terrain Ruggedness Index (TRI)
f.dem = Sys.glob(file.path(path.input, 'srtm*.tif'))
fun.tri = function(file){
  #newname = paste0("tri-", tail(strsplit(file, '/')[[1]], 1))
  gdaldem(mode = "TRI", 
          input_dem = file,
          output_map = file.path(path.output, "tri.tif"))
}
lapply(f.dem, fun.tri)

gdalbuildvrt(gdalfile = Sys.glob(file.path(path.output, 'tri*.tif')),
             output.vrt = file.path(path.output, "VRT_tri.vrt"))


#------------------------------------------------------------
# GDAL: Align rasters by Warping
extent.bbox = extent(bbox_prj)
# Forest Cover
gdalwarp(srcfile = Sys.glob(file.path(path.output, 'VRT_fc.vrt')),
         dstfile = file.path(path.output, "wp_fc.vrt"),
         t_srs = paste0("EPSG:", utm_code),
         te = c(extent.bbox[1], extent.bbox[3], extent.bbox[2], extent.bbox[4]),
         tr = c(as.character(gridSize), as.character(gridSize)),
         r = "sum",
         ot = "UInt16",
         #cutline = file.path(path.input, "gadm_prj.gpkg"),
         #crop_to_cutline = TRUE,
         dstnodata = NA,
         overwrite = TRUE)

# SRTM
gdalwarp(srcfile = Sys.glob(file.path(path.output, 'VRT_srtm.vrt')),
         dstfile = file.path(path.output, "wp_srtm.vrt"),
         t_srs = paste0("EPSG:", utm_code),
         te = c(extent.bbox[1], extent.bbox[3], extent.bbox[2], extent.bbox[4]),
         #te = c(311443.29363404075, 7168029.484370736, 1095443.2936340407, 8678029.484370736),
         tr = c(as.character(gridSize), as.character(gridSize)),
         r = "average",
         ot = "UInt16",
         #cutline = file.path(path.input, "gadm_prj.gpkg"),
         #crop_to_cutline = TRUE,
         dstnodata = NA,
         overwrite = TRUE)

# TRI
gdalwarp(srcfile = Sys.glob(file.path(path.output, 'VRT_tri.vrt')),
         dstfile = file.path(path.output, "wp_tri.vrt"),
         t_srs = paste0("EPSG:", utm_code),
         te = c(extent.bbox[1], extent.bbox[3], extent.bbox[2], extent.bbox[4]),
         #te = c(311443.29363404075, 7168029.484370736, 1095443.2936340407, 8678029.484370736),
         tr = c(as.character(gridSize), as.character(gridSize)),
         r = "average",
         ot = "Float32",
         #cutline = file.path(path.input, "gadm_prj.gpkg"),
         #crop_to_cutline = TRUE,
         dstnodata = NA,
         overwrite = TRUE)

# Travel Time
gdalwarp(srcfile = Sys.glob(file.path(path.output, 'VRT_traveltime.vrt')),
         dstfile = file.path(path.output, "wp_traveltime.vrt"),
         t_srs = paste0("EPSG:", utm_code),
         te = c(extent.bbox[1], extent.bbox[3], extent.bbox[2], extent.bbox[4]),
         #te = c(311443.29363404075, 7168029.484370736, 1095443.2936340407, 8678029.484370736),
         tr = c(as.character(gridSize), as.character(gridSize)),
         r = "average",
         ot = "UInt16",
         #cutline = file.path(path.input, "gadm_prj.gpkg"),
         #crop_to_cutline = TRUE,
         dstnodata = NA,
         overwrite = TRUE)

# Clay Content
gdalwarp(srcfile = Sys.glob(file.path(path.output, 'VRT_clay.vrt')),
         dstfile = file.path(path.output, "wp_clay.vrt"),
         t_srs = paste0("EPSG:", utm_code),
         te = c(extent.bbox[1], extent.bbox[3], extent.bbox[2], extent.bbox[4]),
         #te = c(311443.29363404075, 7168029.484370736, 1095443.2936340407, 8678029.484370736),
         tr = c(as.character(gridSize), as.character(gridSize)),
         r = "average",
         ot = "UInt16",
         #cutline = file.path(path.input, "gadm_prj.gpkg"),
         #crop_to_cutline = TRUE,
         dstnodata = NA,
         overwrite = TRUE)

#------------------------------------------------------------
# Rasters --> Dataframe

lst_aligned_layers = file.path(path.output, "r_wdpaid_test.tif") %>%
  append(file.path(path.output, "r_group_test.tif")) %>%
  append(Sys.glob(file.path(path.output, "wp*[!fc].vrt"))) 

rToDF = function(file){
  r = rast(file)
  df = as.data.frame(r, na.rm=TRUE, xy=TRUE, centroids=TRUE) %>%
    select(-c(centroids))
}
# Load rasters of other covariates as DF, in list
df.covar = lapply(lst_aligned_layers, FUN = rToDF)
# Merge DFs
df.merge.test = Reduce(function(x, y) merge(x, y, by=c("x","y")),
                       df.covar)


# Dataframe of Forest Cover
df.fc = rToDF(file.path(path.output, "wp_fc.vrt"))
# Rename df.fc
fc.names.old = colnames(df.fc)[grepl("wp", colnames(df.fc))]
fc.names.new = lapply(fc.names.old, function(i){
  postfix = as.numeric(strsplit(i, "_")[[1]][3])-1+2000
  newname = paste0("fc_percent_", postfix)
}) %>% unlist()

df.fc1 = df.fc %>% 
  # Rename FC Columns
  rename_with(~ fc.names.new, .cols = fc.names.old) %>%
  # Calculate forest cover percentage in grid cells
  mutate_at(fc.names.new, ~round((.*900)/(gridSize^2)*100, 1)) %>%
  # if value > 100% by bias, set it to 100%
  mutate_at(fc.names.new, ~replace(., .>100, 100))


# Calculate Forest Loss Time Series
# Add new columns: treeloss_tn = treecover_tn - treecover_t(n-1)
for (i in 1:as.numeric(length(fc.names.new)-1)) {
  # Drop the first year
  dropFirst = tail(fc.names.new, -1)
  # Drop the last year
  dropLast = head(fc.names.new, -1)
  # Forest Loss name postfix
  postfix = dropFirst %>% strsplit(., "_")
  # FL Column Name
  fl.colname = paste0("fl_percent_", postfix[[i]][3])
  # Add FL column
  df.fc1[[fl.colname]] = df.fc1[[dropFirst[i]]] - df.fc1[[dropLast[i]]]
}


# Merge DFs
df.merge = merge(df.merge.test, df.fc1, by=c("x","y"))
# DF --> GDF
gdf.merge = df.merge %>%
  st_as_sf(., coords = c("x", "y"),
           crs = utm_code) %>%
  st_transform(crs = 4326)
