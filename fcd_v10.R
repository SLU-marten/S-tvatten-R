library(sp)
library(foreach)
library(doParallel)
library(gdistance)
library(raster)
library(sf)
library(leaflet)
library(leafem)

# Enable parallel processing
(n.cores <- parallel::detectCores() - 1) # Subtract 1 to leave one core free
clust <- parallel::makeCluster(n.cores, type = "PSOCK")
doParallel::registerDoParallel(cl = clust)
foreach::getDoParRegistered()

# Specify input and output folders
path_in <-  "./inputs/"
path_out <- "./outputs/"

# Specify parameters for the model
max.dist <- 10000   # This is the maximum dispersal distance for the species

# Import cost/destination/goal raster
base.cost <- raster(paste0(path_in, "cost_raster_250m.tif"))
base.cost[base.cost == 0] <- NA
plot(base.cost)

# Define extent
ex <- extent(base.cost)

# Create a function for the analyze
fun <- function(r_name, out_name){
  # Import lakes
  lake <- raster(paste0(path_in, r_name))
  #lake <- raster(paste0(path_in,"vattenytorRaster_res250.tif"))
  lake <- crop(lake, ex)
  lake[lake == 0] <- NA
  lake[is.na(base.cost) == TRUE] <- NA
  plot(lake)
  
  # Convert origin and goal rasters to points
  lake.o <- raster::rasterToPoints(lake, spatial = TRUE) # Set the origin points
  lake.g <- raster::rasterToPoints(base.cost, spatial = TRUE) # Set the destination points
  length(lake.o[]) # Number of habitat cells
  length(lake.g[]) # Number cells
  
  # Create transition layer for cost distance
  # In the package gdistance a conductance of 0 is a complete barrier, 0.01 is very costly, and a conductance of 1 has cost is equal to distance
  # So for example if a 250m cell has a conductance of 1 the cost of travelling through the cell is 250
  # If the cell has a conductance of 0.5 the cost of travelling through the cell is 500, etc.
  lake.tr <- gdistance::transition(base.cost, transitionFunction = min, directions = 8) # Transition file
  lake.tr <- gdistance::geoCorrection(lake.tr, type = "c") # Geo-correct transition file
  
  # Focal cost distance function
  fcd <- function(ncells, nhabs, speco, specg, spectr){
    cd.base <- c(1:ncells) * 0
    for(i in 1:nhabs){
      lake.cd <- gdistance::costDistance(spectr, speco[i,], specg)
      cd.temp <- lake.cd[1,]
      cd.temp[cd.temp > max.dist] <- NA
      cd.temp[is.na(cd.temp) == FALSE] <- 1#exp(-a*(cd.temp[is.na(cd.temp) == FALSE]))
      cd.temp[is.na(cd.temp) == TRUE] <- 0
      cd.base <- cd.base+cd.temp
    }
    return(cd.base)
  }
  
  # Run focal cost distance in parallel
  nc <- length(lake.g[])
  nh <- length(lake.o[])
  nchunks <- n.cores # Divides the habitat points into chunks for parallel processing, this should be equal to the number of cores used
  blocks <- seq(from = 0, to = nh, length.out = nchunks + 1)
  blocks <- round(blocks)
  
  pfcd <- foreach(i = 1:nchunks, .combine = '+') %dopar% {
    k <- (blocks[i] + 1):blocks[i + 1]
    fcd(ncells = nc, nhabs = length(k), speco = lake.o[k,], specg = lake.g, spectr = lake.tr)
  }
  
  saveRDS(pfcd, paste0(path_out, out_name,".rds"))
  
  # Rasterize
  r.river <- base.cost # New raster
  points.river <- lake.g # 
  points.river$cost_raster_250m <- pfcd
  r2.river <- raster::rasterize(points.river, r.river, "cost_raster_250m", fun = mean)
  r2.river[r2.river == 0] <- NA
  writeRaster(r2.river, paste0(path_out, out_name, ".tif"))
}

# Run lakes
fun("vattenytorRaster_res250.tif", "lake250")
# Run rivers
fun("vattendrag_no_Overlap_res.tif", "river250")


# Calculate combined
r2.r <- raster("river250.tif")
r2.l <- raster("lake250.tif")

r2.r[is.na(r2.r)] <- 0
r2.l[is.na(r2.l)] <- 0
r3 <- r2.l * 0.047 + r2.r * 0.002
writeRaster(r3, paste0(path_out, "lek_sötvatten.tif"))
r3[r3 == 0] <- NA




# Visualize results

# Shapefiles
survey <- st_read(paste0(path_in,"CoastLand areas.shp"))
s <- nngeo::st_remove_holes(survey) |> st_transform("+init=epsg:4326")

 # Visulalize
pal <- colorNumeric("Spectral", domain = NULL, na.color = "#00000000", reverse = T)
factpal <- colorFactor(palette.colors(9), s$ORIG_FID)

leaflet() |> 
  addProviderTiles(providers$CartoDB.Positron) |>
  addRasterImage(r2, colors = pal, layerId = "values", group = "values") |>
  addImageQuery(r2, type="mousemove", layerId = "values")

leaflet() |> 
  addProviderTiles(providers$CartoDB.Positron) |>
  addRasterImage(r2.river, colors = pal, layerId = "values", group = "values") |>
  addImageQuery(r2.river, type="mousemove", layerId = "values")



l <- leaflet() |> 
  addProviderTiles(providers$CartoDB.Positron) |>
  addRasterImage(r3, colors = pal, layerId = "values", group = "values") |>
  addImageQuery(r3, type="mousemove", layerId = "values") |>
  addPolygons(data = s, 
              label = ~paste0(Fångstomr, ", ", FISKE),
              color = ~factpal(ORIG_FID)
  )

htmlwidgets::saveWidget(l, file = "sötvatten.html")
