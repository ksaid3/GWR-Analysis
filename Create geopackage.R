library(sf)


# Read the CSV file
csv_data <- read.csv("C:/Users/Kalee/Downloads/calenviroscreen-3.0-results-june-2018-update.csv")
csv_data
# Convert to spatial data frame
coordinates <- c("Longitude", "Latitude")
csv_sf <- st_as_sf(csv_data, coords = coordinates)


st_crs(csv_sf) <- st_crs("EPSG:4326")

# Read the shapefile
shapefile_data <- st_read("C:/Users/Kalee/Downloads/calenviroscreen-3.0-shapefile/CES3Results.shp")


st_write(shapefile_data, "Calicountydata.gpkg", driver = "GPKG")
