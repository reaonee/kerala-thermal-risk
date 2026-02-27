# --- MASTER SCRIPT: STEP 1 (SPATIAL FOUNDATION) ---

library(sf)
library(dplyr)
library(spdep)

# 1. Set directory and load the massive India map

india_map <- st_read("Kerala_github/2011_Dist.shp")

# 2. Transform the whole map to UTM Zone 43N first (so the math is flat)
if(st_crs(india_map)$epsg != 32643) {
  india_map <- st_transform(india_map, crs = 32643)
  print("CRS successfully transformed to UTM Zone 43N.")
}

# 3. Use the scissors: Isolate ONLY the rows containing "Kerala"
kerala_only <- india_map[apply(st_drop_geometry(india_map), 1, function(x) any(grepl("Kerala", x, ignore.case = TRUE))), ]

# 4. Clean up the IDs for the hhh4 model so they count from 1 to 14
kerala_only <- kerala_only %>%
  mutate(Region_ID = row_number())

# 5. Verify the foundation visually
#plot(st_geometry(kerala_only), main = "Actual Kerala Spatial Framework (UTM 43N)")

# 6. Print the column names so we know exactly what the district column is called
print(names(kerala_only))




# --- MASTER SCRIPT: STEP 2 (IGNITION COUNTS AND SPATIAL JOIN) ---
library(readr)
library(lubridate)
library(dplyr)
library(sf)

# 1. Load BOTH historical 5-year baselines from the GitHub root folder
modis_history <- read_csv("modis_history.csv", show_col_types = FALSE) %>%
  mutate(confidence = as.character(confidence))

viirs_history <- read_csv("viirs_history.csv", show_col_types = FALSE) %>%
  mutate(confidence = as.character(confidence))

# Stitch them into one massive historical dataset
historical_fires <- bind_rows(modis_history, viirs_history)

# 2. Fetch LIVE data from NASA FIRMS API (Last 7 Days for India)
# The GitHub robot will securely inject your API key here
nasa_key <- Sys.getenv("NASA_FIRMS_KEY")

if(nasa_key == "") {
  print("WARNING: NASA_FIRMS_KEY not found in server. Using only historical data.")
  fire_data_raw <- historical_fires
} else {
  print("Fetching live fire data from NASA FIRMS API...")
  api_url <- paste0("https://firms.modaps.eosdis.nasa.gov/api/country/csv/", nasa_key, "/VIIRS_SNPP_NRT/IND/7")
  
  # Try to download live data, but don't crash if the NASA server is temporarily down
  live_fires <- tryCatch({
    read_csv(api_url, show_col_types = FALSE) %>%
      mutate(confidence = as.character(confidence))
  }, error = function(e) {
    print("API fetch failed. Defaulting to historical data only.")
    return(data.frame())
  })
  
  # 3. Stitch the historical data and the live data together
  fire_data_raw <- bind_rows(historical_fires, live_fires) %>% distinct()
}

# 4. Convert the raw Excel numbers into a spatial "sf" object (WGS84 -> UTM 43N)
fires_sf <- st_as_sf(fire_data_raw, coords = c("longitude", "latitude"), crs = 4326)
fires_sf <- st_transform(fires_sf, crs = 32643)

# 5. The Spatial Join (Point-in-Polygon)
fires_in_kerala <- st_join(fires_sf, kerala_only, join = st_intersects, left = FALSE)

# 6. Print a quick summary of total fires caught inside the borders
print(paste("Total verified fires inside Kerala boundaries:", nrow(fires_in_kerala)))



# --- MASTER SCRIPT: STEP 3 (DEDUPLICATION AND MATRIX BUILDING) ---
library(tidyr)
library(lubridate)

# 1. Extract Date and Coordinates to find the duplicates
fires_clean <- fires_in_kerala %>%
  mutate(
    Date = as.Date(acq_date),
    # Extract X and Y (Easting/Northing) from the UTM geometry
    X_coord = st_coordinates(.)[,1],
    Y_coord = st_coordinates(.)[,2],
    # Round to nearest 1000 meters (1km) to catch MODIS/VIIRS overlap
    X_grid = round(X_coord, -3),
    Y_grid = round(Y_coord, -3)
  )

# 2. Delete the duplicates (Same District, Same 1km Grid, Same Day)
fires_dedup <- fires_clean %>%
  distinct(DISTRICT, Date, X_grid, Y_grid, .keep_all = TRUE)

print(paste("Fires after removing same-day satellite duplicates:", nrow(fires_dedup)))

# 3. Convert Dates to Year-Week format for the hhh4 model
fires_dedup <- fires_dedup %>%
  mutate(
    Year = isoyear(Date),
    Week = isoweek(Date),
    # Create a clean label like "2021-W05"
    YearWeek = paste0(Year, "-W", sprintf("%02d", Week))
  )

# 4. Count the distinct fires per District per Week
fire_counts <- fires_dedup %>%
  st_drop_geometry() %>%
  group_by(DISTRICT, Region_ID, YearWeek) %>%
  summarise(FireCount = n(), .groups = "drop")

# 5. Build the strict 5-year timeline (The Zero-Fill)
# 5 years * 52 weeks = 260 weeks. 260 weeks * 14 districts = 3640 rows.
all_weeks <- expand_grid(
  DISTRICT = unique(kerala_only$DISTRICT),
  Year = 2021:2025,
  Week = 1:52
) %>%
  mutate(YearWeek = paste0(Year, "-W", sprintf("%02d", Week))) %>%
  # Join with the Region IDs so the matrix stays organized
  left_join(st_drop_geometry(kerala_only) %>% select(DISTRICT, Region_ID), by = "DISTRICT")

# 6. Merge the actual fire counts into the blank calendar
final_matrix <- all_weeks %>%
  left_join(fire_counts %>% select(DISTRICT, YearWeek, FireCount), by = c("DISTRICT", "YearWeek")) %>%
  # Replace all the "NA" blanks with actual zeros
  mutate(FireCount = replace_na(FireCount, 0)) %>%
  arrange(Region_ID, Year, Week)

print("Target Matrix built. Preview of the first 15 rows:")
print(head(final_matrix, 15))


# --- MASTER SCRIPT: STEP 4 (SPATIAL ADJACENCY MATRIX) ---

# 1. Create the neighborhood list (calculating which district borders which)
# We use queen = TRUE so even touching at a single corner counts as a shared border
district_neighbors <- poly2nb(kerala_only, queen = TRUE)

# 2. Convert the neighborhood list into a binary mathematical matrix
# 1 means they share a border, 0 means they do not
adjacency_matrix <- nb2mat(district_neighbors, style = "B", zero.policy = TRUE)

# 3. Assign the exact district names to the rows and columns so the hhh4 model can read it
rownames(adjacency_matrix) <- kerala_only$DISTRICT
colnames(adjacency_matrix) <- kerala_only$DISTRICT

print("Spatial Adjacency Matrix successfully built. Here is a preview of the borders:")
print(adjacency_matrix[1:5, 1:5])


# --- MASTER SCRIPT: STEP 5 (BUILDING THE STS ENGINE) ---

library(surveillance)

# 1. Pivot the counts from a long list into a wide mathematical matrix
# Rows = 260 Weeks, Columns = 14 Districts
fire_matrix_wide <- final_matrix %>%
  select(YearWeek, DISTRICT, FireCount) %>%
  pivot_wider(names_from = DISTRICT, values_from = FireCount) %>%
  arrange(YearWeek)

# 2. Strip the text labels to create a pure numeric matrix for the algorithm
# We drop the YearWeek column so only the 14 columns of numbers remain
observed_counts <- as.matrix(fire_matrix_wide[, -1])

# 3. CRITICAL: Ensure the columns of the count matrix perfectly match the adjacency matrix
# If Alappuzha is column 1 in the counts, it MUST be row/col 1 in the adjacency matrix
observed_counts <- observed_counts[, rownames(adjacency_matrix)]

# 4. Compile the actual sts (Surveillance Time Series) object
fire_sts <- sts(
  observed = observed_counts,
  start = c(2021, 1), # Starts Year 2021, Week 1
  frequency = 52,     # 52 weeks in a year
  neighbourhood = adjacency_matrix
)

# 5. Visual Verification of the Engine
# This will generate a massive plot showing the timeline of fires for every single district
#plot(fire_sts, type = observed ~ time | unit, 
     ylab = "Fire Intensity (Count)", 
     main = "Kerala Spatiotemporal Fire Dynamics (2021-2025)")

print("STS Object successfully compiled. The foundation is complete.")



# --- MASTER SCRIPT: STEP 6 (AUTOMATED ENVIRONMENTAL COVARIATES) ---

install.packages('MODISTools')
library(MODISTools)
library(sf)
library(dplyr)

# 1. Find the exact center (centroid) of every district
# We temporarily project back to WGS84 (EPSG:4326) because the NASA API requires standard Lat/Lon
districts_gps <- st_transform(kerala_only, crs = 4326)

# suppressWarnings is used here because st_centroid on curved coordinates throws a standard math warning we can ignore for this API
suppressWarnings({
  centroids <- st_centroid(districts_gps)
})

# 2. Extract the coordinates into a clean table required by the MODISTools API
api_targets <- data.frame(
  site_name = centroids$DISTRICT,
  lat = st_coordinates(centroids)[,2],
  lon = st_coordinates(centroids)[,1]
)

print("Target Coordinates successfully generated:")
print(api_targets)

# 3. Execute the API call for 5 Years of Land Surface Temperature (MOD11A2)
print("Pinging NASA ORNL Server for 5 years of LST data... This will take a minute. Do not interrupt R.")

lst_data <- mt_batch_subset(
  df = api_targets,
  product = "MOD11A2",
  band = "LST_Day_1km",
  start = "2021-01-01",
  end = "2025-12-31",
  internal = TRUE
)

# 4. Clean the downloaded data (DO NOT RE-RUN THE DOWNLOAD STEP)
# MODIS scales temperature data by 0.02. We multiply it to get Kelvin, then subtract 273.15 to get Celsius.
lst_clean <- lst_data %>%
  filter(value != -3000) %>% # Remove NASA error/cloud codes
  mutate(
    Date = as.Date(calendar_date),
    Temp_Celsius = (value * 0.02) - 273.15
  ) %>%
  select(site, Date, Temp_Celsius) %>%
  rename(DISTRICT = site) # Rename immediately to match our foundation

print("Data Cleaned and Converted to Celsius. Preview:")
print(head(lst_clean))




# This will physically save the data as a CSV file into your D:/icti/fire folder
library(readr)
write_csv(lst_clean, "Kerala_LST_2021_2025.csv")

print("Data successfully saved to your root folder as Kerala_LST_2021_2025.csv")




# --- MASTER SCRIPT: STEP 7 (ALIGNING LST TO THE WEEKLY ENGINE - FINAL FIX) ---
library(tidyr)
library(lubridate)

# 1. The Strict Biological Filter
# You correctly identified cloud contamination. We are terminating any reading below 15°C.
lst_weekly <- lst_clean %>%
  filter(Temp_Celsius > 15) %>% 
  mutate(
    Year = isoyear(Date),
    Week = isoweek(Date),
    YearWeek = paste0(Year, "-W", sprintf("%02d", Week))
  ) %>%
  group_by(DISTRICT, YearWeek) %>%
  summarise(Mean_LST = mean(Temp_Celsius, na.rm = TRUE), .groups = "drop")

# 2. Merge this into our perfect 5-year, 260-week calendar
lst_matrix_long <- all_weeks %>%
  left_join(lst_weekly, by = c("DISTRICT", "YearWeek")) %>%
  arrange(DISTRICT, Year, Week) %>%
  group_by(DISTRICT) %>%
  # Fill the gaps left by the NASA 8-day cycle AND our deleted cloud anomalies
  fill(Mean_LST, .direction = "downup") %>%
  ungroup()

# 3. Pivot into the wide mathematical matrix required by the hhh4 model
lst_matrix_wide <- lst_matrix_long %>%
  select(YearWeek, DISTRICT, Mean_LST) %>%
  pivot_wider(names_from = DISTRICT, values_from = Mean_LST) %>%
  arrange(YearWeek)

# 4. Strip text labels to create the pure numeric covariate matrix
lst_covariate <- as.matrix(lst_matrix_wide[, -1])

# 5. CRITICAL: Ensure the columns perfectly match the adjacency matrix
lst_covariate <- lst_covariate[, rownames(adjacency_matrix)]

print("Biologically strict LST Matrix built. The new minimum temperature is:")
print(paste(min(lst_covariate), "Celsius"))





# --- MASTER SCRIPT: STEP 8 (FITTING THE HHH4 MODEL - FIXED PARAMETER) ---
library(surveillance)
# 1. Scale the temperature matrix so the MLE optimization doesn't fail
lst_scaled <- lst_covariate / 10

# 2. Define the absolute mathematical parameters correctly
model_control_lst <- list(
  # The surveillance package strictly requires this to be called 'end', not 'endemic'
  end = list(f = ~ 1 + lst_scaled),
  
  ar = list(f = ~ 1),
  
  ne = list(f = ~ 1, weights = neighbourhood(fire_sts)),
  
  data = list(lst_scaled = lst_scaled)
)

print("Re-running MLE with the correctly named 'end' parameter...")

# 3. Run the forced optimization
fire_model_LST <- hhh4(fire_sts, control = model_control_lst)

# 4. Print the true statistical summary
summary(fire_model_LST)



# --- MASTER SCRIPT: STEP 9 (AUTOMATED VEGETATION COVARIATE - NDVI) ---

# 1. Execute the API call for 5 Years of Vegetation Data (MOD13Q1)
print("Pinging NASA ORNL Server for 5 years of NDVI data... This will take a few minutes.")

ndvi_data <- mt_batch_subset(
  df = api_targets,
  product = "MOD13Q1",
  band = "250m_16_days_NDVI",
  start = "2021-01-01",
  end = "2025-12-31",
  internal = TRUE
)

# 2. Clean the downloaded NDVI data
# NASA scales NDVI by 0.0001. Valid NDVI values range from -1 (water) to 1 (dense green forest).
ndvi_clean <- ndvi_data %>%
  filter(value >= -10000) %>% # Remove strict error codes
  mutate(
    Date = as.Date(calendar_date),
    NDVI_Actual = value * 0.0001
  ) %>%
  select(site, Date, NDVI_Actual) %>%
  rename(DISTRICT = site)

# 3. Convert the 16-day NASA dates into our strict Year-Week format
ndvi_weekly <- ndvi_clean %>%
  mutate(
    Year = isoyear(Date),
    Week = isoweek(Date),
    YearWeek = paste0(Year, "-W", sprintf("%02d", Week))
  ) %>%
  group_by(DISTRICT, YearWeek) %>%
  summarise(Mean_NDVI = mean(NDVI_Actual, na.rm = TRUE), .groups = "drop")

# 4. Merge this into our perfect 5-year, 260-week calendar
ndvi_matrix_long <- all_weeks %>%
  left_join(ndvi_weekly, by = c("DISTRICT", "YearWeek")) %>%
  arrange(DISTRICT, Year, Week) %>%
  group_by(DISTRICT) %>%
  # Fill the gaps left by the 16-day interval
  fill(Mean_NDVI, .direction = "downup") %>%
  ungroup()

# 5. Pivot into the wide mathematical matrix required by the hhh4 model
ndvi_matrix_wide <- ndvi_matrix_long %>%
  select(YearWeek, DISTRICT, Mean_NDVI) %>%
  pivot_wider(names_from = DISTRICT, values_from = Mean_NDVI) %>%
  arrange(YearWeek)

# 6. Strip text labels to create the pure numeric covariate matrix
ndvi_covariate <- as.matrix(ndvi_matrix_wide[, -1])

# 7. CRITICAL: Ensure the columns perfectly match the exact order of the spatial adjacency matrix
ndvi_covariate <- ndvi_covariate[, rownames(adjacency_matrix)]

print("NDVI Covariate Matrix successfully built. Preview of the first 5 weeks:")
print(ndvi_covariate[1:5, 1:5])




# --- MASTER SCRIPT: STEP 10 (THE FINAL MULTIVARIATE HHH4 MODEL) ---

# 1. Define the absolute mathematical parameters incorporating BOTH climate drivers
model_control_final <- list(
  # The endemic baseline is now driven by Heat (lst_scaled) AND Fuel Dryness (ndvi_covariate)
  end = list(f = ~ 1 + lst_scaled + ndvi_covariate),
  
  # Autoregressive risk (fire persisting in the same district)
  ar = list(f = ~ 1),
  
  # Epidemic risk (fire spreading across the physical borders)
  ne = list(f = ~ 1, weights = neighbourhood(fire_sts)),
  
  # Explicitly map BOTH matrices to the algorithm's data environment
  data = list(
    lst_scaled = lst_scaled,
    ndvi_covariate = ndvi_covariate
  )
)


                               

print("Running the Final Multivariate Maximum Likelihood Estimation...")

# 2. Run the forced optimization
fire_model_FINAL <- hhh4(fire_sts, control = model_control_final)

# 3. Print the true statistical summary
print("Final Model Fitting Complete. Here is the absolute mathematical truth:")
summary(fire_model_FINAL)





# --- MASTER SCRIPT: STEP 11 (VISUALIZING THE PREDICTIVE ENGINE - FIXED) ---

# 1. Plot the state-wide predicted vs observed fires (letting the package do its default rendering)
#plot(fire_model_FINAL, type = "fitted", total = TRUE)

# Now we stamp the custom title on top of the finished plot
#title(main = "Kerala Forest Fires: Predicted Risk vs Actual Outbreaks (2021-2025)")

# 2. Plot the exact risk decomposition for a specific high-risk district like Wayanad
#plot(fire_model_FINAL, type = "fitted", units = "Wayanad")

# Stamp the title for the district-specific plot
#title(main = "Wayanad: Environmental Endemic Risk vs Epidemic Spread")




# --- MASTER SCRIPT: STEP 12 (FREEZING THE BRAIN FOR THE DASHBOARD) ---

print("Saving the entire mathematical environment...")
save.image("kerala_fire_model.RData")
print("SUCCESS: kerala_fire_model.RData has been created in your project folder.")


