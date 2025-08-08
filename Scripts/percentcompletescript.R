## **Calcuation of completeness of individual cross-sections**

# Catherine Klesner
# 2025

#R packages are required: 

library(magick)
library(raster)

# Prior to running this script, images need to be processed using an image 
# processing software (Inkscape was used for this research) where the missing 
# components of the cross-section of the vessel are infilled with **RED** 
# following standard archaeological illustration practices. 

------------------ # Function to calculate the color counts --------------------
calculate_color_counts <- function(img_raster) {
  red_channel <- values(img_raster[[1]])
  green_channel <- values(img_raster[[2]])
  blue_channel <- values(img_raster[[3]])
  
  # Define conditions for red, black, and white pixels
  red_condition <- red_channel > (green_channel + 100) & red_channel > (blue_channel + 100)
  black_condition <- red_channel < 50 & green_channel < 50 & blue_channel < 50
  white_condition <- red_channel > 205 & green_channel > 205 & blue_channel > 205
  
  red_count <- sum(red_condition)
  black_count <- sum(black_condition)
  white_count <- sum(white_condition)
  
  # Calculate the proportion of red pixels to red+black pixels
  total_red_black <- red_count + black_count
  if (total_red_black == 0) {
    proportion_red_to_red_black <- 0
  } else {
    proportion_red_to_red_black <- red_count / total_red_black
  }
  
  # Create a dataframe with the results
  result_df <- data.frame(
    red_count = red_count,
    black_count = black_count,
    white_count = white_count,
    proportion_red_to_red_black = proportion_red_to_red_black
  )
  
  # Return the dataframe
  return(result_df)
}

----------------------# Function to process each image--------------------------
process_image <- function(image_file) {
  tryCatch({
    img <- image_read(image_file)
    
    temp_file <- tempfile(fileext = ".png")
    image_write(img, temp_file)
    
    img_raster <- stack(temp_file)
    
    color_counts_and_proportion <- calculate_color_counts(img_raster)
    
    color_counts_and_proportion$file_name <- basename(image_file)
    
    return(color_counts_and_proportion)
  }, error = function(e) {
    cat("Error processing file:", image_file, "- ", conditionMessage(e), "\n")
    return(NULL)
  })
}


# ---------------------------Process cross-sections ---------------------------
 

# Path to cross-sections folder
image_directory <- "./vessel_cross_sections"

# Create a list of image paths
image_files <- list.files(image_directory, pattern = "\\.(jpg|jpeg|tiff)$", full.names = TRUE)

# Apply the processing function to each image
results_list <- lapply(image_files, process_image)

#you may have warnings - ignore them initially

# Combine the results into a single dataframe
combined_results <- do.call(rbind, results_list)

# Write the combined results to a CSV file
write.csv(combined_results, "./Scripts/color_counts_combined.csv", row.names = FALSE)




