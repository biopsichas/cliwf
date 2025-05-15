## This a list to unzip or preprocess large data files, to make it possible 
## to use github.com and git

# Check if farm.csv exists
if (!file.exists(mgt)) {
  message("farmR_input.csv not found.")
  
  # Load required package for extracting .7z
  if (!requireNamespace("archive", quietly = TRUE)) {
    install.packages("archive")
  }
  library(archive)
  # Check if farm.7z exists
  if (file.exists("data/farmR_input.7z")) {
    message("farmR_input.7z found. Extracting...")
    
    # Extract farm.7z
    archive::archive_extract("data/farmR_input.7z", "data")
    
    # Check again if extraction was successful
    if (file.exists(mgt)) {
      message("Extraction successful. farmR_input.csv is now available.")
    } else {
      message("Extraction failed or farmR_input.csv not found inside the archive.")
    }
  } else {
    message("data/farmR_input.7z not found.")
  }
}
