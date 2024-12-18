##Change climate file ending

# Load necessary library
library(fs)
# Find all files in subdirectories ending with .Tmp
files_to_rename <- dir(cli_dir, pattern = "\\.Tmp$", recursive = TRUE, full.names = TRUE)
# Check the files to be renamed
print(files_to_rename)
# Rename files by changing their extension to .tmp
for (file in files_to_rename) {
  new_name <- sub("\\.Tmp$", ".tmp", file)  # Change extension
  file.rename(file, new_name)              # Rename the file
}

# Confirm changes
print("Renaming completed.")