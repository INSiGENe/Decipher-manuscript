import os
import pandas as pd

# Define the base directory
base_dir = 'manuscript_analysis/sepsis/cytosig'

# Define output paths for P values and Z scores
p_value_output_path = 'manuscript_analysis/sepsis/cytosig/0_outputs/p_value'
z_score_output_path = 'manuscript_analysis/sepsis/cytosig/0_outputs/z_score'

# Create the directories if they don't exist
os.makedirs(p_value_output_path, exist_ok=True)
os.makedirs(z_score_output_path, exist_ok=True)

# Function to process files in each folder
def process_files_in_folder(folder_path):
    # Define the filenames
    file_types = [('output.Pvalue', p_value_output_path), ('output.Zscore', z_score_output_path)]

    # Get the folder name as cell type
    cell_type = os.path.basename(folder_path)

    # Loop through each file type and process it
    for file_name, output_path in file_types:
        file_path = os.path.join(folder_path, file_name)

        # Check if the file exists
        if os.path.exists(file_path):
            # Read the file into a pandas DataFrame
            df = pd.read_csv(file_path, sep='\t', index_col=0, header=0)
            
            # Create the output filename
            output_filename = os.path.join(output_path, f"{cell_type}.csv")
            
            # Save the DataFrame as a CSV file
            df.to_csv(output_filename, index=True)

# Loop through each folder in the base directory
for folder_name in os.listdir(base_dir):
    folder_path = os.path.join(base_dir, folder_name)

    # Check if it's a directory
    if os.path.isdir(folder_path):
        process_files_in_folder(folder_path)

print("Files processed and saved successfully.")
