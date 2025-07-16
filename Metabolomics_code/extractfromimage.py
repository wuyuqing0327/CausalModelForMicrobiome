import os
import pytesseract
from PIL import Image
import pandas as pd
import numpy as np

# Specify the folder containing the images
folder_path = r'C:\Users\yuqingw1\Workfolder\Data\Metabolomics\MSC1445_Sample_ID'  # replace with your folder path

# List all PNG files in the folder
image_paths = [os.path.join(folder_path, f) for f in os.listdir(folder_path) if f.endswith('.png')]

# Extract text from each image
extracted_text = []
for path in image_paths:
    image = Image.open(path)
    text = pytesseract.image_to_string(image)
    extracted_text.append(text)
    print(f"Extracted text from {path}")

# Combine extracted text into a single list of lines
all_text_lines = "\n".join(extracted_text).split("\n")

# Parse the lines into a list of lists for DataFrame
data = []
for line in all_text_lines:
    columns = line.split()
    if len(columns) == 3:  # Assuming the table has 3 columns: Index, Sample_ID, Specimen_ID
        data.append(columns)

# Convert to DataFrame
df = pd.DataFrame(data, columns=["Index", "Sample_ID", "Specimen_ID"])

# Save to Excel
excel_path = r"C:\Users\yuqingw1\Workfolder\Data\Metabolomics\MSC1445_Sample_ID_from_images.xlsx"
df.to_excel(excel_path, index=False)

