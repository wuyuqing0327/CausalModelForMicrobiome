import requests
from bs4 import BeautifulSoup
import time
import csv
import pandas as pd
import numpy as np
import xml.etree.ElementTree as ET

data = pd.read_csv(r'C:\Users\yuqingw1\Workfolder\ProcessingData\MSC1445_list_filterID.csv')
# List of HMP IDs to search
kegg_ids = list(data[~pd.isnull(data['KEGG ID'])]['KEGG ID']) # Add more HMP IDs as needed

# Base URL format for KEGG metabolite pages
base_url = 'https://www.kegg.jp/entry/'

# Output file to store results
output_file = r'C:\Users\yuqingw1\Workfolder\ProcessingData\MSC1445_list_KEGGID.csv'


# Initialize a list to store metabolite data
metabolites_data = []

for kegg_id in kegg_ids:
    # Construct the URL to access the XML file
    url = f"{base_url}{kegg_id}"
    response = requests.get(url)

    # Check if the request was successful
    if response.status_code == 200:
        soup = BeautifulSoup(response.content, 'html.parser')


# Extract information
def extract_info(soup):
    data = {}

    # Compound Name
    name_tag = soup.find("th", string="Name")
    if name_tag:
        name = name_tag.find_next("td").text.strip().replace("\n", ", ")
        data["Name"] = name

    # Formula
    formula_tag = soup.find("th", string="Formula")
    if formula_tag:
        formula = formula_tag.find_next("td").text.strip()
        data["Formula"] = formula

    # Exact Mass
    exact_mass_tag = soup.find("th", string="Exact mass")
    if exact_mass_tag:
        exact_mass = exact_mass_tag.find_next("td").text.strip()
        data["Exact Mass"] = exact_mass

    # Molecular Weight
    mol_weight_tag = soup.find("th", string="Mol weight")
    if mol_weight_tag:
        mol_weight = mol_weight_tag.find_next("td").text.strip()
        data["Molecular Weight"] = mol_weight

    # Pathways
    pathways = []
    pathway_tag = soup.find("th", string="Pathway")
    if pathway_tag:
        pathway_entries = pathway_tag.find_next("td").find_all("tr")
        for entry in pathway_entries:
            pathways.append(entry.get_text(separator=": ").strip())
        data["Pathways"] = pathways

    # Modules
    modules = []
    module_tag = soup.find("th", string="Module")
    if module_tag:
        module_entries = module_tag.find_next("td").find_all("tr")
        for entry in module_entries:
            modules.append(entry.get_text(separator=": ").strip())
        data["Modules"] = modules

    # Enzymes
    enzymes = []
    enzyme_tag = soup.find("th", string="Enzyme")
    if enzyme_tag:
        enzyme_entries = enzyme_tag.find_next("td").find_all("a")
        enzymes = [entry.text for entry in enzyme_entries]
        data["Enzymes"] = enzymes

    return data

# Fetch and display data
compound_info = extract_info(soup)
for key, value in compound_info.items():
    print(f"{key}: {value}")