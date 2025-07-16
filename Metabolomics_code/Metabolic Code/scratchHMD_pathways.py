import requests
from bs4 import BeautifulSoup
import time
import csv
import pandas as pd
import numpy as np
import xml.etree.ElementTree as ET


def find_section_by_heading(data, heading):
    """Recursively search for a section by TOCHeading."""
    if isinstance(data, dict):
        if data.get('TOCHeading') == heading:
            return data
        for key, value in data.items():
            result = find_section_by_heading(value, heading)
            if result:
                return result
    elif isinstance(data, list):
        for item in data:
            result = find_section_by_heading(item, heading)
            if result:
                return result
    return None


def extract_pathways(data):
    """Extract pathway names and URLs as a list of dictionaries for each metabolite."""
    # Find the 'Pharmacology and Biochemistry' section
    pharmacology_section = find_section_by_heading(data['Record'], 'Pharmacology and Biochemistry')
    # Find the 'Metabolite Pathways' section within 'Pharmacology and Biochemistry'
    pathways_section = find_section_by_heading(pharmacology_section, 'Metabolite Pathways')
    # Initialize a list to store pathways
    pathways_list = []
    # If the 'Metabolite Pathways' section exists, extract pathway names and URLs
    if pathways_section:
        for item in pathways_section['Information'][0]['Value']['StringWithMarkup']:
            pathway_name = item['String']
            # url = item['Markup'][0]['URL']
            # pathways_list.append({'Pathway Name': pathway_name, 'URL': url})
            pathways_list.append(pathway_name)

    return pathways_list



data = pd.read_csv(r'C:\Users\yuqingw1\Workfolder\ProcessingData\MSC1445_list_filterID.csv')
# List of HMP IDs to search
hmdb_ids = list(data[~pd.isnull(data['HMP ID'])]['HMP ID']) # Add more HMP IDs as needed


# Base URL for HMDB
base_url = "https://hmdb.ca/metabolites/"

# Output file to store results
output_file = r'C:\Users\yuqingw1\Workfolder\ProcessingData\MSC1445_list_HMPID.csv'

# Initialize a list to store metabolite data
metabolites_data = []

for hmdb_id in hmdb_ids:
    print("HMD ID:", hmdb_id)
    # Construct the URL to access the XML file
    url = f"{base_url}{hmdb_id}.xml"
    response = requests.get(url)

    # Check if the request was successful
    if response.status_code == 200:
        soup = BeautifulSoup(response.content, 'html.parser')
        smiles = soup.find('th', text='SMILES').find_next_sibling('td').get_text(strip=True)
        pubchem_url = f'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/{smiles}/cids/TXT'
        response = requests.get(pubchem_url)
        cid = response.text.strip()  # Extract the CID (PubChem ID)
        pathway_url = f'https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/{cid}/JSON'
        response = requests.get(pathway_url)
        if response.status_code == 200:
            data = response.json()
            metabolite_pathways = extract_pathways(data)
        else:
            metabolite_pathways = []
        # Extract relevant information
        metabolite_data = {
            "Creation Date": soup.find(text="Creation Date").find_next("td").text,
            "Update Date": soup.find(text="Update Date").find_next("td").text,
            "Accession": soup.find(text="HMDB ID").find_next("td").text,
            "Status": soup.find(text="Status").find_next("td").text,
            # "Secondary Accessions": [acc.text for acc in
            #                          soup.find(text="Secondary Accession Numbers").find_next("ul").find_all("li")],
            "Name": soup.find(text="Common Name").find_next("td").text,
            "Description": soup.find(text="Description").find_next("td").text,
            "Synonyms": [syn.text for syn in soup.find(text="Synonyms").find_next("ul").find_all("li")],
            "Chemical Formula": soup.find(text="Chemical Formula").find_next("td").text,
            "Average Molecular Weight": soup.find(text="Average Molecular Weight").find_next("td").text,
            "Monoisotopic Molecular Weight": soup.find(text="Monoisotopic Molecular Weight").find_next("td").text,
            # "IUPAC Name": soup.find(text="IUPAC Name").find_next("td").text,
            # "Traditional Name": soup.find(text="Traditional Name").find_next("td").text,
            # "CAS Registry Number": soup.find(text="CAS Registry Number").find_next("td").text,
            "SMILES": soup.find(text="SMILES").find_next("td").text,
            # "InChI": soup.find(text="InChI").find_next("td").text,
            # "InChI Key": soup.find(text="InChI Key").find_next("td").text,
            "Pathways": metabolite_pathways
        }

        # Append the data to the list
        metabolites_data.append(metabolite_data)

    else:
        print(f"Failed to fetch data for {hmdb_id}")


df = pd.DataFrame(metabolites_data)
df.to_csv(output_file, index=False)

print("Data saved to metabolites_identification_data.csv")













basic = pd.read_csv(r'C:\Users\yuqingw1\Workfolder\result\Metobolomics\mummichog\microbiome_001.txt', sep='\t')



