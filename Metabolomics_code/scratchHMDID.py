import requests
from bs4 import BeautifulSoup
import time
import csv
import pandas as pd
import numpy as np
import xml.etree.ElementTree as ET

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


        # Extract relevant information
        metabolite_data = {
            "Version": soup.find(text="Version").find_next("td").text,
            "Creation Date": soup.find(text="Creation Date").find_next("td").text,
            "Update Date": soup.find(text="Update Date").find_next("td").text,
            "Accession": soup.find(text="HMDB ID").find_next("td").text,
            "Status": soup.find(text="Status").find_next("td").text,
            "Secondary Accessions": [acc.text for acc in
                                     soup.find(text="Secondary Accession Numbers").find_next("ul").find_all("li")],
            "Name": soup.find(text="Common Name").find_next("td").text,
            "Description": soup.find(text="Description").find_next("td").text,
            "Synonyms": [syn.text for syn in soup.find(text="Synonyms").find_next("ul").find_all("li")],
            "Chemical Formula": soup.find(text="Chemical Formula").find_next("td").text,
            "Average Molecular Weight": soup.find(text="Average Molecular Weight").find_next("td").text,
            "Monoisotopic Molecular Weight": soup.find(text="Monoisotopic Molecular Weight").find_next("td").text,
            "IUPAC Name": soup.find(text="IUPAC Name").find_next("td").text,
            "Traditional Name": soup.find(text="Traditional Name").find_next("td").text,
            "CAS Registry Number": soup.find(text="CAS Registry Number").find_next("td").text,
            "SMILES": soup.find(text="SMILES").find_next("td").text,
            "InChI": soup.find(text="InChI").find_next("td").text,
            "InChI Key": soup.find(text="InChI Key").find_next("td").text
        }

        # Taxonomy Section
        taxonomy_section = soup.find("tr", id="taxonomy")
        if taxonomy_section:
            metabolite_data["Classification"] = taxonomy_section.find_next("tr").find(
                "td").text.strip() if taxonomy_section.find_next("tr").find("td") else None
        else:
            metabolite_data["Classification"] = None

        # Physical Properties Section
        physical_properties_section = soup.find("tr", id="physical_properties")
        if physical_properties_section:
            state_td = physical_properties_section.find_next("td")
            metabolite_data["State"] = state_td.text.strip() if state_td else "Not Available"
            experimental_properties = physical_properties_section.find_all("tbody")
            if experimental_properties:
                metabolite_data["Experimental Molecular Properties"] = [
                    {
                        "Property": cells[0].text.strip() if len(cells) > 0 else "Not Available",
                        "Value": cells[1].text.strip() if len(cells) > 1 else "Not Available",
                        "Reference": cells[2].text.strip() if len(cells) > 2 else "Not Available"
                    }
                    for row in experimental_properties[0].find_all("tr")
                    if (cells := row.find_all("td"))  # Assign cells to the result of row.find_all("td")
                ]

        # Predicted Molecular Properties Section
        predicted_properties_section = soup.find(text="Predicted Molecular Properties")
        if predicted_properties_section:
            predicted_properties_table = predicted_properties_section.find_next("tbody")
            if predicted_properties_table:
                metabolite_data["Predicted Molecular Properties"] = [
                    {
                        "Property": cells[0].text.strip() if len(cells) > 0 else "Not Available",
                        "Value": cells[1].text.strip() if len(cells) > 1 else "Not Available",
                        "Source": cells[2].text.strip() if len(cells) > 2 else "Not Available"
                    }
                    for row in predicted_properties_table.find_all("tr")
                    if (cells := row.find_all("td"))  # Assign cells to the result of row.find_all("td")
                ]

        # Chromatographic Properties Section
        chromatographic_properties_section = soup.find("tr", id="chromatographic_properties")
        if chromatographic_properties_section:
            chromatographic_properties_table = chromatographic_properties_section.find_next("tbody")
            if chromatographic_properties_table:
                metabolite_data["Predicted Chromatographic Properties"] = [
                    {
                        "Predictor": cells[0].text.strip() if len(cells) > 0 else "Not Available",
                        "Adduct Type": cells[1].text.strip() if len(cells) > 1 else "Not Available",
                        "CCS Value (Å²)": cells[2].text.strip() if len(cells) > 2 else "Not Available",
                        "Reference": cells[3].text.strip() if len(cells) > 3 else "Not Available"
                    }
                    for row in chromatographic_properties_table.find_all("tr")
                    if (cells := row.find_all("td"))  # Assign cells to the result of row.find_all("td")
                ]

        # Biological Properties Section
        biological_properties_section = soup.find("tr", id="biological_properties")
        if biological_properties_section:
            cellular_locations = biological_properties_section.find_next("ul")
            if cellular_locations:
                metabolite_data["Cellular Locations"] = [loc.text.strip() for loc in cellular_locations.find_all("li")]

            biospecimen_locations = biological_properties_section.find_next("ul")
            if biospecimen_locations:
                metabolite_data["Biospecimen Locations"] = [loc.text.strip() for loc in biospecimen_locations.find_all("li")]

            tissue_locations = biological_properties_section.find_next("ul")
            if tissue_locations:
                metabolite_data["Tissue Locations"] = [loc.text.strip() for loc in tissue_locations.find_all("li")]


        # Concentrations Section
        concentrations_section = soup.find("tr", id="concentrations")
        if concentrations_section:
            concentration_table = concentrations_section.find_next("tbody")
            if concentration_table:
                metabolite_data["Normal Concentrations"] = [
                    {
                        "Biospecimen": cells[0].text.strip() if len(cells) > 0 else "Not Available",
                        "Status": cells[1].text.strip() if len(cells) > 1 else "Not Available",
                        "Value": cells[2].text.strip() if len(cells) > 2 else "Not Available",
                        "Age": cells[3].text.strip() if len(cells) > 3 else "Not Available",
                        "Sex": cells[4].text.strip() if len(cells) > 4 else "Not Available",
                        "Condition": cells[5].text.strip() if len(cells) > 5 else "Not Available",
                        "Reference": cells[6].text.strip() if len(cells) > 6 else "Not Available"
                    }
                    for row in concentration_table.find_all("tr")
                    if (cells := row.find_all("td"))  # Assign cells to the result of row.find_all("td")
                ]

        # Associated Disorders and Diseases Section
        associated_disorders_section = soup.find("tr", id="associated-disorders-and-diseases")
        if associated_disorders_section:
            associated_diseases_table = associated_disorders_section.find_next("tbody")
            if associated_diseases_table:
                metabolite_data["Associated Disorders"] = [
                    {
                        "Disease": row.find("th").text.strip() if row.find("th") else "Not Available",
                        "References": [ref.text.strip() for ref in row.find("td").find_all("li")] if row.find(
                            "td") else []
                    }
                    for row in associated_diseases_table.find_all("tr")
                ]

        # External Links Section
        links_section = soup.find("tr", id="links")
        if links_section:
            external_links_table = links_section.find_next("tbody")
            if external_links_table:
                metabolite_data["External Links"] = [
                    {
                        "Type": row.find("th").text.strip() if row.find("th") else "Not Available",
                        "Link": row.find("td").text.strip() if row.find("td") else "Not Available"
                    }
                    for row in external_links_table.find_all("tr")
                ]

        # Enzymes Section
        enzymes_section = soup.find("h2", id="enzymes")
        if enzymes_section:
            enzymes_list = enzymes_section.find_next("div", class_="enzyme").find_all("div", class_="panel panel-default")
            metabolite_data["Enzymes"] = [
                {
                    "Name": enzyme.find("strong").find("a").text.strip(),
                    "Function": enzyme.find("dd").text.strip(),
                    "Gene Name": enzyme.find_all("dd")[1].text.strip(),
                    "Uniprot ID": enzyme.find("dd").find_next("a").text.strip()
                }
                for enzyme in enzymes_list
            ]

        # Append the data to the list
        metabolites_data.append(metabolite_data)

    else:
        print(f"Failed to fetch data for {hmdb_id}")


df = pd.DataFrame(metabolites_data)
df.to_csv(output_file, index=False)

print("Data saved to metabolites_identification_data.csv")
