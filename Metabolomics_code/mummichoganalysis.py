from mummichog1.functional_analysis import PathwayAnalysis, HsaNetwork
import mummichog1.pydata.human_model_mfn as MetabolicModel  # Import the actual model

import mummichog1
print(mummichog1.__file__)

# Define the options with the required keys
options = {
    'input_file': r'C:\Users\yuqingw1\Workfolder\result\Metobolomics\mummichog\microbiome_001.txt',
    'infile': r'C:\Users\yuqingw1\Workfolder\result\Metobolomics\mummichog\microbiome_001.txt',
    'output_directory': r'C:\Users\yuqingw1\Workfolder\result\Metobolomics\mummichogresult',
    'mode': 'positive',
    'p_value_cutoff': 0.05,
    'network': 'human_mfn',
    'workdir': r'C:\Users\yuqingw1\Workfolder\result\Metobolomics',
    'cutoff': 0.05,  # Add the cutoff key with the appropriate value
    'instrument': 10,  # Set instrument to ppm value or a string like 'FTMS'
    'permutation': 1000
}

# Initialize the network with the correct MetabolicModel
metabolic_model = HsaNetwork(MetabolicModel, options)

print(type(metabolic_model))
print(type(metabolic_model.input_mzlist))

if hasattr(metabolic_model, 'network') and metabolic_model.network:
    print("Network initialized successfully.")
else:
    raise ValueError("Failed to initialize the network.")

metabolic_model.network.paradict = metabolic_model.paradict
metabolic_model.network.ref_mzlist = list(metabolic_model.input_mzlist)  # Or another appropriate default value
metabolic_model.input_cpdlist = metabolic_model.total_matched_cpds  # Or another appropriate list
metabolic_model.network.total_matched_cpds = []  # Example: This should be a list of matched compounds

def count_cpd2mz(self, overlap_features):
    count = 0
    for feature in overlap_features:
        if feature in self.cpd_dict:
            count += len(self.cpd_dict[feature])  # Adjust as necessary
    return count

# Assign this function to your HsaNetwork instance
setattr(metabolic_model, 'count_cpd2mz', count_cpd2mz.__get__(metabolic_model, HsaNetwork))


# Perform pathway analysis
pathway_analysis = PathwayAnalysis(metabolic_model.network, metabolic_model)
pathway_analysis.run()
# pathway_analysis.cpd_enrich_test()
# pathway_analysis.run_all_analysis()
pathway_analysis.export_results(options['output_directory'])





