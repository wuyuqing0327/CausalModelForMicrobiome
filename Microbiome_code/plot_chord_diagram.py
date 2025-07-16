import pandas as pd
import holoviews as hv
from holoviews import opts
import hvplot.pandas
from sklearn.preprocessing import MinMaxScaler

hv.extension('bokeh')

path = r'C:\Users\yuqingw1\Workfolder\result'
finalvar_2024 = ['UPDASV323','UPDASV028','UPDASV319','UPDASV019','UPDASV152','UPDASV200','UPDASV018','UPDASV322','UPDASV204','UPDASV093','UPDASV394','UPDASV317','UPDASV182','UPDASV244','UPDASV294','UPDASV110','UPDASV092','UPDASV203','UPDASV144','UPDASV280','UPDASV438','UPDASV428','UPDASV032','UPDASV176','UPDASV040','UPDASV458','UPDASV408','UPDASV047','UPDASV214','UPDASV219','UPDASV442','UPDASV445','UPDASV157','UPDASV486','UPDASV468','UPDASV391','UPDASV229','UPDASV036','UPDASV189','UPDASV338','UPDASV085','UPDASV403','UPDASV215','UPDASV413','UPDASV190','UPDASV002','UPDASV003','UPDASV004','UPDASV011','UPDASV041','UPDASV055','UPDASV056','UPDASV068','UPDASV069','UPDASV079','UPDASV107','UPDASV116','UPDASV132','UPDASV139','UPDASV146','UPDASV148','UPDASV149','UPDASV154','UPDASV156','UPDASV159','UPDASV168','UPDASV171','UPDASV198','UPDASV210','UPDASV224','UPDASV225','UPDASV227','UPDASV228','UPDASV238','UPDASV240','UPDASV246','UPDASV259','UPDASV264','UPDASV270','UPDASV272','UPDASV276','UPDASV283','UPDASV298','UPDASV301','UPDASV302','UPDASV306','UPDASV312','UPDASV318','UPDASV334','UPDASV342','UPDASV346','UPDASV348','UPDASV353','UPDASV355','UPDASV358','UPDASV360','UPDASV362','UPDASV374','UPDASV377','UPDASV382','UPDASV383','UPDASV384','UPDASV411','UPDASV414','UPDASV418','UPDASV426','UPDASV429','UPDASV431','UPDASV432','UPDASV433','UPDASV434','UPDASV435','UPDASV448','UPDASV450','UPDASV459','UPDASV460','UPDASV465','UPDASV490','UPDASV491','UPDASV385','UPDASV170','UPDASV043','UPDASV181','UPDASV213','UPDASV249','UPDASV396','UPDASV281','UPDASV404','UPDASV108','UPDASV057','UPDASV155','UPDASV099','UPDASV288','UPDASV284','UPDASV292','UPDASV286','UPDASV007','UPDASV153','UPDASV304','UPDASV172','UPDASV425','UPDASV255','UPDASV256','UPDASV439','UPDASV274','UPDASV117','UPDASV356','UPDASV380','UPDASV368','UPDASV397','UPDASV482','UPDASV120','UPDASV388','UPDASV050','UPDASV134','UPDASV236','UPDASV347','UPDASV039','UPDASV024','UPDASV437','UPDASV400','UPDASV001','UPDASV042','UPDASV044','UPDASV376','UPDASV066','UPDASV232','UPDASV098','UPDASV094','UPDASV212','UPDASV424','UPDASV010','UPDASV005','UPDASV410','UPDASV341','UPDASV378','UPDASV370','UPDASV105','UPDASV367','UPDASV217','UPDASV297','UPDASV033','UPDASV409','UPDASV035','UPDASV145','UPDASV013','UPDASV275','UPDASV366','UPDASV205','UPDASV014','UPDASV243','UPDASV381','UPDASV130','UPDASV100','UPDASV369','UPDASV054','UPDASV140','UPDASV489','UPDASV142','UPDASV267','UPDASV136','UPDASV354','UPDASV420','UPDASV416','UPDASV084','UPDASV352','UPDASV191','UPDASV359','UPDASV427','UPDASV009','UPDASV311','UPDASV049','UPDASV183','UPDASV021','UPDASV493','UPDASV165','UPDASV234','UPDASV162','UPDASV209','UPDASV278','UPDASV147','UPDASV034','UPDASV315','UPDASV295','UPDASV113','UPDASV096','UPDASV401','UPDASV071','UPDASV080','UPDASV412','UPDASV202','UPDASV492','UPDASV265','UPDASV316','UPDASV123','UPDASV020','UPDASV483','UPDASV336','UPDASV481','UPDASV186','UPDASV330','UPDASV313','UPDASV321','UPDASV070','UPDASV371','UPDASV199','UPDASV150','UPDASV444','UPDASV030','UPDASV257','UPDASV129','UPDASV417','UPDASV379','UPDASV363','UPDASV327','UPDASV485','UPDASV163','UPDASV067','UPDASV314','UPDASV051','UPDASV415','UPDASV103','UPDASV106','UPDASV081','UPDASV310','UPDASV451','UPDASV175','UPDASV329','UPDASV449','UPDASV473','UPDASV365']
# var_insulin_pm25_2025 =

varcols = finalvar_2024
vardf = pd.DataFrame(varcols, columns=['index'])

ghrelin = pd.read_csv(path + "\Biomarker_Ghrelin.csv")
insulin = pd.read_csv(path + "\Biomarker_Insulin.csv")
leptin = pd.read_csv(path + "\Biomarker_Leptin.csv")
resistin = pd.read_csv(path + "\Biomarker_Resistin.csv")
pm25 = pd.read_csv(path + "\PM25_M12result.csv")


data = vardf.merge(ghrelin[['index','coeff','pvalue', 'adjust_Pvalue']], on=['index']).merge(insulin[['index','coeff','pvalue', 'adjust_Pvalue']], on=['index']).\
merge(resistin[['index','coeff','pvalue', 'adjust_Pvalue']], on=['index']).merge(pm25[['index','M12_coeff','M12_pvalue', 'M12_adjpvalue']], on=['index'])
data.columns=['Index', 'Ghrelin_coef','Ghrelin_pvalue', 'Ghrelin_adjp',
              'Insulin_coef','Insulin_pvalue', 'Insulin_adjp',
              'Resistin_coef','Resistin_pvalue', 'Resistin_adjp',
              'PM25_coef','PM25_pvalue', 'PM25_adjp']


# Prepare the data for the MinMaxScaler
coefdata = data[['Index', 'PM25_coef', 'Ghrelin_coef', 'Resistin_coef', 'Insulin_coef']]
coefdata = coefdata.sort_values(['Index']).reset_index(drop=True).loc[:20]
coefdata.columns = ['Index','PM2.5', 'Ghrelin', 'Resistin', 'Insulin']

# Prepare the data for the chord diagram
chord_data = pd.melt(coefdata, id_vars="Index", var_name="Metric", value_name="Coefficient")
chord_data = chord_data[chord_data['Metric'].isin(['PM2.5', 'Ghrelin', 'Resistin', 'Insulin'])]
print("chord_data:", chord_data)

# Ensure the links DataFrame is correctly structured
links = chord_data.rename(columns={'Index': 'source', 'Metric': 'target', 'Coefficient': 'value'})
# Debug: Print links to check structure
print("Links DataFrame:")
print("\nLinks columns:", links.columns)


# Ensure all nodes are included
nodes = pd.DataFrame(list(set(links['source']).union(set(links['target']))), columns=['name'])
# Create a 'group' column to distinguish between source and target nodes in the visualization
nodes['group'] = ['Microbiome' if x.startswith('UPDASV') else 'Air Pollution' if x.startswith('PM') else 'Biomarker' for x in nodes['name']]
# Debug: Print nodes with group to check structure
print("\nNodes DataFrame with group:")
print("\nNodes columns:", nodes.columns)

# Convert nodes DataFrame to Holoviews Dataset
nodes = hv.Dataset(nodes, 'name')
print('nodes:', nodes)

# Check if the source and target nodes are correctly linked
print("\nUnique sources:", links['source'].unique())
print("\nUnique targets:", links['target'].unique())
print("\nUnique node names:", nodes.data['name'].unique())


# Create the chord diagram
chord = hv.Chord((links, nodes)).opts(
    opts.Chord(cmap='Category20', edge_cmap='Category20', labels='name',
               node_color='name', edge_color='value', edge_alpha=0.8,
               node_alpha=0.8, node_size=10, title="Relationships between Microbiomes, PM2.5 and Biomarkers",
               edge_line_width=hv.dim('value')*5))


hv.save(chord, r'C:\Users\yuqingw1\Workfolder\result\chord_diagram.html', fmt='html')
# Save the plot to a PNG file
hv.save(chord, r'C:\Users\yuqingw1\Workfolder\result\chord_diagram.png', fmt='png')




### version2 -----discard
import pandas as pd
import matplotlib.pyplot as plt
import networkx as nx

# Prepare the data
coefdata = data[['Index', 'PM25_coef', 'Ghrelin_coef', 'Resistin_coef', 'Insulin_coef']]
coefdata = coefdata.sort_values(['Index']).reset_index(drop=True).loc[:20]
coefdata.columns = ['Index', 'PM2.5', 'Ghrelin', 'Resistin', 'Insulin']

# Prepare the data for the graph
chord_data = pd.melt(coefdata, id_vars="Index", var_name="Metric", value_name="Coefficient")
chord_data = chord_data[chord_data['Metric'].isin(['PM2.5', 'Ghrelin', 'Resistin', 'Insulin'])]

# Create the graph
G = nx.DiGraph()

# Add nodes and edges
for _, row in chord_data.iterrows():
    G.add_node(row['Index'], type='microbiome')
    G.add_node(row['Metric'], type='biomarker')
    G.add_edge(row['Index'], row['Metric'], weight=row['Coefficient'])

# Set node colors based on type
color_map = []
for node in G:
    if G.nodes[node]['type'] == 'microbiome':
        color_map.append('blue')
    else:
        color_map.append('red')

# Position nodes in a circle
pos = nx.circular_layout(G)

# Adjust positions to have two groups: microbiomes on the left and biomarkers on the right
for node in pos:
    if G.nodes[node]['type'] == 'microbiome':
        pos[node][0] = pos[node][0] - 1.5  # Shift microbiomes left
    else:
        pos[node][0] = pos[node][0] + 1.5  # Shift biomarkers right

# Draw the graph
plt.figure(figsize=(12, 12))

# Draw nodes
nx.draw_networkx_nodes(G, pos, node_color=color_map, node_size=500)

# Draw labels
nx.draw_networkx_labels(G, pos)

# Draw edges
edges = G.edges(data=True)
weights = [abs(edge[2]['weight']) * 5 for edge in edges]  # Adjust the width to represent the coefficient size
nx.draw_networkx_edges(G, pos, edgelist=edges, width=weights)

plt.title("Relationships between Microbiomes, PM2.5 and Biomarkers")
plt.show()