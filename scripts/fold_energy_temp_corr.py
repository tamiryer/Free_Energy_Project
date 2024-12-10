# ###works great! plotting using seabotn. mean (free energy + GC content)
import os
import re
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

# Directories and files
base_path = "/home/projects/zeevid/tamirye"
output_dir = "/home/projects/zeevid/Analyses/2023-Tamir/secondary_structures/10k"
catalog_file = "/home/projects/zeevid/tamirye/Data/Combined_datasets_full_modified.csv"
output_csv = f"{base_path}/merged_results.csv"

# Function to calculate mean free energy and mean GC content from a ViennaRNA output file
def parse_viennarna_file(file_path):
    gc_conts = []
    free_energies = []
    with open(file_path, 'r') as file:
        for line in file:
            if re.search(r'gc_cont=([\d\.]+)', line):
                gc_cont = float(re.findall(r'gc_cont=([\d\.]+)', line)[0])
                gc_conts.append(gc_cont)
            
            if re.search(r'\((\-?\d+\.\d+)\)', line):
                free_energy = float(re.findall(r'\((\-?\d+\.\d+)\)', line)[0])
                free_energies.append(free_energy)

    mean_free_energy = np.mean(free_energies) if free_energies else np.nan
    mean_gc_content = np.mean(gc_conts) if gc_conts else np.nan
    return mean_free_energy, mean_gc_content

# Process each output file and store results
results = []

for filename in os.listdir(output_dir):
    sample_name = filename.split('_')[0]
    file_path = os.path.join(output_dir, filename)
    mean_free_energy, mean_gc_content = parse_viennarna_file(file_path)
    results.append((sample_name, mean_free_energy, mean_gc_content))

# Create a DataFrame from results
results_df = pd.DataFrame(results, columns=["Sample", "MeanFreeEnergy", "MeanGCContent"])

# Load the catalog file
catalog_df = pd.read_csv(catalog_file)

# Split RunID by '|' and explode the DataFrame to create rows for each individual ID
catalog_df['RunID'] = catalog_df['RunID'].astype(str)  # Ensure RunID is a string
catalog_df_expanded = catalog_df.assign(RunID=catalog_df['RunID'].str.split('|')).explode('RunID')

# Filter catalog to remove samples with empty temperature values
catalog_df_expanded = catalog_df_expanded[catalog_df_expanded['Temperature'].notna()]

# Merge results with catalog to get temperatures
merged_df = pd.merge(results_df, catalog_df_expanded[['RunID', 'Temperature']], left_on="Sample", right_on="RunID")

# Determine sample source for coloring
def determine_source(sample):
    if sample.startswith("ERR"):
        return "Tara"
    elif sample.startswith("B"):
        return "NEON"
    elif sample.startswith("SRR"):
        return "GEOTRACES"
    else:
        return "Unknown"

merged_df['Source'] = merged_df['Sample'].apply(determine_source)

# Calculate correlation between free energy and temperature
correlation_free_energy = merged_df["MeanFreeEnergy"].corr(merged_df["Temperature"])
print(f"Correlation between mean free energy and temperature: {correlation_free_energy}")

# Calculate correlation between GC content and temperature
correlation_gc_content = merged_df["MeanGCContent"].corr(merged_df["Temperature"])
print(f"Correlation between mean GC content and temperature: {correlation_gc_content}")

# Save the merged data to a CSV file
merged_df.to_csv(output_csv, index=False)

# Count the number of samples
num_samples = len(merged_df)

# Plot mean free energy vs temperature using Seaborn
plt.figure(figsize=(10, 6))
sns.scatterplot(data=merged_df, x="Temperature", y="MeanFreeEnergy", hue="Source", palette="Set1")
plt.title(f"Mean Free Energy vs Temperature\nCorrelation: {correlation_free_energy:.2f}")
plt.xlabel("Temperature")
plt.ylabel("Mean Free Energy")
plt.legend(title="Source")
plt.annotate(f"Number of samples: {num_samples}", xy=(0.05, 0.95), xycoords='axes fraction', fontsize=10)
plt.savefig(f"{base_path}/mean_free_energy_vs_temperature.png")
plt.close()

# Plot mean GC content vs temperature using Seaborn
plt.figure(figsize=(10, 6))
sns.scatterplot(data=merged_df, x="Temperature", y="MeanGCContent", hue="Source", palette="Set1")
plt.title(f"Mean GC Content vs Temperature\nCorrelation: {correlation_gc_content:.2f}")
plt.xlabel("Temperature")
plt.ylabel("Mean GC Content")
plt.legend(title="Source")
plt.annotate(f"Number of samples: {num_samples}", xy=(0.05, 0.95), xycoords='axes fraction', fontsize=10)
plt.savefig(f"{base_path}/mean_gc_content_vs_temperature.pdf")
plt.close()

print(f"Results and plots have been saved.")

# Find samples in output_dir not represented in scatter plot
samples_with_temperature = set(merged_df['Sample'])
all_samples = set([filename.split('_')[0] for filename in os.listdir(output_dir)])
samples_without_temperature = all_samples - samples_with_temperature

# Output the list of samples without temperature
samples_without_temperature_file = f"{base_path}/samples_without_temperature.txt"
with open(samples_without_temperature_file, 'w') as file:
    for sample in samples_without_temperature:
        file.write(f"{sample}\n")

print(f"Samples without temperature have been saved to {samples_without_temperature_file}.")

#########################################################################################

# ###works great! plotting using plotly. mean (free energy + GC content)

# import os
# import re
# import pandas as pd
# import numpy as np
# import plotly.express as px

# # Directories and files
# base_path = "/home/projects/zeevid/tamirye"
# output_dir = "/home/projects/zeevid/Analyses/2023-Tamir/secondary_structures/Marine_normalized_len"
# catalog_file = "/home/projects/zeevid/tamirye/Data/Combined_datasets_full_modified.csv"
# output_csv = f"{base_path}/100bp_len_secondary_structure.csv"

# # Function to calculate mean free energy and mean GC content from a ViennaRNA output file
# def parse_viennarna_file(file_path):
#     gc_conts = []
#     free_energies = []
#     with open(file_path, 'r') as file:
#         for line in file:
#             if re.search(r'gc_cont=([\d\.]+)', line):
#                 gc_cont = float(re.findall(r'gc_cont=([\d\.]+)', line)[0])
#                 gc_conts.append(gc_cont)
            
#             if re.search(r'\((\-?\d+\.\d+)\)', line):
#                 free_energy = float(re.findall(r'\((\-?\d+\.\d+)\)', line)[0])
#                 free_energies.append(free_energy)

#     mean_free_energy = np.mean(free_energies) if free_energies else np.nan
#     mean_gc_content = np.mean(gc_conts) if gc_conts else np.nan
#     return mean_free_energy, mean_gc_content

# # Process each output file and store results
# results = []
# processed_samples = set()

# for filename in os.listdir(output_dir):
#     sample_name = filename.split('_')[0]
#     file_path = os.path.join(output_dir, filename)
#     mean_free_energy, mean_gc_content = parse_viennarna_file(file_path)
#     results.append((filename, mean_free_energy, mean_gc_content))

# # Create a DataFrame from results
# results_df = pd.DataFrame(results, columns=["Sample", "MeanFreeEnergy", "MeanGCContent"])

# # Load the catalog file
# catalog_df = pd.read_csv(catalog_file)

# # Split RunID by '|' and explode the DataFrame to create rows for each individual ID
# catalog_df['RunID'] = catalog_df['RunID'].astype(str)  # Ensure RunID is a string
# catalog_df_expanded = catalog_df.assign(RunID=catalog_df['RunID'].str.split('|')).explode('RunID')

# # Filter catalog to remove samples with empty temperature values
# catalog_df_expanded = catalog_df_expanded[catalog_df_expanded['Temperature'].notna()]

# # Merge results with catalog to get temperatures
# merged_df = pd.merge(results_df, catalog_df_expanded[['RunID', 'Temperature']], left_on="Sample", right_on="RunID")

# # Calculate correlation between free energy and temperature
# correlation_free_energy = merged_df["MeanFreeEnergy"].corr(merged_df["Temperature"])
# print(f"Correlation between mean free energy and temperature: {correlation_free_energy}")

# # Calculate correlation between GC content and temperature
# correlation_gc_content = merged_df["MeanGCContent"].corr(merged_df["Temperature"])
# print(f"Correlation between mean GC content and temperature: {correlation_gc_content}")

# # Save the merged data to a CSV file
# merged_df.to_csv(output_csv, index=False)

# # Count the number of samples
# num_samples = len(merged_df)

# # Plot mean free energy vs temperature using Plotly
# fig_free_energy = px.scatter(
#     merged_df, x="Temperature", y="MeanFreeEnergy",
#     title=f"Mean Free Energy vs Temperature\nCorrelation: {correlation_free_energy:.2f}",
#     labels={"Temperature": "Temperature", "MeanFreeEnergy": "Mean Free Energy"},
#     hover_data={"Sample": True}
# )
# fig_free_energy.update_traces(marker=dict(size=8))
# fig_free_energy.add_annotation(
#     text=f"Number of samples: {num_samples}",
#     xref="paper", yref="paper",
#     x=0.05, y=0.95, showarrow=False, 
#     font=dict(size=10)
# )
# fig_free_energy.write_html(f"{base_path}/mean_Marine_100bp_len_free_energy_vs_temperature.html")

# # Plot mean GC content vs temperature using Plotly
# fig_gc_content = px.scatter(
#     merged_df, x="Temperature", y="MeanGCContent",
#     title=f"Mean GC Content vs Temperature\nCorrelation: {correlation_gc_content:.2f}",
#     labels={"Temperature": "Temperature", "MeanGCContent": "Mean GC Content"},
#     hover_data={"Sample": True}
# )
# fig_gc_content.update_traces(marker=dict(size=10))
# fig_gc_content.add_annotation(
#     text=f"Number of samples: {num_samples}",
#     xref="paper", yref="paper",
#     x=0.05, y=0.95, showarrow=False,
#     font=dict(size=18)
# )
# fig_gc_content.write_html(f"{base_path}/mean_gc_content_vs_temperature.html")

# print(f"Results and plots have been saved.")

# # Find samples in output_dir not represented in scatter plot
# samples_with_temperature = set(merged_df['Sample'])
# all_samples = set([filename.split('_')[0] for filename in os.listdir(output_dir)])
# samples_without_temperature = all_samples - samples_with_temperature

# # Output the list of samples without temperature
# samples_without_temperature_file = f"{base_path}/samples_without_temperature.txt"
# with open(samples_without_temperature_file, 'w') as file:
#     for sample in samples_without_temperature:
#         file.write(f"{sample}\n")

# print(f"Samples without temperature have been saved to {samples_without_temperature_file}.")


###works great! plotting using plotly. median (free energy + GC content)

# import os
# import re
# import pandas as pd
# import numpy as np
# import plotly.express as px

# # Directories and files
# base_path = "/home/projects/zeevid/tamirye"
# output_dir = "/home/projects/zeevid/Analyses/2023-Tamir/secondary_structures/100k"
# catalog_file = "/home/projects/zeevid/tamirye/Data/Combined_datasets_full_modified.csv"
# output_csv = f"{base_path}/merged_results.csv"

# # Function to calculate median free energy and median GC content from a ViennaRNA output file
# def parse_viennarna_file(file_path):
#     gc_conts = []
#     free_energies = []
#     with open(file_path, 'r') as file:
#         for line in file:
#             if re.search(r'gc_cont=([\d\.]+)', line):
#                 gc_cont = float(re.findall(r'gc_cont=([\d\.]+)', line)[0])
#                 gc_conts.append(gc_cont)
            
#             if re.search(r'\((\-?\d+\.\d+)\)', line):
#                 free_energy = float(re.findall(r'\((\-?\d+\.\d+)\)', line)[0])
#                 free_energies.append(free_energy)

#     median_free_energy = np.median(free_energies) if free_energies else np.nan
#     median_gc_content = np.median(gc_conts) if gc_conts else np.nan
#     return median_free_energy, median_gc_content

# # Process each output file and store results
# results = []
# processed_samples = set()

# for filename in os.listdir(output_dir):
#     sample_name = filename.split('_')[0]
#     file_path = os.path.join(output_dir, filename)
#     median_free_energy, median_gc_content = parse_viennarna_file(file_path)
#     results.append((filename, median_free_energy, median_gc_content))

# # Create a DataFrame from results
# results_df = pd.DataFrame(results, columns=["Sample", "MedianFreeEnergy", "MedianGCContent"])

# # Load the catalog file
# catalog_df = pd.read_csv(catalog_file)

# # Split RunID by '|' and explode the DataFrame to create rows for each individual ID
# catalog_df['RunID'] = catalog_df['RunID'].astype(str)  # Ensure RunID is a string
# catalog_df_expanded = catalog_df.assign(RunID=catalog_df['RunID'].str.split('|')).explode('RunID')

# # Filter catalog to remove samples with empty temperature values
# catalog_df_expanded = catalog_df_expanded[catalog_df_expanded['Temperature'].notna()]

# # Merge results with catalog to get temperatures
# merged_df = pd.merge(results_df, catalog_df_expanded[['RunID', 'Temperature']], left_on="Sample", right_on="RunID")

# # Calculate correlation between free energy and temperature
# correlation_free_energy = merged_df["MedianFreeEnergy"].corr(merged_df["Temperature"])
# print(f"Correlation between median free energy and temperature: {correlation_free_energy}")

# # Calculate correlation between GC content and temperature
# correlation_gc_content = merged_df["MedianGCContent"].corr(merged_df["Temperature"])
# print(f"Correlation between median GC content and temperature: {correlation_gc_content}")

# # Save the merged data to a CSV file
# merged_df.to_csv(output_csv, index=False)

# # Count the number of samples
# num_samples = len(merged_df)

# # Plot median free energy vs temperature using Plotly
# fig_free_energy = px.scatter(
#     merged_df, x="Temperature", y="MedianFreeEnergy",
#     title=f"Median Free Energy vs Temperature\nCorrelation: {correlation_free_energy:.2f}",
#     labels={"Temperature": "Temperature", "MedianFreeEnergy": "Median Free Energy"},
#     hover_data={"Sample": True}
# )
# fig_free_energy.update_traces(marker=dict(size=10))
# fig_free_energy.add_annotation(
#     text=f"Number of samples: {num_samples}",
#     xref="paper", yref="paper",
#     x=0.05, y=0.95, showarrow=False, 
#     font=dict(size=10)
# )
# fig_free_energy.write_html(f"{base_path}/median_100k_free_energy_vs_temperature.html")

# # Plot median GC content vs temperature using Plotly
# fig_gc_content = px.scatter(
#     merged_df, x="Temperature", y="MedianGCContent",
#     title=f"Median GC Content vs Temperature\nCorrelation: {correlation_gc_content:.2f}",
#     labels={"Temperature": "Temperature", "MedianGCContent": "Median GC Content"},
#     hover_data={"Sample": True}
# )
# fig_gc_content.update_traces(marker=dict(size=10))
# fig_gc_content.add_annotation(
#     text=f"Number of samples: {num_samples}",
#     xref="paper", yref="paper",
#     x=0.05, y=0.95, showarrow=False,
#     font=dict(size=18)
# )
# fig_gc_content.write_html(f"{base_path}/median_gc_content_vs_temperature.html")

# print(f"Results and plots have been saved.")

# # Find samples in output_dir not represented in scatter plot
# samples_with_temperature = set(merged_df['Sample'])
# all_samples = set([filename.split('_')[0] for filename in os.listdir(output_dir)])
# samples_without_temperature = all_samples - samples_with_temperature

# # Output the list of samples without temperature
# samples_without_temperature_file = f"{base_path}/samples_without_temperature.txt"
# with open(samples_without_temperature_file, 'w') as file:
#     for sample in samples_without_temperature:
#         file.write(f"{sample}\n")

# print(f"Samples without temperature have been saved to {samples_without_temperature_file}.")


# # # Plot mean GC content vs temperature using Plotly
# # fig_gc_content = px.scatter(
# #     merged_df, x="Temperature", y="MedianGCContent",
# #     title=f"Median GC Content vs Temperature\nCorrelation: {correlation_gc_content:.2f}",
# #     labels={"Temperature": "Temperature", "MedianGCContent": "Median GC Content"},
# #     hover_data={"Sample": True}
# # )
# # fig_gc_content.update_traces(marker=dict(size=10))
# # fig_gc_content.add_annotation(
# #     text=f"Number of samples: {num_samples}",
# #     xref="paper", yref="paper",
# #     x=0.05, y=0.95, showarrow=False,
# #     font=dict(size=18)
# # )
# # fig_gc_content.write_html(f"{base_path}/median_gc_content_vs_temperature.html")

# print(f"Results and plots have been saved.")

# # Find samples in output_dir not represented in scatter plot
# samples_with_temperature = set(merged_df['Sample'])
# all_samples = set([filename.split('_')[0] for filename in os.listdir(output_dir)])
# samples_without_temperature = all_samples - samples_with_temperature

# # Output the list of samples without temperature
# samples_without_temperature_file = f"{base_path}/samples_without_temperature.txt"
# with open(samples_without_temperature_file, 'w') as file:
#     for sample in samples_without_temperature:
#         file.write(f"{sample}\n")

# print(f"Samples without temperature have been saved to {samples_without_temperature_file}.")



####################################### 
###Median calculation, no CG content, plotly

# import os
# import re
# import pandas as pd
# import numpy as np
# import plotly.express as px

# # Directories and files
# base_path = "/home/projects/zeevid/tamirye"
# output_dir = "/home/projects/zeevid/Analyses/2023-Tamir/secondary_structures/100k_Marine_normalized_len"
# catalog_file = "/home/projects/zeevid/tamirye/Data/Combined_datasets_full_modified.csv"
# output_csv = f"{base_path}/merged_results.csv"

# # Function to calculate median free energy from a ViennaRNA output file
# def parse_viennarna_file(file_path):
#     free_energies = []
#     with open(file_path, 'r') as file:
#         for line in file:
#             if re.search(r'\((\-?\d+\.\d+)\)', line):
#                 free_energy = float(re.findall(r'\((\-?\d+\.\d+)\)', line)[0])
#                 free_energies.append(free_energy)

#     median_free_energy = np.median(free_energies) if free_energies else np.nan
#     return median_free_energy

# # Process each output file and store results
# results = []
# processed_samples = set()

# for filename in os.listdir(output_dir):
#     sample_name = filename.split('_')[0]
#     file_path = os.path.join(output_dir, filename)
#     median_free_energy = parse_viennarna_file(file_path)
#     results.append((filename, median_free_energy))

# # Create a DataFrame from results
# results_df = pd.DataFrame(results, columns=["Sample", "MedianFreeEnergy"])

# # Load the catalog file
# catalog_df = pd.read_csv(catalog_file)

# # Split RunID by '|' and explode the DataFrame to create rows for each individual ID
# catalog_df['RunID'] = catalog_df['RunID'].astype(str)  # Ensure RunID is a string
# catalog_df_expanded = catalog_df.assign(RunID=catalog_df['RunID'].str.split('|')).explode('RunID')

# # Filter catalog to remove samples with empty temperature values
# catalog_df_expanded = catalog_df_expanded[catalog_df_expanded['Temperature'].notna()]

# # Merge results with catalog to get temperatures
# merged_df = pd.merge(results_df, catalog_df_expanded[['RunID', 'Temperature']], left_on="Sample", right_on="RunID")

# # Calculate correlation between free energy and temperature
# correlation_free_energy = merged_df["MedianFreeEnergy"].corr(merged_df["Temperature"])
# print(f"Correlation between median free energy and temperature: {correlation_free_energy}")

# # Save the merged data to a CSV file
# merged_df.to_csv(output_csv, index=False)

# # Count the number of samples
# num_samples = len(merged_df)

# # Plot median free energy vs temperature using Plotly
# fig_free_energy = px.scatter(
#     merged_df, x="Temperature", y="MedianFreeEnergy",
#     title=f"Marine, 10k reads: Median Free Energy vs Temperature\nCorrelation: {correlation_free_energy:.2f}",
#     labels={"Temperature": "Temperature", "MedianFreeEnergy": "Median Free Energy"},
#     hover_data={"Sample": True}
# )
# fig_free_energy.update_traces(marker=dict(size=8))
# fig_free_energy.add_annotation(
#     text=f"Number of samples: {num_samples}",
#     xref="paper", yref="paper",
#     x=0.05, y=0.95, showarrow=False, 
#     font=dict(size=10)
# )
# fig_free_energy.write_html(f"{base_path}/100k_median_Marine_100bp_len_free_energy_vs_temperature.html")

# print(f"Results and plots have been saved.")

# # Find samples in output_dir not represented in scatter plot
# samples_with_temperature = set(merged_df['Sample'])
# all_samples = set([filename.split('_')[0] for filename in os.listdir(output_dir)])
# samples_without_temperature = all_samples - samples_with_temperature

# # Output the list of samples without temperature
# samples_without_temperature_file = f"{base_path}/samples_without_temperature.txt"
# with open(samples_without_temperature_file, 'w') as file:
#     for sample in samples_without_temperature:
#         file.write(f"{sample}\n")

# print(f"Samples without temperature have been saved to {samples_without_temperature_file}.")

