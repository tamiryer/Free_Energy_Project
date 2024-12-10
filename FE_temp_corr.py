import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import pearsonr

# File paths
rank_summary_path = "/home/projects/zeevid/Analyses/2023-Tamir/secondary_structures/plots/Marine/median/same_gc_10_mod/summary_marine.csv"
catalog_path = "/home/projects/zeevid/tamirye/Data/United_all_Ocean-with_fasta_no_dup.csv"
output_plot_path = "/home/projects/zeevid/Analyses/2023-Tamir/secondary_structures/plots/Marine/mean/same_gc_10_mod/correlation_plot_mean_rank.pdf"  # Path to save the plot

# Read the rank summary data
rank_summary_df = pd.read_csv(rank_summary_path)
rank_summary_df["File Name"] = rank_summary_df["File Name"].str.replace('.txt', '', regex=False)

# Read the catalog data (assuming it has columns "RunID" and "Temperature_degC")
catalog_df = pd.read_csv(catalog_path)

# Merge the rank summary with the catalog data
merged_df = pd.merge(catalog_df, rank_summary_df, left_on="RunID", right_on="File Name")

# Remove rows with missing Temperature_degC values
merged_df = merged_df.dropna(subset=["Temperature_degC"])

# Calculate the correlation between Temperature_degC and Median Rank
correlation, p_value = pearsonr(merged_df["Temperature_degC"], merged_df["Mean Rank"])
print(f"Pearson correlation: {correlation:.2f}, p-value: {p_value:.2e}")

# Plotting the correlation with a regression line
plt.figure(figsize=(10, 6))
sns.regplot(x="Temperature_degC", y="Mean Rank", data=merged_df, ci=None, scatter_kws={"s": 50}, line_kws={"color": "red"})
plt.xlabel("Temperature (°C)")
plt.ylabel("Mean Rank")
plt.title(f"Marine - Mean Rank, Correlation: {correlation:.2f}") #, p-value: {p_value:.2e}")
plt.grid(True)

# Save the plot as a PDF file
plt.savefig(output_plot_path)
#plt.show()
