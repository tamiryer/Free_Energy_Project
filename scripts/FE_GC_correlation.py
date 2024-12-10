# import pandas as pd
# import matplotlib.pyplot as plt

# # Load the CSV data into a pandas DataFrame
# csv_file = '/home/projects/zeevid/Analyses/2023-Tamir/secondary_structures/plots/Marine/same_gc_10_mod/rank_summary_Marine.csv' 
# output_plot_path = '/home/projects/zeevid/Analyses/2023-Tamir/secondary_structures/plots/Marine/same_gc_10_mod/FE_GC_correlation.pdf'

# data = pd.read_csv(csv_file)

# # Create a scatter plot
# plt.figure(figsize=(10, 6))
# plt.scatter(data['Median GC Rank'], data['Median FE Rank'], color='blue', edgecolor='black', s=10, alpha=0.6)

# # Add labels and title
# plt.xlabel('Median GC Rank')
# plt.ylabel('Median FE Rank')
# plt.title('Scatter Plot of GC Content vs Free Energy')

# # Display the plot
# plt.grid(True)
# plt.savefig(output_plot_path)
# #plt.show()


# import pandas as pd
# import numpy as np
# import matplotlib.pyplot as plt
# import seaborn as sns
# from scipy.stats import gaussian_kde

# # Load the CSV data into a pandas DataFrame
# csv_file = '/home/projects/zeevid/Analyses/2023-Tamir/secondary_structures/plots/Marine/mean/diff_gc_10_mod/10_mod_diff_gc_merged_summary_Marine_mean.csv' 
# output_plot_path = '/home/projects/zeevid/Analyses/2023-Tamir/secondary_structures/plots/Marine/mean/diff_gc_10_mod/FE_GC_correlation_density_merged.pdf'
# data = pd.read_csv(csv_file)

# # Perform Gaussian Kernel Density Estimation (KDE) to get density of points
# x = data['Mean GC Rank']
# y = data['Mean FE Rank']
# xy = np.vstack([x, y])
# z = gaussian_kde(xy)(xy)  # Calculate point density

# # Sort points by density for proper plotting order (densest on top)
# idx = z.argsort()
# x, y, z = x[idx], y[idx], z[idx]

# # Create a scatter plot with density-based transparency
# plt.figure(figsize=(10, 6))
# plt.scatter(x, y, c=z, cmap='viridis', s=50, alpha=0.6, edgecolor='black')

# # Add labels and title
# plt.xlabel('Mean GC Rank')
# plt.ylabel('Mean FE Rank')
# plt.title('Marine - GC Content vs Free Energy per sample')

# # Add colorbar to show density scale
# plt.colorbar(label='Density')

# # Display the plot
# plt.set_axisbelow(True)
# plt.grid(True)
# plt.savefig(output_plot_path)

# #plt.show()

###works great - for thesis 
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats

# Load the CSV data into a pandas DataFrame
csv_file = '/home/projects/zeevid/Analyses/2023-Tamir/secondary_structures/plots/Soil/mean/diff_gc_10_mod_merged/diff_gc_10_mod_summary__mean.csv'
data = pd.read_csv(csv_file)

# Create figure and axis with a specific size
plt.figure(figsize=(10, 8))

# Create the scatter plot
sns.scatterplot(data=data, 
                x='Mean GC Rank', 
                y='Mean FE Rank',
                alpha=0.6)

# Calculate and add regression line
slope, intercept, r_value, p_value, std_err = stats.linregress(data['Mean GC Rank'], 
                                                              data['Mean FE Rank'])
r_squared = r_value**2
x_range = np.linspace(data['Mean GC Rank'].min(), data['Mean GC Rank'].max(), 100)
y_range = slope * x_range + intercept

# Plot regression line
plt.plot(x_range, y_range, color='red', linestyle='--', 
         label=f'y = {slope:.2f}x + {intercept:.2f}\nR² = {r_squared:.2f} n={len(data)}')

# Customize the plot
plt.title('Soil - GC Content vs Free Energy per Sample', 
          fontsize=14, pad=20)
plt.xlabel('Mean GC Rank', fontsize=12)
plt.ylabel('Mean FE Rank', fontsize=12)

# Add grid for better readability
plt.grid(True, linestyle='--', alpha=0.7)

# Add legend
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)

# Tight layout to prevent label cutoff
plt.tight_layout()

# Save the plot
plt.savefig('/home/projects/zeevid/Analyses/2023-Tamir/secondary_structures/plots/Soil/mean/diff_gc_10_mod_merged/FE_GC_correlation_scatter_merged.png', 
            dpi=300, 
            bbox_inches='tight')

# Clear the current figure
plt.close()


# ###plotly version of the scatter plot above

# import pandas as pd
# import numpy as np
# import plotly.express as px
# from scipy import stats

# # Load the CSV data into a pandas DataFrame
# csv_file = '/home/projects/zeevid/Analyses/2023-Tamir/secondary_structures/plots/Soil/mean/diff_gc_10_mod_merged/diff_gc_10_mod_summary_soil_mean.csv'
# data = pd.read_csv(csv_file)
# num_samples = len(data)
# # Calculate linear regression
# slope, intercept, r_value, p_value, std_err = stats.linregress(data['Mean GC Rank'], data['Mean FE Rank'])
# r_squared = r_value**2

# # Create the regression line data
# x_range = np.linspace(data['Mean GC Rank'].min(), data['Mean GC Rank'].max(), 100)
# y_range = slope * x_range + intercept

# # Create a Plotly scatter plot
# fig = px.scatter(data, 
#                  x='Mean GC Rank', 
#                  y='Mean FE Rank', 
#                  hover_name='Sample Name',  # Assuming the column with sample names is 'Sample Name'
#                  labels={'Mean GC Rank': 'Mean GC Rank', 'Mean FE Rank': 'Mean FE Rank'},
#                  title='Soil - GC Content vs Free Energy per sample')

# # Add the regression line to the plot
# fig.add_scatter(x=x_range, y=y_range, mode='lines', name=f'y = {slope:.2f}x + {intercept:.2f}\nR² = {r_squared:.2f} n={num_samples}', line=dict(color='red', dash='dash'))

# # Show the interactive plot
# #fig.show()

# # Optionally save the plot as an HTML file
# output_plot_path = '/home/projects/zeevid/Analyses/2023-Tamir/secondary_structures/plots/Soil/mean/diff_gc_10_mod_merged/FE_GC_correlation_scatter_merged.html'
# fig.write_html(output_plot_path)
