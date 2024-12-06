import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import sys
import os

def plot_distributions(csv_files, output_folder):
    # Initialize a figure for multiple plots
    plt.figure(figsize=(15, 9))
    
    # Loop through the CSV files and plot each TF's distance distribution
    for csv_file in csv_files:
        # Extract the TF name from the CSV file name
        tf_name = os.path.basename(csv_file).split('_')[0]

        # Load the distances for the TF
        data = pd.read_csv(csv_file)

        # Plot the distribution for both regions (first and second) using seaborn
        sns.kdeplot(
            data=data[data['RegionSet'] == 'first']['Distance'],
            label=f"{tf_name} (Domains)",
            linestyle='-'
        )
        

    # Customize the plot
    plt.title("Distance Distributions Between TF Binding Sites and BEDPE Regions")
    plt.xlabel("Distance")
    plt.ylabel("Density")
    plt.legend()
    plt.tight_layout()

    # Save the plot
    plt.savefig(os.path.join(output_folder, "distance_distributions.png"))
    plt.show()

if __name__ == "__main__":
    # The first argument is a space-separated list of CSV files
    csv_files = sys.argv[1:-1]
    output_folder = sys.argv[-1]

    # Plot distributions
    plot_distributions(csv_files, output_folder)

