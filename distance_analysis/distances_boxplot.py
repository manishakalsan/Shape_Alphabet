import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import sys
import os

def plot_boxplots(csv_files, tf_preferences_file, output_folder):
    # Load the TF preferences file
    tf_preferences = pd.read_csv(tf_preferences_file)

    # Initialize a list to collect all TF data
    all_data = []

    # Loop through the CSV files and collect TF data
    for csv_file in csv_files:
        # Extract the TF name from the CSV file name
        tf_name = os.path.basename(csv_file).split('_')[0]

        # Load the distances for the TF
        data = pd.read_csv(csv_file)

        # Add TF name to the data
        data['TF'] = tf_name

        # Append to the all_data list
        all_data.append(data)

    # Concatenate all data into a single DataFrame
    all_distances = pd.concat(all_data, ignore_index=True)

    # Merge with the TF preferences to add class information
    merged_data = all_distances.merge(tf_preferences[['TF', 'Class']], on='TF', how='left')

    # Create the boxplot, grouping by TF class
    plt.figure(figsize=(12, 8))
    sns.boxplot(x='Class', y='Distance', data=merged_data, palette='Set2')

    # Customize the plot
    plt.title("Boxplot of TF Distance Distributions by Class")
    plt.xlabel("TF Class")
    plt.ylabel("Distance")
    plt.xticks(rotation=45)
    plt.tight_layout()

    # Save the plot
    plt.savefig(os.path.join(output_folder, "tf_class_boxplot.png"))
    plt.show()

if __name__ == "__main__":
    # The first argument is a space-separated list of CSV files
    csv_files = sys.argv[1:-2]
    tf_preferences_file = sys.argv[-2]
    output_folder = sys.argv[-1]

    # Plot boxplots
    plot_boxplots(csv_files, tf_preferences_file, output_folder)

