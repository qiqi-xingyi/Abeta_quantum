# --*-- conding:utf-8 --*--
# @Time : 11/18/24 9:43 AM
# @Author : Yuqi Zhang
# @Email : yzhan135@kent.edu
# @File : plt_box.py

import matplotlib.pyplot as plt

# Function to read data from the text file and extract values
def read_data_from_txt(file_path):
    quantum_data = {"affinity": [], "rmsd_lower": [], "rmsd_upper": []}
    af2_data = {"affinity": [], "rmsd_lower": [], "rmsd_upper": []}
    is_quantum = True  # Flag to differentiate between Quantum and AF2 results

    with open(file_path, 'r') as file:
        for line in file:
            if line.strip().startswith("Quantum Results:"):
                is_quantum = True
            elif line.strip().startswith("AF2 Results:"):
                is_quantum = False
            elif "Affinity" in line:
                # Extract Affinity, RMSD Lower Bound, and RMSD Upper Bound values
                affinity = float(line.split("Affinity =")[1].split(",")[0].strip())
                rmsd_lower = float(line.split("RMSD Lower Bound =")[1].split(",")[0].strip())
                rmsd_upper = float(line.split("RMSD Upper Bound =")[1].strip())

                # Append to the appropriate dictionary
                if is_quantum:
                    quantum_data["affinity"].append(affinity)
                    quantum_data["rmsd_lower"].append(rmsd_lower)
                    quantum_data["rmsd_upper"].append(rmsd_upper)
                else:
                    af2_data["affinity"].append(affinity)
                    af2_data["rmsd_lower"].append(rmsd_lower)
                    af2_data["rmsd_upper"].append(rmsd_upper)

    return quantum_data, af2_data

if __name__ == '__main__':

    # Path to the input text file
    file_path = 'docking_output_2/summary_results.txt'  # Update this with your actual file path

    # Read data from the text file
    quantum_data, af2_data = read_data_from_txt(file_path)

    # List of data types and corresponding colors
    data_types = ["affinity", "rmsd_lower", "rmsd_upper"]
    titles = ["Affinity Scores", "RMSD Lower Bound", "RMSD Upper Bound"]
    colors = [
        dict(facecolor='lightblue', edgecolor='darkblue'),
        dict(facecolor='lightgreen', edgecolor='darkgreen'),
        dict(facecolor='lightcoral', edgecolor='darkred'),
    ]

    for i, data_type in enumerate(data_types):
        # Create boxplots for each data type
        fig, ax = plt.subplots(figsize=(6, 7))
        ax.boxplot(
            [quantum_data[data_type], af2_data[data_type]],
            labels=['Quantum Results', 'AF3 Results'],
            patch_artist=True,
            boxprops=colors[i],  # Apply custom colors for boxes
            medianprops=dict(color='red', linewidth=2),  # Median style
            whiskerprops=dict(color='black', linewidth=1, linestyle='-'),  # Whisker style
            capprops=dict(color='black', linewidth=1.5),  # Cap style
            flierprops=dict(marker='o', color='orange', alpha=0.5),  # Flier style
            widths=0.3
        )

        # Customization for better visualization
        ax.set_title(f'Comparison of {titles[i]}')
        ax.set_ylabel(f'{titles[i]}')
        ax.grid(axis='y', linestyle='--', alpha=0.7)

        # Save each plot
        plt.savefig(f'./img2/{data_type}_box.png', dpi=300, bbox_inches='tight')

        # Show the plot
        plt.tight_layout()
        plt.show()

