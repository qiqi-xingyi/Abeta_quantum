# --*-- conding:utf-8 --*--
# @Time : 11/19/24 2:09 PM
# @Author : Yuqi Zhang
# @Email : yzhan135@kent.edu
# @File : plt_state.py

import matplotlib.pyplot as plt
import numpy as np

# Load data from the file
file_path = './6mu3_L/Prob_distribution/prob_distribution.txt'
x_labels = []
y_labels = []
z_values = []

if __name__ == '__main__':


    with open(file_path, 'r') as file:
        for line in file:
            state, count = line.strip().split(':')
            state = state.strip()
            count = int(count.strip())
            x_labels.append(state[:4])  # First 4 qubits as x-axis labels
            y_labels.append(state[4:])  # Last 3 qubits as y-axis labels
            z_values.append(count)

    # Convert binary labels to integers for plotting
    x_indices = [int(x, 2) for x in x_labels]
    y_indices = [int(y, 2) for y in y_labels]

    # Get unique x and y indices for the grid
    x_unique = sorted(set(x_indices))
    y_unique = sorted(set(y_indices))

    # Create a 2D grid and populate it with z-values
    z_grid = np.zeros((len(y_unique), len(x_unique)))

    for x, y, z in zip(x_indices, y_indices, z_values):
        x_idx = x_unique.index(x)
        y_idx = y_unique.index(y)
        z_grid[y_idx, x_idx] = z

    # Plot the 3D heatmap
    fig = plt.figure(figsize=(12, 8))
    ax = fig.add_subplot(111, projection='3d')

    # Create a meshgrid for x, y coordinates
    x, y = np.meshgrid(x_unique, y_unique)

    # Flatten the grid for bar3d
    x = x.flatten()
    y = y.flatten()
    z = np.zeros_like(x)
    dx = dy = 0.8
    dz = z_grid.flatten()

    # Use a colormap to indicate frequency
    cmap = plt.cm.viridis
    colors = cmap(dz / dz.max())  # Normalize to [0, 1]

    # Plot the 3D bar chart
    ax.bar3d(x, y, z, dx, dy, dz, color=colors, shade=True)

    # Customize the plot
    ax.set_title('Measurement Distribution by Qubit States', fontsize=14)
    ax.set_xlabel('First 4 Qubits', fontsize=13, labelpad=20)
    ax.set_ylabel('Last 3 Qubits', fontsize=13,labelpad=20)
    ax.set_zlabel('Frequency', fontsize=13,labelpad=20)

    # Adjust tick labels
    ax.set_xticks(x_unique)
    ax.set_xticklabels([bin(x)[2:].zfill(4) for x in x_unique], rotation=45)
    ax.set_yticks(y_unique)
    ax.set_yticklabels([bin(y)[2:].zfill(3) for y in y_unique])

    plt.savefig('./img/qubit_distribution.png', dpi=300, bbox_inches='tight')
    plt.tight_layout()
    plt.show()
