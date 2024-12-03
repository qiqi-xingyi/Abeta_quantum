# --*-- conding:utf-8 --*--
# @Time : 11/19/24 2:07 PM
# @Author : Yuqi Zhang
# @Email : yzhan135@kent.edu
# @File : plt_line.py

import matplotlib.pyplot as plt

if __name__ == '__main__':

    # Load energy data from the file
    file_path = './4zb8_B/System_Enegry/energy_list.txt'  # Replace with the correct path if necessary
    with open(file_path, 'r') as file:
        energy_data = [float(line.strip()) for line in file]

    # Plot the energy data
    plt.figure(figsize=(10, 6))
    plt.plot(range(1, len(energy_data) + 1), energy_data, marker='o', linestyle='-', color='blue', label='Energy')

    # Add labels and title
    plt.title('Energy Variation During Training', fontsize=16)
    plt.xlabel('Training Step', fontsize=14)
    plt.ylabel('Energy', fontsize=14)

    # Add grid and legend
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.legend(fontsize=12)

    # Save the plot
    plt.savefig('./img2/energy_variation_2.png', dpi=300, bbox_inches='tight')

    # Show the plot
    plt.tight_layout()
    plt.show()
