# --*-- conding:utf-8 --*--
# @Time : 12/1/24 4:57 PM
# @Author : Yuqi Zhang
# @Email : yzhan135@kent.edu
# @File : integration_Abeta.py

import os
import numpy as np

def read_xyz_file(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()
        num_atoms = int(lines[0].strip())
        atoms = []
        coords = []
        for line in lines[2:2+num_atoms]:
            parts = line.strip().split()
            atoms.append(parts[0])
            coords.append([float(parts[1]), float(parts[2]), float(parts[3])])
        return atoms, np.array(coords)

if __name__ == '__main__':
    protein_name = 'Abeta_A'
    output_dir = f"Post_Processing/process_data/{protein_name}"
    window_size = 7  # Sliding window size

    # Get all .xyz files for sliding windows, sorted by window index
    xyz_files = sorted([f for f in os.listdir(output_dir) if f.endswith('.xyz') and 'window' in f],
                       key=lambda x: int(x.split('_window_')[1].split('.xyz')[0]))
    total_windows = len(xyz_files)

    # Read prediction results for all sliding windows
    window_atoms_list = []
    window_coords_list = []
    for xyz_file in xyz_files:
        file_path = os.path.join(output_dir, xyz_file)
        atoms, coords = read_xyz_file(file_path)
        window_atoms_list.append(atoms)
        window_coords_list.append(coords)

    # Total length of the protein
    protein_length = len(window_atoms_list[0]) + total_windows - 1

    # Initialize the list of direction vector sets, length is protein length minus one
    direction_vectors_list = [[] for _ in range(protein_length - 1)]

    # Traverse each sliding window and extract direction vectors
    for idx, (atoms, coords) in enumerate(zip(window_atoms_list, window_coords_list)):
        # Starting position of the current sliding window in the global sequence
        start_pos = idx
        for i in range(len(atoms) - 1):
            global_pos = start_pos + i
            # Calculate the direction vector
            direction = coords[i+1] - coords[i]
            # Add it to the corresponding position's direction vector set
            direction_vectors_list[global_pos].append(direction)

    # Calculate the average direction vector for each position
    average_directions = []
    for idx, vectors in enumerate(direction_vectors_list):
        if vectors:
            # Compute the average direction vector
            avg_direction = np.mean(vectors, axis=0)
            average_directions.append(avg_direction)
        else:
            # Use a zero vector if no direction vectors are available
            average_directions.append(np.zeros(3))

    # Reconstruct global coordinates
    global_atoms = []
    global_coords = []

    # Assume the first residue is at the origin
    current_coord = np.array([0.0, 0.0, 0.0])
    global_coords.append(current_coord)
    global_atoms.append(window_atoms_list[0][0])  # Name of the first residue

    for i in range(len(average_directions)):
        current_coord = current_coord + average_directions[i]
        global_coords.append(current_coord)
        # Get the name of the (i+1)-th residue
        residue_name = window_atoms_list[0][i+1] if i+1 < len(window_atoms_list[0]) else window_atoms_list[-1][-1]
        global_atoms.append(residue_name)

    # Convert the coordinate list to an array
    global_coords = np.array(global_coords)

    # Save the complete protein structure
    final_xyz_file = os.path.join(output_dir, f'{protein_name}_final.xyz')
    num_atoms = len(global_atoms)
    with open(final_xyz_file, 'w') as file:
        file.write(f"{num_atoms}\n\n")
        for atom, coord in zip(global_atoms, global_coords):
            file.write(f"{atom} {coord[0]} {coord[1]} {coord[2]}\n")
    print(f"The complete protein structure has been saved as {final_xyz_file}")


