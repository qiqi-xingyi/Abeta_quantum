# --*-- conding:utf-8 --*--
# @Time : 11/25/24 1:22 PM
# @Author : Yuqi Zhang
# @Email : yzhan135@kent.edu
# @File : Abeta_Quantum.py

import os
from Protein_Folding import Peptide
from Protein_Folding.interactions.miyazawa_jernigan_interaction import MiyazawaJerniganInteraction
from Protein_Folding.penalty_parameters import PenaltyParameters
from Protein_Folding.protein_folding_problem import ProteinFoldingProblem
from qiskit_ibm_runtime import QiskitRuntimeService
from Qiskit_VQE import VQE, StateCalculator

# Configurable parameters: window size and starting window index
window_size = 12
# Option 1: Specify the 1-based starting window index directly
start_window = 7

# Option 2: Calculate start_window from a residue position
# Uncomment the following lines and comment out Option 1 if needed:
# start_pos = 16
# start_window = max(start_pos - window_size + 1, 1)

main_chain_residue_seq = "DAEFRHDSGYEVHHQKLVFFAEDVGSNKGAIIGLMVGGVVIA"
protein_name = 'Abeta_A'

def main():
    # Initialize IBM Quantum Runtime service
    service = QiskitRuntimeService(
        channel='ibm_quantum',
        instance='ibm-q-ccf/qradle-catalyzer/qradle-catalyzer',
        token='YOUR_TOKEN_HERE'
    )

    # Create output directory
    output_dir = os.path.join("Abeta_42", protein_name)
    os.makedirs(output_dir, exist_ok=True)

    # Set up interaction and penalty parameters
    mj_interaction = MiyazawaJerniganInteraction()
    penalty_terms = PenaltyParameters(10, 10, 10)

    total_windows = len(main_chain_residue_seq) - window_size + 1
    print(f"Total windows available: {total_windows}. Starting from window {start_window}.")

    for idx in range(start_window - 1, total_windows):
        window_idx = idx + 1
        window_seq = main_chain_residue_seq[idx:idx + window_size]
        print(f"Processing window {window_idx}: positions {idx+1} to {idx+window_size}, sequence: {window_seq}")

        peptide = Peptide(window_seq, [''] * window_size)
        problem = ProteinFoldingProblem(peptide, mj_interaction, penalty_terms)
        hamiltonian = problem.qubit_op()

        qubit_num = hamiltonian.num_qubits + 2
        print(f"Window {window_idx} requires {qubit_num} qubits")

        vqe_inst = VQE(service=service, hamiltonian=hamiltonian, min_qubit_num=qubit_num, maxiter=30)
        energy_list, res, ansatz = vqe_inst.run_vqe()

        # Save energy convergence history
        energy_file = os.path.join(output_dir, f'energy_list_window_{window_idx}.txt')
        with open(energy_file, 'w') as f:
            for e in energy_list:
                f.write(f"{e}\n")

        # Calculate and save probability distribution
        prob_dist = StateCalculator(service, qubit_num, ansatz).get_probability_distribution(res)
        prob_file = os.path.join(output_dir, f'prob_distribution_window_{window_idx}.txt')
        with open(prob_file, 'w') as f:
            for state, p in prob_dist.items():
                f.write(f"{state}: {p}\n")
        print(f"VQE result for window {window_idx}: {prob_dist}")

        # Interpret and save structure as XYZ
        result = problem.interpret(prob_dist)
        result.save_xyz_file(name=f'{protein_name}_window_{window_idx}', path=output_dir)
        print(f"Saved structure .xyz for window {window_idx}")

if __name__ == '__main__':
    main()

