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
from Qiskit_VQE import VQE
from Qiskit_VQE import StateCalculator

main_chain_residue_seq = "DAEFRHDSGYEVHHQKLVFFAEDVGSNKGAIIGLMVGGVVIA" #Abeta-42
side_chain_residue_sequences = ['' for _ in range(len(main_chain_residue_seq))]
protein_name = 'Abeta_A'

window_size = 7

if __name__ == '__main__':

    service = QiskitRuntimeService(
        channel='ibm_quantum',
        instance='ibm-q-ccf/qradle-catalyzer/qradle-catalyzer',
        token='token'
    )

    output_dir = f"Post_Processing/process_data/{protein_name}"
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    mj_interaction = MiyazawaJerniganInteraction()

    penalty_terms = PenaltyParameters(10, 10, 10)


    for i in range(len(main_chain_residue_seq) - window_size + 1):

        window_main_chain = main_chain_residue_seq[i:i+window_size]
        window_side_chain = ['' for _ in range(len(window_main_chain))]

        print(f'Processing sliding window {i+1}, positions {i+1} to {i+window_size}, sequence: {window_main_chain}')

        peptide = Peptide(window_main_chain, window_side_chain)

        protein_folding_problem = ProteinFoldingProblem(peptide, mj_interaction, penalty_terms)

        hamiltonian = protein_folding_problem.qubit_op()

        qubits_num = hamiltonian.num_qubits + 2
        print(f'Sliding window {i+1} requires {qubits_num} qubits')

        vqe_instance = VQE(service=service, hamiltonian=hamiltonian, min_qubit_num=qubits_num, maxiter=30)

        energy_list, res, ansatz = vqe_instance.run_vqe()

        energy_file = os.path.join(output_dir, f'energy_list_window_{i+1}.txt')
        with open(energy_file, 'w') as file:
            for item in energy_list:
                file.write(str(item) + '\n')

        state_calculator = StateCalculator(service, qubits_num, ansatz)
        prob_distribution = state_calculator.get_probability_distribution(res)

        print(f'VQE result for sliding window {i+1}: {prob_distribution}')

        prob_file = os.path.join(output_dir, f'prob_distribution_window_{i+1}.txt')
        with open(prob_file, 'w') as file:
            for key, value in prob_distribution.items():
                file.write(f'{key}: {value}\n')

        protein_result = protein_folding_problem.interpret(prob_distribution)

        protein_result.save_xyz_file(name=f'{protein_name}_window_{i+1}', path=output_dir)
        print(f"Protein structure for sliding window {i+1} has been saved as .xyz file")
