# --*-- conding:utf-8 --*--
# @Time : 1/3/25 9:10 AM
# @Author : Yuqi Zhang
# @Email : yzhan135@kent.edu
# @File : Supplementary_structure_case1.py

import os
from Protein_Folding import Peptide
from Protein_Folding.interactions.miyazawa_jernigan_interaction import MiyazawaJerniganInteraction
from Protein_Folding.penalty_parameters import PenaltyParameters
from Protein_Folding.protein_folding_problem import ProteinFoldingProblem
from qiskit_ibm_runtime import QiskitRuntimeService
from Qiskit_VQE import VQE5
from Qiskit_VQE import StateCalculator

main_chain_residue_seq = "YFASGQPYRYER"
side_chain_residue_sequences = ['' for _ in range(len(main_chain_residue_seq))]
protein_name = '4zb8_B'

if __name__ == '__main__':

    ########## Create a Quantum Service ############

    service = QiskitRuntimeService(
        channel='ibm_quantum',
        instance='ibm-q-ccf/qradle-catalyzer/qradle-catalyzer',
        token='token'
    )

    char_count = len(main_chain_residue_seq)
    print(f'Num of Acid:{char_count}')

    side_site = len(side_chain_residue_sequences)
    print(side_chain_residue_sequences)
    print(f'Num of Side cite:{side_site}')

    # create Peptide
    peptide = Peptide(main_chain_residue_seq, side_chain_residue_sequences)

    # Interaction definition (e.g. Miyazawa-Jernigan)
    mj_interaction = MiyazawaJerniganInteraction()

    # Penalty Parameters Definition
    penalty_terms = PenaltyParameters(10, 10, 10)

    # Create Protein Folding case
    protein_folding_problem = ProteinFoldingProblem(peptide, mj_interaction, penalty_terms)

    # create quantum Op
    hamiltonian = protein_folding_problem.qubit_op()

    # print('Operator',hamiltonian)


    # ansatz = EfficientSU2(hamiltonian.num_qubits)
    qubits_num = hamiltonian.num_qubits + 2
    print(f'Num of qubits:{qubits_num}')
    vqe_instance = VQE5(service=service, hamiltonian=hamiltonian, min_qubit_num=qubits_num, maxiter=150)

    # Run the VQE algorithm
    energy_list, res, ansatz, top_5_results = vqe_instance.run_vqe()


    with open('./QC_Status_Analysis/System_Enegry/energy_list_4zb8_B.txt', 'w') as file:
        for item in energy_list:
            file.write(str(item) + '\n')

    state_calculator = StateCalculator(service, qubits_num, ansatz)
    prob_distribution = state_calculator.get_probability_distribution(res)

    protein_result = protein_folding_problem.interpret(prob_distribution)

    # save to .xyz
    output_dir = f"Post_Processing/create_structure/{protein_name}"

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    protein_result.save_xyz_file(name=protein_name, path=output_dir)
    print("Protein structure saved as .xyz file")

    for rank, (energy_val, best_params) in enumerate(top_5_results, start=1):
        print(f"Top {rank} best energy = {energy_val}")

        prob_distribution_best = state_calculator.get_probability_distribution(best_params)

        protein_result_best = protein_folding_problem.interpret(prob_distribution_best)

        protein_result_best.save_xyz_file(
            name=f"{protein_name}_top_{rank}",
            path=output_dir
        )

        print(f"Protein structure for top {rank} best result has been saved.")

