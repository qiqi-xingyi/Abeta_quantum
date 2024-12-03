# --*-- conding:utf-8 --*--
# @Time : 10/23/24 PM3:16
# @Author : Yuqi Zhang
# @Email : yzhan135@kent.edu
# @File : Get_result_from_job.py

from qiskit_ibm_runtime import QiskitRuntimeService
from Protein_Folding import Peptide, ProteinFoldingProblem
from Protein_Folding.interactions.miyazawa_jernigan_interaction import MiyazawaJerniganInteraction
from Protein_Folding.penalty_parameters import PenaltyParameters
from qiskit import result

main_chain_residue_seq = "ELFDRIEPDIGMPDER"
side_chain_residue_sequences = ['', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '']

# 创建 Peptide 对象
peptide = Peptide(main_chain_residue_seq, side_chain_residue_sequences)

# 创建相互作用和惩罚参数
mj_interaction = MiyazawaJerniganInteraction()
penalty_terms = PenaltyParameters(15, 15, 15)

# 创建 ProteinFoldingProblem 实例
protein_folding_problem = ProteinFoldingProblem(peptide, mj_interaction, penalty_terms)

if __name__ == '__main__':
    # 初始化 IBM Quantum 服务
    service = QiskitRuntimeService(
        channel='ibm_quantum',
        instance='ibm-q-ccf/qradle-catalyzer/qradle-catalyzer',
        token='98da9815dd1fbbe8d3010882e9a317f9495f2d61652ec33f19429c2136da25975a0728843211b0b389d731778c600c27e30b5edfeee39c318793a925668dbfae'
    )

    job = service.job('cws3cc9ehebg008j1zbg')
    job_result = job.result()

    for idx, pub_result in enumerate(job_result):
        print(f"Expectation values for pub {idx}: {pub_result.data.evs}")






