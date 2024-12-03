# --*-- conding:utf-8 --*--
# @Time : 12/1/24 4:14 PM
# @Author : Yuqi Zhang
# @Email : yzhan135@kent.edu
# @File : Read_PDB.py

from Bio.PDB import PDBParser


def extract_amino_acid_sequence(pdb_file, chain_id):
    # Mapping of three-letter codes to single-letter codes
    aa_map = {
        'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
        'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
        'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
        'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'
    }

    # Create a PDB parser
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", pdb_file)

    sequence = []
    for model in structure:  # Loop over all models
        for chain in model:  # Loop over all chains
            if chain.id == chain_id:  # Check if the chain matches the target chain
                for residue in chain:  # Loop over all residues in the chain
                    # Check if the residue is a standard amino acid
                    if residue.resname in aa_map:
                        sequence.append(aa_map[residue.resname])
                break  # Exit loop after processing the target chain

    return ''.join(sequence)

if __name__ == '__main__':


    # Specify the PDB file path and chain ID
    pdb_file_path = "2nao.pdb"  # Replace with your PDB file path
    target_chain_id = "A"  # Specify the chain you want to extract

    # Get the amino acid sequence for the specified chain
    sequence = extract_amino_acid_sequence(pdb_file_path, target_chain_id)
    print(sequence)

