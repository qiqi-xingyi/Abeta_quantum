# --*-- conding:utf-8 --*--
# @Time : 11/14/24 12:48 AM
# @Author : Yuqi Zhang
# @Email : yzhan135@kent.edu
# @File : create_docking_file.py


from files_tool import DockingFilePreparer


if __name__ == "__main__":

    preparer = DockingFilePreparer("process_data/6mu3_L/PDBbind_data/6mu3/6mu3_pocket.pdb")
    preparer.prepare_pdbqt()


