#/Users/maghoi/opt/anaconda3/envs/py310/bin/python
import sys
import numpy as np
import biotite.structure.io as strucio
import biotite.structure.io.pdb as pdb

def get_n_residues_chains(A):
    n_residues = (A.atom_name == 'CA').sum()
    chains = list(np.unique(A.chain_id))

    return n_residues, chains

def extract_and_save_chains(input_pdb_path, output_pdb_path, output_chains):

    # Load atom array
    pdb_file = pdb.PDBFile.read(input_pdb_path)
    A_in = pdb_file.get_structure(model=1, extra_fields=["b_factor"])
    n_residues, chains = get_n_residues_chains(A_in)
    print(f"Read {n_residues} residues, chains {chains} from {input_pdb_path}")

    # Extract the chains
    mask = np.full(len(A_in), False)
    for c in output_chains:
        m = A_in.chain_id == c
        mask += m

    # Only extract selected chains
    A_out = A_in[mask >= 1]

    # Save the new structure to a PDB file
    n_residues, chains = get_n_residues_chains(A_out)
    strucio.save_structure(output_pdb_path, A_out)
    print(f"Saving {n_residues} residues, chains {chains} to {output_pdb_path}")

if __name__ == "__main__":
    if len(sys.argv) < 4:
        print("Usage: python script.py input.pdb output.pdb chain1 chain2 ...")
        sys.exit(1)

    input_pdb_path = sys.argv[1]
    output_pdb_path = sys.argv[2]
    chains = sys.argv[3:]

    assert output_pdb_path[-4:] == ".pdb"

    extract_and_save_chains(input_pdb_path, output_pdb_path, chains)
