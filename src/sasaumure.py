import os
import pickle
import sys
from pprint import pprint

import freesasa
import pdb_parser
import pdb_analyzer


def compute_sasa(pdb_file):
    """
    Computes the Solvent Accessible Surface Area (SASA) for each residue in a protein structure.

    Args:
    - pdb_file (str): Path to the PDB file to analyze.

    Returns:
    - d_pdb (dict): Updated dictionary with SASA information for each residue.
    """

    # Extract the base name of the PDB file without the extension
    rootname = os.path.basename(pdb_file).split(".")[0]

    # Define the output file name for the serialized data
    outfile = os.path.splitext(pdb_file)[0] + ".pkl"
    print(f"Processing PDB file: {pdb_file} -> Output: {outfile}")

    # Parse the PDB file to get a dictionary representation of the structure
    d_pdb = pdb_parser.pdb_parser_optimized(pdb_file)

    # Use freesasa to calculate SASA for the structure
    strpdb = freesasa.Structure(pdb_file)
    outASA = freesasa.calc(strpdb)
    sasa = outASA.residueAreas()

    # Get the chains present in both the parsed structure and the SASA results
    pdb_chains = set(d_pdb['chains'])
    sasa_chains = set(sasa.keys())
    common_chains = sasa_chains & pdb_chains

    # Warn if there's a mismatch in the chains between the parsed structure and SASA results
    if common_chains != sasa_chains or common_chains != pdb_chains:
        print(f"WARNING: Not all chains have SASA ({', '.join(common_chains)} out of {', '.join(pdb_chains)})")

    # Update the parsed structure dictionary with SASA information for each residue
    for ch in common_chains:
        for res in d_pdb[ch]:
            if res in sasa[ch]:
                d_pdb[ch][res]['rSASA'] = sasa[ch][res].relativeTotal  # Relative SASA
                d_pdb[ch][res]['SASA'] = sasa[ch][res].total  # Absolute SASA

    # Serialize the updated dictionary to a file for later use
    pickle.dump(d_pdb, open(outfile, 'wb'))
    print("SASA calculation successful!")


if __name__ == "__main__":
    pass
