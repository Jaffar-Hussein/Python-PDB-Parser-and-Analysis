import argparse
from pprint import pprint
import pdb_analyzer
import pdb_parser
import atom_distance as ad
import visualizations as vis
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import freesasa
import pickle

from src.sasaumure import compute_sasa


def get_contact_card(residues_one, residues_two, mode):
    border_residues = []
    for threshold in range(1, 11):
        contact_card_chains = ad.residue_residue(residues_one, residues_two, threshold=threshold, mode=mode)
        border_residues.append(contact_card_chains)
    return border_residues


def calculate_interface_residues(contact_data):
    """Calculate unique residues at the interface for each threshold."""
    unique_residues_counts = []
    for threshold_data in contact_data:
        flattened_list = [item for sublist in threshold_data for item in sublist]
        unique_items = set(flattened_list)
        unique_residues_counts.append(len(unique_items))
    return unique_residues_counts


def main():
    pickle_file = '../data/1brs.pkl'
    # Parsing of the pdb file
    pdb_file = "../data/1brs.pdb"
    pdb_file_two = "../data/1FFW_AB_c.pdb"
    compute_sasa(pdb_file)

    with open(pickle_file, 'rb') as file:
        protein_structure = pickle.load(file)

    # Load the PDB data using your custom parser
    pdb_data = pdb_parser.pdb_parser_optimized(pdb_file)

    # Initialize a dictionary to hold the count of residues below each threshold
    data = {key: 0 for key in range(1, 11)}

    # Iterate over each threshold
    for threshold in data.keys():
        # Iterate over each chain and residue in the protein structure
        for chain in pdb_data["chains"]:
            for residue in pdb_data[chain]:
                # Access the SASA value for the current residue
                sasa = protein_structure[chain][residue]['SASA']
                # Increment the count for the current threshold if the SASA is below the threshold
                if sasa < threshold:
                    data[threshold] += 1
    sasa = data.values()
    # Parsing of the pdb file

    pdb_data = pdb_parser.pdb_parser_optimized(pdb_file)
    # print(f"{calculate_sasa(pdb_file) = }")
    pdf_data_2 = pdb_parser.pdb_parser_optimized(pdb_file_two)
    # pprint(pdb_data)
    name1 = pdb_file.split("/")[-1].split(".")[0]
    name2 = pdb_file_two.split("/")[-1].split(".")[0]
    # Analysis of the pdb file

    chain_count = pdb_analyzer.chain_count(pdb_data)
    print(f"Chain count: {chain_count}")
    residue_count = pdb_analyzer.residue_count(pdb_data)
    print(f"Residue count: {residue_count}")
    aa_count = pdb_analyzer.total_aa(pdb_data, "LYS")
    print(f"{aa_count[0]} count: {aa_count[1]}")
    aa_per_chain = pdb_analyzer.total_aa_per_chain(pdb_data, "LYS")
    pprint(aa_per_chain)

    # Lists all the residues in the pdb dictionary
    residues = pdb_analyzer.residue_list(pdb_data)
    residues_two = pdb_analyzer.residue_list(pdb_data)

    # Returns one residue
    residue_1 = pdb_data["A"]["8"]
    residue_2 = pdb_data["A"]["12"]

    # Manipulation of the pdb file
    print(f"shortest: {ad.calculate_distance("atom", residue_1, residue_2)}")
    print(f"centroid: {ad.calculate_distance("centroid", residue_1, residue_2)}")
    # Inter chain contact card
    chain_one_residues = pdb_analyzer.residues_in_chain(pdb_data, "A")
    chain_two_residues = pdb_analyzer.residues_in_chain(pdb_data, "D")

    # Calculate contacts
    contact_atom = get_contact_card(chain_one_residues, chain_two_residues, mode="atom")
    contact_centroid = get_contact_card(chain_one_residues, chain_two_residues, mode="centroid")

    # Prepare DataFrame
    center = calculate_interface_residues(contact_centroid)
    atom = calculate_interface_residues(contact_atom)
    diff = np.array(list(sasa)) - np.array(center)
    data = {
        'Threshold': list(range(1, 11)),
        'Atom_Mode': calculate_interface_residues(contact_atom),
        'Centroid_Mode': calculate_interface_residues(contact_centroid),
        'SASA_Mode': sasa,
        'SASA - Centroid': diff
    }
    df = pd.DataFrame(data)

    vis.plot_interface_residues(df, "threshold.png")
    thresholds = list(range(1, 11))
    # Prepare the data for plotting
    data = {
        'Threshold': thresholds,
        'Difference': diff
    }

    # Convert the data to a DataFrame
    df = pd.DataFrame(data)
    vis.plot_interface_difference(df, "Difference_Threshold.png")

    # vis.density_mapper(contact_atom_chains, "atom")
    # vis.density_mapper(contact_centroid_chains, "centroid")

    # generate heat map
    # vis.heat_mapper(contact_card_chains, "chain_a", "chain_b")
    # Get the contact card for residues on the pdb file
    # contact_card = ad.residue_residue(residues, residues_two)
    # generate heat map
    # vis.heat_mapper(contact_card, name1, name2)


if __name__ == "__main__":
    main()
