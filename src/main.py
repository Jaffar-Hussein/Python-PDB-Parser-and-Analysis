from pprint import pprint
import argparse
import pdb_analyzer
import pdb_parser
import atom_distance as ad
import visualizations as vis


def main():
    THRESHOLD: int = 5
    # Parsing of the pdb file
    pdb_file = "../data/1brs.pdb"
    pdb_file_two = "../data/1FFW_AB_c.pdb"

    # Parsing of the pdb file

    pdb_data = pdb_parser.pdb_parser_optimized(pdb_file)
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
    residues_two = pdb_analyzer.residue_list(pdf_data_2)

    # Returns one residue
    residue_1 = pdb_data["A"]["8"]
    residue_2 = pdb_data["A"]["12"]

    # Manipulation of the pdb file
    print(f"shortest: {ad.calculate_distance("atom", residue_1, residue_2)}")
    print(f"centroid: {ad.calculate_distance("centroid", residue_1, residue_2)}")
    # Inter chain contact card
    chain_one_residues = pdb_analyzer.residues_in_chain(pdb_data, "A")
    chain_two_residues = pdb_analyzer.residues_in_chain(pdb_data, "B")
    pprint(chain_one_residues)

    # Get the contact card for residues on different pdb chains
    border_residues = []
    for threshold in range(4, 11):
        contact_card_chains = ad.residue_residue(chain_one_residues, chain_two_residues, threshold=threshold)
        border_residues.append(contact_card_chains)

    vis.density_mapper(border_residues)
    # generate heat map
    # vis.heat_mapper(contact_card_chains, "chain_a", "chain_b")
    # Get the contact card for residues on the pdb file
    contact_card = ad.residue_residue(residues, residues_two)
    # generate heat map
    vis.heat_mapper(contact_card, name1, name2)

    # move to output file
    with open("../data/contact_card.txt", "w") as f:
        for row in contact_card:
            f.write(" ".join([str(i) for i in row]) + "\n")


if __name__ == "__main__":
    main()
