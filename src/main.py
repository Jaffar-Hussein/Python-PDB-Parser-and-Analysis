from pprint import pprint

import pdb_parser, pdb_analyzer
from src.atom_distance import calculate_distance


def main():
    pdb_file = "../data/1brs.pdb"

    # Parsing of the pdb file

    pdb_data = pdb_parser.pdb_parser_optimized(pdb_file)
    # pprint(pdb_data)

    # Analysis of the pdb file

    chain_count = pdb_analyzer.chain_count(pdb_data)
    print(f"Chain count: {chain_count}")
    residue_count = pdb_analyzer.residue_count(pdb_data)
    print(f"Residue count: {residue_count}")
    aa_count = pdb_analyzer.total_aa(pdb_data, "LYS")
    print(f"{aa_count[0]} count: {aa_count[1]}")
    aa_per_chain = pdb_analyzer.total_aa_per_chain(pdb_data, "LYS")
    pprint(aa_per_chain)

    # Manipulation of the pdb file
    print(f"shortest: {calculate_distance(pdb_data, "atom", [["A", "8"], ["A", "12"]])}")
    print(f"centroid: {calculate_distance(pdb_data, "centroid", [["A", "8"], ["A", "12"]])}")


if __name__ == "__main__":
    main()
