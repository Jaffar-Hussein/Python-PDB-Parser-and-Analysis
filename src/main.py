from pprint import pprint

import pdb_parser


def main():
    pdb_file = "../data/1brs.pdb"
    pdb_data = pdb_parser.pdb_parser_optimized(pdb_file)
    pprint(pdb_data)


if __name__ == "__main__":
    main()
