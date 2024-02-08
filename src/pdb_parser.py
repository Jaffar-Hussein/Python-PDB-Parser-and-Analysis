def pdb_parser_optimized(file: str) -> dict:
    """
    This function reads a PDB file line by line, and for each line that starts with "ATOM", it splits the line into
    fields, extracts the necessary information, and stores it in a dictionary.

    :param file: The name of the PDB file to parse.
    :return: A dictionary containing the parsed PDB data.
    """
    pdb_dict = {"chains": []}

    with open(file, "r") as f:
        for line in f:
            if line.startswith("ATOM"):
                fields = line.split()

                chain = fields[4]

                chain_dict = pdb_dict.setdefault(chain, {})

                current_residue = fields[5]
                residue_dict = chain_dict.setdefault(
                    current_residue, {"atomlist": [], "resname": fields[3]}
                )

                atom_type = fields[2]
                atom_dict = residue_dict.setdefault(
                    atom_type, {"id": fields[1], "bfactor": fields[10]}
                )

                atom_dict["x"] = float(fields[6])
                atom_dict["y"] = float(fields[7])
                atom_dict["z"] = float(fields[8])

                if chain not in pdb_dict["chains"]:
                    pdb_dict["chains"].append(chain)
                if atom_type not in residue_dict["atomlist"]:
                    residue_dict["atomlist"].append(atom_type)
    return pdb_dict
