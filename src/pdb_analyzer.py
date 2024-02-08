from pprint import pprint


def chain_count(pdb_data: dict) -> int:
    """
    Gives the number of chains in the pdb
    :param pdb_data: The dictionary holding the pdb data
    :return: The number of chains
    """
    return len(pdb_data["chains"])


def residue_count(pdb_data: dict) -> int:
    """
    Gives the number of residues in the pdb file
    :param pdb_data: The dictionary holding the pdb data
    :return: The number of residues
    """
    total = 0
    for chain in pdb_data["chains"]:
        total += len(pdb_data[chain])
    return total


def total_aa(pdb_data: dict, aa: str) -> list:
    """
    Gives the number of amino acids in the pdb
    :param pdb_data: The dictionary holding the pdb data
    :param aa: The aa to count
    :return: [aa, total]
    """
    total = 0
    for chain in pdb_data["chains"]:
        for residue in pdb_data[chain]:
            if pdb_data[chain][residue]["resname"] == aa:
                total += 1
    return [aa, total]


def total_aa_per_chain(pdb_data: dict, aa: str) -> dict:
    """
    Gives the number of amino acids per chain
    :param pdb_data: The dictionary holding the pdb data
    :param aa: The aa to count
    :return: dictionary with the number of aa per chain
    """
    response = {aa: {}}
    for chain in pdb_data["chains"]:
        total = 0
        for residue in pdb_data[chain]:
            if pdb_data[chain][residue]["resname"] == aa:
                total += 1
        response[aa][chain] = total
    return response
