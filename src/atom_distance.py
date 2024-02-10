def euclidean_distance(coordinate_1: list, coordinate_2: list):
    """
    This function takes two lists as input, each representing a point in 3D space with x, y, and z coordinates.
    It calculates the Euclidean distance between these two points using the formula:
    sqrt((x2-x1)^2 + (y2-y1)^2 + (z2-z1)^2)
    :param coordinate_1: The first point, represented as a list [x,y,z].
    :param coordinate_2: The second point, represented as a list [x,y,z].
    :return: float: The Euclidean distance between the two points.
    """
    return (
        ((coordinate_1[0] - coordinate_2[0]) ** 2)
        + ((coordinate_1[1] - coordinate_2[1]) ** 2)
        + ((coordinate_1[2] - coordinate_2[2]) ** 2)
    ) ** 0.5


def coordinate_extractor(residue: dict, atom: str) -> list:
    """
    This function takes a residue and an atom name as input. The residue is represented as a dictionary where each
    atom is a key and its value is another dictionary containing the atom's x, y, and z coordinates. The function
    returns the x, y, and z coordinates of the specified atom as a list.
    :param residue: A dictionary representing a residue. The keys are atom names and the values are dictionaries
    containing the atom's x, y, and z coordinates.
    :param atom: The name of the atom whose coordinates are to be
    extracted.
    :return: A list representing the x, y, and z coordinates of the specified atom.
    """
    n = residue[atom]
    return [n["x"], n["y"], n["z"]]


def centroid_searcher(residue: dict) -> list:
    """
    Calculate the centroid of a residue.
    This function calculates the centroid (geometric center) of a residue in a protein structure. The residue is
    represented as a dictionary where each atom is a key and its value is another dictionary containing the atom's x,
    y, and z coordinates. The centroid is calculated by averaging the x, y, and z coordinates of all atoms in the
    residue.
    :param residue: A dictionary representing a residue. The keys are atom names and the values are dictionaries
    containing the atom's x, y, and z coordinates.
    :return: A list representing the x, y, and z coordinates of the
    centroid.
    """
    total_x = total_y = total_z = 0.0
    atom_count = len(residue["atomlist"])
    for atom in residue["atomlist"]:
        total_x += residue[atom]["x"]
        total_y += residue[atom]["y"]
        total_z += residue[atom]["z"]
    return [total_x / atom_count, total_y / atom_count, total_z / atom_count]


def calculate_distance(pdb_data: dict, mode: str, residues: list) -> float:
    """
    Calculate the distance between two residues based on the specified mode.
    This function calculates the distance between two residues in a protein structure. The residues are specified as
    a list of tuples, where each tuple contains the chain identifier and the residue number. The mode parameter
    determines how the distance is calculated: - "atom": The function calculates all pairwise distances between atoms
    of the two residues and returns the minimum distance. - "centroid": The function calculates the distance between
    the centroids of the two residues.

    :param pdb_data: A dictionary representing the protein structure. The keys are chain identifiers and the values
    are dictionaries, where the keys are residue numbers and the values are dictionaries representing the residue.
    :param mode: A string specifying the mode of distance calculation. It can be either "atom" or "centroid".
    :param residues: A list of two tuples, where each tuple contains the chain identifier and the residue number of a
    residue.
    :return: The calculated distance between the two residues.
    """
    first_residue = pdb_data[residues[0][0]][residues[0][1]]
    second_residue = pdb_data[residues[1][0]][residues[1][1]]
    if mode == "atom":
        distances = []
        for atom_one in first_residue["atomlist"]:
            for atom_two in second_residue["atomlist"]:
                distance = euclidean_distance(
                    coordinate_extractor(first_residue, atom_one),
                    coordinate_extractor(second_residue, atom_two),
                )
                distances.append(distance)
        return min(distances)
    elif mode == "centroid":
        centroid_1 = centroid_searcher(first_residue)
        centroid_2 = centroid_searcher(second_residue)
        return euclidean_distance(centroid_1, centroid_2)
