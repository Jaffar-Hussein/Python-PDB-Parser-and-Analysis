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


def calculate_distance(mode: str, first_residue: list, second_residue: list) -> float:
    """
    This function calculates the minimum distance between all pairs of atoms in two given residues or the distance between their centroids.

    If the mode is "atom", the function iterates over each atom in the first residue and each atom in the second residue, and calculates the
    Euclidean distance between them using the 'euclidean_distance' and 'coordinate_extractor' functions. The minimum of these distances is returned.

    If the mode is "centroid", the function calculates the centroids of the two residues using the 'centroid_searcher' function, and then calculates
    the Euclidean distance between these centroids using the 'euclidean_distance' function. This distance is returned.

    :param mode: A string that specifies the mode of distance calculation. It can be "atom" or "centroid".
    :param first_residue: A list of dictionaries, each representing an atom in the first residue. The keys are atom names and the values are dictionaries
    containing the atom's x, y, and z coordinates.
    :param second_residue: A list of dictionaries, each representing an atom in the second residue. The keys are atom names and the values are dictionaries
    containing the atom's x, y, and z coordinates.

    :return: A float representing the minimum distance between all pairs of atoms in the two residues (if mode is "atom"), or the distance between the centroids
    of the two residues (if mode is "centroid").
    """
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


def residue_residue(residue1: list, residue2: list) -> list:
    """
    This function calculates the minimum distance between all pairs of residues in two given sets of residues.

    The function iterates over each residue in the first set and each residue in the second set, and calculates the
    minimum distance between them using the 'calculate_distance' function. The distances are stored in a 2D list (contact_list),
    where the element at index [i][j] represents the minimum distance between the i-th residue in the first set and the j-th
    residue in the second set.

    :param residue1: A dictionary representing a residue. The keys are atom names and the values are dictionaries
    containing the atom's x, y, and z coordinates.
    :param residue2: A dictionary representing a residue. The keys are atom names and the values are dictionaries
    containing the atom's x, y, and z coordinates.

    :return: A 2D list (contact_list) where the element at index [i][j] represents the minimum distance between the i-th
    residue in the first set and the j-th residue in the second set.
    """
    contact_list = []
    for residue_one in residue1:
        distances = []
        for residue_two in residue2:
            distance = calculate_distance("atom", residue_one, residue_two)
            distances.append(int(distance))
        contact_list.append(distances)
    return contact_list
