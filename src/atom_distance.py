from pprint import pprint
import numpy as np


def euclidean_distance(coordinate_1: list, coordinate_2: list):
    """
    This function takes two lists as input, each representing a point in 3D space with x, y, and z coordinates.
    It calculates the Euclidean distance between these two points using the formula:
    sqrt((x2-x1)^2 + (y2-y1)^2 + (z2-z1)^2)
    :param coordinate_1: The first point, represented as a list [x,y,z].
    :param coordinate_2: The second point, represented as a list [x,y,z].
    :return: float: The Euclidean distance between the two points.
    """
    return np.linalg.norm(np.subtract(coordinate_1, coordinate_2))


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
    return [residue[atom]['x'], residue[atom]['y'], residue[atom]['z']]


def centroid_searcher(residue: dict):
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
    coordinates = [coordinate_extractor(residue, atom) for atom in residue['atomlist']]
    return np.mean(coordinates, axis=0)


def calculate_distance(mode: str, first_residue: dict, second_residue: dict) -> float:
    """
    Calculate the minimum distance between all pairs of atoms or the distance between centroids of two residues.
    """
    if mode not in ["atom", "centroid"]:
        raise ValueError("Mode must be 'atom' or 'centroid'.")

    if mode == "atom":
        # Extract and store coordinates for all atoms in each residue to prevent repetitive calculations
        coords1 = [np.array(coordinate_extractor(first_residue, atom)) for atom in first_residue["atomlist"]]
        coords2 = [np.array(coordinate_extractor(second_residue, atom)) for atom in second_residue["atomlist"]]

        # Calculate all pairwise distances and return the minimum
        distances = [np.linalg.norm(c1 - c2) for c1 in coords1 for c2 in coords2]
        return min(distances)

    elif mode == "centroid":
        # Use NumPy for efficient centroid calculation
        centroid_1 = np.mean(
            np.array([coordinate_extractor(first_residue, atom) for atom in first_residue["atomlist"]]), axis=0)
        centroid_2 = np.mean(
            np.array([coordinate_extractor(second_residue, atom) for atom in second_residue["atomlist"]]), axis=0)

        # Calculate and return the distance between centroids
        return np.linalg.norm(centroid_1 - centroid_2)


def residue_residue(residue1: list, residue2: list, threshold: float = float('inf'), mode: str = "atom") -> list:
    """
    Refactored function description...
    """
    contact_list = []

    if mode == "atom":
        for residue_one in residue1:
            distances = []
            for residue_two in residue2:
                if residue_one == residue_two:
                    # Skip redundant calculations for the same residue
                    continue
                distance = calculate_distance(mode, residue_one, residue_two)
                if distance < threshold:
                    distances.append(distance)
            # Append the list of distances for the current residue pair
            contact_list.append(distances)

    elif mode == "centroid":
        # Using numpy to calculate centroid distances
        centroids1 = np.array([centroid_searcher(r) for r in residue1])
        centroids2 = np.array([centroid_searcher(r) for r in residue2])
        # Calculate the distances between all centroids
        distances = np.linalg.norm(centroids1[:, np.newaxis, :] - centroids2[np.newaxis, :, :], axis=2)
        # Filter the distances by the threshold
        for dist_row in distances:
            contact_list.append(dist_row[dist_row < threshold].tolist())

    return contact_list
