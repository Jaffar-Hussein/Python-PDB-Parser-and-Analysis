# TODO: Takes in a contact card 2D and gives back a heat map
import matplotlib.pyplot as plt
import seaborn as sns


def heat_mapper(data: list, residue_one: str = "PDB One", residue_two: str = "PDB Two") -> None:
    """
    This function generates a heatmap based on the provided data.

    Parameters:
    data (list): A 2D list representing the contact card data.
    residue_one (str, optional): The name of the first residue. Defaults to "PDB One".
    residue_two (str, optional): The name of the second residue. Defaults to "PDB Two".

    Returns:
    None: The function saves the generated heatmap as a .png file and does not return any value.
    """
    # Create a new figure with specific size and resolution
    plt.figure(figsize=(12, 10), dpi=300)

    # Generate a heatmap using seaborn
    sns.heatmap(data, cmap='Spectral_r', yticklabels=20, xticklabels=20, annot=False, cbar=True,
                cbar_kws={'label': 'Distance (Å)'})

    # Set the title of the heatmap
    plt.title(f"Residue Heatmap of {residue_one} and {residue_two}", fontsize=12)

    # Invert the y-axis
    plt.gca().invert_yaxis()

    # Set the labels for the x and y axes
    plt.xlabel(f"{residue_one.title()} Residues")
    plt.ylabel(f"{residue_two.title()} Residues")

    # Save the heatmap as a .png file
    plt.savefig(f"../data/{residue_one[:3]}_{residue_two[:3]}_residue_heatmap.png", dpi=300)

    # Close the figure to free up memory
    plt.close()
