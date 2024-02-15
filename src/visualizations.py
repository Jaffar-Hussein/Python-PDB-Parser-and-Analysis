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


def density_mapper(distances):
    """
    This function generates a series of Kernel Density Estimation (KDE) plots for a given list of distances.
    Each plot corresponds to a different threshold value.

    Parameters:
    distances (list): A list of lists, where each sublist represents the distances for a specific threshold.

    Returns:
    None: The function saves the generated plots as a .png file and does not return any value.
    """
    # Set the style of the plots
    sns.set(style="ticks")

    # Generate a color palette for the plots
    palette = sns.color_palette("tab10", len(distances))

    # Create a new figure with specific size
    fig, axes = plt.subplots(3, 3, figsize=(15, 10))

    # For each distance list in the distances list
    for i, ax in enumerate(axes.flatten()[:7]):
        # Generate a KDE plot for the current distance list
        sns.kdeplot(distances[i], fill=True, ax=ax, color=palette[i % 10], label=f'Threshold {i + 4}')

        # Set the labels for the x and y axes
        ax.set_xlabel('Distance (Å)', fontsize=12)
        ax.set_ylabel('Density', fontsize=12)

        # Set the parameters for the ticks on both axes
        ax.tick_params(axis='both', which='major', labelsize=10)

        # Add a legend to the plot
        ax.legend(fontsize=10, loc='upper left')

    # For the remaining subplots that are not used, turn off the axis
    for ax in axes.flatten()[7:]:
        ax.axis('off')

    # Adjust the spacing between the subplots
    plt.subplots_adjust(hspace=0.4, wspace=0.3)

    # Adjust the layout of the subplots to fit the figure area
    plt.tight_layout()

    # Save the figure as a .png file
    plt.savefig('../data/density_plot.png', dpi=300)

    # Close the figure to free up memory
    plt.close()

