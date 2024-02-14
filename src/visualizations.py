# TODO: Takes in a contact card 2D and gives back a heat map
import matplotlib.pyplot as plt
import seaborn as sns


def heat_mapper(data: list, residue_one: str = "PDB One", residue_two: str = "PDB Two") -> None:
    plt.figure(figsize=(12, 10), dpi=300)
    ax = sns.heatmap(data, cmap='Spectral_r', yticklabels=20, xticklabels=20, annot=False, cbar=True,
                     cbar_kws={'label': 'Distance (Å)'})
    plt.title(f"Residue Heatmap of {residue_one} and {residue_two}", fontsize=12)
    plt.gca().invert_yaxis()
    plt.xlabel(f"{residue_one.title()} Residues")
    plt.ylabel(f"{residue_two.title()} Residues")
    plt.savefig(f"../data/{residue_one[:3]}_{residue_two[:3]}_residue_heatmap.png", dpi=300)
    plt.close()
