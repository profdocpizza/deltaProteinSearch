import os
import numpy as np
import pandas as pd
from cath_parsing import get_cath_domains
from pymol import cmd
import ampal

# External utility imports (ensure these functions are available in your PYTHONPATH)
from dp_utils.helix_assembly.utils import apply_transform_to_assembly
from dp_utils.permutation_data import read_flip_permutations_and_rearrangements
from dp_utils.run_pipeline import generate_helix_assembly, select_path_choice
from isambard.specifications.deltaprot import DeltaProt
from deltaprot_search import get_decompressed_pdb, load_domain_as_ampal_assembly


# These functions must be defined or imported as appropriate:
# get_cath_domains, get_decompressed_pdb


def visualize_with_pymol(
    full_pdb,
    domain_pdb,
    deltaprot_pdb,
    pdb_code,
    domain_string,
    orientation_code,
    save_dir,
):
    """Load structures into PyMOL, set colors, and save a session."""
    cmd.load(full_pdb, pdb_code)
    cmd.load(domain_pdb, domain_string)
    cmd.load(deltaprot_pdb, orientation_code)

    # Color the objects:
    cmd.color("palecyan", f"{domain_string} and polymer.protein")
    cmd.color("grey70", f"{pdb_code} and polymer.protein")
    cmd.color("violet", f"{orientation_code} and polymer.protein")
    cmd.set("cartoon_transparency", 0.1, "polymer.protein")
    cmd.orient(orientation_code)

    # Save the PyMOL session
    pymol_sessions_dir = os.path.join(save_dir, "pymol_sessions")
    if not os.path.exists(pymol_sessions_dir):
        os.makedirs(pymol_sessions_dir)
    session_file = os.path.join(pymol_sessions_dir, f"{orientation_code}_session.pse")
    cmd.save(session_file)
    cmd.delete("all")


def process_orientation(orientation_code, df, cath_domains, all_paths, save_dir):
    """Process a single orientation, generate assemblies, save PDBs, and visualize in PyMOL."""
    # Select the best alignment row for the given orientation
    best_alignment_row = df[df["orientation_code"] == orientation_code].iloc[0]
    domain_string = best_alignment_row["filename"]
    pdb_code, domain = domain_string.split("_")[:2]

    # Get CATH domain info (assuming a unique match)
    cath_domain_info = [
        i for i in cath_domains if (i["pdb_code"] == pdb_code and i["domain"] == domain)
    ][0]

    # Load assemblies for the domain and full PDB
    cath_domain_assembly = load_domain_as_ampal_assembly(cath_domain_info)
    pdb_file = get_decompressed_pdb(pdb_code, biounit=False)
    cath_pdb_assembly = ampal.load_pdb(pdb_file, path=False)
    if isinstance(cath_pdb_assembly, ampal.assembly.AmpalContainer):
        cath_pdb_assembly = cath_pdb_assembly[0]

    # Generate the delta prot assembly
    path_choice = select_path_choice(all_paths, orientation_code, "Best", 3)
    deltaprot_assembly = generate_helix_assembly(
        path_choice,
        deltahedron_edge_length=best_alignment_row["total_scale"],
        residues_per_helix=10,
        single_chain=False,
    )
    apply_transform_to_assembly(
        deltaprot_assembly,
        rotation_matrix=np.array(best_alignment_row["rotation_matrix"]),
        translation_vector=np.array(best_alignment_row["translation_matrix"]),
        translate_first=True,
    )

    # Translate the assemblies to align centers
    deltaprot_center = np.mean(
        [i["center"] for i in best_alignment_row["transformed_orientation_data"]],
        axis=0,
    )
    deltaprot_assembly.translate(deltaprot_center - deltaprot_assembly.centre_of_mass)
    translation_to_zero = -deltaprot_assembly.centre_of_mass
    deltaprot_assembly.translate(translation_to_zero)
    cath_pdb_assembly.translate(translation_to_zero)
    cath_domain_assembly.translate(translation_to_zero)

    # Save the assemblies as PDB files
    orientation_save_dir = os.path.join(save_dir, orientation_code)
    if not os.path.exists(orientation_save_dir):
        os.makedirs(orientation_save_dir)
    full_pdb_path = os.path.join(
        orientation_save_dir, f"{orientation_code}_{domain_string}_full.pdb"
    )
    domain_pdb_path = os.path.join(
        orientation_save_dir, f"{orientation_code}_{domain_string}_domain.pdb"
    )
    deltaprot_pdb_path = os.path.join(orientation_save_dir, f"{orientation_code}.pdb")

    with open(full_pdb_path, "w") as f:
        f.write(cath_pdb_assembly.pdb)
    with open(domain_pdb_path, "w") as f:
        f.write(cath_domain_assembly.pdb)
    with open(deltaprot_pdb_path, "w") as f:
        f.write(deltaprot_assembly.pdb)

    # Visualize the structures using PyMOL
    visualize_with_pymol(
        full_pdb_path,
        domain_pdb_path,
        deltaprot_pdb_path,
        pdb_code,
        domain_string,
        orientation_code,
        save_dir,
    )


def main():
    # Define paths and directories
    pickle_path = (
        "/home/tadas/code/deltaProteinSearch/outputs/filtered_cath_deltaprots.pkl"
    )
    save_dir = "/home/tadas/code/deltaProteinSearch/outputs/best_hits_per_orientation"
    if not os.path.exists(save_dir):
        os.makedirs(save_dir)

    # Load dataframe and perform filtering/sorting in-line
    df = pd.read_pickle(pickle_path)
    df.sort_values("total_cost", inplace=True)
    df = df[df["total_cost"] < 1.5]
    df = df[(df["total_scale"] > 9.0) & (df["total_scale"] < 14.0)]

    # Get CATH domains (assumed to be implemented elsewhere)
    cath_domains = get_cath_domains(C=1)

    # Read permutation and rearrangement paths
    all_paths = read_flip_permutations_and_rearrangements(read_json=False)

    # Process each orientation code
    for orientation_code in DeltaProt.orientation_codes:
        process_orientation(orientation_code, df, cath_domains, all_paths, save_dir)


if __name__ == "__main__":
    main()
