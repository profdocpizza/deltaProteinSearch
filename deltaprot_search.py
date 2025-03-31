import gzip
import json
import os
import ampal
import multiprocessing
from dpFinder.cluster_helix_data import cluster_helix_axes
from dpFinder.deltaprot_finder_utils import (
    calculate_iterations,
    file_size_too_large,
    find_dp_in_helix_axes,
    get_helix_axes,
)
from dpFinder.find_deltaprots import dict_to_dataframe
import pandas as pd
from tqdm import tqdm

from cath_parsing import get_cath_domains
import re

CONFIG = {
    "cath_class": 1,  # Mainly alpha
    "max_file_size_MB": 1,
    "max_helices_in_domain": 6,
}
OUTPUT_DIR = "/home/tadas/code/deltaProteinSearch/outputs"


def get_decompressed_pdb(pdb_code, biounit=True):
    compressed_file = locate_pdb(pdb_code, biounit)
    return gzip.open(compressed_file, "rb").read().decode()


def locate_pdb(pdb_code, biounit):
    if biounit:
        return f"/mnt/scratch/datasets/biounit/{pdb_code[1:-1]}/{pdb_code}.pdb1.gz"
    else:
        return f"/mnt/scratch/datasets/pdb/{pdb_code[1:-1]}/pdb{pdb_code}.ent.gz"


############ MODIFIED DP FINDER FUNCTIONS ################
def find_dp_in_assembly(assembly, known_orientation=None):

    try:
        helix_axes = get_helix_axes(assembly=assembly)
    except Exception as e:
        print(f"Error processing {assembly.id}: {e}")
        return {
            "helix_axes": None,
            "MF_orientations_results": None,
        }

    if len(helix_axes) == 0:
        print(f"No helix axes for {assembly.id} identified")
        return {
            "helix_axes": helix_axes,
            "MF_orientations_results": None,
        }
    
    helix_axes_clusters = cluster_helix_axes(helix_axes)
    if not helix_axes_clusters:
        print(f"No valid clusters found for {assembly.id} identified")
        return {
            "helix_axes": helix_axes,
            "MF_orientations_results": None,
        }

    MF_orientations_results = find_dp_in_helix_axes(
        helix_axes, helix_axes_clusters, known_orientation=known_orientation, name=assembly.id
    )

    pdb_data = {
        "helix_axes": helix_axes,
        "MF_orientations_results": MF_orientations_results,
    }
    return pdb_data


def process_single_assembly(assembly):
    print(f"Processing {assembly}...")
    dp_finder_results = find_dp_in_assembly(assembly)
    return assembly.id, dp_finder_results


def get_domain_name(domain):
    return domain["pdb_code"] + "_" + domain["domain"]



def load_domain_as_ampal_assembly(domain_row):
    try:
        pdb_file = get_decompressed_pdb(domain_row["pdb_code"], biounit=False)
        domain_name = get_domain_name(domain_row)
        structure = ampal.load_pdb(pdb_file, path=False)
        if structure.__class__ == ampal.assembly.AmpalContainer:
            structure = structure[0]
        domain_structures = []
        chains_info = domain_row["chain_info_str"].split("/")
        for chain in chains_info:

            match = re.fullmatch(r"([A-Za-z])(\d+)-(\d+)", chain)
            if not match:
                print(
                    f"Error parsing {chain} in {load_domain_as_ampal_assembly.__name__}"
                )
                continue
            chain_id = match.group(1)
            start = int(match.group(2))
            end = int(match.group(3))
            selection = structure[chain_id][start:end]
            domain_structures.append(selection)
        return ampal.Assembly(domain_structures, assembly_id=domain_name)
    except FileNotFoundError:
        print(
            f"File not found for {domain_row['pdb_code']} in {load_domain_as_ampal_assembly.__name__}"
        )
        return None
    except Exception as e:
        print(
            f"Error processing {domain_row['pdb_code']} in {load_domain_as_ampal_assembly.__name__}: {e}"
        )
        return None


def process_domain(domain_row):
    print(domain_row)
    """Process a single domain and return (domain_name, dp_finder_results)."""
    compressed_file = locate_pdb(domain_row["pdb_code"], biounit=False)
    if not os.path.exists(compressed_file):
        print(f"File not found: {compressed_file}")
        return None
    if file_size_too_large(compressed_file, threshold_mb=CONFIG["max_file_size_MB"]):
        print(f"File too large: {compressed_file}")
        return None
    
    assembly = load_domain_as_ampal_assembly(domain_row)
    if assembly is None:
        return None

    dp_finder_results = find_dp_in_assembly(assembly)
    return (get_domain_name(domain_row), dp_finder_results)


def run_serial(domains_df):
    results = {}
    for index, domain_row in tqdm(domains_df.iterrows(), total=len(domains_df)):
        processed = process_domain(domain_row)
        if processed is not None:
            domain_name, dp_data = processed
            results[domain_name] = dp_data
    return results


def run_parallel(domains_df, num_cpu):
    results = {}
    # Convert the DataFrame into a list of dictionaries
    domains = domains_df.to_dict('records')
    with multiprocessing.Pool(processes=num_cpu) as pool:
        for processed in tqdm(
            pool.imap_unordered(process_domain, domains), total=len(domains)
        ):
            if processed is not None:
                domain_name, dp_data = processed
                results[domain_name] = dp_data
    return results


if __name__ == "__main__":

    pkl_file = os.path.join(OUTPUT_DIR, "cath_deltaprots.pkl")

    # Get domains from CATH
    cath_domains = get_cath_domains()

#     good_deltaprots = pd.read_csv(
#         "/home/tadas/code/deltaProteinSearch/outputs/pdb_with_uniprot_details_async.csv"
#     )
#     good_deltaprots = good_deltaprots[good_deltaprots["total_cost"]<1]
#     good_deltaprots = good_deltaprots[good_deltaprots['orientation_code'].astype(str).str.contains("6", na=False)]
#     # 
#     # only keep rows in cath_domains where pdb_code value exists in good_deltaprots pdb_code values
#     cath_domains = cath_domains.merge(
#     good_deltaprots[['pdb_code', 'domain']], 
#     on=['pdb_code', 'domain'], 
#     how='inner'
# )
#     print(f"remaining len {len(cath_domains)}")


    # Determine number of CPUs to use (can be parameterized as needed)
    num_cpu = 25  # multiprocessing.cpu_count()
    # For debugging, you might want to set: num_cpu = 1
    if num_cpu == 1:
        print("Running in serial mode...")
        cath_deltaprots = run_serial(cath_domains)
    else:
        print(f"Running in parallel mode with {num_cpu} CPUs...")
        cath_deltaprots = run_parallel(cath_domains, num_cpu)
    # Save results as JSON and pickle dataframe.
    json_file = os.path.join(OUTPUT_DIR, "cath_deltaprots.json")
    with open(json_file, "w") as f:
        json.dump(cath_deltaprots, f)
    print(f"Results saved to JSON: {json_file}")
    df = dict_to_dataframe(cath_deltaprots)
    df.to_pickle(pkl_file)
    print(f"Dataframe saved as pickle: {pkl_file}")

    df = pd.read_pickle(pkl_file)


# cath_domains = get_cath_domains(C=1)


# cath_deltaprots = {}
# for domain in tqdm(cath_domains[:50]):
#     compressed_file = locate_pdb(domain["pdb_code"], biounit=False)
#     if not os.path.exists(compressed_file):
#         print(f"File not found: {compressed_file}")
#         continue
#     if file_size_too_large(compressed_file, threshold_mb=CONFIG["max_file_size_MB"]):
#         print(f"File too large: {compressed_file}")
#         continue

#     assembly = load_domain_as_ampal_assembly(domain)
#     if assembly is None:
#         continue

#     dp_finder_results = find_dp_in_assembly(assembly)
#     cath_deltaprots[get_domain_name(domain)] = dp_finder_results


# # save dict as json in /home/tadas/code/deltaProteinSearch/outputs
# output_dir = "/home/tadas/code/deltaProteinSearch/outputs"
# with open(os.path.join(output_dir, "cath_deltaprots.json"), "w") as f:
#     json.dump(cath_deltaprots, f)

# df = dict_to_dataframe(cath_deltaprots)
# df.to_pickle(os.path.join(output_dir, "cath_deltaprots.pkl"))
