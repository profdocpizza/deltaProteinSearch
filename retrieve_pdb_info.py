import requests
import pandas as pd
import time
from ratelimit import limits, sleep_and_retry
from tqdm import tqdm

import requests
import pandas as pd
import time
from ratelimit import limits, sleep_and_retry

from cath_parsing import get_cath_domains

# API endpoints
PDBE_API_URL = "https://www.ebi.ac.uk/pdbe/api/mappings/uniprot/"
UNIPROT_API_URL = "https://rest.uniprot.org/uniprotkb/"
REQUESTS_PER_MINUTE = 600  # Conservative rate limit
PERIOD = 60  # 1 minute


# Rate-limited function to fetch UniProt accession from PDB code
@sleep_and_retry
@limits(calls=REQUESTS_PER_MINUTE, period=PERIOD)
def fetch_uniprot_accession(pdb_code):
    try:
        response = requests.get(f"{PDBE_API_URL}{pdb_code.lower()}", timeout=10)
        response.raise_for_status()
        data = response.json()
        uniprot_data = data.get(pdb_code.lower(), {}).get("UniProt", {})
        if uniprot_data:
            return list(uniprot_data.keys())[0]  # Return first UniProt ID
        return None
    except requests.exceptions.RequestException as e:
        print(f"Error fetching UniProt for {pdb_code}: {e}")
        return None


# Rate-limited function to fetch keywords and EC number from UniProt accession
@sleep_and_retry
@limits(calls=REQUESTS_PER_MINUTE, period=PERIOD)
def fetch_uniprot_details(uniprot_id):
    try:
        response = requests.get(f"{UNIPROT_API_URL}{uniprot_id}", timeout=10)
        response.raise_for_status()
        data = response.json()
        keywords = [kw["name"] for kw in data.get("keywords", [])]
        ec_numbers = [
            ec["value"]
            for ec in data.get("proteinDescription", {})
            .get("recommendedName", {})
            .get("ecNumbers", [])
        ]
        return keywords, ec_numbers[0] if ec_numbers else None
    except requests.exceptions.RequestException as e:
        print(f"Error fetching details for {uniprot_id}: {e}")
        return [], None


# def main():

#     cath_domains = get_cath_domains(C=1)

#     # Load subset of data (first 1
#     df = pd.read_pickle(
#         "/home/tadas/code/deltaProteinSearch/outputs/filtered_cath_deltaprots.pkl"
#     )
#     df = df[df["orientation_code"].notna()]
#     df.sort_values("total_cost", inplace=True)
#     df = df[df["total_cost"] < 1.5]
#     df = df[(df["total_scale"] > 9.0) & (df["total_scale"] < 13.0)]
#     df["pdb_code"] = df["filename"].apply(lambda x: x.split("_")[0])
#     df["domain"] = df["filename"].apply(lambda x: x.split("_")[1])
#     df = pd.merge(df, cath_domains, how="left", on=["domain", "pdb_code"])
#     df["ec_number"] = None

#     # Process each row and fetch EC number
#     start_time = time.time()
#     for index, row in df.iterrows():
#         pdb_code = row["pdb_code"]
#         print(f"Processing {pdb_code}...")

#         # Step 1: Get UniProt accession
#         uniprot_id = fetch_uniprot_accession(pdb_code)
#         df.at[index, "uniprot_id"] = uniprot_id

#         # Step 2: Get keywords and EC number if UniProt ID exists
#         if uniprot_id:
#             keywords, ec_number = fetch_uniprot_details(uniprot_id)
#             df.at[index, "keywords"] = "; ".join(keywords) if keywords else None
#             df.at[index, "ec_number"] = ec_number
#         print(uniprot_id, ec_number, keywords)

#         time.sleep(0.001)  # Reduced delay for faster testing, adjust as needed

#     # Calculate elapsed time
#     elapsed_time = time.time() - start_time

#     # Display the updated dataframe
#     print("\nUpdated DataFrame:")
#     print(df)

#     # Save to CSV for inspection
#     df.to_csv(
#         "/home/tadas/code/deltaProteinSearch/outputs/filtered_cath_deltaprots_anotated.csv",
#         index=False,
#     )
#     print(f"\nProcessing completed in {elapsed_time:.2f} seconds.")


# if __name__ == "__main__":
#     main()
import asyncio
import aiohttp
import pandas as pd
import time
from tqdm import tqdm
from aiohttp import ClientSession

# API endpoints
RCSB_API_URL = "https://data.rcsb.org/rest/v1/core/entry/"
PDBE_API_URL = "https://www.ebi.ac.uk/pdbe/api/mappings/uniprot/"
UNIPROT_API_URL = "https://rest.uniprot.org/uniprotkb/"

# Throttling control
REQUESTS_PER_SECOND = 50  # Adjust this to control rate (e.g., 10 req/s = 600 req/min)


async def fetch_uniprot_accession(session, pdb_code, semaphore):
    """Fetch UniProt accession asynchronously with throttling."""
    async with semaphore:  # Limit concurrency
        try:
            async with session.get(
                f"{PDBE_API_URL}{pdb_code.lower()}", timeout=10
            ) as response:
                response.raise_for_status()
                data = await response.json()
                uniprot_data = data.get(pdb_code.lower(), {}).get("UniProt", {})
                return list(uniprot_data.keys())[0] if uniprot_data else None
        except (aiohttp.ClientError, asyncio.TimeoutError) as e:
            print(f"Error fetching UniProt for {pdb_code}: {e}")
            return None


async def fetch_uniprot_details(session, uniprot_id, semaphore):
    """Fetch keywords and EC number asynchronously with throttling."""
    async with semaphore:
        try:
            async with session.get(
                f"{UNIPROT_API_URL}{uniprot_id}", timeout=10
            ) as response:
                response.raise_for_status()
                data = await response.json()
                keywords = [kw["name"] for kw in data.get("keywords", [])]
                ec_numbers = [
                    ec["value"]
                    for ec in data.get("proteinDescription", {})
                    .get("recommendedName", {})
                    .get("ecNumbers", [])
                ]
                return keywords, ec_numbers[0] if ec_numbers else None
        except (aiohttp.ClientError, asyncio.TimeoutError) as e:
            print(f"Error fetching details for {uniprot_id}: {e}")
            return [], None


async def fetch_ec_from_rcsb(pdb_code, session, semaphore):
    """Fetch EC number directly from RCSB PDB API asynchronously."""
    async with semaphore:  # Limit concurrency
        try:
            async with session.get(
                f"{RCSB_API_URL}{pdb_code.lower()}", timeout=10
            ) as response:
                response.raise_for_status()
                data = await response.json()
                # Extract EC number from enzyme classification if available
                enzyme = data.get("rcsb_entry_info", {}).get(
                    "enzyme_classification", []
                )
                return enzyme[0] if enzyme else None
        except (aiohttp.ClientError, asyncio.TimeoutError) as e:
            print(f"Error fetching EC for {pdb_code}: {e}")
            return None


async def process_row(session, row, semaphore):
    """Process a single row asynchronously."""
    pdb_code = row["pdb_code"]
    print(f"Processing {pdb_code}...")

    # Step 1: Fetch UniProt ID
    uniprot_id = await fetch_uniprot_accession(session, pdb_code, semaphore)

    # Step 2: Fetch keywords and EC number if UniProt ID exists
    keywords, ec_number = [], None
    if uniprot_id:
        keywords, ec_number = await fetch_uniprot_details(
            session, uniprot_id, semaphore
        )

    return uniprot_id, "; ".join(keywords) if keywords else None, ec_number


async def main():
    cath_domains = get_cath_domains(C=1)

    # Load subset of data (first 1
    df = pd.read_pickle(
        "/home/tadas/code/deltaProteinSearch/outputs/filtered_cath_deltaprots.pkl"
    )
    df = df[df["orientation_code"].notna()]
    df.sort_values("total_cost", inplace=True)
    df = df[df["total_cost"] < 1.5]
    df = df[(df["total_scale"] > 9.0) & (df["total_scale"] < 13.0)]
    df["pdb_code"] = df["filename"].apply(lambda x: x.split("_")[0])
    df["domain"] = df["filename"].apply(lambda x: x.split("_")[1])
    df = pd.merge(df, cath_domains, how="left", on=["domain", "pdb_code"])
    df["ec_number"] = None

    # Throttling semaphore based on requests per second
    semaphore = asyncio.Semaphore(REQUESTS_PER_SECOND)

    # Start timing
    start_time = time.time()

    # Create async session and process rows
    async with ClientSession() as session:
        tasks = [process_row(session, row, semaphore) for _, row in df.iterrows()]
        results = await asyncio.gather(*tasks, return_exceptions=True)

    # Update DataFrame with results
    for index, result in enumerate(results):
        if isinstance(result, tuple):  # Successful result
            uniprot_id, keywords, ec_number = result
            df.at[index, "uniprot_id"] = uniprot_id
            df.at[index, "keywords"] = keywords
            df.at[index, "ec_number"] = ec_number
        else:  # Exception occurred
            print(f"Error at index {index}: {result}")

    # Calculate elapsed time
    elapsed_time = time.time() - start_time

    # Display results
    print("\nUpdated DataFrame:")
    print(df)

    # Save to CSV
    df.to_csv("pdb_with_uniprot_details_async.csv", index=False)
    print(f"\nProcessing completed in {elapsed_time:.2f} seconds.")


if __name__ == "__main__":
    asyncio.run(main())
