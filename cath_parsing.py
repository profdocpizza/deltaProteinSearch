import pandas as pd


def get_cath_string(C=None, A=None, T=None, H=None):
    # join CATH with . if not none
    cath_string = ""
    if C is not None:
        cath_string += str(C)
    if A is not None:
        cath_string += "." + str(A)
    if T is not None:
        cath_string += "." + str(T)
    if H is not None:
        cath_string += "." + str(H)
    return cath_string


def get_cath_domains(C=None, A=None, T=None, H=None):
    cath_string = get_cath_string(C, A, T, H)

    # on alanine server
    file_path = "/mnt/scratch/cath/cath-b-newest-all"
    with open(file_path, "r") as file:
        domains = file.read().split("\n")
    domains_df = parse_domains_to_df(domains)
    if cath_string is not None:
        # return filtered where cath value starts with cath string
        domains_df = domains_df[domains_df["cath"].str.startswith(cath_string)]
    return domains_df


def parse_domains_to_df(domains):
    rows = []
    for line in domains:
        if not line.strip():
            continue
        parts = line.split()
        pdb_code = parts[0][:4]
        domain = parts[0][4:]
        version = parts[1]
        cath = parts[2]
        chain_info_str = ""
        chains = parts[3].split(",")
        chain_parts = []
        for chain in chains:
            interval, chain_id = chain.split(":")
            try:
                start, end = parse_interval(interval)
            except Exception as e:
                print(
                    f"Something is wrong with {pdb_code+domain} interval: {interval}, chain {chain_id}",
                    e,
                )
                continue
            chain_parts.append(f"{chain_id}{start}-{end}")
        chain_info_str = "/".join(chain_parts)
        cath_categories = cath.split(".")
        rows.append(
            {
                "pdb_code": pdb_code,
                "domain": domain,
                "version": version,
                "cath": cath,
                "chain_info_str": chain_info_str,
                "C": int(cath_categories[0]),
                "A": int(cath_categories[1]),
                "T": int(cath_categories[2]),
                "H": int(cath_categories[3]),
            }
        )
    return pd.DataFrame(rows)


def parse_interval(s):
    # Check if the first character is '-'
    if s[0] == "-":
        start = "-" + s[1 : s.index("-", 1)]  # The first number is negative
    else:
        start = s[0 : s.index("-")]  # The first number is positive

    # Check if there are two '-' symbols indicating the second number is negative
    if "--" in s:
        end = "-" + s[s.rindex("-") + 1 :]  # The second number is negative
    else:
        end = s[s.rindex("-") + 1 :]  # The second number is positive

    return start, end
