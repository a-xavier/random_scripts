import requests
import pandas as pd
from io import StringIO

def _string_api_request(input, score_threshold: float) -> list[str]:
    """
    Make a request to the STRING API.
    
    Parameters:
    input (str): Input gene symbol.
    
    Returns:
    list[str]: List of network names.
    """
    string_api_url = "https://version-12-0.string-db.org/api"
    output_format = "tsv"
    method = "interaction_partners"

    request_url = "/".join([string_api_url, output_format, method])
    params = {
        "identifiers": "%0d".join(input),
        "species": 9606,
        "limit": 500,
        "caller_identity": "cairns_lab"
    }
    
    response = requests.post(request_url, data=params)
    
    try:
        df = pd.read_csv(StringIO(response.text), sep="\t")
        # df = df[df.score >= 0.9]
        # Only keep row if dscore or escore is greater than 0.9
        df = df[(df["dscore"] >= score_threshold) | (df["escore"] >= score_threshold)]
        return list(set(df["preferredName_A"]).union(set(df["preferredName_B"])))
    except Exception as e:
        raise("Error in STRING API request:", e)


def gene_to_network(gene: str, score_threshold: float = 0.9) -> list[str]:
    """
    Does first and second shell network analysis using STRING API.
    """
    print("--------------------")
    print("Getting gene network...")
    print(f"Gene: {gene}")
    print("--------------------")

    first_shell = _string_api_request([gene], score_threshold)
    print(f"First shell: {first_shell}")
    if not first_shell:
        raise ValueError(f"No first shell found for gene {gene}. Please check the gene symbol or try a different one.")
    second_shell = _string_api_request(first_shell, score_threshold)
    # Sort second shell
    second_shell = sorted(second_shell)
    print(f"Second shell: {second_shell}")

    return second_shell

def get_gene_boundaries(genes: list[str], assembly: str = "hg38", upstream_buffer: int = 5000, downstream_buffer: int = 1000) -> pd.DataFrame:
    """
    Get gene boundaries for a list of genes.
    
    Parameters:
    genes (list[str]): List of gene symbols.
    assembly (str): Genome assembly version ('hg19' or 'hg38').
    
    Returns:
    pd.DataFrame: DataFrame containing gene boundaries.
    """
    print("--------------------")
    print("Getting gene boundaries...")
    print(f"Genes: {genes}")
    print("--------------------")

    if not isinstance(genes, list):
        if isinstance(genes, str):
            genes = [genes]
        else:
            raise ValueError("genes must be a list of gene symbols or a single gene symbol as a string.")
    
    if assembly == "hg19":
        server = "https://grch37.rest.ensembl.org"
    elif assembly == "hg38":
        server = "https://rest.ensembl.org"
    else:
        raise ValueError("Invalid assembly. Must be 'hg19' or 'hg38'.")
    
    ext = "/lookup/symbol/homo_sapiens"
    headers={ "Content-Type" : "application/json", "Accept" : "application/json"}

    data = '{{ "symbols" : ["{}" ] }}'.format('", "'.join(genes))


    r = requests.post(server+ext, headers=headers, data=data)

    if not r.ok:
      r.raise_for_status()
    
    decoded = r.json()

    final_dict = {}

    for item in decoded:
        assembly = decoded[item]["assembly_name"]
        start = decoded[item]["start"]
        end = decoded[item]["end"]
        chrom = decoded[item]["seq_region_name"]
        symbol = decoded[item]["display_name"]
        strand = decoded[item]["strand"]
        final_dict[item] = {"start_BP":start,
                      "end_BP": end,
                      "CHR": chrom,
                      "Gene": symbol,
                      "strand": strand,
                      "assembly": assembly}

    # Adjust based on upstream and downsteam buffer 

    for gene, gene_dict in final_dict.items():
        if gene_dict["strand"] == 1:
            gene_dict["start_BP"] = gene_dict["start_BP"] - upstream_buffer
            gene_dict["end_BP"] = gene_dict["end_BP"] + downstream_buffer
        else:
            gene_dict["start_BP"] = gene_dict["start_BP"] - downstream_buffer
            gene_dict["end_BP"] = gene_dict["end_BP"] + upstream_buffer
      
    final_df = pd.DataFrame.from_dict(final_dict).transpose()

    # replace chrosmome names with numbers
    final_df["CHR"] = final_df["CHR"].astype(str)
    final_df["CHR"] = final_df["CHR"].replace({"X": 23, "Y": 24, "MT": 25, "M": 25})
    # Convert the column to integers
    final_df["CHR"] = pd.to_numeric(final_df["CHR"], errors="coerce")
    # Drop rows where conversion failed (NaN values)
    final_df = final_df.dropna(subset=["CHR"])
    # Convert the column to integers again after cleaning
    final_df["CHR"] = final_df["CHR"].astype(int)
    # Make sure start_BP and end_BP are ints
    final_df["start_BP"] = final_df["start_BP"].astype(int)
    final_df["end_BP"] = final_df["end_BP"].astype(int)

    return(final_df)