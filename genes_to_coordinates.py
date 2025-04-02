import requests
import pandas as pd

def genes_to_coordinates(genes: list[str], assembly: str) -> pd.DataFrame:
    """
    This function takes a list of genes and an assembly as input and returns a pandas DataFrame containing the coordinates of those genes.
    """

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
        
    final_df = pd.DataFrame.from_dict(final_dict).transpose()

    return(final_df)