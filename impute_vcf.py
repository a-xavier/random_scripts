#!/usr/bin/env python3

# Script to impute missing genotypes in a VCF file
# Uses Beagle 5.5 
# Split vcf with SnpSift
# Needs bcftools to merge / filter 

# INPUT VCF NEEDS TO BE GRCH37 ** WITHOUT ** chr prefix

from subprocess import run, CalledProcessError
import argparse
import os
from glob import glob
from sys import exit
# DEBUG
DEBUG = True

# Constants to store path for beagle and snpsift
BEAGLE = "/Users/alex/Library/CloudStorage/OneDrive-TheUniversityofNewcastle/Cairns_Group/NF1_Array_Project/beagle.jar"
BEAGLE_PANEL_FOLDER="/Users/alex/Library/CloudStorage/OneDrive-TheUniversityofNewcastle/Cairns_Group/NF1_Array_Project/1000kg_VCFs"
BEAGLE_MAP_FOLDER="/Users/alex/Library/CloudStorage/OneDrive-TheUniversityofNewcastle/Cairns_Group/NF1_Array_Project/1000kg_VCFs/plink.GRCh37.map"
SNPSIFT = "/Users/alex/Library/CloudStorage/OneDrive-TheUniversityofNewcastle/Cairns_Group/NF1_Array_Project/VCFs/snpEff/SnpSift.jar"
# bcftools is in the PATH
# bgzip is in the PATH
# tabix is in the PATH
# Java is in the PATH


def debug_print(message: str, end = "\n", flush = False) -> None:
    if DEBUG:
        print(message, end = end, flush = flush)

def get_arguments():
    # Set up arguments (only 2 for input and output)
    parser = argparse.ArgumentParser(description="Impute missing genotypes in a VCF file")
    parser.add_argument("-i", "--input", help="Input VCF file",required=True)
    parser.add_argument("-o", "--output", help="Output VCF file", required=True)
    parser.add_argument("-d", "--debug", help="Debug mode", action="store_true", default=False)
    parser.add_argument("-t", "--threshold", help="Threshold for filtering DR2 values", type=float, default=0.8)

    args = parser.parse_args()
    return args

def check_requirements():
    # Check exit code when running 
    debug_print("---------------------")
    debug_print("Checking requirements")
    debug_print("---------------------")

    # Beagle
    try:
        run(["java", "-jar", BEAGLE], check=True, capture_output=True)
        debug_print("{:<10} {:<5}".format("Beagle", "✅"))
    except CalledProcessError:
        print("Beagle not found. Please check the path")
        debug_print("{:<10} {:<5}".format("Beagle", "❌"))
        exit(1)

    # SnpSift (will exit 1 even if just runninf with help)
    output = run(["java", "-jar", SNPSIFT], capture_output=True)
    if "SnpSift version" not in output.stderr.decode():
        print("SnpSift not found. Please check the path")
        debug_print("{:<10} {:<5}".format("SnpSift", "❌"))
        exit(1)
    else:
        debug_print("{:<10} {:<5}".format("SnpSift", "✅"))

    # Use which for the rest 
    # bcftools
    try:
        run(["which", "bcftools"], check=True, capture_output=True)
        debug_print("{:<10} {:<5}".format("bcftools", "✅"))
    except CalledProcessError:
        print("bcftools not found. Please check the path")
        debug_print("{:<10} {:<5}".format("bcftools", "❌"))
        exit(1)

    # bgzip
    try:
        run(["which", "bgzip"], check=True, capture_output=True)
        debug_print("{:<10} {:<5}".format("bgzip", "✅"))
    except CalledProcessError:
        print("bgzip not found. Please check the path")
        debug_print("{:<10} {:<5}".format("bgzip", "❌"))
        exit(1)
    
    # tabix
    try:
        run(["which", "tabix"], check=True, capture_output=True)
        debug_print("{:<10} {:<5}".format("tabix", "✅"))
    except CalledProcessError:
        print("tabix not found. Please check the path")
        debug_print("{:<10} {:<5}".format("tabix", "❌"))
        exit(1)

def is_there_chr_prefix(input_vcf : str) -> bool:
    debug_print("------------------------------")
    debug_print("Checking for chr prefix in VCF")
    debug_print("------------------------------")

    with open(input_vcf, "r") as f:
        for line in f:
            if line.startswith("##"):
                if line.startswith("##contig=<ID="):
                    if "chr" in line:
                        debug_print("{:<10} {:<5}".format("NO chr Prefix", "❌"))
                        return True
                    else:
                        debug_print("{:<10} {:<5}".format("NO chr Prefix", "✅"))
                        return False
            else:
                break
    exit("No contig line found in VCF")

def convert_to_plain_vcf(input_vcf: str) -> str:

    debug_print("--------------------------------")
    debug_print("Converting INPUT to plain format")
    debug_print("--------------------------------")

    # Save old name to remove tmp file if needed
    if input_vcf.endswith(".vcf.gz"):
        old_vcf_name = input_vcf
        new_vcf_name = input_vcf.replace(".vcf.gz", ".tmp.vcf")
        # Run bgzip
        run(["bgzip", "-d", old_vcf_name, "-o", new_vcf_name ], check=True)
      
        debug_print("Decompressed VCF ✅")

    elif input_vcf.endswith(".bcf"):
        old_vcf_name = input_vcf
        new_vcf_name = input_vcf.replace(".bcf", ".tmp.vcf")
        run(["bcftools", "view", "-O", "v", old_vcf_name, "-o", new_vcf_name])
        debug_print("Converted BCF to VCF ✅")
    elif input_vcf.endswith(".vcf"):
        old_vcf_name = None
        new_vcf_name = input_vcf
        debug_print("VCF already in plain format ✅")
    else :
        exit("❌ Input VCF file is not in a supported format: Must be .vcf, .vcf.gz or .bcf (based on extension) ❌")
    return (new_vcf_name, old_vcf_name)

def split_VCF(input_vcf: str) -> None:

    debug_print("--------------------------")
    debug_print("Split VCF into chromosomes")
    debug_print("--------------------------")

    run(["java", "-jar", SNPSIFT, "split", input_vcf], check=True)

    # Delete chromosomes that are not autosomes 
    # Find all files with the format input_vcf.*.vcf
    # if wildcard is not a number within 1-22, delete it

    # Grab all that fit the format
    files = glob(input_vcf.replace(".vcf", ".*.vcf"))
    for file in files:
        # Get the number
        number = file.split(".")[-2]
        if not number in [str(i) for i in range(1, 23)]:
            os.remove(file)
            debug_print("Deleted non autosomal file {}".format(file))

def impute_single_chromosome(input_vcf: str, chromosome_number: int) -> None:
    filename = input_vcf.replace(".vcf", ".{}.vcf".format(chromosome_number))
    output_filename = filename.replace(".vcf", ".imputed")
    debug_print("Imputing chromosome {}".format(chromosome_number), end="", flush=True)
    try:
        output = run(["java", "-jar", BEAGLE,
                      "ref={}/chr{}.1kg.phase3.v5a.b37.bref3".format(BEAGLE_PANEL_FOLDER, chromosome_number),
                      "map={}/plink.chr{}.GRCh37.map".format(BEAGLE_MAP_FOLDER, chromosome_number),
                      "gt={}".format(filename),
                      "out={}".format(output_filename),
                      ], check=True, capture_output=True)
    except CalledProcessError as e:
        print("Error running beagle")
        print(e.stderr.decode())
        exit(1)
    debug_print(" ✅")


def index_vcf(vcf_file: str, index: int) -> None:
    debug_print("Indexing imputed VCF file for chromosome {}".format(index), end="", flush=True)
    vcf_file = vcf_file.replace(".vcf", ".{}.imputed.vcf.gz".format(index))
    run(["tabix", vcf_file], check=True)
    debug_print(" ✅")


def filter_vcf(vcf_file: str, index: int, threshold: float) -> None:
    debug_print("Filtering VCF file for chromosome {}".format(index), end="", flush=True)
    vcf_file = vcf_file.replace(".vcf", ".{}.imputed.vcf.gz".format(index))
    output_file = vcf_file.replace(".vcf.gz", ".filtered.vcf.gz")
    run(["bcftools", "view", "-i", "DR2>{}".format(threshold), vcf_file, "-o", output_file], check=True)
    debug_print(" ✅")


def impute(input_vcf: str, output_vcf:str,  threshold: float) -> None:
    # Create big string to store all the logs then delete them then make on big log file
    beagle_logs = ""

    # Don't multithread because beagle will use all cores already
    for i in range(1, 23):
        impute_single_chromosome(input_vcf, i)
        index_vcf(input_vcf, i)
        filter_vcf(input_vcf, i, threshold)
        # Cleanup temporary files for the current chromosome
        os.remove(input_vcf.replace(".vcf", ".{}.vcf".format(i)))
        debug_print("Deleted pre-imputation chromosome {} VCF file".format(i))
        os.remove(input_vcf.replace(".vcf", ".{}.imputed.vcf.gz".format(i)))
        os.remove(input_vcf.replace(".vcf", ".{}.imputed.vcf.gz.tbi".format(i))) # Also delete index
        debug_print("Deleted imputed chromosome {} VCF file".format(i))
        # Read log file 
        with open(input_vcf.replace(".vcf", ".{}.imputed.log".format(i)), "r") as f:
            beagle_logs += f.read()
        os.remove(input_vcf.replace(".vcf", ".{}.imputed.log".format(i)))
    
    # Write all logs to one file
    with open(input_vcf.replace(".vcf", ".beagle_imputation.log"), "w") as f:
        f.write(beagle_logs)

    merge_all_files(input_vcf, output_vcf)


def merge_all_files(input_vcf: str, output_vcf: str) -> None:
    debug_print("Merging all chromosomes together ", end="", flush=True)
    all_vcfs = [input_vcf.replace(".vcf", ".{}.imputed.filtered.vcf.gz".format(x)) for x in range(1, 23)]
    # Merge using BCFTools
    run(["bcftools", "concat", *all_vcfs, "-o", output_vcf, "-O", "z"], check=True)


def process_vcf(input_vcf: str, output_vcf: str, threshold: float) -> None:
    # Step 1 Split the VCF by chromosome
    split_VCF(input_vcf)

    # Step 2 Impute each chromosome
    impute(input_vcf, output_vcf, threshold=threshold)

    # Step 4 Merge all the files
    merge_all_files(input_vcf, output_vcf)


def cleanup(input_vcf: str, tmp_vcf_to_delete: str = None) -> None:
    debug_print("---------------------------")
    debug_print("Cleaning up temporary files")
    debug_print("---------------------------")

    # Delete the temporary VCF file if there
    if tmp_vcf_to_delete:
        os.remove(tmp_vcf_to_delete)
        debug_print("Deleted temporary VCF file")

    # Delete filtered files (after merge)
    for i in range(1, 23):
        os.remove(input_vcf.replace(".vcf", ".{}.imputed.filtered.vcf.gz".format(i)))
        debug_print("Deleted filtered chromosome {} VCF file".format(i))


if __name__ == "__main__":
    check_requirements()

    args = get_arguments()

    input_vcf = args.input
    output_vcf = args.output
    threshold = args.threshold

    input_vcf, original_vcf_name = convert_to_plain_vcf(input_vcf)

    if is_there_chr_prefix(input_vcf):
        print("There is a chr prefix in the VCF file")

    print()
    print("===================")
    print("Starting imputation")
    print("===================")
    print()

    process_vcf(input_vcf=input_vcf,
                output_vcf=output_vcf,
                threshold=threshold)

    cleanup(input_vcf, original_vcf_name)

