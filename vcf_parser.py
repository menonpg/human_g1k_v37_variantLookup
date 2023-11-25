import requests
import csv
import pickle
import cyvcf2

def parse_vcf_entry(entry):
    """
    Parses a single VCF entry to extract essential information.

    Args:
        entry (cyvcf2.VCFRecord): A record from a VCF file.

    Returns:
        dict: Essential information extracted from the VCF entry.
    """
    return {
        "chrom": entry.CHROM,
        "pos": entry.POS,
        "ref": entry.REF,
        "alt": entry.ALT[0]  # Assuming a single alternate allele
    }

def query_ensembl(variant):
    """
    Queries the Ensembl VEP endpoint using SPDI notation for a given variant.

    Args:
        variant (dict): A dictionary containing variant information.

    Returns:
        dict: The response from the Ensembl VEP API.
    """
    url = "http://grch37.rest.ensembl.org/vep/human/region"
    headers = {"Content-Type": "application/json"}
    adjusted_position = variant['pos'] + 1  # Adjusting for 1-based indexing
    spdi_notation = f"{variant['chrom']}:{adjusted_position}:{variant['ref']}:{variant['alt']}"

    payload = {"variants": [spdi_notation]}

    try:
        response = requests.post(url, json=payload, headers=headers)
        return response.json() if response.status_code == 200 else {}
    except Exception as e:
        print(f"An error occurred during the API call: {e}")
        return {}

def write_to_csv(data, filename):
    """
    Writes processed data to a CSV file.

    Args:
        data (list of dicts): Processed variant data.
        filename (str): Path to the output CSV file.
    """
    fieldnames = ["chrom", "pos", "ref", "alt", "depth", "alt_reads", "percent_alt_reads", "percent_ref_reads", "gene", "variant_effect", "minor_allele", "minor_allele_frequency", "somatic", "id"]
    with open(filename, 'w', newline='') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        for row in data:
            writer.writerow(row)

""" 
all_data list format:
# all_data = [
#     # ... (your data here) ...
# ]
"""
def extract_csv_fields(entry):
    """
    Extracts necessary fields from a VCF entry for CSV output.

    Args:
        entry (dict): A dictionary containing information about a VCF entry.

    Returns:
        list: A list of values corresponding to the CSV columns.
    """
    # Extracting basic data from the entry
    chrom = entry.get('chrom', '')
    pos = entry.get('pos', '')
    ref = entry.get('ref', '')
    alt = entry.get('alt', '')

    # Extracting additional data if available
    depth = entry.get('INFO').get('DP', '') if 'INFO' in entry else ''
    alt_reads = entry.get('INFO').get('AO', '') if 'INFO' in entry else ''
    percent_alt_reads = str(round(alt_reads / depth * 100, 2)) if depth and alt_reads else ''
    percent_ref_reads = str(round((depth - alt_reads) / depth * 100, 2)) if depth and alt_reads else ''

    gene, variant_effect = '', ''
    if 'transcript_consequences' in entry and entry['transcript_consequences']:
        gene = entry['transcript_consequences'][0].get('gene_symbol', '')
        variant_effect = ', '.join(entry['transcript_consequences'][0].get('consequence_terms', []))

    minor_allele = entry.get('INFO').get('MA', '') if 'INFO' in entry else ''
    minor_allele_frequency = entry.get('INFO').get('MAF', '') if 'INFO' in entry else ''

    id = entry.get('id', '')
    somatic = '1' if 'COSV' in id else ''

    return [chrom, pos, ref, alt, depth, alt_reads, percent_alt_reads, percent_ref_reads, gene, variant_effect, minor_allele, minor_allele_frequency, somatic, id]


def main(vcf_file, output_csv, output_pkl):
    """
    Main function to process VCF file, query Ensembl API, and write to CSV and pickle files.

    Args:
        vcf_file (str): Path to the VCF file.
        output_csv (str): Path to the output CSV file.
        output_pkl (str): Path to the output pickle file.
    """
    vcf_reader = cyvcf2.VCF(vcf_file)
    all_data = []

    for record in vcf_reader:
        vcf_data = parse_vcf_entry(record)
        ensembl_data = query_ensembl(vcf_data)
        combined_data = {**vcf_data, **ensembl_data}
        all_data.append(combined_data)

    # Writing data to a pickle file for serialization
    with open(output_pkl, 'wb') as pkl_file:
        pickle.dump(all_data, pkl_file)

    # Writing data to a CSV file
    # write_to_csv(all_data, output_csv)
        
    with open(output_csv, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(['chrom', 'pos', 'ref', 'alt', 'depth', 'alt_reads', 'percent_alt_reads', 'percent_ref_reads', 'gene', 'variant_effect', 'minor_allele', 'minor_allele_frequency', 'somatic', 'id'])

        for entry in all_data:
            writer.writerow(extract_csv_fields(entry))

    print(f"Data written to {output_csv}")

if __name__ == "__main__":
    main("path/to/vcf_file.vcf", "path/to/output.csv", "path/to/output.pkl")
    
    
