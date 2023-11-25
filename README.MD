# VCF Variant Processor

This Python project processes VCF (Variant Call Format) files to extract variant information, queries the Ensembl VEP (Variant Effect Predictor) API, and saves the processed data in both CSV and pickle formats.

## Files

- `vcf_parser.py`: Main script for processing VCF files and querying the Ensembl VEP API.
- `test_vcf_parser.py`: Unit tests for the functions in `vcf_parser.py`.

## Functions

- `parse_vcf_entry(entry)`: Parses a single entry from a VCF file.
- `query_ensembl(variant)`: Queries the Ensembl VEP API for variant information using SPDI notation.
- `write_to_csv(data, filename)`: Writes the processed variant data to a CSV file.
- `main(vcf_file, output_csv, output_pkl)`: Main function that orchestrates the VCF processing and data output.

## Usage

Run `vcf_parser.py` with the path to a VCF file, and paths for the output CSV and pickle files:

```bash
python vcf_parser.py /path/to/vcf_file.vcf /path/to/output.csv /path/to/output.pkl
```

## Testing
Run test_vcf_parser.py to perform unit tests:
```bash
python test_vcf_parser.py
```

## Dependencies
- cyvcf2: For parsing VCF files.
- requests: For making API requests.
- csv: For CSV file operations.
- pickle: For data serialization.
