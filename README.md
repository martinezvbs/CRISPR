# CRISPR Data Analysis Pipeline
This pipeline is designed to analyze CRISPR data from Next Generation Sequencing (NGS) experiments. 
It is designed to be run on both high performance computing cluster or personal computer. The pipeline is written in Python and uses Snakemake to manage the workflow. [In progress]

### [The Augert Lab](https://www.augertlab.org/)
#### Yale University, Department of Pathology
* PI: [Dr. Arnaud Augert](mailto:arnaud.augert@yale.edu)
* PhD Student: [Danny Gallant](mailto:danny.gallant@yale.edu)
* Postgraduate: [Juan M. Martinez-Villalobos](mailto:juan.martinezvillalobos@yale.edu) 

### Pre-requisites
* Snakemake 5.10.0 or higher
* Python 3.6 or higher
* R 3.6 or higher

### Installation
1. Clone the repository
```bash
git clone https://github.com/martinezvbs/CRISPR.git
```

### Usage
1. Construct a sample sheet with each line corresponding to a separate barcode.
2. Run the pipeline
```bash
python python3 count_barcodes.py -i CRISPR_library.csv -f ORF_Library_R1_001.fastq -o File.csv -no-g
```

### Output
The pipeline will generate the following files:
* `File.csv` - A CSV file containing 
    * Unique barcode name 
    * Unique barcode sequence
    * Counts of each barcode
    * Gene length of the ORF 
    * RefSeq ID of the ORF
* `statistics_file.txt` - A TXT file containing following statistics:
    * Total number of reads
    * Number of perfect barcode matches:
    * Number of nonperfect barcode matches:
    * Number of reads processed:
    * Percentage of barcodes that matched perfectly:
    * Percentage of undetected barcodes:
    * Skew ratio of top 10% to bottom 10%:
* `CRISPR-scatter.tiff` - A scatter plot of the CRISPR data
* `CRISPR-perfect-matches.csv` - A CSV file containing the perfect matches
* `CRISPR-nonperfect-matches.csv` - A CSV file containing the nonperfect matches

### Contact
For questions or comments, please contact [Juan M. Martinez-Villalobos](mailto:juan.martinezvillalobos@yale.edu)

### Acknowledgements
Part of the code was adapted from
Joung, J., Konermann, S., Gootenberg, J. et al. _Genome-scale CRISPR-Cas9 knockout and transcriptional activation screening_. **Nat Protoc** 12, 828–863 (2017). https://doi.org/10.1038/nprot.2017.016