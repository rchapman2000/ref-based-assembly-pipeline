# Nextflow Reference Based Assembly Pipeline

This pipeline automates the process of creating a reference-based assembly from Illumina NGS Data. It supports the processing of both typical and amplicon sequencing approaches (depending on the options provided).

## Technical Considerations
Some aspects of this analysis are important to consider and understand when interpreting the output.

### Variant Calling
This pipeline calls majority variants. Thus, the variant allele must be applied in the final consensus so long as it is more prevalent than the reference allele.

As well, the pipeline uses an in-house script to filter multiallelic sites from the variant call files. A multiallelic site is one where multiple alternative alleles are present. See the "CCA,GCA,ACA,ACG,ATA" in the following example.
```
K03455.1        2911    .       ATG     CCA,GCA,ACA,ACG,ATA     10482.7 .       AB=0,0,0,0,0;ABP=0,0,0,0,0;AC=0,0,2,0,0;AF=0,0,1,0,0;AN=2;AO=1,2,310,1,2;CIGAR=3X,3X,1M2X,1M1X1M,2M1X;DP=316;DPB=316;DPRA=0,0,0,0,0;EPP=5.18177,3.0103,5.8122,5.18177,3.0103;EPPR=0;GTI=0;LEN=3,3,2,1,1;MEANALT=5,5,5,5,5;MQM=60,60,59.9032,60,60;MQMR=0;NS=1;NUMALT=5;ODDS=413.347;PAIRED=1,1,1,1,1;PAIREDR=0;PAO=0,0,0,0,0;PQA=0,0,0,0,0;PQR=0;PRO=0;QA=27,76,11724,39,72;QR=0;RO=0;RPL=1,2,165,1,2;RPP=5.18177,7.35324,5.8122,5.18177,7.35324;RPPR=0;RPR=0,0,145,0,0;RUN=1,1,1,1,1;SAF=0,1,178,1,1;SAP=5.18177,3.0103,17.8324,5.18177,3.0103;SAR=1,1,132,0,1;SRF=0;SRP=0;SRR=0;TYPE=complex,complex,mnp,snp,snp      GT:DP:AD:RO:QR:AO:QA:GL 3/3:316:0,1,2,310,1,2:0:0:1,2,310,1,2:27,76,11724,39,72:-1054.34,-1052.22,-1051.91,-1048.11,-1045.98,-1047.51,-93.3193,-91.2078,-87.0746,0,-1051.14,-1049.01,-1044.9,-90.1068,-1050.84,-1048.47,-1046.34,-1042.24,-87.4429,-1045.26,-1047.87
```

This script uses the AO field (Alternative allele appearances) to identify the **most prevalent alternative allele** (ACA in the previous example has an AO value of 310 and is thus the most prevalent). Also note that the **genotype field is removed** by this script as these values are not important for downstream processing.
## Installation

To install this pipeline, enter the following commands:
```
# Clone the repository
git clone https://github.com/rchapman2000/ref-based-assembly-pipeline.git

# Create a conda environment using the provided environment.yml file
conda env create -f environment.yml

# Activate the conda environment
conda activate RefBasedAssembly
```
### Updating the Pipeline
If you already have the pipeline installed, you can update it using the following commands:
```
# Navigate to your installation directory
cd ref-based-assembly-pipeline

# Use git pull to get the latest update
git pull

# Activate the conda environment and use the environment.yml file to download updates
conda activate RefBasedAssembly
conda env update --file environment.yml --prune
```

## Usage
To run the pipeline, use the following command:
```
# You must either be in the same directory as the main.nf file or reference the file location.
nextflow run main.nf [options] --input INPUT_DIR --output OUTPUT_DIR --reference REFERENCE_FASTA
```

### Options
The pipeline also supports the following optional arguments:

| Option | Type | Description |
|---|---|---|
| --xgen | *File* | This pipeline supports the IDT xGen amplicon panels. If one of their 'masterfiles' (provided with the technology) is supplied to this option, the pipeline will perform primer clipping using the PrimerClip software. |
| --host_reference | *FASTA File* | If provided a host fasta file, the pipeline will perform host read removal via alignment (Supplying this option will require a bowtie2 index to be built). |
| --host_bt2_index | *Directory* | Alternative to providing a fasta file, a pre-built bowtie2 index can be provided for host removal. The index must be in its own directory. |
| --noPicard | *None* | Supplying this option disables the use of picard (which is used by default). Supplying a primer file will do this automatically. [Default = use picard unless primers supplied] |
| --minCov | *int* | The minimum depth of coverage, below which a a position will be masked. [Default = 20] |
| --minLen | *int* | The minimum length of a read (in base pairs) to keep post trimming. [Default = 75] |
| --minBQ | *int* |  The minimum base call quality for a site to be considered in variant calling and depth-masking [Default = 10] |
| --minMapQ | *int* | The minimum mapping quality for a site to be considered in variant calling and depth-masking [Default = 0] |
| --threads | *int* | The number of threads that can be use to run pipeline tools in parallel. [Default = 1] |

To view the list of options from the commandline, use the following command:
```
nextflow rn main.nf --help
```

