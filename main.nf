#!/usr/bin/env/ nextflow

nextflow.enable.dsl=2

def helpMessage() {
    log.info"""
Reference Based Assembly Pipeline

Takes an input of paired-end fastq files from illumina
sequencers, processes the reads, aligns to a reference, performs indel
realignment, calls variants, and then generates a consensus genome sequence.
As well, the pipeline steps are variable depending on the specified
parameters to meet difference analysis needs. These include a primer
clipping step for amplicon sequencing, and a host read removal step.

USAGE: nextflow run main.nf [options] --input INPUT_DIR --output OUTPUT_DIR --ref REFERENCE_FASTA --abraLoc PATH_TO_ABRA2_JAR [--swift PRIMER_MASTER_FILE | --picardLoc PATH_TO_PICARD_JAR]
OPTIONS:

--input INPUT_DIR - [Required] A directory containing paired-end fastq files

--output OUTPUT_DIR - [Required] A directory to place output files (If not existing, pipeline will create)

--reference REFERENCE_FASTA - [Required] A reference genome to align reads to.

--abraLoc PATH_TO_ABRA2_JAR - [Required] The path to the abra2 jar file used for realignment.

ONE OF THE FOLLOWING IS REQUIRED:
    --swift PRIMER_MASTER_FILE - tell the pipeline to perform primer clipping based on the supplied masterfile
    
    --picardLoc PATH_TO_PICARD_JAR  - The path to the picard jar file to be used for read deduplication

OPTIONAL:
    --host_reference HOST_REF_FASTA - a fasta file containing a host reference sequence. Supplying this option will require a bowtie2 index to be built

    --host_bt2_index INDEX_DIRECTORY - To save time, an existing bowtie2 index can be supplied. Must be in its own directory

    --minCov INT - The minimum coverage below which a position will be masked [Default = 20]
    
    --threads INT - the number of threads that can be use to run pipeline tools in parallel
    """
}

// Function that checks whether a directory ends in a trailing slash.
// This is useful for directory variables that are not parsed into
// file objects in the pipeline (such as the output directory).
def checkDirectoryEnding (fileName) {
    // Grabs the last character in the directory name.
    lastChar = fileName.substring(fileName.length() - 1, fileName.length())

    // Checks whether the last character is slash
    if (lastChar != "/") {
        // If it is not a slash, add that to the directory name.
        fileName = fileName + "/"
    }
    
    // Return the directory name.
    return fileName
}

// If the help parameter is supplied, link display the help message
// and quit the pipeline
params.help = false
if (params.help){
    helpMessage()
    exit 0
}

// Defines input parameters. Setting to false by default
// allows us to check that these have been set by the user.
params.input = false
params.reference = false
params.output = false
params.swift = false
params.threads = 1
params.picardLoc = false
params.abraLoc = false
params.minCov = 20
params.host_bt2_index = false
params.host_reference = false

// Imports the adapters file present in the installation directory.
adapters = file("${baseDir}/adapters.fa")

println "Input Directory: ${params.input}"

// Inports modules
include { IndexReference } from './modules.nf'
// The same module cannot be used more than once,
// thus it is aliased to be used multiple times.
include { IndexReference as IndexHostReference } from './modules.nf'
include { QCReport } from './modules.nf'
include { QCReport as QCReport_Trimmed } from './modules.nf'
include { Trimming } from './modules.nf'
include { HostReadRemoval } from './modules.nf'
include { Bowtie2Alignment } from './modules.nf'
include { SwiftPrimerClip } from './modules.nf'
include { MarkDuplicates } from './modules.nf'
include { Realignment } from './modules.nf'
include { CallVariants } from './modules.nf'
include { GenerateConsensus } from "./modules.nf"

// Checks the input parameter
if (params.input == false) {
    // If the parameter is not set, notify the user and exit.
    println "ERROR: No input directory provided. Pipeline requires an input directory."
    exit(1)
}
else if (!(file(params.input).isDirectory())) {
    // If the input directory is not set, notify the user and exit.
    println "ERROR: ${params.input} is not an existing directory."
    exit(1)
}

// Create a channel for hte input files.
inputFiles_ch = Channel
    // Pull from pairs of files (illumina fastq files denoted by having R1 or R2 in
    // the file name).
    .fromFilePairs("${params.input}*_R{1,2}*.fastq.gz")
    // The .fromFilePairs() function spits out a list where the first 
    // item is the base file name, and the second is a list of the files.
    // This command creates a tuple with the base file name and two files.
    .map { it -> [it[0], it[1][0], it[1][1]]}


// Checks the output parameter.
outDir = ''
if (params.output == false) {
    // If the parameter is not set, notify the user and exit.
    println "ERROR: No output directory provided. Pipeline requires an output directory."
}
else {
    // If the parameter is set, ensure that the directory provided ends
    // in a trailing slash (to keep things consistent throughout) the
    // pipeline code.
    outDir = checkDirectoryEnding(params.output)
}

// Checks the reference parameter. For this, we cannot use an
// input channel like was used for the input files. Using an input channel
// will cause Nextflow to only iterate once as the reference 
// channel would only only have 1 file in it. Thus, we manually parse
// the reference file into a tuple.
refData = ''
if (params.reference == false) {
    // If the parameter is not set, notify the user and exit.
    println "ERROR: no reference file proivded. Pipeline requires a reference file."
    exit(1)
}
else if (!(file(params.reference).exists())) {
    // If the reference file provided does not exist, notify the user and exit.
    println "ERROR: ${params.reference} does not exist."
    exit(1)
}
else {
    // Process the reference file to be supplied to the index step.
    
    // Parse the file provided into a file object.
    ref = file(params.reference)

    // Grab the basename of the file.
    refName = ref.getBaseName()

    // Place the file basename and file object into
    // a tuple.
    refData = tuple(refName, ref)
}


// Process the swift primers parameter.
primerfile = ''
if (params.swift != false) {
    // If the parameter is provided, check if the file exists.
    if (file(params.swift).isFile())
    {
        // If the file does exist, parse it into
        // a file object.
        primerfile = file(params.swift)
    }
    else {
        // If the file does not exist, notify the user and exit.
        println "ERROR: ${params.swift} does not exist."
        exit(1)
    }
}

// Parses the host options (--host_reference and --host_bt2_index).
// Again, we cannot use a channel like we did for the input files
// because then Nextfow would only run other modules once.
// Thus, we need to manually create a tuple of input data to pass to the indexing
// step and the alignment step.
hostRefData = ''
hostRefIdxData = ''
if (params.host_reference != false && params.host_bt2_index != false) {
    // If both options are supplied, notify the user and exit.
    println "ERROR: you have specified both a host fasta file and bowtie2 index. Please only supply one."
    exit(1)
}
else {
    // If the user supplied the --host_reference option
    if (params.host_reference != false) {
        if (!(file(params.host_reference).exists())) {
            // If the file supplied does not exist, notify the user and exit.
            println "ERROR: ${params.host_reference} does not exist."
            exit(1)
        }
        else {
            // Parse the file into a file object
            hostRef = file(params.host_reference)
            // Use the getBaseName() function to 
            // get the reference name. This will be
            // used to name the bowtie2 index.
            hostRefName = hostRef.getBaseName()
            // Place these both into a tuple.
            hostRefData = tuple(hostRefName, hostRef)
        }
    }
    // If the user supplied the --host_bt2_index
    else if (params.host_bt2_index != false) {
        if (!(file(params.host_bt2_index).exists())) {
            // If the index provided does not exist, notify the user and exit.
            println "Error: ${params.host_bt2_idx} does not exist."
            exit(1)
        }
        else {
            // Parse the directory into a file object
            hostRefDir = file(params.host_bt2_index)
            println hostRefDir
            // Grab a list of file objects from the directory
            // ending in .bt2
            indexFiles = file("${hostRefDir}/*.bt2")
            if (indexFiles.size() == 0){
                // If there are no file in the directory ending in bt2, notify the user and exit
                println "Index Directory provided (${params.host_bt2_index}) does not contain any bt2 files"
                exit(1)
            }
            else {
                // Use the getSimpleName() function to grab the base name
                // of the index files (getSimpleName() removes anything following
                // the last . in a file name.)
                hostRefName = indexFiles[0].getSimpleName()
                println hostRefName
                // Place the index dir and name into a tuple.
                hostRefIdxData = tuple(hostRefDir, hostRefName)
            }
        }
    }
}

workflow {

    // Index the reference file provided.
    IndexReference( refData, params.output, params.threads )
    
    // Use FASTQC to perform an initial QC check on the reads
    QCReport( inputFiles_ch, params.output, "FASTQC-Pre-Processing", params.threads )

    // Perform adapter and quality trimming with trimmomatic.
    Trimming( inputFiles_ch, params.output, adapters )

    // Use FASTQC to perform a QC check on the trimmed reads
    QCReport_Trimmed( Trimming.out[0], params.output, "FASTQC-Trimmed", params.threads )
    
    // If the user supplied a host bowtie2 index
    if (params.host_bt2_index != false) {
        
        // Perform host removal by aligning to the host reference using bowtie2
        HostReadRemoval( Trimming.out[0], params.output, hostRefIdxData, params.threads)

        // Align the non-host reads to the reference using bowtie2
        Bowtie2Alignment( HostReadRemoval.out[0], params.output, IndexReference.out, params.threads )
    }
    // If the user supplied a host fasta file
    else if (params.host_reference != false) {

        // Index the host reference file.
        IndexHostReference( hostRefData, params.output, params.threads )
        
        // Perform host removal by aligning to the host reference using bowtie2
        HostReadRemoval( Trimming.out[0], params.output, IndexHostReference.out, params.threads )
        
        // Align the non-host reads to the reference using bowtie2
        Bowtie2Alignment( HostReadRemoval.out[0], params.output, IndexReference.out, params.threads )
    }
    else {
        // Align the reads to the reference using bowtie2
        Bowtie2Alignment( Trimming.out[0], params.output, IndexReference.out, params.threads )
    }

    // If the user supplied the swift option
    if (params.swift) {
        // Perform realignment to improve indel quality
        Realignment( Bowtie2Alignment.out[0], params.output, params.abraLoc, refData, params.threads ) 
        
        // Perform primer clipping using primerclip
        SwiftPrimerClip( Realignment.out, primerfile, params.output, params.threads )

        // Call and filter variants
        CallVariants( SwiftPrimerClip.out, baseDir, params.output, refData, params.minCov )
    }
    // The user did not supply the swift option
    else {
        // Mark duplicates using picard
        MarkDuplicates( Bowtie2Alignment.out[0], params.output, params.picardLoc )

        // Perform Realignment to improve indel quality
        Realignment( MarkDuplicates.out, params.output, params.abraLoc, refData, params.threads ) 

        // Call and filter variants
        CallVariants( Realignment.out, baseDir, params.output, refData, params.minCov )
    }

    // Generate a consensus from the alignment and variant calling.
    GenerateConsensus( CallVariants.out[0], baseDir, params.output, refData, params.minCov )
}
