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

USAGE: nextflow run main.nf [options] --input INPUT_DIR --output OUTPUT_DIR --ref REFERENCE_FASTA
OPTIONS:

--input INPUT_DIR - [Required] A directory containing paired-end fastq files

--output OUTPUT_DIR - [Required] A directory to place output files (If not existing, pipeline will create)

--reference REFERENCE_FASTA - [Required] A reference genome to align reads to.


OPTIONAL:
    
    --primers PRIMER_BED_FILE - Supply a .bed file to perform primer clipping using iVAR (Cannot be used with --xgen option)

    --xgen PRIMER_MASTER_FILE - Supply a masterfile to perform primer clipping using PrimerClip (for IDT xGen Assays) (Cannot be used with --primers option)

    --host_reference HOST_REF_FASTA - a fasta file containing a host reference sequence. Supplying this option will require a bowtie2 index to be built

    --host_bt2_index INDEX_DIRECTORY - To save time, an existing bowtie2 index can be supplied. Must be in its own directory

    --localAlignment - Supplying this option enables bowtie2's local alignment mode (will be used for host removal if that option is supplied as well) [Default = end-to-end mode]

    --noPicard - Supplying this option disables the use of picard (which is used by default). Supplying a primer file will do this automatically. [Default = use picard unless primers supplied]

    --minCov INT - The minimum coverage below which a position will be masked [Default = 20]

    --minBQ INT - The minimum base call quality for a site to be considered in variant calling and depth-masking [Default = 10]

    --minMapQ INT - The minimum mapping quality for a site to be considered in variant calling and depth-masking [Default = 0]

    --minLen INT - the minimum length of a read to keep post trimming [Default = 75bp]

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

// Function creates the header for the summary file based on the parameters
// supplied by the user. Because the pipeline dynamically changes based on
// what the user specifies, the summary file must also be alter to reflect
// the analysis. This will also make incorporating new modules easier.
def createSummaryHeader (hostRef, hostIdx, primers, xgen, noPicard) {
    
    // The header will always start with the sample.
    FinalHeader = 'Sample,'

    // The summary sheet will also always contain the Raw and trimmed read counts.
    FinalHeader = FinalHeader + "Raw Reads,Trimmed Reads,"

    // Next, the user may supply host removal, which would be the next useful
    // statistic to know. Thus, if a host reference or bowtie2 index is supplied,
    // add this field to the header.
    if (hostRef != false || hostIdx != false) {
        FinalHeader = FinalHeader + "Non-Host Reads,"
    }

    // Next, the pipeline will also perform alignment. Thus, this the next
    // useful statistic.
    FinalHeader = FinalHeader + "Mapped Reads,"

    // After alignment, the user may supply a set of primers to be clipped. Thus,
    // if the primers or xgen parameters have been supplied, a field for the number of
    // clipped mapped reads will be added.
    if (primers != false || xgen != false) {
        FinalHeader = FinalHeader + "Clipped Mapped Reads,"
    }
    // If the user has supplied the noPicard option, then do not add anything to the summary header.
    else if (noPicard != false) {
        FinalHeader = FinalHeader
    }
    // If not primer clipping has been done, the pipeline will instead perform deduplication.
    // Thus, a field for the number of deduplicated reads will be added to the
    // summary file.
    else {
        FinalHeader = FinalHeader + "Deduped Mapped Reads,"
    }

    // Finally, the pipeline will always report the number of SNPs, indels, masking positions,
    // and coverage.
    FinalHeader = FinalHeader + "Average Read Depth,SNPs,Indels,Masked Positions,Coverage"

    return FinalHeader
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
params.xgen = false
params.primers = false
params.threads = 1
params.localAlignment = false
params.minCov = 20
params.minLen = 75
params.trimQualThreshold = 20
params.minBQ = 10
params.minMapQ = 0
params.noPicard = false
params.host_bt2_index = false
params.host_reference = false

// Imports the adapters file present in the installation directory.
adapters = file("${baseDir}/data/adapters.fa")

// Inports modules
include { Setup } from './modules.nf'
include { IndexReference } from './modules.nf'
// The same module cannot be used more than once,
// thus it is aliased to be used multiple times.
include { IndexReference as IndexHostReference } from './modules.nf'
include { QCReport } from './modules.nf'
include { QCReport as QCReport_Trimmed } from './modules.nf'
include { Trimming } from './modules.nf'
include { HostReadRemoval } from './modules.nf'
include { Bowtie2Alignment } from './modules.nf'
include { XGenPrimerClip } from './modules.nf'
//include { BedPrimerClip } from './modules.nf'
include { MarkDuplicates } from './modules.nf'
include { Realignment } from './modules.nf'
include { CallVariants } from './modules.nf'
include { GenerateConsensus } from "./modules.nf"
include { WriteSummary } from "./modules.nf"

// Checks the input parameter
inDir = ''
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
else {
    // If the parameter is set, convert the value provided to a file type
    // to get the absolute path, and then convert back to a string to be
    // used in the pipeline.
    inDir = file(params.input).toString()
    println "Input Directory: ${inDir}"
}

// Create a channel for hte input files.
inputFiles_ch = Channel
    // Pull from pairs of files (illumina fastq files denoted by having R1 or R2 in
    // the file name).
    .fromFilePairs("${inDir}/*_R{1,2}*.fastq*")
    // The .fromFilePairs() function spits out a list where the first 
    // item is the base file name, and the second is a list of the files.
    // This command creates a tuple with the base file name and two files.
    .map { it -> [it[0], it[1][0], it[1][1]]}


// Checks the output parameter.
outDir = ''
if (params.output == false) {
    // If the parameter is not set, notify the user and exit.
    println "ERROR: No output directory provided. Pipeline requires an output directory."
    exit(1)
}
else {
    // If the parameter is set, convert the value provided to a file type
    // to get the absolute path, and then convert back to a string to be
    // used in the pipeline.
    outDir = file(params.output).toString()
    println(outDir)
}
// Checks the reference parameter. For this, we cannot use an
// input channel like was used for the input files. Using an input channel
// will cause Nextflow to only iterate once as the reference 
// channel would only only have 1 file in it. Thus, we manually parse
// the reference file into a tuple.
refData = ''
refName = ''
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


// Process the --xgen, --primers, and/or --noPicard parameters.
primerfile = ''
primerFileName = 'NONE'
picardUsed = 'True'
if (params.xgen != false && params.primers != false) {
    // Both the --xgen and --primers options cannot be provided. If they are,
    // notify the user and exit.
    println "ERROR: Both --primers and --xgen were supplied. Only one can be supplied per run. Please adjust the parameters."
    exit(1)
}
else if ((params.xgen != false && params.noPicard != false) || (params.primers != false && params.noPicard != false)) {
    // If a primer file is supplied to either --xgen or --primers
    // then the --noPicard option is not needed. Notify the user
    // and exit.
    println "ERROR: When the --xgen or --primers options are supplied. The --noPicard option is not needed. Please remove this option."
    exit(1)
}
else if (params.xgen != false) {
    // If the --xgen parameter is provided, check if the file exists.
    if (file(params.xgen).isFile()) {
        // If the file does exist, parse it into
        // a file object.
        primerfile = file(params.xgen)
        primerFileName = primerfile.getName()
        
        // Set the picardUsed variable to false (for analysis parameters file)
        picardUsed = 'False'
    }
    else {
        // If the file does not exist, notify the user and exit.
        println "ERROR: ${params.xgen} does not exist."
        exit(1)
    }
}
else if (params.primers != false) {
    // If the --primers paramter is provided, check if the file exists.
    if (file(params.primers).isFile()) {
        println "UNDER DEVELOPMENT - Exiting"
        exit(0)
        // If the file does exist, parse it into
        // a file object.
        primerfile = file(params.primers)
        primerFileName = primerfile.getName()

        // Set the picardUsed variable to false (for analysis parameters file)
        picardUsed = 'False'

        // iVar requires a bed file containing primers and coordiantes. Check to ensure
        // that the provided file ends in .bed, and if not notify the user and exit.
        if (primerfile.getExtension() != 'bed') {
            println "ERROR: The primer file ${params.primers} is not a .bed file. iVar requires a bed file containing primer sequence and positions. Please supply a .bed file."
            exit(1)
        }
    }
    else {
        // If the file does not exist, notify the user and exit.
        println "ERROR: ${params.primers} does not exist."
        exit(1)
    }
}
else if (params.noPicard != false) {
    // If the noPicard flag was provided, set the picardUsed variable to false
    // (for analysis parameters file)
    picardUsed = 'False'

}

// Parses the --localAlignment Parameter
alignmentParam = "--end-to-end"
alignmentModeSummary = "End-to-End Alignment Mode"
if (params.localAlignment != false) {
    alignmentParam = "--local"
    alignmentModeSummary = "Local Alignment Mode"
}


// Parses the host options (--host_reference and --host_bt2_index).
// Again, we cannot use a channel like we did for the input files
// because then Nextfow would only run other modules once.
// Thus, we need to manually create a tuple of input data to pass to the indexing
// step and the alignment step.
hostRefData = ''
hostRefIdxData = ''
hostRefName = 'NONE'
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
            println "Error: ${params.host_bt2_index} does not exist."
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

// Creates the summary header based on the options provided.
summaryHeader = createSummaryHeader(params.host_reference, params.host_bt2_index, params.primers, params.xgen, params.noPicard)

// To keep track of the summary, a string will be passed between modules that report
// metrics. At the end of the module, the metric will be added to the string and
// returned as output to the next module.
workflow {

    // Creates a parameters and summary file.
    Setup( refName, params.minLen, params.trimQualThreshold, params.minCov, params.minBQ, params.minMapQ, picardUsed, alignmentModeSummary, primerFileName, hostRefName, summaryHeader, outDir )

    // Index the reference file provided.
    IndexReference( refData, outDir, params.threads )
    
    // Use FASTQC to perform an initial QC check on the reads
    QCReport( inputFiles_ch, outDir, "FASTQC-Pre-Processing", params.threads )

    // Perform adapter and quality trimming with trimmomatic.
    Trimming( inputFiles_ch, outDir, adapters, params.minLen, params.trimQualThreshold )

    // Use FASTQC to perform a QC check on the trimmed reads
    QCReport_Trimmed( Trimming.out[0], outDir, "FASTQC-Trimmed", params.threads )
    
    // If the user supplied a host bowtie2 index
    if (params.host_bt2_index != false) {
        
        // Perform host removal by aligning to the host reference using bowtie2
        HostReadRemoval( Trimming.out[0], outDir, hostRefIdxData, alignmentParam, params.threads, Trimming.out[2])

        // Align the non-host reads to the reference using bowtie2
        Bowtie2Alignment( HostReadRemoval.out[0], outDir, IndexReference.out, alignmentParam, params.threads, HostReadRemoval.out[2] )
    }
    // If the user supplied a host fasta file
    else if (params.host_reference != false) {

        // Index the host reference file.
        IndexHostReference( hostRefData, outDir, params.threads )
        
        // Perform host removal by aligning to the host reference using bowtie2
        HostReadRemoval( Trimming.out[0], outDir, IndexHostReference.out, alignmentParam, params.threads, Trimming.out[2] )
        
        // Align the non-host reads to the reference using bowtie2
        Bowtie2Alignment( HostReadRemoval.out[0], outDir, IndexReference.out, alignmentParam, params.threads, HostReadRemoval.out[2])
    }
    else {
        // Align the reads to the reference using bowtie2
        Bowtie2Alignment( Trimming.out[0], outDir, IndexReference.out, alignmentParam, params.threads, Trimming.out[2] )
    }

    // If the user supplied the xgen option
    if (params.xgen) {
        // Perform realignment to improve indel quality
        //Realignment( Bowtie2Alignment.out[0], outDir, refData, params.threads, Bowtie2Alignment.out[2])
        
        // Perform primer clipping using primerclip
        XGenPrimerClip( Bowtie2Alignment.out[0], primerfile, outDir, params.threads, Bowtie2Alignment.out[1] )

        // Call and filter variants
        CallVariants( XGenPrimerClip.out[0], baseDir, outDir, refData, params.minCov, params.minBQ, params.minMapQ, XGenPrimerClip.out[1] )
    }
    /*

    UNDER DEVELOPMENT

    else if (params.primers) {
        // Perform realignment to improve indel quality
        Realignment( Bowtie2Alignment.out[0], outDir, refData, params.threads ) 
        
        // Perform primer clipping (UNDER DEVELOPMENT)
        BedPrimerClip( Realignment.out, primerfile, outDir, params.threads )

        // Call and filter variants
        CallVariants( IVarPrimerClip.out, baseDir, outDir, refData, params.minCov )
    }
    */
    // The user did not supply the xgen option
    else if (params.noPicard) {
        // Perform Realignment to improve indel quality
        //Realignment( Bowtie2Alignment.out[0], outDir, refData, params.threads, Bowtie2Alignment.out[2])

        // Call and filter variants
        CallVariants( Bowtie2Alignment.out[0], baseDir, outDir, refData, params.minCov, params.minBQ, params.minMapQ, Bowtie2Alignment.out[1] )
    }
    else {
        // Mark duplicates using picard
        MarkDuplicates( Bowtie2Alignment.out[0], outDir, Bowtie2Alignment.out[1] )

        // Perform Realignment to improve indel quality
        //Realignment( MarkDuplicates.out[0], outDir, refData, params.threads, MarkDuplicates.out[1] ) 

        // Call and filter variants
        CallVariants( MarkDuplicates.out[0], baseDir, outDir, refData, params.minCov, params.minBQ, params.minMapQ, MarkDuplicates.out[1] )
    }

    // Generate a consensus from the alignment and variant calling.
    GenerateConsensus( CallVariants.out[0], baseDir, outDir, refData, params.minCov, params.minBQ, params.minMapQ, CallVariants.out[2] )

    WriteSummary(GenerateConsensus.out[1], outDir)
}
