// Builds a bowtie2 index for a provided reference file
process IndexReference {
    input:
        // Tuple contains the reference name and reference file.
        tuple val(refName), file(ref)
        // The name of the output directory
        val outDir
        // The number of threads provided
        val threads

    output:
        // Tuple containing the reference index directory and the index file name.
        tuple file("${refName}-idx/"), val(refName)

    publishDir "${outDir}", mode: 'copy'

    /*
    Creates an index Directory and then runs bowtie2-build to create the index
    */
    script:
    """
    #!/bin/bash
    
    mkdir ${refName}-idx/

    bowtie2-build --threads ${threads} ${ref} ${refName}-idx/${refName}
    """
}


// Creates a fastqc report for a set of paired-end reads provided.
process QCReport {
    input:
        // Tuple contains the file basename as well as the paired-end read files
        tuple val(base), file(F1), file(F2)
        // The name of the output directory
        val outDir
        // The name of the directory to place the fastqc output into (allows
        // for this command to be used multiple times and to separate the output
        // i.e. pre-processed-reads vs trimmed-reads.)
        val dirName
        // The number of threads provided.
        val threads

    output:
        // The directory cotaining the fastqc files produced.
        file("${base}")
    
    publishDir "${outDir}${dirName}/", mode: 'copy'

    // Creates a Directory and then runs fastqc on a set of reads.
    script:
    """
    #!/bin/bash

    mkdir ${base}

    fastqc ${F1} ${F2} -o ${base} --threads ${threads}
    """
}


// Performs quality and adapter trimming on a set of paired-end reads.
process Trimming {
    input: 
        // Tuple cotains the file basename as well as the paired-end read files.
        tuple val(base), file(R1), file(R2)
        // The name of the output directory
        val outDir
        // The adapter file in fasta format
        file adapters

    output:
        // Tuple containing the file basename and the trimmed forward and reverse reads
        tuple val(base), file("${base}_1.trimmed.fq.gz"), file("${base}_2.trimmed.fq.gz")
        // Tuple containing the unpaired read files.
        tuple file("${base}_1.unpaired.fq.gz"), file("${base}_2.unpaired.fq.gz")

    publishDir "${outDir}${base}-Processed-Reads", mode: 'copy'

    script:
    /*
    Uses trimmomatic to trim the read files.
    Trimmomatic performs trimming steps in the order provided in the 
    command line. Our steps:
    1. ILLUMINACLIP: Removes illumina adapters
    2. LEADING: Begins at the start of the read and trims bases with quality less than
        5 until it hits a base above that threshold.
    3. TRAILING: Begins at the end of the read and trims bases with quality less than 5
        until it hits a base above that threshold.
    4. SLIDINGWINDOW: Begins at the 5' end and scans in 4 base windows, trimming when it hits
        an average quality less than 20.
    5. MINLEN: Removes reads less than 75 bp long.

    Prinseq (used in the following process) cannot take gzipped reads as input. Thus,
    only the unpaired reads are gzipped to save space.
    */
    """
    #!/bin/bash

    trimmomatic PE ${R1} ${R2} ${base}_1.trimmed.fq ${base}_1.unpaired.fq \
    ${base}_2.trimmed.fq ${base}_2.unpaired.fq ILLUMINACLIP:${adapters}:2:30:10:1:true \
    LEADING:5 TRAILING:5 SLIDINGWINDOW:4:20 MINLEN:75
    
    gzip ${base}_1.trimmed.fq ${base}_1.unpaired.fq ${base}_2.trimmed.fq ${base}_2.unpaired.fq
    """
        
}

// Aligns reads to a reference gneome using bowtie2
process Bowtie2Alignment {
    input:
        // Tuple contains the file basename and paired-end reads
        tuple val(base), file(R1), file(R2)
        // The output directory name
        val outDir
        // Tuple contains the bowtie2 index directory and the name of the reference used
        tuple file(refDir), val(refName)
        // The number of threads provided.
        val threads

    output:
        // Tuple contains the file basename and the alignment in a sorted bam file
        tuple val(base), file("${base}.bam")
        // A directory containing the alignment in sam format.
        file "${base}-align.sam"
    
    publishDir "${outDir}", mode: 'copy'

    script:
    /*
    Aligns the reads to the reference genome using bowtie2. ocal alignment is used (--local)
    to ensure the reads are aligne without gaps.

    The alignment is then converted to bam format and sorted using samtools.
    */
    """
    #!/bin/bash

    bowtie2 --threads ${threads} -x ${refDir}/${refName} -1 ${R1} -2 ${R2} --local -S ${base}-align.sam

    samtools view -b ${base}-align.sam | samtools sort > ${base}.bam
    """
}

// Removes host reads by aligning to a host genome using bowtie2.
process HostReadRemoval {
    input:
        // Tuple contains the file basename and paired-end reads
        tuple val(base), file(R1), file(R2)
        // The output directory name
        val outDir
        // Tuple contains the bt2 index directory and basename of the index files.
        tuple file(refDir), val(refName)
        // The number of threads provided.
        val threads

    output:
        // Tuple contains the file basename and the paired-end read files with host reads removed.
        tuple val(base), file("${base}_host_removed_1.fq.gz"), file("${base}_host_removed_2.fq.gz")
        // A directory containing the alignment to the host file.
        file "host-reads/${base}-host.sam"
    
    publishDir "${outDir}${base}-Processed-Reads", mode: 'copy'

    script:
    /*
    Aligns the reads to a host reference genome using bowtie2. Local alignment is used (--local)
    to ensure the reads are aligne without gaps. The unaligned reads are sent back ot into a
    fastq files (--un-conc).

    The unaligned read files are then renamed to give them the typical paired end
    read name scheme.

    Finally, the reads are gzipped to preserve space.
    */
    """
    #!/bin/bash
    mkdir host-reads
    bowtie2 --threads ${threads} -x ${refDir}/${refName} -1 ${R1} -2 ${R2} --local -S host-reads/${base}-host.sam --un-conc ${base}_host_removed

    mv ${base}_host_removed.1 ${base}_host_removed_1.fq
    mv ${base}_host_removed.2 ${base}_host_removed_2.fq

    gzip ${base}_host_removed_1.fq ${base}_host_removed_2.fq
    """
}


// Soft clips primers from an alignment using primerclip
process SwiftPrimerClip {
    input:
        // Tuple contains the file basename and alignment bam file
        tuple val(base), file(bam)
        // The primer masterfile for input into primerclip.
        file primersFile
        // The name of the output directory
        val outDir
        // THe number of threads provided
        val threads

    output:
        // Tuple contains the file basename and the clipped and sorted bam alignment file.
        tuple val(base), file("${base}-clipped-sorted.bam")

    publishDir "${outDir}", mode: 'copy'

    script:
    /*
    Primerclips takes a namesorted alignment in sam format. Thus, the first
    step uses samtools to convert the alignment bam to that format.

    Then, primerclip is run to softclip primers from the alignment.

    Finally, the primerclip output is converted back to bam format and sorted using
    samtools.
    */
    """
    #!/bin/bash

    samtools view -bS ${bam} | samtools sort -@ ${threads} -n -O sam > ${base}-sort.sam

    primerclip ${primersFile} ${base}-sort.sam ${base}-clip.sam

    samtools view -b ${base}-clip.sam | samtools sort -@ ${threads} > ${base}-clipped-sorted.bam
    """
}


// Uses Picard to mark duplicate reads in an alignment.
process MarkDuplicates {
    input:
        // Tuple contains the file basename and alignment bam file
        tuple val(base), file(bam)
        // The output directory name
        val outDir
        // The path to the picard jar file
        val picardLoc

    output:
        // Tuple contains the file basename and alignment bam file with
        // duplicates removed.
        tuple val(base), file("${base}-align-nodups.bam")

    publishDir "${outDir}", mode: 'copy'

    script:
    /*
    Uses picard to mark duplicates. Includes the following options:
        M - a metrics file (for later pipeline development)
        ASSUME_SORT_ORDER=coordinate - tells the tool to assume the file is
                                       sorted by coordinate
        REMOVE_DUPLICATES=true - does to write the duplicates in the output file
    */
    """
    #!/bin/bash

    java -jar ${picardLoc} MarkDuplicates I=${bam} \
    O=${base}-align-nodups.bam M=${base}-dup-metrics.txt ASSUME_SORT_ORDER=coordinate \
    REMOVE_DUPLICATES=true
    """
}

// Uses Abra2 to perform a realignment to improve indel quality.
process Realignment {
    input:
        // Tuple contains the file basename and the alignment bam file
        tuple val(base), file(bam)
        // The output directory name
        val outDir
        // The path to the abra2 jar file
        val abraLoc
        // Tuple contains the reference file name and reference fasta file
        tuple val(refName), file(ref)
        // The number of threads provided
        val threads

    output:
        // Tuple contains the file basename and the realigned and sorted bam file
        tuple val(base), file("${base}-realigned-sorted.bam")
    
    publishDir "${outDir}", mode: 'copy'

    script:
    /*
    Uses samtools to index the bam file

    Runs Abra2 to realign the reads and improve indel quality

    Sorts the realigned bam file using samtools.
    */
    """
    #!/bin/bash

    samtools index ${bam}

    java -Xmx8g -jar ${abraLoc} --in ${bam} --out ${base}-realigned.bam --ref ${ref} --threads ${threads}

    samtools sort -@ ${threads} ${base}-realigned.bam > ${base}-realigned-sorted.bam
    """
}


// Performs variant calling and filtering.
process CallVariants {
    input:
        // Tuple contains the file basename and alignment bam file
        tuple val(base), file(bam)
        // The script base directory name (to call python scripts)
        val baseDir
        // The output directory name
        val outDir
        // Tuple contains the reference name and reference fasta file
        tuple val(refName), file(ref)
        // The minimum coverage cutoff
        val minCov

    output:
        // Tuple contains the file basename, alignment bamfile, filtered snp vcf, and filtered indel vcf
        tuple val(base), file(bam), file("${base}-snps-filtered.vcf"), file("${base}-indels-filtered.vcf")
        // Tuple contains the unfiltered vcf and multiallelic filtered vcf files.
        tuple file("${base}.vcf"), file("${base}-biallelic.vcf")

    publishDir "${outDir}", mode: 'copy'

    script:
    /*
    Beings by indexing the bam file using samtools.

    Next, variants are called using freebayes. Freebayes can produce multiallelic
    variants (variants were multiple alternative alleles are provided). Thus, a 
    custom python script is used to select only the most prevalent allele at these 
    multiallelic positions. The INFO fields are corrected to account for this
    filtering step. (NOTE: this script also removes genotype information as
    this data is not used for downstreadm processing).

    VCFtools is used to split indels and snps apart, and each set of variatns are
    compressed with bgzip and tab indexed. (Required for filtering)

    BCFtools is then used to remove variants at position with less than the
    minimum coverage and that are not majority. The alternative and reference
    allele frequencies are calculated by dividing the occurrances (AO = alternative occurrance
    and RO = reference occurance) by the depth.
    */
    """
    #!/bin/bash

    samtools index ${bam}
    freebayes -f ${ref} ${bam} > ${base}.vcf

    python3 ${baseDir}/scripts/fix_multi_allelic.py -i ${base}.vcf -o ${base}-biallelic.vcf

    vcftools --keep-only-indels --vcf ${base}-biallelic.vcf --recode --recode-INFO-all --stdout > ${base}-indels.vcf
    bgzip ${base}-indels.vcf
    tabix ${base}-indels.vcf.gz

    vcftools --remove-indels --vcf ${base}-biallelic.vcf --recode --recode-INFO-all --stdout > ${base}-snps.vcf
    bgzip ${base}-snps.vcf
    tabix ${base}-snps.vcf.gz

    bcftools view -i "(INFO/DP >= ${minCov}) && ((INFO/AO / INFO/DP) > (INFO/RO / INFO/DP))" ${base}-indels.vcf.gz > ${base}-indels-filtered.vcf
    bcftools view -i "(INFO/DP >= ${minCov}) && ((INFO/AO / INFO/DP) > (INFO/RO / INFO/DP))" ${base}-snps.vcf.gz > ${base}-snps-filtered.vcf
    """
}

// Compiles a consensus by masking the reference and applying the variants.
process GenerateConsensus {
    input:
        // Tuple contains the file basename, the alignment bam, the snp vcf file
        // and the indel vcf file
        tuple val(base), file(bam), file(snps), file(indels)
        // The name of the base directory
        val baseDir
        // The name of the base directory
        val outDir
        // Tuple contains the reference file name and reference file
        tuple val(refName), file(ref)
        // The minimum coverage threshold
        val minCov

    output:
        // Tuple contains the consensus fasta and the sites that were masked in 
        // a bed file.
        tuple file("${base}-consensus.fasta"), file("${base}-mask-sites.bed")

    publishDir "${outDir}", mode: 'copy'

    script:
    /*
    The script first computes sites to mask by identifying sites that have less than
    the minimum coverage provided. However, there are formatting issues with the pileup
    format that make this difficult. Samtools mpileup's output has the depth in 0-based
    format, which makes it impossible to distinguish between a site with 0 and 1 coverage.

    Thus, the pipeline instead makes use of bedtools subtract. It first creates a pileup for only sites with
    coverage, uses an in-house script to filter sites with less than the minimum coverage
    and converts these into a bed file.

    Next, the pipeline creates a pileup containing every site, and converts this into a bed file using
    the in-house script. 

    Finally, the sites we want to keep (those above the minimum coverage threshold) are substracted
    from the bed file with every site, to give us the sites we want to mask.

    The bedtools maskfasta command is then used to mask the reference at these positions.

    Then the variants are applied to the mask fasta. The reason this is done after masking is 
    because the pileup (and therefore masking) positions do not account for indels, which would
    shift the genomic coordinates (we would end up masking things we did not want to).

    Finally, the fasta is wrapped to make it visually appealing. 
    */
    """
    #!/bin/bash

    samtools mpileup --no-BAQ -d 100000 -x -A -f ${ref} ${bam} > ${base}.pileup
    python3 ${baseDir}/scripts/pileup_to_bed.py -i ${base}.pileup -o passed-sites.bed --minCov ${minCov}

    samtools mpileup --no-BAQ -d 100000 -x -A -a -f ${ref} ${bam} > all-sites.pileup
    python3 ${baseDir}/scripts/pileup_to_bed.py -i all-sites.pileup -o all-sites.bed 

    bedtools subtract -a all-sites.bed -b passed-sites.bed > ${base}-mask-sites.bed

    bedtools maskfasta -fi ${ref} -bed ${base}-mask-sites.bed -fo masked.fasta

    bgzip ${snps}
    tabix ${snps}.gz

    bgzip ${indels}
    tabix ${indels}.gz

    bcftools consensus -f masked.fasta ${snps}.gz > with-snps.fasta

    bcftools consensus -f with-snps.fasta ${indels}.gz > with-indels-snps.fasta
    
    bioawk -c fastx '{ gsub(/\\n/,"",seq); print ">${base}"; print \$seq }' with-indels-snps.fasta > ${base}-consensus.fasta
    """
}
