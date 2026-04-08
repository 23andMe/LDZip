#! /usr/bin/env nextflow

nextflow.enable.dsl=2
include { validateParameters; paramsSummaryLog } from 'plugin/nf-schema'
validateParameters()

CHROMS = params.chroms.tokenize(',')*.trim()

process getChromosomeBounds {
    tag { "chr${chr}" }
    memory { 8.GB * task.attempt }

    input:
        tuple val(chr), path(pgen), path(pvar), path(psam)

    output:
        tuple val(chr), env(min_pos), env(max_pos), path(pgen), path(pvar), path(psam), emit: bounds

    script:
    """
	${PLINK2} \\
		--pfile ${pgen.baseName} \\
		--chr ${chr} \\
		--make-just-pvar \\
		--threads 1 \\
		--out plink.chr${chr}
	min_pos=\$(awk '!/^#/ {print \$2; exit}' plink.chr${chr}.pvar)
	max_pos=\$(awk '!/^#/ {pos=\$2} END {print pos}' plink.chr${chr}.pvar)
    """

    stub:
    """
    min_pos=60343
    max_pos=227860
    """
}

process vcfToPgen {
    tag { "chr${chr}" }
    memory { 8.GB * task.attempt }
    publishDir "${params.outdir}/logs/${task.process}/", mode: 'copy', pattern: ".command.log", overwrite: true, saveAs: {"${task.tag}.log"}
    publishDir "${params.outdir}/pgen/", mode: 'link', overwrite: true, pattern: "converted.chr${chr}.*", enabled: params.stage_pgen

    input:
        tuple val(chr), path(vcf)

    output:
        tuple val(chr), path("converted.chr${chr}.pgen"), path("converted.chr${chr}.pvar"), path("converted.chr${chr}.psam"), emit: pfiles
        path(".command.log"), emit: log

    script:
    def vcf_file = params.vcf_template.replace('{CHR}', chr)
    """
    ${PLINK2} \\
        --vcf ${vcf_file}.vcf.gz \\
        --make-pgen \\
        --threads ${params.ld_threads} \\
        --out converted.chr${chr}
    """

    stub:
    """
    touch converted.chr${chr}.pgen
    touch converted.chr${chr}.pvar
    touch converted.chr${chr}.psam
    """
}

process ldPlink {
    tag { "chr${chr}-chunk${chunk_id}" }
    memory { 8.GB * task.attempt }
    publishDir "${params.outdir}/logs/${task.process}/", mode: 'copy', pattern: ".command.log", overwrite: true, saveAs: {"${task.tag}.log"}
    publishDir "${params.outdir}/plinkLD/", mode: 'link', overwrite: true, pattern: "plink.chr${chr}_${chunk_id}.*", enabled: params.stage_plink

    input:
        tuple val(chr), val(chunk_id), val(start_bp), val(end_bp), path(pgen), path(pvar), path(psam)

    output:
        tuple val(chr), val(chunk_id), path("plink.chr${chr}_${chunk_id}.vcor"), path("plink.chr${chr}_${chunk_id}.pvar"), emit: vcor
        path(".command.log"), emit: log

    script:
    def pfile_base = pgen.baseName  // Get basename without .pgen extension
    def subset_snps = params.extract ? "--extract ${params.extract}" : ""
    def exclude_snps = params.exclude ? "--exclude ${params.exclude}" : ""
    def subset_samples = params.keep ? "--keep ${params.keep}" : ""
    def plink_ld_command = params.ld_command ? "${params.ld_command}" : "--r-phased ref-based cols=id,ref,alt,dprime"
    def chunk_filter = (end_bp > 0) ? "--chr ${chr} --from-bp ${start_bp} --to-bp ${end_bp}" : "--chr ${chr}"
    """
    ${PLINK2} \\
        --pfile ${pfile_base} ${subset_snps} ${exclude_snps} ${chunk_filter} \\
        --rm-dup exclude-all \\
        --make-just-pvar \\
        --threads ${params.ld_threads} \\
        --out plink.chr${chr}_${chunk_id} || touch plink.chr${chr}_${chunk_id}.pvar

    if [ -s plink.chr${chr}_${chunk_id}.pvar ]; then
        cat plink.chr${chr}_${chunk_id}.pvar | grep -v '#' | cut -f 3 > ids
        ${PLINK2} \\
            --pfile ${pfile_base} --extract ids ${subset_samples} \\
            --ld-window-kb ${params.ld_window_kb} \\
            --ld-window-r2 ${params.ld_window_r2} \\
            ${plink_ld_command} \\
            --threads ${params.ld_threads} \\
            --memory ${params.ld_threads * (task.memory.toMega()-1024)} \\
            --out plink.chr${chr}_${chunk_id}
    else
        touch plink.chr${chr}_${chunk_id}.vcor
    fi
    """

    stub:
    """
    touch plink.chr${chr}_${chunk_id}.{vcor,pvar}
    """
}

process compressLD {
    tag { "chr${chr}-chunk${chunk_id}" }
    memory { 8.GB * task.attempt }
    publishDir "${params.outdir}/logs/${task.process}/", mode: 'copy', pattern: ".command.log", overwrite: true, saveAs: {"${task.tag}.log"}
    publishDir "${params.outdir}/chunks/", mode: 'link', overwrite: true, pattern: "chr${chr}_${chunk_id}.ldzip.*", enabled: params.stage_chunk

    input:
    tuple val(chr), val(chunk_id), path(vcor), path(snp_file)

    output:
        tuple val(chr), val(chunk_id), path("chr${chr}_${chunk_id}.ldzip.*"), emit: ldzip
        path(".command.log"), emit: log

    script:
    """
    if [ -s ${vcor} ]; then
        ${LDZIP} compress plinkTabular \\
          --ld_file ${vcor} \\
          --snp_file ${snp_file} \\
          --output_prefix chr${chr}_${chunk_id}.ldzip \\
          --min ${params.min} \\
          --bits ${params.bits} --min_col ${params.min_col}
    else
        touch chr${chr}_${chunk_id}.ldzip.{i.bin,io.bin,io.index,x.PHASED_R.bin,x.DPRIME.bin,p.bin,meta.json,vars.txt}
    fi
    """

    stub:
    """
    touch chr${chr}_${chunk_id}.ldzip.{i.bin,io.bin,io.index,x.PHASED_R.bin,x.DPRIME.bin,p.bin,meta.json,vars.txt}
    """
}

process concatChromosome {
    tag { "chr${chr}" }
    memory { 8.GB * task.attempt }
    publishDir "${params.outdir}/chr/", mode: 'link', overwrite: true, saveAs: { filename -> filename.replace('concat_chr', "${params.prefix}.chr${chr}.ldzip") }, enabled: params.stage_chr
    publishDir "${params.outdir}/logs/${task.process}/", mode: 'copy', pattern: ".command.log", overwrite: true, saveAs: {"${task.tag}.log"}

    input:
        tuple val(chr), path("*")

    output:
        path("concat_chr${chr}.*"), emit: data
        path(".command.log"), emit: log

    script:
    """
    ${LDZIP} concat \\
    --inputs \$(ls -1 *.i.bin 2>/dev/null \\
              | sed "s/\\.i\\.bin\$//" \\
              | sort -V -u) \\
    --output_prefix concat_chr${chr} \\
    --overlapping
    """

    stub:
    """
    touch concat_chr${chr}.{i.bin,io.bin,io.index,x.PHASED_R.bin,x.DPRIME.bin,p.bin,meta.json,vars.txt}
    """
}

process concatGenome {
    tag "genome"
    memory { 8.GB * task.attempt }
    publishDir "${params.outdir}/whole_genome/", mode: 'link', overwrite: true, saveAs: { filename -> filename.replace('concat', params.prefix) }
    publishDir "${params.outdir}/logs/${task.process}/", mode: 'copy', pattern: ".command.log", overwrite: true, saveAs: {"concat_genome.log"}

    input:
        path("*")

    output:
        path("concat.*"),                           emit: data
        path("concat.vars.txt"),                    emit: variants   
        path(".command.log"),           emit: log

    script:
    """
    ${LDZIP} concat \\
    --inputs \$(ls -1 *.i.bin 2>/dev/null \\
              | sed "s/\\.i\\.bin\$//" \\
              | sort -V -u) \\
    --output_prefix concat
    """
    stub:
    """
    touch concat.{i.bin,io.bin,io.index,x.PHASED_R.bin,x.DPRIME.bin,p.bin,meta.json,vars.txt}
    """
}

process indexVariants {
    tag "sqlite"
    memory { 8.GB * (4 ** (task.attempt - 1)) }
    publishDir "${params.outdir}/whole_genome/", pattern: "*.sqlite", mode: 'link', overwrite: true, saveAs: { filename -> filename.replace('concat', params.prefix) }
    publishDir "${params.outdir}/logs/${task.process}/", mode: 'copy', pattern: ".command.log", overwrite: true, saveAs: {"${task.tag}.log"}

    input:
        path("*")

    output:
        path("concat.sqlite"),      emit: sqlite
        path(".command.log"),           emit: log

    script:
    """
    #!/usr/bin/env Rscript

    library(LDZipMatrix)
    ld = LDZipMatrix("concat")
    buildIndex(ld)
    """


    stub:
    """
    touch concat.sqlite
    """
}


workflow {

        date = new Date().format( 'yyyyMMdd' )
    log.info ""
    log.info "__________________________________________________"
    log.info "--------------------------------------------------"
    log.info "*         Whole Genome LDZip Compression         *"
    log.info "__________________________________________________"
    log.info "--------------------------------------------------"
    log.info ""
    log.info " Temp Directory = ${workDir}"
    log.info ""

    // Validate that either vcf_template or pfile_template is provided
    if (params.vcf_template && params.pfile_template) {
        error "Please provide either --vcf_template OR --pfile_template, not both"
    }
    if (!params.vcf_template && !params.pfile_template) {
        error "Please provide either --vcf_template or --pfile_template"
    }

    // Set chunk_size and overlap_size (in bp) based on ld_window_kb if not provided
    def chunk_size = params.chunk_size_kb ? (params.chunk_size_kb * 1000) : (params.ld_window_kb * 2 * 1000)
    def overlap_size = params.overlap_size_kb ? (params.overlap_size_kb * 1000) : (params.ld_window_kb * 1000)

    if (params.vcf_template) {
        vcf_files = Channel.from(CHROMS)
                           .map { chr ->
                               def vcf = params.vcf_template.replace('{CHR}', chr)
                               tuple(chr, file(vcf))
                           }

        pgen_files = vcfToPgen(vcf_files).pfiles
    } else {
        pgen_files = Channel.from(CHROMS)
                           .map { chr ->
                               def base = params.pfile_template.replace('{CHR}', chr)
                               tuple(chr,
                                     file("${base}.pgen"),
                                     file("${base}.pvar"),
                                     file("${base}.psam"))
                           }
    }

    // Get chromosome bounds
    chr_bounds = getChromosomeBounds(pgen_files)

    // Calculate number of chunks per chromosome for grouping
    def chunks_per_chr = [:]

    chunked_pgen = chr_bounds.bounds.flatMap { chr, min_pos_str, max_pos_str, pgen, pvar, psam ->
        def min_pos = min_pos_str.toInteger()
        def max_pos = max_pos_str.toInteger()

        def step = chunk_size - overlap_size
        def chunks = (min_pos..max_pos).step(step).withIndex(1).collect { start_bp, chunk_id ->
            def end_bp = Math.min(start_bp + chunk_size - 1, max_pos)
            [chr, chunk_id, start_bp, end_bp, pgen, pvar, psam]
        }

        chunks_per_chr[chr.toString()] = chunks.size()
        return chunks
    }

    plink_ld_files = ldPlink(chunked_pgen).vcor
    compressed_ld_files = compressLD(plink_ld_files).ldzip

    // Group by chromosome with groupKey for immediate processing when complete
    chr_grouped = compressed_ld_files
        .map { chr, chunk_id, files -> tuple(groupKey(chr, chunks_per_chr[chr.toString()]), files) }
        .groupTuple()
        .map { chr, file_lists ->
            def filtered = file_lists.findAll { chunk_files ->
                chunk_files.find { it.name.endsWith('.vars.txt') }?.size() > 0
            }
            tuple(chr, filtered.flatten())
        }
        .filter { chr, files -> files.size() > 0 }

    // Concat chunks within each chromosome with --overlapping
    chr_concat = concatChromosome(chr_grouped).data

    // Concat all chromosomes without --overlapping
    whole_genome = concatGenome(chr_concat.collect())

    indexVariants(whole_genome.data)
}



workflow.onComplete = {
    summary = """
    Pipeline execution summary
    ---------------------------
    Completed at: ${workflow.complete}
    Duration    : ${workflow.duration}
    Success     : ${workflow.success}
    workDir     : ${workflow.workDir}
    exit status : ${workflow.exitStatus}
    Command Line: ${workflow.commandLine}
    """

    println summary
    def outlog = new File("nf_log.txt")
    outlog.newWriter().withWriter {
        outlog << summary
   }

}


workflow.onError = {
    println "Oops .. something went wrong"
}

