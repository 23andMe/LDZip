#! /usr/bin/env nextflow

nextflow.enable.dsl=2
include { validateParameters; paramsSummaryLog } from 'plugin/nf-schema'
validateParameters()

CHROMS = params.chroms.tokenize(',')*.trim()

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
    tag { "chr${chr}" }
    memory { 8.GB * task.attempt }
    publishDir "${params.outdir}/logs/${task.process}/", mode: 'copy', pattern: ".command.log", overwrite: true, saveAs: {"${task.tag}.log"}
    publishDir "${params.outdir}/plinkLD/", mode: 'link', overwrite: true, pattern: "plink.chr${chr}.*", enabled: params.stage_plink

    input:
        tuple val(chr), path(pgen), path(pvar), path(psam)

    output:
        tuple val(chr), path("plink.chr${chr}.vcor"), path("plink.chr${chr}.pvar"), emit: vcor
        path(".command.log"), emit: log

    script:
    def pfile_base = pgen.baseName  // Get basename without .pgen extension
    def subset_snps = params.extract ? "--extract ${params.extract}" : ""
    def exclude_snps = params.exclude ? "--exclude ${params.exclude}" : ""
    def subset_samples = params.keep ? "--keep ${params.keep}" : ""
    def plink_ld_command = params.ld_command ? "${params.ld_command}" : "--r-phased ref-based cols=id,ref,alt,dprime"
    """
    ${PLINK2} \\
        --pfile ${pfile_base} ${subset_snps} ${exclude_snps} \\
        --rm-dup exclude-all \\
        --chr ${chr} \\
        --make-just-pvar \\
        --threads ${params.ld_threads} \\
        --out plink.chr${chr}
    cat plink.chr${chr}.pvar | grep -v '#' | cut -f 3 > ids
    ${PLINK2} \\
        --pfile ${pfile_base} --extract ids ${subset_samples} \\
        --ld-window-kb ${params.ld_window_kb} \\
        --ld-window-r2 ${params.ld_window_r2} \\
        ${plink_ld_command} \\
        --threads ${params.ld_threads} \\
        --memory ${params.ld_threads * (task.memory.toMega()-1024)} \\
        --out plink.chr${chr}
    """

    stub:
    """
    touch plink.chr${chr}.vcor
    touch plink.chr${chr}.pvar
    """
}

process compressLD {
    tag { "chr${chr}" }
    memory { 8.GB * task.attempt }
    publishDir "${params.outdir}/logs/${task.process}/", mode: 'copy', pattern: ".command.log", overwrite: true, saveAs: {"${task.tag}.log"}
    publishDir "${params.outdir}/chr/", mode: 'link', overwrite: true, pattern: "chr${chr}.ldzip.*", saveAs: { fname -> "${params.prefix}.${fname}" }

    input:
    tuple val(chr), path(vcor), path(snp_file)

    output:
        path("chr${chr}.ldzip.*"),              emit: ldzip
        path(".command.log"),           emit: log

    script:
    """
    ${LDZIP} compress plinkTabular \
      --ld_file ${vcor} \
      --snp_file ${snp_file} \
      --output_prefix chr${chr}.ldzip \
      --min ${params.min} \
      --bits ${params.bits} --min_col ${params.min_col}
    """

    stub:
    """
    touch chr${chr}.ldzip.i.bin
    touch chr${chr}.ldzip.io.bin
    touch chr${chr}.ldzip.x.PHASED_R.bin
    touch chr${chr}.ldzip.x.DPRIME.bin
    touch chr${chr}.ldzip.p.bin
    touch chr${chr}.ldzip.meta.json
    touch chr${chr}.ldzip.vars.txt
    """
}

process concatLD {
    tag "concat"
    memory { 8.GB * task.attempt }
    publishDir "${params.outdir}/whole_genome/", mode: 'link', overwrite: true, saveAs: { filename -> filename.replace('concat', params.prefix) }
    publishDir "${params.outdir}/logs/${task.process}/", mode: 'copy', pattern: ".command.log", overwrite: true, saveAs: {"${task.tag}.log"}

    input:
        path("*")

    output:
        path("concat.*"),                           emit: data
        path("concat.vars.txt"),                    emit: variants   
        path(".command.log"),           emit: log

    script:
    """
    ${LDZIP} concat \
    --inputs \$(ls -1 *.i.bin 2>/dev/null \
              | sed "s/\\.i\\.bin\$//" \
              | sort -V -u) \
    --output_prefix concat
    """
   
    stub:
    """
    touch concat.i.bin
    touch concat.io.bin
    touch concat.x.PHASED_R.bin
    touch concat.x.DPRIME.bin
    touch concat.p.bin
    touch concat.meta.json
    touch concat.vars.txt
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

    plink_ld_files      = ldPlink(pgen_files).vcor
    compressed_ld_files = compressLD(plink_ld_files).ldzip
    whole_genome = concatLD(compressed_ld_files.collect())
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

