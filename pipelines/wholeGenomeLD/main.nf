#! /usr/bin/env nextflow

nextflow.enable.dsl=2
include { validateParameters; paramsSummaryLog } from 'plugin/nf-schema'

process getChromosomeBounds {
    tag { "chr${chr}" }
    memory { 8.GB * task.attempt }

    input:
        tuple val(chr), path(pgen), path(pvar), path(psam)

    output:
        tuple val(chr), env('min_pos'), env('max_pos'), path(pgen), path(pvar), path(psam), emit: bounds

    script:
    """
    \$PLINK2 \\
        --pfile ${pgen.baseName} \\
        --chr ${chr} \\
        --make-just-pvar \\
        --threads 1 \\
        --out plink.chr${chr}
    export min_pos=\$(awk '!/^#/ {print \$2; exit}' plink.chr${chr}.pvar)
    export max_pos=\$(awk '!/^#/ {pos=\$2} END {print pos}' plink.chr${chr}.pvar)
    """

    stub:
    """
    min_pos=60343
    max_pos=227860
    """
}

process vcfToPgen {
    tag { "chr${chr}" }
    cpus { params.ld_threads }
    memory { 2.GB * task.cpus * task.attempt }
    publishDir { "${params.outdir}/logs/${task.process}/" }, mode: 'copy', pattern: ".command.log", overwrite: true, saveAs: {"${task.tag}.log"}
    publishDir "${params.outdir}/pgen/", mode: 'link', overwrite: true, pattern: "converted.chr*.*", enabled: params.stage_pgen

    input:
        tuple val(chr), path(vcf)

    output:
        tuple val(chr), path("converted.chr${chr}.pgen"), path("converted.chr${chr}.pvar"), path("converted.chr${chr}.psam"), emit: pfiles
        path(".command.log"), emit: log

    script:
    def vcf_file = params.vcf_template.replace('{CHR}', chr)
    """
    \$PLINK2 \\
        --vcf ${vcf_file}.vcf.gz \\
        --make-pgen \\
        --threads ${task.cpus} \\
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
    cpus { params.ld_threads }
    memory { 8.GB * task.cpus * task.attempt }
    publishDir { "${params.outdir}/logs/${task.process}/" }, mode: 'copy', pattern: ".command.log", overwrite: true, saveAs: {"${task.tag}.log"}
    publishDir "${params.outdir}/plinkLD/", mode: 'link', overwrite: true, pattern: "plink.chr*.*", enabled: params.stage_plink

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
    def ld_filter = params.ld_filter ? "${params.ld_filter}" : ""
    def plink_ld_command = params.ld_command ? "${params.ld_command}" : "--r-phased ref-based cols=id,ref,alt,dprime"
    def chunk_filter = (end_bp > 0) ? "--chr ${chr} --from-bp ${start_bp} --to-bp ${end_bp}" : "--chr ${chr}"
    """
    \$PLINK2 \\
        --pfile ${pfile_base} ${subset_snps} ${exclude_snps} ${chunk_filter} ${ld_filter} --force-intersect \\
        --rm-dup exclude-all \\
        --make-just-pvar \\
        --threads ${task.cpus} \\
        --out plink.chr${chr}_${chunk_id} || touch plink.chr${chr}_${chunk_id}.pvar

    if [ -s plink.chr${chr}_${chunk_id}.pvar ]; then
        cat plink.chr${chr}_${chunk_id}.pvar | grep -v '#' | cut -f 3 > ids
        \$PLINK2 \\
            --pfile ${pfile_base} --extract ids ${subset_samples} \\
            --ld-window-kb ${params.ld_window_kb} \\
            --ld-window-r2 ${params.ld_window_r2} \\
            ${plink_ld_command} \\
            --threads ${task.cpus} \\
            --memory ${(task.memory.toMega()-1024)} \\
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
    publishDir { "${params.outdir}/logs/${task.process}/" }, mode: 'copy', pattern: ".command.log", overwrite: true, saveAs: {"${task.tag}.log"}
    publishDir "${params.outdir}/chunks/", mode: 'link', overwrite: true, pattern: "chr*.ldzip.*", enabled: params.stage_chunk

    input:
        tuple val(chr), val(chunk_id), val(type), path(ld_file), path(snp_file)

    output:
        tuple val(chr), val(chunk_id), path("chr${chr}_${chunk_id}.ldzip.*"), emit: ldzip
        path(".command.log"), emit: log

    script:
    def compress_cmd = type == "tabular" ? "plinkTabular" : "plinkSquare"
    def extra_args = type == "tabular" ? "--min_col ${params.min_col}" : "--type ${params.npz_ld_type}"
    """
    if [ -s ${ld_file} ]; then
        \$LDZIP compress ${compress_cmd} \\
          --ld_file ${ld_file} \\
          --snp_file ${snp_file} \\
          --output_prefix chr${chr}_${chunk_id}.ldzip \\
          --min ${params.min} \\
          --bits ${params.bits} \\
          ${extra_args}
    else
        touch chr${chr}_${chunk_id}.ldzip.{i.bin,i.bin.index,x.PHASED_R.bin,x.PHASED_R.bin.index,x.DPRIME.bin,x.DPRIME.bin.index,p.bin,meta.json,vars.txt}
    fi
    """

    stub:
    """
    touch chr${chr}_${chunk_id}.ldzip.{i.bin,i.bin.index,x.PHASED_R.bin,x.PHASED_R.bin.index,x.DPRIME.bin,x.DPRIME.bin.index,p.bin,meta.json,vars.txt}
    """
}

process convertNpzToBinary {
    tag { "chr${chr}-chunk${chunk_id}" }
    memory { env('CONSTRAIN_MEMORY') ? 8.GB : 16.GB * task.attempt }
    publishDir { "${params.outdir}/logs/${task.process}/" }, mode: 'copy', pattern: ".command.log", overwrite: true, saveAs: {"${task.tag}.log"}
    publishDir "${params.outdir}/plinkLD/", mode: 'link', overwrite: true, pattern: "plink.chr*.*", enabled: params.stage_binary

    input:
        tuple val(chr), val(chunk_id), path(npz_file), path(gz_file)

    output:
        tuple val(chr), val(chunk_id), path("chr${chr}_${chunk_id}.bin"), path("chr${chr}_${chunk_id}.bin.vars"), emit: binary
        path(".command.log"), emit: log

    script:
    """
    python3 ${projectDir}/../../scripts/npz_to_binary_matrix.py \\
        ${npz_file} \\
        ${gz_file} \\
        chr${chr}_${chunk_id}.bin
    """

    stub:
    """
    touch chr${chr}_${chunk_id}.bin
    touch chr${chr}_${chunk_id}.bin.vars
    """
}

process concatPairwise {
    tag { "chr${chr}-pair${pair_id}" }
    memory { 8.GB * task.attempt }
    publishDir { "${params.outdir}/logs/${task.process}/" }, mode: 'copy', pattern: ".command.log", overwrite: true, saveAs: {"${task.tag}.log"}

    input:
        tuple val(chr), val(pair_id), path("*")

    output:
        tuple val(chr), path("pair_${chr}_${pair_id}.*"), emit: paired
        path(".command.log"), emit: log

    script:
    """
    \$LDZIP concat \\
    --inputs \$(ls -1 *.i.bin 2>/dev/null \\
              | sed "s/\\.i\\.bin\$//" \\
              | sort -V -u) \\
    --output_prefix pair_${chr}_${pair_id}
    """

    stub:
    """
    touch pair_${chr}_${pair_id}.{i.bin,i.bin.index,x.PHASED_R.bin,x.PHASED_R.bin.index,p.bin,meta.json,vars.txt}
    """
}


process concatChromosome {
    tag { "chr${chr}" }
    memory { 8.GB * task.attempt }
    publishDir "${params.outdir}/chr/", mode: 'link', overwrite: true, saveAs: { filename -> filename.replace('concat_chr', "${params.prefix}.chr${chr}.ldzip") }, enabled: params.stage_chr
    publishDir { "${params.outdir}/logs/${task.process}/" }, mode: 'copy', pattern: ".command.log", overwrite: true, saveAs: {"${task.tag}.log"}

    input:
        tuple val(chr), path("*")

    output:
        path("concat_chr${chr}.*"), emit: data
        path(".command.log"), emit: log

    script:
    """
    inputs=\$(for f in *.i.bin; do [ -s "\$f" ] && echo "\${f%.i.bin}"; done | sort -V -u)
    \$LDZIP concat --inputs \$inputs --output_prefix concat_chr${chr}
    """

    stub:
    """
    touch concat_chr${chr}.{i.bin,i.bin.index,x.PHASED_R.bin,x.PHASED_R.bin.index,x.DPRIME.bin,x.DPRIME.bin.index,p.bin,meta.json,vars.txt}
    """
}

process concatGenome {
    tag "genome"
    memory { 8.GB * task.attempt }
    publishDir "${params.outdir}/whole_genome/", mode: 'link', overwrite: true, saveAs: { filename -> filename.replace('concat', params.prefix) }
    publishDir { "${params.outdir}/logs/${task.process}/" }, mode: 'copy', pattern: ".command.log", overwrite: true, saveAs: {"concat_genome.log"}

    input:
        path("*")

    output:
        path("concat.*"),                           emit: data
        path("concat.vars.txt"),                    emit: variants   
        path(".command.log"),           emit: log

    script:
    def chrom_order = params.chroms.tokenize(',')*.trim().collect { "concat_chr${it}" }.join(' ')
    """
    \$LDZIP concat \\
        --inputs ${chrom_order} \\
        --naive \\
        --output_prefix concat
    """
    stub:
    """
    touch concat.{i.bin,i.bin.index,x.PHASED_R.bin,x.PHASED_R.bin.index,x.DPRIME.bin,x.DPRIME.bin.index,p.bin,meta.json,vars.txt}
    """
}

process indexVariants {
    tag "sqlite"
    memory { env('CONSTRAIN_MEMORY') ? 8.GB : 8.GB * (4 ** (task.attempt - 1)) }
    publishDir "${params.outdir}/whole_genome/", pattern: "*.sqlite", mode: 'link', overwrite: true, saveAs: { filename -> filename.replace('concat', params.prefix) }
    publishDir { "${params.outdir}/logs/${task.process}/" }, mode: 'copy', pattern: ".command.log", overwrite: true, saveAs: {"${task.tag}.log"}

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

    workflow.onComplete = {
        def summary = """
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
        outlog.newWriter().withWriter { w ->
            w << summary
        }

    }

    workflow.onError = {
        println "Oops .. something went wrong"
    }

    validateParameters()

    CHROMS = params.chroms.tokenize(',')*.trim()

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

    // Validate that only one input type is provided
    def input_count = [params.vcf_template, params.pfile_template, params.npz_template].count { it -> it != null }
    if (input_count > 1) {
        error "Please provide only ONE of: --vcf_template, --pfile_template, or --npz_template"
    }
    if (input_count == 0) {
        error "Please provide one of: --vcf_template, --pfile_template, or --npz_template"
    }

    // Track chunks per chromosome for groupKey
    def chunks_per_chr = [:]

    // Branch based on input type
    if (params.npz_template) {
        // NPZ workflow: use template with {CHR} and {CHUNK} placeholders
        // Find all NPZ files by globbing the template pattern
        def npz_glob = params.npz_template.replace('{CHR}', '*').replace('{CHUNK}', '*')

        // Create regex pattern from template for parsing
        // TODO: generalize this to work with any template format
        def template_basename = new File(params.npz_template).name.replaceAll(/\.npz$/, '')
        def template_pattern = template_basename
            .replaceAll(/\{CHR\}/, '([^_]+)')
            .replaceAll(/\{CHUNK\}/, '(.+)')

        npz_files = channel.fromPath(npz_glob, followLinks: true, type: 'file', checkIfExists: true)
            .map { npz_file ->
                // Extract chr and chunk from filename using template pattern
                def basename = npz_file.name.replaceAll(/\.npz$/, '')
                def matcher = basename =~ template_pattern

                if (!matcher.matches()) {
                    return null
                }

                def chr = matcher.group(1)
                def chunk_str = matcher.group(2)

                // Try to extract starting position for sorting (if chunk is position-based)
                def start_pos_match = chunk_str =~ /^(\d+)/
                def sort_key = start_pos_match ? start_pos_match.group(1).toLong() : chunk_str

                def gz_file = file(npz_file.toString().replace('.npz', '.gz'))
                return tuple(chr, sort_key, chunk_str, npz_file, gz_file)
            }
            .filter { it -> it != null && it[0] in CHROMS }
            .groupTuple(by: 0)  // Group by chromosome
            .flatMap { chr, sort_keys, chunk_strs, npz_files_list, gz_files_list ->
                // Sort by sort_key within chromosome
                def sorted = [sort_keys, chunk_strs, npz_files_list, gz_files_list].transpose()
                    .sort { it -> it[0] }

                // Assign sequential chunk IDs and track count
                def chunk_list = sorted.withIndex().collect { item, idx ->
                    def chunk_id = idx + 1
                    [chr, chunk_id, item[2], item[3]]
                }
                chunks_per_chr[chr.toString()] = chunk_list.size()
                return chunk_list
            }

        // Convert NPZ to binary
        binary_files = convertNpzToBinary(npz_files).binary
            .map { chr, chunk_id, bin_file, vars_file -> tuple(chr, chunk_id, 'square', bin_file, vars_file) }

        // Compress from binary
        compressed_ld_files = compressLD(binary_files).ldzip

    } else {
        // PLINK workflow (VCF or PGEN input)
        // Set chunk_size and overlap_size (in bp) based on ld_window_kb if not provided
        def chunk_size = params.chunk_size_kb ? (params.chunk_size_kb * 1000) : (params.ld_window_kb * 2 * 1000)
        def overlap_size = params.overlap_size_kb ? (params.overlap_size_kb * 1000) : (params.ld_window_kb * 1000)

        if (params.vcf_template) {
            vcf_files = channel.from(CHROMS)
                               .map { chr ->
                                   def vcf = params.vcf_template.replace('{CHR}', chr)
                                   tuple(chr, file(vcf))
                               }

            pgen_files = vcfToPgen(vcf_files).pfiles
        } else {
            pgen_files = channel.from(CHROMS)
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
            .map { chr, chunk_id, vcor, snp_file -> tuple(chr, chunk_id, 'tabular', vcor, snp_file) }
        compressed_ld_files = compressLD(plink_ld_files).ldzip
    }

    // Optional pairwise concatenation stage
    if (params.concat_pairwise) {
        // Stage 1: Pairwise concatenation
        pairwise_input = compressed_ld_files
            .map { chr, chunk_id, files ->
                def pair_id = (chunk_id - 1).intdiv(2)
                tuple(chr, pair_id, chunk_id, files)
            }
            .groupTuple(by: [0, 1], size: 2, remainder: true)
            .map { chr, pair_id, chunk_ids, file_lists ->
                def sorted = [chunk_ids, file_lists].transpose().sort { it -> it[0] }
                def all_files = sorted.collect { it -> it[1] }.flatten()
                tuple(chr, pair_id, all_files)
            }

        // Run pairwise concat
        paired_results = concatPairwise(pairwise_input).paired

        // Group by chromosome for stage 2 with groupKey (number of pairs per chr)
        chr_input = paired_results
            .map { chr, files ->
                def n_pairs = (chunks_per_chr[chr.toString()] + 1).intdiv(2)
                tuple(groupKey(chr, n_pairs), files)
            }
            .groupTuple()
            .map { chr, file_lists -> tuple(chr, file_lists.flatten()) }
    } else {
        // Skip pairwise, go directly to chromosome concatenation with groupKey
        chr_input = compressed_ld_files
            .map { chr, chunk_id, files ->
                tuple(groupKey(chr, chunks_per_chr[chr.toString()]), files)
            }
            .groupTuple()
            .map { chr, file_lists -> tuple(chr, file_lists.flatten()) }
    }

    // Concat all chunks within each chromosome
    chr_concat = concatChromosome(chr_input).data

    // Concat all chromosomes without --overlapping
    whole_genome = concatGenome(chr_concat.collect())

    indexVariants(whole_genome.data)
}
