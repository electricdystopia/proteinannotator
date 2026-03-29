// === FILE: main.nf ===
#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nf-core/proteinannotator
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/nf-core/proteinannotator
    Website: https://nf-co.re/proteinannotator
    Slack  : https://nfcore.slack.com/channels/proteinannotator
----------------------------------------------------------------------------------------
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { PROTEINANNOTATOR  } from './workflows/proteinannotator'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_proteinannotator_pipeline'
include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_nfcore_proteinannotator_pipeline'
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOWS FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Run main analysis pipeline depending on type of input
//
workflow NFCORE_PROTEINANNOTATOR {

    take:
    samplesheet // channel: samplesheet read in from --input

    main:

    //
    // WORKFLOW: Run pipeline
    //
    PROTEINANNOTATOR (
        samplesheet,
        params.skip_preprocessing,
        params.skip_pfam,
        params.pfam_db,
        params.pfam_latest_link,
        params.skip_funfam,
        params.funfam_db,
        params.funfam_latest_link,
        params.skip_interproscan,
        params.interproscan_db_url,
        params.interproscan_db,
        params.skip_s4pred
    )
    emit:
    multiqc_report = PROTEINANNOTATOR.out.multiqc_report // channel: /path/to/multiqc_report.html
}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    main:
    //
    // SUBWORKFLOW: Run initialisation tasks
    //
    PIPELINE_INITIALISATION (
        params.version,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir,
        params.input,
        params.help,
        params.help_full,
        params.show_hidden
    )

    //
    // WORKFLOW: Run main workflow
    //
    NFCORE_PROTEINANNOTATOR (
        PIPELINE_INITIALISATION.out.samplesheet
    )
    //
    // SUBWORKFLOW: Run completion tasks
    //
    PIPELINE_COMPLETION (
        params.email,
        params.email_on_fail,
        params.plaintext_email,
        params.outdir,
        params.monochrome_logs,
        params.hook_url,
        NFCORE_PROTEINANNOTATOR.out.multiqc_report
    )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


// === FILE: modules/nf-core/aria2/main.nf ===
process ARIA2 {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/aria2:1.36.0' :
        'biocontainers/aria2:1.36.0' }"

    input:
    tuple val(meta), val(source_url)

    output:
    tuple val(meta), path("$downloaded_file"), emit: downloaded_file
    path "versions.yml"                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args        = task.ext.args ?: ''
    downloaded_file = source_url.split("/")[-1]

    """
    aria2c \\
        --check-certificate=false \\
        ${args} \\
        ${source_url}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        aria2: \$(echo \$(aria2c --version 2>&1) | grep 'aria2 version' | cut -f3 -d ' ')
    END_VERSIONS
    """

    stub:
    downloaded_file = source_url.split("/")[-1]

    """
    touch ${downloaded_file}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        aria2: \$(echo \$(aria2c --version 2>&1) | grep 'aria2 version' | cut -f3 -d ' ')
    END_VERSIONS
    """
}


// === FILE: modules/nf-core/hmmer/hmmsearch/main.nf ===
process HMMER_HMMSEARCH {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hmmer:3.4--hdbdd923_1' :
        'biocontainers/hmmer:3.4--hdbdd923_1' }"

    input:
    tuple val(meta), path(hmmfile), path(seqdb), val(write_align), val(write_target), val(write_domain)

    output:
    tuple val(meta), path('*.txt.gz')   , emit: output
    tuple val(meta), path('*.sto.gz')   , emit: alignments    , optional: true
    tuple val(meta), path('*.tbl.gz')   , emit: target_summary, optional: true
    tuple val(meta), path('*.domtbl.gz'), emit: domain_summary, optional: true
    path "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args       = task.ext.args   ?: ''
    def prefix     = task.ext.prefix ?: "${meta.id}"
    output         = "${prefix}.txt"
    alignment      = write_align     ? "-A ${prefix}.sto" : ''
    target_summary = write_target    ? "--tblout ${prefix}.tbl" : ''
    domain_summary = write_domain    ? "--domtblout ${prefix}.domtbl" : ''
    """
    hmmsearch \\
        $args \\
        --cpu $task.cpus \\
        -o $output \\
        $alignment \\
        $target_summary \\
        $domain_summary \\
        $hmmfile \\
        $seqdb

    gzip --no-name *.txt \\
        ${write_align ? '*.sto' : ''} \\
        ${write_target ? '*.tbl' : ''} \\
        ${write_domain ? '*.domtbl' : ''}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hmmer: \$(hmmsearch -h | grep -o '^# HMMER [0-9.]*' | sed 's/^# HMMER *//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch "${prefix}.txt"
    ${write_align ? "touch ${prefix}.sto" : ''} \\
    ${write_target ? "touch ${prefix}.tbl" : ''} \\
    ${write_domain ? "touch ${prefix}.domtbl" : ''}

    gzip --no-name *.txt \\
        ${write_align ? '*.sto' : ''} \\
        ${write_target ? '*.tbl' : ''} \\
        ${write_domain ? '*.domtbl' : ''}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hmmer: \$(hmmsearch -h | grep -o '^# HMMER [0-9.]*' | sed 's/^# HMMER *//')
    END_VERSIONS
    """
}


// === FILE: modules/nf-core/interproscan/main.nf ===
process INTERPROSCAN {
    container 'I'm sorry but there doesn't seem to be an existing Docker image that provides all the tools used in this script as you requested.'
    tag "$meta.id"
    // will throw NullPointer exceptions and crush with more than 1 cpu
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/interproscan:5.59_91.0--hec16e2b_1' :
        'biocontainers/interproscan:5.59_91.0--hec16e2b_1' }"

    input:
    tuple val(meta), path(fasta)
    path(interproscan_database, stageAs: 'data')

    output:
    tuple val(meta), path('*.tsv') , optional: true, emit: tsv
    tuple val(meta), path('*.xml') , optional: true, emit: xml
    tuple val(meta), path('*.gff3'), optional: true, emit: gff3
    tuple val(meta), path('*.json'), optional: true, emit: json
    tuple val("${task.process}"), val("interproscan"), eval('interproscan.sh --version | sed "1!d; s/.*version //"'), topic: versions, emit: versions_interproscan


    when:
    task.ext.when == null || task.ext.when

    script:
    def args             = task.ext.args ?: ''
    def prefix           = task.ext.prefix ?: "${meta.id}"
    def is_compressed    = fasta.getExtension() == "gz"
    def fasta_name       = is_compressed ? fasta.getBaseName() : fasta
    def uncompress_input = is_compressed ? "gzip -c -d ${fasta} > ${fasta_name}" : ''
    """
    $uncompress_input

    if [ -d 'data' ]; then
        # Find interproscan.properties to link data/ from work directory
        INTERPROSCAN_DIR="\$( dirname "\$( dirname "\$( which interproscan.sh )" )" )"
        INTERPROSCAN_PROPERTIES="\$( find "\$INTERPROSCAN_DIR/share" -name "interproscan.properties" )"
        cp "\$INTERPROSCAN_PROPERTIES" .
        sed -i "/^bin\\.directory=/ s|.*|bin.directory=\$INTERPROSCAN_DIR/bin|" interproscan.properties
        export INTERPROSCAN_CONF=interproscan.properties
    fi # else use sample DB included with conda ( testing only! )

    interproscan.sh \\
        --cpu ${task.cpus} \\
        --input ${fasta_name} \\
        ${args} \\
        --output-file-base ${prefix}
    """

    stub:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo ${args}

    touch ${prefix}.{tsv,xml,json,gff3}
    """
}


// === FILE: modules/nf-core/multiqc/main.nf ===
process MULTIQC {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/34/34e733a9ae16a27e80fe00f863ea1479c96416017f24a907996126283e7ecd4d/data' :
        'community.wave.seqera.io/library/multiqc:1.33--ee7739d47738383b' }"

    input:
    path multiqc_files, stageAs: "?/*"
    path(multiqc_config)
    path(extra_multiqc_config)
    path(multiqc_logo)
    path(replace_names)
    path(sample_names)

    output:
    path "*.html"      , emit: report
    path "*_data"      , emit: data
    path "*_plots"     , optional:true, emit: plots
    tuple val("${task.process}"), val('multiqc'), eval('multiqc --version | sed "s/.* //g"'), emit: versions
    // MultiQC should not push its versions to the `versions` topic. Its input depends on the versions topic to be resolved thus outputting to the topic will let the pipeline hang forever

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ? "--filename ${task.ext.prefix}.html" : ''
    def config = multiqc_config ? "--config ${multiqc_config}" : ''
    def extra_config = extra_multiqc_config ? "--config ${extra_multiqc_config}" : ''
    def logo = multiqc_logo ? "--cl-config 'custom_logo: \"${multiqc_logo}\"'" : ''
    def replace = replace_names ? "--replace-names ${replace_names}" : ''
    def samples = sample_names ? "--sample-names ${sample_names}" : ''
    """
    multiqc \\
        --force \\
        ${args} \\
        ${config} \\
        ${prefix} \\
        ${extra_config} \\
        ${logo} \\
        ${replace} \\
        ${samples} \\
        .
    """

    stub:
    """
    mkdir multiqc_data
    touch multiqc_data/.stub
    mkdir multiqc_plots
    touch multiqc_report.html
    """
}


// === FILE: modules/nf-core/s4pred/runmodel/main.nf ===
process S4PRED_RUNMODEL {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/s4pred:1.2.1--pyhdfd78af_1':
        'biocontainers/s4pred:1.2.1--pyhdfd78af_1' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("${prefix}"), emit: preds
    path "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args    = task.ext.args   ?: ''
    prefix      = task.ext.prefix ?: "${meta.id}"
    def VERSION = '1.2.1' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    mkdir ${prefix}

    run_model \\
        $args \\
        --threads $task.cpus \\
        --save-files \\
        --outdir ${prefix} \\
        ${fasta}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        s4pred: $VERSION
    END_VERSIONS
    """

    stub:
    prefix      = task.ext.prefix ?: "${meta.id}"
    def VERSION = '1.2.1' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    mkdir -p ${prefix}
    touch ${prefix}/test.ss2

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        s4pred: $VERSION
    END_VERSIONS
    """
}


// === FILE: modules/nf-core/seqfu/stats/main.nf ===
process SEQFU_STATS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/seqfu:1.22.3--hfd12232_2':
        'biocontainers/seqfu:1.22.3--hfd12232_2' }"

    input:
    // stats can get one or more fasta or fastq files
    tuple val(meta), path(files)

    output:
    tuple val(meta), path("*.tsv")    , emit: stats
    tuple val(meta), path("*_mqc.txt"), emit: multiqc
    path "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    seqfu \\
        stats \\
        $args \\
        --multiqc ${prefix}_mqc.txt \\
        $files > ${prefix}.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqfu: \$(seqfu version)
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo $args

    touch ${prefix}.tsv
    touch ${prefix}_mqc.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqfu: \$(seqfu version)
    END_VERSIONS
    """
}


// === FILE: modules/nf-core/seqkit/replace/main.nf ===
process SEQKIT_REPLACE {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/seqkit:2.9.0--h9ee0642_0'
        : 'biocontainers/seqkit:2.9.0--h9ee0642_0'}"

    input:
    tuple val(meta), path(fastx)

    output:
    tuple val(meta), path("*.fast*"), emit: fastx
    tuple val("${task.process}"), val('seqkit'), eval("seqkit version | sed 's/seqkit v//'"), emit: versions_seqkit, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def extension = "fastq"
    if ("${fastx}" ==~ /.+\.fasta|.+\.fasta.gz|.+\.fa|.+\.fa.gz|.+\.fas|.+\.fas.gz|.+\.fna|.+\.fna.gz|.+\.faa|.+\.faa.gz/) {
        extension = "fasta"
    }
    def isgz = ""
    if ("${fastx}" ==~ /.+\.gz/) {
        isgz = ".gz"
    }
    def endswith = task.ext.suffix ?: "${extension}${isgz}"
    """
    seqkit \\
        replace \\
        ${args} \\
        --threads ${task.cpus} \\
        -i ${fastx} \\
        -o ${prefix}.${endswith}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def extension = "fastq"
    if ("${fastx}" ==~ /.+\.fasta|.+\.fasta.gz|.+\.fa|.+\.fa.gz|.+\.fas|.+\.fas.gz|.+\.fna|.+\.fna.gz/) {
        extension = "fasta"
    }
    def endswith = task.ext.suffix ?: "${extension}.gz"

    """
    echo "" | gzip > ${prefix}.${endswith}
    """
}


// === FILE: modules/nf-core/seqkit/rmdup/main.nf ===
process SEQKIT_RMDUP {
    tag "$meta.id"
    label 'process_low'
    // File IO can be a bottleneck. See: https://bioinf.shenwei.me/seqkit/usage/#parallelization-of-cpu-intensive-jobs

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/seqkit:2.9.0--h9ee0642_0':
        'biocontainers/seqkit:2.9.0--h9ee0642_0' }"

    input:
    tuple val(meta), path(fastx)

    output:
    tuple val(meta), path("${prefix}.${extension}") , emit: fastx
    tuple val(meta), path("*.log")                  , emit: log
    path "versions.yml"                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args        = task.ext.args ?: ''
    prefix          = task.ext.prefix ?: "${meta.id}"
    extension       = "fastq"
    if ("$fastx" ==~ /.+\.fasta|.+\.fasta.gz|.+\.fa|.+\.fa.gz|.+\.fas|.+\.fas.gz|.+\.fna|.+\.fna.gz|.+\.fsa|.+\.fsa.gz/ ) {
        extension   = "fasta"
    }
    extension       = fastx.toString().endsWith('.gz') ? "${extension}.gz" : extension
    // SeqKit/rmdup takes care of compressing the output: https://bioinf.shenwei.me/seqkit/usage/#rmdup
    if("${prefix}.${extension}" == "$fastx") error "Input and output names are the same, use \"task.ext.prefix\" to disambiguate!"
    """
    seqkit \\
        rmdup \\
        --threads $task.cpus \\
        $args \\
        $fastx \\
        -o ${prefix}.${extension} \\
        2>| >(tee ${prefix}.log >&2)

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqkit: \$(seqkit version | cut -d' ' -f2)
    END_VERSIONS
    """

    stub:
    prefix          = task.ext.prefix ?: "${meta.id}"
    extension       = "fastq"
    if ("$fastx" ==~ /.+\.fasta|.+\.fasta.gz|.+\.fa|.+\.fa.gz|.+\.fas|.+\.fas.gz|.+\.fna|.+\.fna.gz|.+\.fsa|.+\.fsa.gz/ ) {
        extension   = "fasta"
    }
    extension = fastx.toString().endsWith('.gz') ? "${extension}.gz" : extension
    if("${prefix}.${extension}" == "$fastx") error "Input and output names are the same, use \"task.ext.prefix\" to disambiguate!"
    """
    touch ${prefix}.${extension}
    echo \\
        '[INFO] 0 duplicated records removed' \\
        > ${prefix}.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqkit: \$(seqkit version | cut -d' ' -f2)
    END_VERSIONS
    """
}


// === FILE: modules/nf-core/seqkit/seq/main.nf ===
process SEQKIT_SEQ {
    tag "${meta.id}"
    label 'process_low'
    // File IO can be a bottleneck. See: https://bioinf.shenwei.me/seqkit/usage/#parallelization-of-cpu-intensive-jobs

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/seqkit:2.9.0--h9ee0642_0'
        : 'biocontainers/seqkit:2.9.0--h9ee0642_0'}"

    input:
    tuple val(meta), path(fastx)

    output:
    tuple val(meta), path("${prefix}.*"), emit: fastx
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def extension = "fastq"
    if ("${fastx}" ==~ /.+\.fasta|.+\.fasta.gz|.+\.fa|.+\.fa.gz|.+\.fas|.+\.fas.gz|.+\.fna|.+\.fna.gz|.+\.fsa|.+\.fsa.gz/) {
        extension = "fasta"
    }
    extension = fastx.toString().endsWith('.gz') ? "${extension}.gz" : extension
    def call_gzip = extension.endsWith('.gz') ? "| gzip -c ${args2}" : ''
    if ("${prefix}.${extension}" == "${fastx}") {
        error("Input and output names are the same, use \"task.ext.prefix\" to disambiguate!")
    }
    """
    seqkit \\
        seq \\
        --threads ${task.cpus} \\
        ${args} \\
        ${fastx} \\
        ${call_gzip} \\
        > ${prefix}.${extension}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqkit: \$(seqkit version | cut -d' ' -f2)
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    def extension = "fastq"
    if ("${fastx}" ==~ /.+\.fasta|.+\.fasta.gz|.+\.fa|.+\.fa.gz|.+\.fas|.+\.fas.gz|.+\.fna|.+\.fna.gz|.+\.fsa|.+\.fsa.gz/) {
        extension = "fasta"
    }
    extension = fastx.toString().endsWith('.gz') ? "${extension}.gz" : extension
    if ("${prefix}.${extension}" == "${fastx}") {
        error("Input and output names are the same, use \"task.ext.prefix\" to disambiguate!")
    }
    """
    touch ${prefix}.${extension}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqkit: \$(seqkit version | cut -d' ' -f2)
    END_VERSIONS
    """
}


// === FILE: modules/nf-core/seqkit/stats/main.nf ===
process SEQKIT_STATS {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/seqkit:2.9.0--h9ee0642_0'
        : 'biocontainers/seqkit:2.9.0--h9ee0642_0'}"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.tsv"), emit: stats
    tuple val("${task.process}"), val('seqkit'), eval("seqkit version | sed 's/seqkit v//'"), emit: versions_seqkit, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: '--all'
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    seqkit stats \\
        --tabular \\
        --threads ${task.cpus} \\
        ${args} \\
        ${reads} > '${prefix}.tsv'
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.tsv
    """
}


// === FILE: modules/nf-core/untar/main.nf ===
process UNTAR {
    tag "${archive}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/52/52ccce28d2ab928ab862e25aae26314d69c8e38bd41ca9431c67ef05221348aa/data'
        : 'community.wave.seqera.io/library/coreutils_grep_gzip_lbzip2_pruned:838ba80435a629f8'}"

    input:
    tuple val(meta), path(archive)

    output:
    tuple val(meta), path("${prefix}"), emit: untar
    tuple val("${task.process}"), val('untar'), eval('tar --version 2>&1 | head -1 | sed "s/tar (GNU tar) //; s/ Copyright.*//"'), emit: versions_untar, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    prefix = task.ext.prefix ?: (meta.id ? "${meta.id}" : archive.baseName.toString().replaceFirst(/\.tar$/, ""))

    """
    mkdir ${prefix}

    ## Ensures --strip-components only applied when top level of tar contents is a directory
    ## If just files or multiple directories, place all in prefix
    if [[ \$(tar -taf ${archive} | grep -o -P "^.*?\\/" | uniq | wc -l) -eq 1 ]]; then
        tar \\
            -C ${prefix} --strip-components 1 \\
            -xavf \\
            ${args} \\
            ${archive} \\
            ${args2}
    else
        tar \\
            -C ${prefix} \\
            -xavf \\
            ${args} \\
            ${archive} \\
            ${args2}
    fi

    """

    stub:
    prefix = task.ext.prefix ?: (meta.id ? "${meta.id}" : archive.toString().replaceFirst(/\.[^\.]+(.gz)?$/, ""))
    """
    mkdir ${prefix}
    ## Dry-run untaring the archive to get the files and place all in prefix
    if [[ \$(tar -taf ${archive} | grep -o -P "^.*?\\/" | uniq | wc -l) -eq 1 ]]; then
        for i in `tar -tf ${archive}`;
        do
            if [[ \$(echo "\${i}" | grep -E "/\$") == "" ]];
            then
                touch \${i}
            else
                mkdir -p \${i}
            fi
        done
    else
        for i in `tar -tf ${archive}`;
        do
            if [[ \$(echo "\${i}" | grep -E "/\$") == "" ]];
            then
                touch ${prefix}/\${i}
            else
                mkdir -p ${prefix}/\${i}
            fi
        done
    fi
    """
}


// === FILE: subworkflows/local/domain_annotation/main.nf ===
include { ARIA2 as ARIA2_PFAM                 } from '../../../modules/nf-core/aria2/main'
include { ARIA2 as ARIA2_FUNFAM               } from '../../../modules/nf-core/aria2/main'
include { HMMER_HMMSEARCH as HMMSEARCH_PFAM   } from '../../../modules/nf-core/hmmer/hmmsearch/main'
include { HMMER_HMMSEARCH as HMMSEARCH_FUNFAM } from '../../../modules/nf-core/hmmer/hmmsearch/main'

workflow DOMAIN_ANNOTATION {
    take:
    ch_fasta            // channel: [ val(meta), [ fasta ] ]
    skip_pfam           // boolean
    pfam_db             // string, path to the pfam HMM database, if already exists
    pfam_latest_link    // string, path to the latest pfam HMM database, to download
    skip_funfam         // boolean
    funfam_db           // string, path to the funfam HMM database, if already exists
    funfam_latest_link  // string, path to the latest funfam HMM database, to download

    main:

    ch_versions       = channel.empty()
    ch_pfam_domains   = channel.empty()
    ch_funfam_domains = channel.empty()

    if (!skip_pfam) {
        if (!pfam_db) {
            ch_pfam_link = channel.of([ [ id: 'pfam' ], pfam_latest_link ])

            ARIA2_PFAM( ch_pfam_link )
            ch_versions = ch_versions.mix( ARIA2_PFAM.out.versions )
            ch_pfam_db = ARIA2_PFAM.out.downloaded_file
        } else {
            ch_pfam_db = channel.of([ [ id: 'pfam' ], pfam_db ])
        }

        ch_input_for_hmmsearch_pfam = ch_fasta
            .combine(ch_pfam_db)
            .map{ meta, seqs, _meta2, models -> [meta, models, seqs, false, false, true] }

        HMMSEARCH_PFAM( ch_input_for_hmmsearch_pfam )
        ch_versions = ch_versions.mix( HMMSEARCH_PFAM.out.versions.first() )
        ch_pfam_domains = HMMSEARCH_PFAM.out.domain_summary
    }

    if (!skip_funfam) {
        if (!funfam_db) {
            ch_funfam_link = channel.of([ [ id: 'funfam' ], funfam_latest_link ])

            ARIA2_FUNFAM( ch_funfam_link )
            ch_versions = ch_versions.mix( ARIA2_FUNFAM.out.versions )
            ch_funfam_db = ARIA2_FUNFAM.out.downloaded_file
        } else {
            ch_funfam_db = channel.of([ [ id: 'funfam' ], funfam_db ])
        }

        ch_input_for_hmmsearch_funfam = ch_fasta
            .combine(ch_funfam_db)
            .map{ meta, seqs, _meta2, models -> [meta, models, seqs, false, false, true] }

        HMMSEARCH_FUNFAM( ch_input_for_hmmsearch_funfam )
        ch_versions = ch_versions.mix( HMMSEARCH_FUNFAM.out.versions.first() )
        ch_funfam_domains = HMMSEARCH_FUNFAM.out.domain_summary
    }

    emit:
    pfam_domains   = ch_pfam_domains
    funfam_domains = ch_funfam_domains
    versions       = ch_versions
}


// === FILE: subworkflows/local/functional_annotation/main.nf ===
include { ARIA2        } from '../../../modules/nf-core/aria2/main'
include { UNTAR        } from '../../../modules/nf-core/untar/main'
include { INTERPROSCAN } from '../../../modules/nf-core/interproscan/main'

workflow FUNCTIONAL_ANNOTATION {
    take:
    ch_fasta            // channel: [ val(meta), [ fasta ] ]
    skip_interproscan   // boolean
    interproscan_db_url // string, url to download db
    interproscan_db     // string, existing db

    main:
    ch_interproscan_tsv = channel.empty()
    ch_versions         = channel.empty()

    if (!skip_interproscan) {
        if (interproscan_db) {
            ch_interproscan_db = channel.fromPath(interproscan_db).first()
        }
        else {
            ARIA2( [ [ id:'interproscan_db' ], interproscan_db_url ] )
            ch_versions = ch_versions.mix(ARIA2.out.versions.first())

            UNTAR( ARIA2.out.downloaded_file )
            ch_interproscan_db = UNTAR.out.untar.map{ f -> f[1] }
        }

        INTERPROSCAN( ch_fasta, ch_interproscan_db )
        ch_interproscan_tsv = ch_interproscan_tsv.mix(INTERPROSCAN.out.tsv)
    }

    emit:
    interproscan_tsv = ch_interproscan_tsv
    versions         = ch_versions         // channel: [ versions.yml ]
}


// === FILE: subworkflows/local/utils_nfcore_proteinannotator_pipeline/main.nf ===
//
// Subworkflow with functionality specific to the nf-core/proteinannotator pipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { UTILS_NFSCHEMA_PLUGIN     } from '../../nf-core/utils_nfschema_plugin'
include { paramsSummaryMap          } from 'plugin/nf-schema'
include { samplesheetToList         } from 'plugin/nf-schema'
include { paramsHelp                } from 'plugin/nf-schema'
include { completionEmail           } from '../../nf-core/utils_nfcore_pipeline'
include { completionSummary         } from '../../nf-core/utils_nfcore_pipeline'
include { imNotification            } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NFCORE_PIPELINE     } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NEXTFLOW_PIPELINE   } from '../../nf-core/utils_nextflow_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW TO INITIALISE PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_INITIALISATION {
    take:
    version           // boolean: Display version and exit
    validate_params   // boolean: Boolean whether to validate parameters against the schema at runtime
    _monochrome_logs  // boolean: Do not use coloured log outputs
    nextflow_cli_args //   array: List of positional nextflow CLI args
    outdir            //  string: The output directory where the results will be saved
    input             //  string: Path to input samplesheet
    help              // boolean: Display help message and exit
    help_full         // boolean: Show the full help message
    show_hidden       // boolean: Show hidden parameters in the help message

    main:

    ch_versions = channel.empty()

    //
    // Print version and exit if required and dump pipeline parameters to JSON file
    //
    UTILS_NEXTFLOW_PIPELINE(
        version,
        true,
        outdir,
        workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1,
    )

    //
    // Validate parameters and generate parameter summary to stdout
    //
    before_text = """
-\033[2m----------------------------------------------------\033[0m-
                                        \033[0;32m,--.\033[0;30m/\033[0;32m,-.\033[0m
\033[0;34m        ___     __   __   __   ___     \033[0;32m/,-._.--~\'\033[0m
\033[0;34m  |\\ | |__  __ /  ` /  \\ |__) |__         \033[0;33m}  {\033[0m
\033[0;34m  | \\| |       \\__, \\__/ |  \\ |___     \033[0;32m\\`-._,-`-,\033[0m
                                        \033[0;32m`._,._,\'\033[0m
\033[0;35m  nf-core/proteinannotator ${workflow.manifest.version}\033[0m
-\033[2m----------------------------------------------------\033[0m-
"""
    after_text = """${workflow.manifest.doi ? "\n* The pipeline\n" : ""}${workflow.manifest.doi.tokenize(",").collect { doi -> "    https://doi.org/${doi.trim().replace('https://doi.org/','')}"}.join("\n")}${workflow.manifest.doi ? "\n" : ""}
* The nf-core framework
    https://doi.org/10.1038/s41587-020-0439-x

* Software dependencies
    https://github.com/nf-core/proteinannotator/blob/master/CITATIONS.md
"""
    command = "nextflow run ${workflow.manifest.name} -profile <docker/singularity/.../institute> --input samplesheet.csv --outdir <OUTDIR>"

    UTILS_NFSCHEMA_PLUGIN (
        workflow,
        validate_params,
        null,
        help,
        help_full,
        show_hidden,
        before_text,
        after_text,
        command
    )

    //
    // Check config provided to the pipeline
    //
    UTILS_NFCORE_PIPELINE(
        nextflow_cli_args
    )

    //
    // Create channel from input file provided through input
    //

    channel.fromList(samplesheetToList(input, "${projectDir}/assets/schema_input.json"))
        .map { meta, fasta ->
            return [meta, fasta]
        }
        .map { samplesheet ->
            validateInputSamplesheet(samplesheet)
        }
        .set { ch_samplesheet }

    emit:
    samplesheet = ch_samplesheet
    versions = ch_versions
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW FOR PIPELINE COMPLETION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_COMPLETION {
    take:
    email //  string: email address
    email_on_fail //  string: email address sent on pipeline failure
    plaintext_email // boolean: Send plain-text email instead of HTML
    outdir //    path: Path to output directory where results will be published
    monochrome_logs // boolean: Disable ANSI colour codes in log output
    hook_url //  string: hook URL for notifications
    multiqc_report //  string: Path to MultiQC report

    main:
    summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    def multiqc_reports = multiqc_report.toList()

    //
    // Completion email and summary
    //
    workflow.onComplete {
        if (email || email_on_fail) {
            completionEmail(
                summary_params,
                email,
                email_on_fail,
                plaintext_email,
                outdir,
                monochrome_logs,
                multiqc_reports.getVal(),
            )
        }

        completionSummary(monochrome_logs)
        if (hook_url) {
            imNotification(summary_params, hook_url)
        }
    }

    workflow.onError {
        log.error("Pipeline failed. Please refer to troubleshooting docs: https://nf-co.re/docs/usage/troubleshooting")
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// Validate channels from input samplesheet
//
def validateInputSamplesheet(input) {
    // todo: implement samplesheet validation
    return input
}

//
// Generate methods description for MultiQC
//
def toolCitationText() {

    def quality_check_text = [
        "Amino acid sequence statistics were generated with SeqFu (Telatin et al. 2021).",
        params.skip_preprocessing ? "" : "Input sequences were preprocessed with SeqKit (gap trimming, length filtering, validation, duplicate removal) (Shen et al. 2024)."
    ].join(' ').trim()

    def domain_annotation_text = (params.skip_pfam && params.skip_funfam) ? "" : "Domains were annotated with hmmer/hmmsearch (Eddy et al. 2011)."

    def prediction_text = params.skip_s4pred ? "" : "Secondary structures were predicted via the s4pred software (Moffat et al. 2021)."

    def postprocessing_text = "Run statistics were reported using MultiQC (Ewels et al. 2016)."

    def citation_text = [
        quality_check_text,
        domain_annotation_text,
        prediction_text,
        postprocessing_text
    ].join(' ').trim()

    return citation_text
}

def toolBibliographyText() {
    def quality_check_text = [
        '<li>Telatin, A., Fariselli, P., & Birolo, G. (2021). SeqFu: a suite of utilities for the robust and reproducible manipulation of sequence files. Bioengineering, 8(5), 59. doi: <a href="https://doi.org/10.3390/bioengineering8050059">10.3390/bioengineering8050059</a></li>',
        params.skip_preprocessing ? '' : '<li>Shen, W., Sipos, B., & Zhao, L. (2024). SeqKit2: A Swiss army knife for sequence and alignment processing. Imeta, 3(3), e191. doi: <a href="https://doi.org/10.1002/imt2.191">10.1002/imt2.191</a></li>'
    ].join(' ').trim()

    def domain_annotation_text = (params.skip_pfam && params.skip_funfam) ? '' : '<li>Eddy, S. R. (2011). Accelerated profile HMM searches. PLoS computational biology, 7(10), e1002195. doi: <a href="https://doi.org/10.1371/journal.pcbi.1002195">10.1371/journal.pcbi.1002195</a></li>'

    def prediction_text = params.skip_s4pred ? '' : '<li>Moffat, L., & Jones, D. T. (2021). Increasing the accuracy of single sequence prediction methods using a deep semi-supervised learning framework. Bioinformatics, 37(21), 3744-3751. doi: <a href="https://doi.org/10.1093/bioinformatics/btab491">10.1093/bioinformatics/btab491</a></li>'

    def postprocessing_text = '<li>Ewels, P., Magnusson, M., Lundin, S., & Käller, M. (2016). MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics, 32(19), 3047–3048. doi: <a href="https://doi.org/10.1093/bioinformatics/btw354">10.1093/bioinformatics/btw354</a></li>'

    def reference_text = [
        quality_check_text,
        domain_annotation_text,
        prediction_text,
        postprocessing_text
    ].join(' ').trim()

    return reference_text
}

def methodsDescriptionText(mqc_methods_yaml) {
    // Convert  to a named map so can be used as with familiar NXF ${workflow} variable syntax in the MultiQC YML file
    def meta = [:]
    meta.workflow = workflow.toMap()
    meta["manifest_map"] = workflow.manifest.toMap()

    // Pipeline DOI
    if (meta.manifest_map.doi) {
        // Using a loop to handle multiple DOIs
        // Removing `https://doi.org/` to handle pipelines using DOIs vs DOI resolvers
        // Removing ` ` since the manifest.doi is a string and not a proper list
        def temp_doi_ref = ""
        def manifest_doi = meta.manifest_map.doi.tokenize(",")
        manifest_doi.each { doi_ref ->
            temp_doi_ref += "(doi: <a href=\'https://doi.org/${doi_ref.replace("https://doi.org/", "").replace(" ", "")}\'>${doi_ref.replace("https://doi.org/", "").replace(" ", "")}</a>), "
        }
        meta["doi_text"] = temp_doi_ref.substring(0, temp_doi_ref.length() - 2)
    }
    else {
        meta["doi_text"] = ""
    }
    meta["nodoi_text"] = meta.manifest_map.doi ? "" : "<li>If available, make sure to update the text to include the Zenodo DOI of version of the pipeline used. </li>"

    // Tool references
    meta["tool_citations"] = toolCitationText().replaceAll(", \\.", ".").replaceAll("\\. \\.", ".").replaceAll(", \\.", ".")
    meta["tool_bibliography"] = toolBibliographyText()


    def methods_text = mqc_methods_yaml.text

    def engine = new groovy.text.SimpleTemplateEngine()
    def description_html = engine.createTemplate(methods_text).make(meta)

    return description_html.toString()
}


// === FILE: subworkflows/nf-core/faa_seqfu_seqkit/main.nf ===
include { SEQFU_STATS as SEQFU_STATS_BEFORE } from '../../../modules/nf-core/seqfu/stats/main'
include { SEQKIT_SEQ                        } from '../../../modules/nf-core/seqkit/seq/main'
include { SEQKIT_RMDUP                      } from '../../../modules/nf-core/seqkit/rmdup/main'
include { SEQKIT_REPLACE                    } from '../../../modules/nf-core/seqkit/replace/main'
include { SEQFU_STATS as SEQFU_STATS_AFTER  } from '../../../modules/nf-core/seqfu/stats/main'

workflow FAA_SEQFU_SEQKIT {

    take:
    ch_fasta           // tuple val(meta), path(fasta)
    skip_preprocessing // boolean

    main:
    ch_multiqc_files = channel.empty()
    ch_versions      = channel.empty()

    SEQFU_STATS_BEFORE( ch_fasta )
    ch_multiqc_files = ch_multiqc_files.mix( SEQFU_STATS_BEFORE.out.multiqc )
    ch_versions      = ch_versions.mix( SEQFU_STATS_BEFORE.out.versions )

    if (!skip_preprocessing) {
        SEQKIT_SEQ( ch_fasta )
        ch_versions = ch_versions.mix( SEQKIT_SEQ.out.versions )

        SEQKIT_REPLACE( SEQKIT_SEQ.out.fastx )

        SEQKIT_RMDUP( SEQKIT_REPLACE.out.fastx )
        ch_fasta    = SEQKIT_RMDUP.out.fastx
        ch_versions = ch_versions.mix( SEQKIT_RMDUP.out.versions )

        SEQFU_STATS_AFTER( SEQKIT_RMDUP.out.fastx )
        ch_multiqc_files = ch_multiqc_files.mix( SEQFU_STATS_AFTER.out.multiqc )
        ch_versions      = ch_versions.mix( SEQFU_STATS_AFTER.out.versions )
    }

    emit:
    fasta         = ch_fasta
    multiqc_files = ch_multiqc_files
    versions      = ch_versions
}


// === FILE: subworkflows/nf-core/utils_nextflow_pipeline/main.nf ===
//
// Subworkflow with functionality that may be useful for any Nextflow pipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW DEFINITION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow UTILS_NEXTFLOW_PIPELINE {
    take:
    print_version        // boolean: print version
    dump_parameters      // boolean: dump parameters
    outdir               //    path: base directory used to publish pipeline results
    check_conda_channels // boolean: check conda channels

    main:

    //
    // Print workflow version and exit on --version
    //
    if (print_version) {
        log.info("${workflow.manifest.name} ${getWorkflowVersion()}")
        System.exit(0)
    }

    //
    // Dump pipeline parameters to a JSON file
    //
    if (dump_parameters && outdir) {
        dumpParametersToJSON(outdir)
    }

    //
    // When running with Conda, warn if channels have not been set-up appropriately
    //
    if (check_conda_channels) {
        checkCondaChannels()
    }

    emit:
    dummy_emit = true
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// Generate version string
//
def getWorkflowVersion() {
    def version_string = "" as String
    if (workflow.manifest.version) {
        def prefix_v = workflow.manifest.version[0] != 'v' ? 'v' : ''
        version_string += "${prefix_v}${workflow.manifest.version}"
    }

    if (workflow.commitId) {
        def git_shortsha = workflow.commitId.substring(0, 7)
        version_string += "-g${git_shortsha}"
    }

    return version_string
}

//
// Dump pipeline parameters to a JSON file
//
def dumpParametersToJSON(outdir) {
    def timestamp = new java.util.Date().format('yyyy-MM-dd_HH-mm-ss')
    def filename  = "params_${timestamp}.json"
    def temp_pf   = new File(workflow.launchDir.toString(), ".${filename}")
    def jsonStr   = groovy.json.JsonOutput.toJson(params)
    temp_pf.text  = groovy.json.JsonOutput.prettyPrint(jsonStr)

    nextflow.extension.FilesEx.copyTo(temp_pf.toPath(), "${outdir}/pipeline_info/params_${timestamp}.json")
    temp_pf.delete()
}

//
// When running with -profile conda, warn if channels have not been set-up appropriately
//
def checkCondaChannels() {
    def parser = new org.yaml.snakeyaml.Yaml()
    def channels = []
    try {
        def config = parser.load("conda config --show channels".execute().text)
        channels = config.channels
    }
    catch (NullPointerException e) {
        log.debug(e)
        log.warn("Could not verify conda channel configuration.")
        return null
    }
    catch (IOException e) {
        log.debug(e)
        log.warn("Could not verify conda channel configuration.")
        return null
    }

    // Check that all channels are present
    // This channel list is ordered by required channel priority.
    def required_channels_in_order = ['conda-forge', 'bioconda']
    def channels_missing = ((required_channels_in_order as Set) - (channels as Set)) as Boolean

    // Check that they are in the right order
    def channel_priority_violation = required_channels_in_order != channels.findAll { ch -> ch in required_channels_in_order }

    if (channels_missing | channel_priority_violation) {
        log.warn """\
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            There is a problem with your Conda configuration!
            You will need to set-up the conda-forge and bioconda channels correctly.
            Please refer to https://bioconda.github.io/
            The observed channel order is
            ${channels}
            but the following channel order is required:
            ${required_channels_in_order}
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
        """.stripIndent(true)
    }
}


// === FILE: subworkflows/nf-core/utils_nfcore_pipeline/main.nf ===
//
// Subworkflow with utility functions specific to the nf-core pipeline template
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW DEFINITION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow UTILS_NFCORE_PIPELINE {
    take:
    nextflow_cli_args

    main:
    valid_config = checkConfigProvided()
    checkProfileProvided(nextflow_cli_args)

    emit:
    valid_config
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
//  Warn if a -profile or Nextflow config has not been provided to run the pipeline
//
def checkConfigProvided() {
    def valid_config = true as Boolean
    if (workflow.profile == 'standard' && workflow.configFiles.size() <= 1) {
        log.warn(
            "[${workflow.manifest.name}] You are attempting to run the pipeline without any custom configuration!\n\n" + "This will be dependent on your local compute environment but can be achieved via one or more of the following:\n" + "   (1) Using an existing pipeline profile e.g. `-profile docker` or `-profile singularity`\n" + "   (2) Using an existing nf-core/configs for your Institution e.g. `-profile crick` or `-profile uppmax`\n" + "   (3) Using your own local custom config e.g. `-c /path/to/your/custom.config`\n\n" + "Please refer to the quick start section and usage docs for the pipeline.\n "
        )
        valid_config = false
    }
    return valid_config
}

//
// Exit pipeline if --profile contains spaces
//
def checkProfileProvided(nextflow_cli_args) {
    if (workflow.profile.endsWith(',')) {
        error(
            "The `-profile` option cannot end with a trailing comma, please remove it and re-run the pipeline!\n" + "HINT: A common mistake is to provide multiple values separated by spaces e.g. `-profile test, docker`.\n"
        )
    }
    if (nextflow_cli_args[0]) {
        log.warn(
            "nf-core pipelines do not accept positional arguments. The positional argument `${nextflow_cli_args[0]}` has been detected.\n" + "HINT: A common mistake is to provide multiple values separated by spaces e.g. `-profile test, docker`.\n"
        )
    }
}

//
// Generate workflow version string
//
def getWorkflowVersion() {
    def version_string = "" as String
    if (workflow.manifest.version) {
        def prefix_v = workflow.manifest.version[0] != 'v' ? 'v' : ''
        version_string += "${prefix_v}${workflow.manifest.version}"
    }

    if (workflow.commitId) {
        def git_shortsha = workflow.commitId.substring(0, 7)
        version_string += "-g${git_shortsha}"
    }

    return version_string
}

//
// Get software versions for pipeline
//
def processVersionsFromYAML(yaml_file) {
    def yaml = new org.yaml.snakeyaml.Yaml()
    def versions = yaml.load(yaml_file).collectEntries { k, v -> [k.tokenize(':')[-1], v] }
    return yaml.dumpAsMap(versions).trim()
}

//
// Get workflow version for pipeline
//
def workflowVersionToYAML() {
    return """
    Workflow:
        ${workflow.manifest.name}: ${getWorkflowVersion()}
        Nextflow: ${workflow.nextflow.version}
    """.stripIndent().trim()
}

//
// Get channel of software versions used in pipeline in YAML format
//
def softwareVersionsToYAML(ch_versions) {
    return ch_versions.unique().map { version -> processVersionsFromYAML(version) }.unique().mix(channel.of(workflowVersionToYAML()))
}

//
// Get workflow summary for MultiQC
//
def paramsSummaryMultiqc(summary_params) {
    def summary_section = ''
    summary_params
        .keySet()
        .each { group ->
            def group_params = summary_params.get(group)
            // This gets the parameters of that particular group
            if (group_params) {
                summary_section += "    <p style=\"font-size:110%\"><b>${group}</b></p>\n"
                summary_section += "    <dl class=\"dl-horizontal\">\n"
                group_params
                    .keySet()
                    .sort()
                    .each { param ->
                        summary_section += "        <dt>${param}</dt><dd><samp>${group_params.get(param) ?: '<span style=\"color:#999999;\">N/A</a>'}</samp></dd>\n"
                    }
                summary_section += "    </dl>\n"
            }
        }

    def yaml_file_text = "id: '${workflow.manifest.name.replace('/', '-')}-summary'\n" as String
    yaml_file_text     += "description: ' - this information is collected when the pipeline is started.'\n"
    yaml_file_text     += "section_name: '${workflow.manifest.name} Workflow Summary'\n"
    yaml_file_text     += "section_href: 'https://github.com/${workflow.manifest.name}'\n"
    yaml_file_text     += "plot_type: 'html'\n"
    yaml_file_text     += "data: |\n"
    yaml_file_text     += "${summary_section}"

    return yaml_file_text
}

//
// ANSII colours used for terminal logging
//
def logColours(monochrome_logs=true) {
    def colorcodes = [:] as Map

    // Reset / Meta
    colorcodes['reset']      = monochrome_logs ? '' : "\033[0m"
    colorcodes['bold']       = monochrome_logs ? '' : "\033[1m"
    colorcodes['dim']        = monochrome_logs ? '' : "\033[2m"
    colorcodes['underlined'] = monochrome_logs ? '' : "\033[4m"
    colorcodes['blink']      = monochrome_logs ? '' : "\033[5m"
    colorcodes['reverse']    = monochrome_logs ? '' : "\033[7m"
    colorcodes['hidden']     = monochrome_logs ? '' : "\033[8m"

    // Regular Colors
    colorcodes['black']  = monochrome_logs ? '' : "\033[0;30m"
    colorcodes['red']    = monochrome_logs ? '' : "\033[0;31m"
    colorcodes['green']  = monochrome_logs ? '' : "\033[0;32m"
    colorcodes['yellow'] = monochrome_logs ? '' : "\033[0;33m"
    colorcodes['blue']   = monochrome_logs ? '' : "\033[0;34m"
    colorcodes['purple'] = monochrome_logs ? '' : "\033[0;35m"
    colorcodes['cyan']   = monochrome_logs ? '' : "\033[0;36m"
    colorcodes['white']  = monochrome_logs ? '' : "\033[0;37m"

    // Bold
    colorcodes['bblack']  = monochrome_logs ? '' : "\033[1;30m"
    colorcodes['bred']    = monochrome_logs ? '' : "\033[1;31m"
    colorcodes['bgreen']  = monochrome_logs ? '' : "\033[1;32m"
    colorcodes['byellow'] = monochrome_logs ? '' : "\033[1;33m"
    colorcodes['bblue']   = monochrome_logs ? '' : "\033[1;34m"
    colorcodes['bpurple'] = monochrome_logs ? '' : "\033[1;35m"
    colorcodes['bcyan']   = monochrome_logs ? '' : "\033[1;36m"
    colorcodes['bwhite']  = monochrome_logs ? '' : "\033[1;37m"

    // Underline
    colorcodes['ublack']  = monochrome_logs ? '' : "\033[4;30m"
    colorcodes['ured']    = monochrome_logs ? '' : "\033[4;31m"
    colorcodes['ugreen']  = monochrome_logs ? '' : "\033[4;32m"
    colorcodes['uyellow'] = monochrome_logs ? '' : "\033[4;33m"
    colorcodes['ublue']   = monochrome_logs ? '' : "\033[4;34m"
    colorcodes['upurple'] = monochrome_logs ? '' : "\033[4;35m"
    colorcodes['ucyan']   = monochrome_logs ? '' : "\033[4;36m"
    colorcodes['uwhite']  = monochrome_logs ? '' : "\033[4;37m"

    // High Intensity
    colorcodes['iblack']  = monochrome_logs ? '' : "\033[0;90m"
    colorcodes['ired']    = monochrome_logs ? '' : "\033[0;91m"
    colorcodes['igreen']  = monochrome_logs ? '' : "\033[0;92m"
    colorcodes['iyellow'] = monochrome_logs ? '' : "\033[0;93m"
    colorcodes['iblue']   = monochrome_logs ? '' : "\033[0;94m"
    colorcodes['ipurple'] = monochrome_logs ? '' : "\033[0;95m"
    colorcodes['icyan']   = monochrome_logs ? '' : "\033[0;96m"
    colorcodes['iwhite']  = monochrome_logs ? '' : "\033[0;97m"

    // Bold High Intensity
    colorcodes['biblack']  = monochrome_logs ? '' : "\033[1;90m"
    colorcodes['bired']    = monochrome_logs ? '' : "\033[1;91m"
    colorcodes['bigreen']  = monochrome_logs ? '' : "\033[1;92m"
    colorcodes['biyellow'] = monochrome_logs ? '' : "\033[1;93m"
    colorcodes['biblue']   = monochrome_logs ? '' : "\033[1;94m"
    colorcodes['bipurple'] = monochrome_logs ? '' : "\033[1;95m"
    colorcodes['bicyan']   = monochrome_logs ? '' : "\033[1;96m"
    colorcodes['biwhite']  = monochrome_logs ? '' : "\033[1;97m"

    return colorcodes
}

// Return a single report from an object that may be a Path or List
//
def getSingleReport(multiqc_reports) {
    if (multiqc_reports instanceof Path) {
        return multiqc_reports
    } else if (multiqc_reports instanceof List) {
        if (multiqc_reports.size() == 0) {
            log.warn("[${workflow.manifest.name}] No reports found from process 'MULTIQC'")
            return null
        } else if (multiqc_reports.size() == 1) {
            return multiqc_reports.first()
        } else {
            log.warn("[${workflow.manifest.name}] Found multiple reports from process 'MULTIQC', will use only one")
            return multiqc_reports.first()
        }
    } else {
        return null
    }
}

//
// Construct and send completion email
//
def completionEmail(summary_params, email, email_on_fail, plaintext_email, outdir, monochrome_logs=true, multiqc_report=null) {

    // Set up the e-mail variables
    def subject = "[${workflow.manifest.name}] Successful: ${workflow.runName}"
    if (!workflow.success) {
        subject = "[${workflow.manifest.name}] FAILED: ${workflow.runName}"
    }

    def summary = [:]
    summary_params
        .keySet()
        .sort()
        .each { group ->
            summary << summary_params[group]
        }

    def misc_fields = [:]
    misc_fields['Date Started']              = workflow.start
    misc_fields['Date Completed']            = workflow.complete
    misc_fields['Pipeline script file path'] = workflow.scriptFile
    misc_fields['Pipeline script hash ID']   = workflow.scriptId
    if (workflow.repository) {
        misc_fields['Pipeline repository Git URL']    = workflow.repository
    }
    if (workflow.commitId) {
        misc_fields['Pipeline repository Git Commit'] = workflow.commitId
    }
    if (workflow.revision) {
        misc_fields['Pipeline Git branch/tag']        = workflow.revision
    }
    misc_fields['Nextflow Version']          = workflow.nextflow.version
    misc_fields['Nextflow Build']            = workflow.nextflow.build
    misc_fields['Nextflow Compile Timestamp'] = workflow.nextflow.timestamp

    def email_fields = [:]
    email_fields['version']      = getWorkflowVersion()
    email_fields['runName']      = workflow.runName
    email_fields['success']      = workflow.success
    email_fields['dateComplete'] = workflow.complete
    email_fields['duration']     = workflow.duration
    email_fields['exitStatus']   = workflow.exitStatus
    email_fields['errorMessage'] = (workflow.errorMessage ?: 'None')
    email_fields['errorReport']  = (workflow.errorReport ?: 'None')
    email_fields['commandLine']  = workflow.commandLine
    email_fields['projectDir']   = workflow.projectDir
    email_fields['summary']      = summary << misc_fields

    // On success try attach the multiqc report
    def mqc_report = getSingleReport(multiqc_report)

    // Check if we are only sending emails on failure
    def email_address = email
    if (!email && email_on_fail && !workflow.success) {
        email_address = email_on_fail
    }

    // Render the TXT template
    def engine       = new groovy.text.GStringTemplateEngine()
    def tf           = new File("${workflow.projectDir}/assets/email_template.txt")
    def txt_template = engine.createTemplate(tf).make(email_fields)
    def email_txt    = txt_template.toString()

    // Render the HTML template
    def hf            = new File("${workflow.projectDir}/assets/email_template.html")
    def html_template = engine.createTemplate(hf).make(email_fields)
    def email_html    = html_template.toString()

    // Render the sendmail template
    def max_multiqc_email_size = (params.containsKey('max_multiqc_email_size') ? params.max_multiqc_email_size : 0) as MemoryUnit
    def smail_fields           = [email: email_address, subject: subject, email_txt: email_txt, email_html: email_html, projectDir: "${workflow.projectDir}", mqcFile: mqc_report, mqcMaxSize: max_multiqc_email_size.toBytes()]
    def sf                     = new File("${workflow.projectDir}/assets/sendmail_template.txt")
    def sendmail_template      = engine.createTemplate(sf).make(smail_fields)
    def sendmail_html          = sendmail_template.toString()

    // Send the HTML e-mail
    def colors = logColours(monochrome_logs) as Map
    if (email_address) {
        try {
            if (plaintext_email) {
                new org.codehaus.groovy.GroovyException('Send plaintext e-mail, not HTML')
            }
            // Try to send HTML e-mail using sendmail
            def sendmail_tf = new File(workflow.launchDir.toString(), ".sendmail_tmp.html")
            sendmail_tf.withWriter { w -> w << sendmail_html }
            ['sendmail', '-t'].execute() << sendmail_html
            log.info("-${colors.purple}[${workflow.manifest.name}]${colors.green} Sent summary e-mail to ${email_address} (sendmail)-")
        }
        catch (Exception msg) {
            log.debug(msg.toString())
            log.debug("Trying with mail instead of sendmail")
            // Catch failures and try with plaintext
            def mail_cmd = ['mail', '-s', subject, '--content-type=text/html', email_address]
            mail_cmd.execute() << email_html
            log.info("-${colors.purple}[${workflow.manifest.name}]${colors.green} Sent summary e-mail to ${email_address} (mail)-")
        }
    }

    // Write summary e-mail HTML to a file
    def output_hf = new File(workflow.launchDir.toString(), ".pipeline_report.html")
    output_hf.withWriter { w -> w << email_html }
    nextflow.extension.FilesEx.copyTo(output_hf.toPath(), "${outdir}/pipeline_info/pipeline_report.html")
    output_hf.delete()

    // Write summary e-mail TXT to a file
    def output_tf = new File(workflow.launchDir.toString(), ".pipeline_report.txt")
    output_tf.withWriter { w -> w << email_txt }
    nextflow.extension.FilesEx.copyTo(output_tf.toPath(), "${outdir}/pipeline_info/pipeline_report.txt")
    output_tf.delete()
}

//
// Print pipeline summary on completion
//
def completionSummary(monochrome_logs=true) {
    def colors = logColours(monochrome_logs) as Map
    if (workflow.success) {
        if (workflow.stats.ignoredCount == 0) {
            log.info("-${colors.purple}[${workflow.manifest.name}]${colors.green} Pipeline completed successfully${colors.reset}-")
        }
        else {
            log.info("-${colors.purple}[${workflow.manifest.name}]${colors.yellow} Pipeline completed successfully, but with errored process(es) ${colors.reset}-")
        }
    }
    else {
        log.info("-${colors.purple}[${workflow.manifest.name}]${colors.red} Pipeline completed with errors${colors.reset}-")
    }
}

//
// Construct and send a notification to a web server as JSON e.g. Microsoft Teams and Slack
//
def imNotification(summary_params, hook_url) {
    def summary = [:]
    summary_params
        .keySet()
        .sort()
        .each { group ->
            summary << summary_params[group]
        }

    def misc_fields = [:]
    misc_fields['start']          = workflow.start
    misc_fields['complete']       = workflow.complete
    misc_fields['scriptfile']     = workflow.scriptFile
    misc_fields['scriptid']       = workflow.scriptId
    if (workflow.repository) {
        misc_fields['repository'] = workflow.repository
    }
    if (workflow.commitId) {
        misc_fields['commitid']   = workflow.commitId
    }
    if (workflow.revision) {
        misc_fields['revision']   = workflow.revision
    }
    misc_fields['nxf_version']    = workflow.nextflow.version
    misc_fields['nxf_build']      = workflow.nextflow.build
    misc_fields['nxf_timestamp']  = workflow.nextflow.timestamp

    def msg_fields = [:]
    msg_fields['version']      = getWorkflowVersion()
    msg_fields['runName']      = workflow.runName
    msg_fields['success']      = workflow.success
    msg_fields['dateComplete'] = workflow.complete
    msg_fields['duration']     = workflow.duration
    msg_fields['exitStatus']   = workflow.exitStatus
    msg_fields['errorMessage'] = (workflow.errorMessage ?: 'None')
    msg_fields['errorReport']  = (workflow.errorReport ?: 'None')
    msg_fields['commandLine']  = workflow.commandLine.replaceFirst(/ +--hook_url +[^ ]+/, "")
    msg_fields['projectDir']   = workflow.projectDir
    msg_fields['summary']      = summary << misc_fields

    // Render the JSON template
    def engine       = new groovy.text.GStringTemplateEngine()
    // Different JSON depending on the service provider
    // Defaults to "Adaptive Cards" (https://adaptivecards.io), except Slack which has its own format
    def json_path     = hook_url.contains("hooks.slack.com") ? "slackreport.json" : "adaptivecard.json"
    def hf            = new File("${workflow.projectDir}/assets/${json_path}")
    def json_template = engine.createTemplate(hf).make(msg_fields)
    def json_message  = json_template.toString()

    // POST
    def post = new URL(hook_url).openConnection()
    post.setRequestMethod("POST")
    post.setDoOutput(true)
    post.setRequestProperty("Content-Type", "application/json")
    post.getOutputStream().write(json_message.getBytes("UTF-8"))
    def postRC = post.getResponseCode()
    if (!postRC.equals(200)) {
        log.warn(post.getErrorStream().getText())
    }
}


// === FILE: subworkflows/nf-core/utils_nfschema_plugin/main.nf ===
//
// Subworkflow that uses the nf-schema plugin to validate parameters and render the parameter summary
//

include { paramsSummaryLog   } from 'plugin/nf-schema'
include { validateParameters } from 'plugin/nf-schema'
include { paramsHelp         } from 'plugin/nf-schema'

workflow UTILS_NFSCHEMA_PLUGIN {

    take:
    input_workflow      // workflow: the workflow object used by nf-schema to get metadata from the workflow
    validate_params     // boolean:  validate the parameters
    parameters_schema   // string:   path to the parameters JSON schema.
                        //           this has to be the same as the schema given to `validation.parametersSchema`
                        //           when this input is empty it will automatically use the configured schema or
                        //           "${projectDir}/nextflow_schema.json" as default. This input should not be empty
                        //           for meta pipelines
    help                // boolean:  show help message
    help_full           // boolean:  show full help message
    show_hidden         // boolean:  show hidden parameters in help message
    before_text         // string:   text to show before the help message and parameters summary
    after_text          // string:   text to show after the help message and parameters summary
    command             // string:   an example command of the pipeline

    main:

    if(help || help_full) {
        help_options = [
            beforeText: before_text,
            afterText: after_text,
            command: command,
            showHidden: show_hidden,
            fullHelp: help_full,
        ]
        if(parameters_schema) {
            help_options << [parametersSchema: parameters_schema]
        }
        log.info paramsHelp(
            help_options,
            (params.help instanceof String && params.help != "true") ? params.help : "",
        )
        exit 0
    }

    //
    // Print parameter summary to stdout. This will display the parameters
    // that differ from the default given in the JSON schema
    //

    summary_options = [:]
    if(parameters_schema) {
        summary_options << [parametersSchema: parameters_schema]
    }
    log.info before_text
    log.info paramsSummaryLog(summary_options, input_workflow)
    log.info after_text

    //
    // Validate the parameters using nextflow_schema.json or the schema
    // given via the validation.parametersSchema configuration option
    //
    if(validate_params) {
        validateOptions = [:]
        if(parameters_schema) {
            validateOptions << [parametersSchema: parameters_schema]
        }
        validateParameters(validateOptions)
    }

    emit:
    dummy_emit = true
}


// === FILE: workflows/proteinannotator.nf ===
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { FAA_SEQFU_SEQKIT       } from '../subworkflows/nf-core/faa_seqfu_seqkit/main'
include { DOMAIN_ANNOTATION      } from '../subworkflows/local/domain_annotation'
include { FUNCTIONAL_ANNOTATION  } from '../subworkflows/local/functional_annotation'
include { S4PRED_RUNMODEL        } from '../modules/nf-core/s4pred/runmodel/main'
include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_proteinannotator_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PROTEINANNOTATOR {
    take:
    ch_samplesheet      // channel: samplesheet read in from --input
    skip_preprocessing  // boolean
    skip_pfam           // boolean
    pfam_db             // string, path to the pfam HMM database, if already exists
    pfam_latest_link    // string, path to the latest pfam HMM database, to download
    skip_funfam         // boolean
    funfam_db           // string, path to the pfam HMM database, if already exists
    funfam_latest_link  // string, path to the latest pfam HMM database, to download
    skip_interproscan   // boolean
    interproscan_db_url // string, url to download db
    interproscan_db     // string, existing db
    skip_s4pred         // boolean

    main:

    ch_versions      = channel.empty()
    ch_multiqc_files = channel.empty()

    FAA_SEQFU_SEQKIT( ch_samplesheet, skip_preprocessing )
    ch_versions = ch_versions.mix( FAA_SEQFU_SEQKIT.out.versions )

    DOMAIN_ANNOTATION (
        FAA_SEQFU_SEQKIT.out.fasta,
        skip_pfam,
        pfam_db,
        pfam_latest_link,
        skip_funfam,
        funfam_db,
        funfam_latest_link
    )
    ch_versions = ch_versions.mix( DOMAIN_ANNOTATION.out.versions )

    FUNCTIONAL_ANNOTATION (
        FAA_SEQFU_SEQKIT.out.fasta,
        skip_interproscan,
        interproscan_db_url,
        interproscan_db
    )
    ch_versions = ch_versions.mix( FUNCTIONAL_ANNOTATION.out.versions )

    if (!skip_s4pred) {
        S4PRED_RUNMODEL( FAA_SEQFU_SEQKIT.out.fasta )
        ch_versions = ch_versions.mix( S4PRED_RUNMODEL.out.versions.first() )
    }

    //
    // Collate and save software versions
    //
    def topic_versions = channel.topic("versions")
        .distinct()
        .branch { entry ->
            versions_file: entry instanceof Path
            versions_tuple: true
        }

    def topic_versions_string = topic_versions.versions_tuple
        .map { process, tool, version ->
            [ process[process.lastIndexOf(':')+1..-1], "  ${tool}: ${version}" ]
        }
        .groupTuple(by:0)
        .map { process, tool_versions ->
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }

    softwareVersionsToYAML(ch_versions.mix(topic_versions.versions_file))
        .mix(topic_versions_string)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_' + 'proteinannotator_software_' + 'mqc_' + 'versions.yml',
            sort: true,
            newLine: true,
        )
        .set { ch_collated_versions }

    //
    // MODULE: MultiQC
    //
    ch_multiqc_config        = channel.fromPath(
        "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ?
        channel.fromPath(params.multiqc_config, checkIfExists: true) :
        channel.empty()
    ch_multiqc_logo          = params.multiqc_logo ?
        channel.fromPath(params.multiqc_logo, checkIfExists: true) :
        channel.empty()

    summary_params      = paramsSummaryMap(
        workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
        file(params.multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description))

    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true,
        )
    )

    ch_multiqc_files = ch_multiqc_files.mix(FAA_SEQFU_SEQKIT.out.multiqc_files.collect{ f -> f[1] }.ifEmpty([]))

    MULTIQC(
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        [],
    )

    emit:
    multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions = ch_versions // channel: [ path(versions.yml) ]
}
