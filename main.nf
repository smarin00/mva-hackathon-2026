#!/usr/bin/env nextflow
/*
 * MVA Hackathon 2026 — singleton WGS pipeline
 * Nextflow >= 26.04 (strict syntax). All defaults live in nextflow.config.
 *
 *   Full run (FASTQ -> CRAM -> calls -> annotation -> ranking):
 *     nextflow run main.nf -params-file params.yaml
 *
 *   Vendor-VCF only (fast: Exomiser + VEP on the provided VCF, no alignment):
 *     nextflow run main.nf -params-file params.yaml --mode vcf
 *
 *   Validate wiring without running anything:
 *     nextflow lint main.nf nextflow.config
 *     nextflow run main.nf -params-file params.yaml -preview
 *
 * Tool placement:
 *   native (host)  : fastp, bwa, samtools, bcftools, mosdepth, delly, ExpansionHunter,
 *                    Exomiser (java), whatshap / multiqc / aneuploidy (uv venv)
 *   containers     : DeepVariant, GATK (Mutect2), VEP 116 (local .sif)
 */


// ===========================================================================
// Alignment
// ===========================================================================

/*
 * fastp -> bwa mem -> fixmate -> sort, streamed in one pipe per lane, so no
 * trimmed FASTQ is ever written to disk.
 */
process ALIGN_LANE {
    tag "${meta.lane}"
    label 'high'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.lane}.bam"), emit: bam
    path "${meta.lane}.fastp.json",              emit: json

    script:
    def rg = "@RG\\tID:${params.sample_id}.${meta.lane}\\tSM:${params.sample_id}\\tPL:ILLUMINA\\tLB:${params.sample_id}\\tPU:${meta.lane}"
    """
    fastp -i ${reads[0]} -I ${reads[1]} --stdout \\
        --detect_adapter_for_pe -w 3 \\
        -j ${meta.lane}.fastp.json -h ${meta.lane}.fastp.html \\
    | bwa mem -p -t ${task.cpus} -R '${rg}' ${params.fasta} - \\
    | samtools fixmate -m -u - - \\
    | samtools sort -@ 4 -m 1500M -T sort_${meta.lane} -o ${meta.lane}.bam -
    """
}

/*
 * Merge lanes, mark duplicates, write CRAM directly (half the size of BAM).
 */
process MARKDUP_CRAM {
    tag "${params.sample_id}"
    label 'high'
    publishDir "${params.outdir}/alignment", mode: 'link'

    input:
    path bams

    output:
    tuple path("${params.sample_id}.cram"), path("${params.sample_id}.cram.crai"), emit: cram
    path "${params.sample_id}.markdup.txt",                                         emit: metrics

    script:
    """
    samtools merge -@ ${task.cpus} -u - ${bams} \\
    | samtools markdup -@ ${task.cpus} -s -f ${params.sample_id}.markdup.txt \\
        --reference ${params.fasta} -O cram - ${params.sample_id}.cram
    samtools index -@ ${task.cpus} ${params.sample_id}.cram
    """
}

process ALIGN_QC {
    tag "${params.sample_id}"
    label 'medium'
    publishDir "${params.outdir}/qc", mode: 'link'

    input:
    tuple path(cram), path(crai)

    output:
    path "${params.sample_id}.{flagstat,stats}.txt", emit: reports

    script:
    """
    samtools flagstat -@ ${task.cpus} ${cram} > ${params.sample_id}.flagstat.txt
    samtools stats -@ ${task.cpus} --reference ${params.fasta} ${cram} > ${params.sample_id}.stats.txt
    """
}

process MOSDEPTH {
    tag "${params.sample_id}"
    label 'medium'
    publishDir "${params.outdir}/qc/mosdepth", mode: 'link'

    input:
    tuple path(cram), path(crai)

    output:
    path "${params.sample_id}.regions.bed.gz",       emit: bins
    path "${params.sample_id}.mosdepth.summary.txt", emit: summary
    path "${params.sample_id}.*",                    emit: all

    script:
    """
    mosdepth -t 4 -n --fast-mode -f ${params.fasta} \\
        --by ${params.bin_size} ${params.sample_id} ${cram}
    """
}


// ===========================================================================
// Variant calling
// ===========================================================================

process DEEPVARIANT {
    tag "${params.sample_id}"
    label 'high'
    publishDir "${params.outdir}/variants", mode: 'link'

    input:
    tuple path(cram), path(crai)

    output:
    tuple path("${params.sample_id}.dv.vcf.gz"), path("${params.sample_id}.dv.vcf.gz.tbi"), emit: vcf

    script:
    // Restrict CALLING to primary chromosomes. The analysis-set reference carries
    // ~2,385 hs38d1 decoy, unplaced and unlocalised contigs, and make_examples
    // sweeps every one of them whether or not reads are present -- that sweep, not
    // the read volume, dominated the 1 h 58 m smoke test on empty input.
    //
    // This does not discard the decoys' benefit. Alignment in ALIGN_LANE used the
    // full reference, so decoys still act as a sink for reads that would otherwise
    // mismap onto the primary assembly. --regions changes only which windows are
    // genotyped, not what bwa had available to align against.
    //
    // Excluded and not missed: chrEBV and the unplaced scaffolds carry no gene
    // annotation in the VEP cache or the Exomiser database, so variants there
    // could never reach a ranked candidate list. chrM is kept for completeness,
    // but the WGS model is not built for heteroplasmy -- chrM calls are not
    // interpreted.
    """
    awk '\$1 ~ /^chr([0-9]+|X|Y|M)\$/ { print \$1"\\t0\\t"\$2 }' \\
        ${params.fasta}.fai > primary.bed
    echo "DeepVariant regions: \$(wc -l < primary.bed) primary contigs" >&2

    mkdir -p dv_tmp
    /opt/deepvariant/bin/run_deepvariant \\
        --model_type=WGS \\
        --ref=${params.fasta} \\
        --reads=${cram} \\
        --regions primary.bed \\
        --sample_name=${params.sample_id} \\
        --output_vcf=${params.sample_id}.dv.vcf.gz \\
        --num_shards=${task.cpus} \\
        --intermediate_results_dir=dv_tmp
    rm -rf dv_tmp
    """
}

/*
 * PASS only, primary contigs only, split multiallelics, left-align. Used for
 * both the vendor VCF and the DeepVariant callset; `label` keeps their outputs
 * apart.
 *
 * The primary-contig restriction is applied HERE, not inherited from the caller.
 * DeepVariant is now called on primary contigs only; the vendor VCF was called
 * genome-wide. Without this filter the two callsets would not be comparable, and
 * any "the vendor found something DeepVariant missed" claim would be confounded
 * by region, not by caller. Filtering both to the same intervals makes the
 * comparison honest.
 */
process NORMALIZE {
    tag "${label}"
    label 'low'
    publishDir "${params.outdir}/variants", mode: 'link'

    input:
    tuple val(label), path(vcf), path(tbi)

    output:
    tuple val(label), path("${label}.norm.vcf.gz"), path("${label}.norm.vcf.gz.tbi"), emit: vcf

    script:
    """
    awk '\$1 ~ /^chr([0-9]+|X|Y|M)\$/ { print \$1"\\t0\\t"\$2 }' \\
        ${params.fasta}.fai > primary.bed

    n_in=\$(bcftools index -n ${vcf})
    bcftools view -f PASS -R primary.bed -O u ${vcf} \\
    | bcftools norm -m -any -f ${params.fasta} -c w -O z -o ${label}.norm.vcf.gz -
    tabix -p vcf ${label}.norm.vcf.gz
    n_out=\$(bcftools index -n ${label}.norm.vcf.gz)

    # Audit trail: how much the PASS + primary-contig filter removed.
    echo "${label}: \${n_in} records in, \${n_out} after PASS + primary-contig filter" >&2
    """
}

/*
 * Mosaic-aware calling, restricted to the candidate panel. Germline callers
 * assume ~0.5 allele fraction; a mosaic allele can sit well below that.
 * Genome-wide Mutect2 would take days on this box; the panel takes minutes.
 */
process MUTECT2_PANEL {
    tag "${params.sample_id}"
    label 'medium'
    publishDir "${params.outdir}/variants", mode: 'link'

    input:
    tuple path(cram), path(crai)
    path panel

    output:
    tuple path("${params.sample_id}.mosaic.vcf.gz"), path("${params.sample_id}.mosaic.vcf.gz.tbi"), emit: vcf

    script:
    """
    gatk --java-options "-Xmx8g" Mutect2 \\
        -R ${params.fasta} -I ${cram} -L ${panel} \\
        --native-pair-hmm-threads ${task.cpus} \\
        -O raw.vcf.gz
    gatk --java-options "-Xmx4g" FilterMutectCalls \\
        -R ${params.fasta} -V raw.vcf.gz \\
        -O ${params.sample_id}.mosaic.vcf.gz
    """
}

process DELLY {
    tag "${params.sample_id}"
    label 'sv'
    publishDir "${params.outdir}/variants", mode: 'link'

    input:
    tuple path(cram), path(crai)

    output:
    tuple path("${params.sample_id}.delly.vcf.gz"), path("${params.sample_id}.delly.vcf.gz.tbi"), emit: vcf

    script:
    // delly 2.x replaced the `call` subcommand with `sr` (short-read discovery)
    """
    delly sr -g ${params.fasta} -o delly.bcf ${cram}
    bcftools view -O z -o ${params.sample_id}.delly.vcf.gz delly.bcf
    tabix -p vcf ${params.sample_id}.delly.vcf.gz
    """
}

process EXPANSIONHUNTER {
    tag "${params.sample_id}"
    label 'medium'
    publishDir "${params.outdir}/variants", mode: 'link'

    input:
    tuple path(cram), path(crai)
    path catalog

    output:
    path "${params.sample_id}.eh.*", emit: all

    script:
    // ExpansionHunter 5.0.0 bundles an htslib that cannot decode CRAM written
    // by recent samtools ("Failure to decode slice"). It only inspects catalog
    // loci, so extract those regions (and their off-target regions) to a small
    // BAM first. The bundled GRCh38 catalog names contigs without "chr"; they
    // are renamed to match the reference, or every locus silently gets no reads.
    """
    python3 - <<'PY'
    import json
    contigs = {l.split("\\t")[0] for l in open("${params.fasta}.fai")}

    def fix(region):
        chrom, span = region.rsplit(":", 1)
        for c in (chrom, "chr" + chrom, "chrM" if chrom in ("MT", "M") else None):
            if c in contigs:
                return f"{c}:{span}"
        raise SystemExit(f"catalog contig {chrom!r} not in the reference")

    catalog, regions = json.load(open("${catalog}")), set()
    for locus in catalog:
        for key in ("ReferenceRegion", "OfftargetRegions"):
            r = locus.get(key)
            if not r:
                continue
            fixed = [fix(x) for x in r] if isinstance(r, list) else fix(r)
            locus[key] = fixed
            for x in (fixed if isinstance(fixed, list) else [fixed]):
                chrom, span = x.rsplit(":", 1)
                start, end = span.split("-")
                regions.add((chrom, max(0, int(start) - 2000), int(end) + 2000))
    json.dump(catalog, open("catalog.ref.json", "w"), indent=1)
    with open("catalog_regions.bed", "w") as fh:
        for chrom, start, end in sorted(regions):
            fh.write(f"{chrom}\\t{start}\\t{end}\\n")
    PY

    samtools view -b -M -T ${params.fasta} -L catalog_regions.bed \\
        -o eh_input.bam ${cram}
    samtools index eh_input.bam
    n=\$(samtools view -c eh_input.bam)
    echo "reads extracted for ExpansionHunter: \$n" >&2
    [ "\$n" -gt 0 ] || { echo "ERROR: no reads in catalog regions (contig names?)" >&2; exit 1; }

    # Sex decides how X-linked repeats (FMR1, AR, AFF2) are genotyped. 'auto':
    # chrY vs chrX depth per Mb from the CRAM index (male ~0.5-1, female ~0).
    sex=${params.sex}
    if [ "\$sex" = auto ]; then
        r=\$(samtools idxstats ${cram} | awk '\$1=="chrX"||\$1=="X"{x=\$3/\$2} \$1=="chrY"||\$1=="Y"{y=\$3/\$2} END{print (x>0) ? y/x : 0}')
        sex=\$(awk -v r=\$r 'BEGIN{print (r>0.3) ? "male" : (r<0.1) ? "female" : "unclear"}')
        echo "sex from chrY/chrX depth ratio \$r: \$sex" >&2
        if [ "\$sex" = unclear ]; then
            echo "WARNING: sex unclear from depth; set sex: in the case file. Using female." >&2
            sex=female
        fi
    fi

    ExpansionHunter \\
        --reads eh_input.bam \\
        --reference ${params.fasta} \\
        --variant-catalog catalog.ref.json \\
        --sex \$sex \\
        --threads ${task.cpus} \\
        --output-prefix ${params.sample_id}.eh

    rm -f eh_input.bam eh_input.bam.bai
    """
}


// ===========================================================================
// Aneuploidy (the differentiating analysis)
// ===========================================================================

/*
 * Mosaic aneuploidy profiling: GC-corrected per-bin depth plus B-allele
 * frequency at heterozygous sites, summarised per chromosome as an estimated
 * fraction of cells carrying a gain or loss. Reports its own detection limit.
 *
 * Replaces CNVpytor, which targets focal CNV calling rather than
 * chromosome-wide mosaic shifts, and whose 1.3.2 resource downloader is broken.
 */
process ANEUPLOIDY {
    tag "${params.sample_id}"
    label 'medium'
    publishDir "${params.outdir}/aneuploidy", mode: 'copy'

    input:
    path bins
    tuple path(vcf), path(tbi)

    output:
    path "aneuploidy_karyotype.pdf",   emit: figure
    path "aneuploidy_karyotype.png",   emit: figure_png
    path "mosaic_fractions.tsv",       emit: fractions
    path "mosaic_fractions.tsv.limit", emit: limit
    path "aneuploidy.log",             emit: log

    script:
    """
    aneuploidy_plot.py \\
        --bins ${bins} \\
        --fasta ${params.fasta} \\
        --vcf ${vcf} \\
        --bin-size ${params.bin_size} \\
        --out-figure aneuploidy_karyotype.pdf \\
        --out-table mosaic_fractions.tsv 2> aneuploidy.log
    cat aneuploidy.log >&2
    """
}


// ===========================================================================
// Annotation, ranking, phasing
// ===========================================================================

/*
 * VEP 116 from the pinned local image. --everything carries gnomAD v4
 * frequencies, ClinVar significance and MANE. Plugins switch on as soon as
 * their data paths are set in params.
 */
process VEP {
    tag "${label}"
    label 'high'
    publishDir "${params.outdir}/annotation", mode: 'link'

    input:
    tuple val(label), path(vcf), path(tbi)

    output:
    tuple val(label), path("${label}.vep.vcf.gz"), path("${label}.vep.vcf.gz.tbi"), emit: vcf

    script:
    def plugins = []
    if (params.vep_plugins_dir) {
        plugins.add("--dir_plugins ${params.vep_plugins_dir}")
    }
    if (params.spliceai_snv && params.spliceai_indel) {
        plugins.add("--plugin SpliceAI,snv=${params.spliceai_snv},indel=${params.spliceai_indel}")
    }
    if (params.alphamissense) {
        plugins.add("--plugin AlphaMissense,file=${params.alphamissense}")
    }
    if (params.clinvar) {
        plugins.add("--custom file=${params.clinvar},short_name=ClinVar,format=vcf,type=exact,coords=0,fields=CLNSIG%CLNREVSTAT%CLNDN")
    }
    """
    vep \\
        --input_file ${vcf} \\
        --output_file ${label}.vep.vcf.gz \\
        --format vcf --vcf --compress_output bgzip \\
        --force_overwrite --no_stats \\
        --offline --cache --dir_cache ${params.vep_cache} \\
        --cache_version ${params.vep_cache_version} \\
        --assembly GRCh38 --fasta ${params.fasta} \\
        --everything \\
        --fork ${task.cpus} --buffer_size 5000 \\
        ${plugins.join(' ')}
    tabix -p vcf ${label}.vep.vcf.gz
    """
}

/*
 * Exomiser 15.1.0, native jar. Memory: ~16 GB for a whole genome in full mode,
 * so the heap is capped at params.exomiser_xmx with G1GC, and the 'exomiser'
 * label reserves enough of the 30 GB executor budget that nothing heavy runs
 * alongside it.
 *
 * Coding and splice-region variants only, which keeps Exomiser's memory in
 * bounds. Deep intronic and regulatory variants are not ranked; see
 * docs/decisions.md for why that was acceptable for this case.
 */
process EXOMISER {
    tag "${label}"
    label 'exomiser'
    publishDir "${params.outdir}/ranking", mode: 'link'

    input:
    tuple val(label), path(vcf), path(tbi)

    output:
    path "*.variants.tsv",                   emit: variants
    path "*.genes.tsv",                      emit: genes
    path "${params.sample_id}_${label}.*",   emit: all

    script:
    def hpo = params.hpo_ids.toString().tokenize(',').collect { id -> "'${id.trim()}'" }.join(', ')
    """
    PROBAND=\$(bcftools query -l ${vcf} | head -1)

    cat > analysis.yml <<EOF
    ---
    analysis:
      genomeAssembly: hg38
      vcf: ${vcf}
      ped:
      proband: \$PROBAND
      hpoIds: [${hpo}]
      inheritanceModes: {
        AUTOSOMAL_DOMINANT: 0.1,
        AUTOSOMAL_RECESSIVE_HOM_ALT: 0.1,
        AUTOSOMAL_RECESSIVE_COMP_HET: 2.0,
        X_DOMINANT: 0.1,
        X_RECESSIVE_HOM_ALT: 0.1,
        X_RECESSIVE_COMP_HET: 2.0,
        MITOCHONDRIAL: 0.2
      }
      analysisMode: PASS_ONLY
      frequencySources: [
        UK10K,
        GNOMAD_E_AFR, GNOMAD_E_AMR, GNOMAD_E_ASJ, GNOMAD_E_EAS, GNOMAD_E_FIN, GNOMAD_E_NFE, GNOMAD_E_SAS,
        GNOMAD_G_AFR, GNOMAD_G_AMR, GNOMAD_G_ASJ, GNOMAD_G_EAS, GNOMAD_G_FIN, GNOMAD_G_NFE, GNOMAD_G_SAS
      ]
      pathogenicitySources: [REVEL, MVP, ALPHA_MISSENSE, SPLICE_AI]
      steps: [
        failedVariantFilter: {},
        variantEffectFilter: {
          remove: [
            FIVE_PRIME_UTR_EXON_VARIANT, FIVE_PRIME_UTR_INTRON_VARIANT,
            THREE_PRIME_UTR_EXON_VARIANT, THREE_PRIME_UTR_INTRON_VARIANT,
            NON_CODING_TRANSCRIPT_EXON_VARIANT, NON_CODING_TRANSCRIPT_INTRON_VARIANT,
            CODING_TRANSCRIPT_INTRON_VARIANT,
            UPSTREAM_GENE_VARIANT, DOWNSTREAM_GENE_VARIANT,
            INTERGENIC_VARIANT, REGULATORY_REGION_VARIANT
          ]
        },
        frequencyFilter: { maxFrequency: 2.0 },
        pathogenicityFilter: { keepNonPathogenic: true },
        inheritanceFilter: {},
        omimPrioritiser: {},
        hiPhivePrioritiser: {}
      ]
    outputOptions:
      outputContributingVariantsOnly: false
      numGenes: 0
      outputDirectory: .
      outputFileName: ${params.sample_id}_${label}
      outputFormats: [HTML, JSON, TSV_GENE, TSV_VARIANT]
    EOF

    java -Xms4g -Xmx${params.exomiser_xmx}g -XX:+UseG1GC \\
        -Dspring.config.location=${params.exomiser_dir}/application.properties \\
        -jar ${params.exomiser_dir}/${params.exomiser_jar} \\
        analyse --analysis analysis.yml
    """
}

/*
 * Read-backed phasing over the candidate panel. With no parents, this is the
 * only evidence for cis vs trans in a compound-het pair.
 */
process WHATSHAP_PANEL {
    tag "${label}"
    label 'low'
    publishDir "${params.outdir}/phasing", mode: 'link'

    input:
    tuple val(label), path(vcf), path(tbi)
    tuple path(cram), path(crai)
    path panel

    output:
    path "${label}.phased.vcf.gz", emit: vcf
    path "${label}.phase_blocks.tsv", emit: blocks

    script:
    """
    bcftools view -R ${panel} -O z -o panel.vcf.gz ${vcf}
    tabix -p vcf panel.vcf.gz
    whatshap phase \\
        --reference ${params.fasta} \\
        --ignore-read-groups \\
        -o ${label}.phased.vcf.gz \\
        panel.vcf.gz ${cram}
    whatshap stats --tsv ${label}.phase_blocks.tsv ${label}.phased.vcf.gz
    """
}


// ===========================================================================
// Outputs
// ===========================================================================

process MAKE_SUBMISSION {
    label 'low'
    publishDir "${params.outdir}/submission", mode: 'copy'

    input:
    path ranked
    path phase
    tuple path(vcf), path(tbi)

    output:
    path "submission.csv", emit: csv

    script:
    def phase_arg = phase ? "--phase ${phase}" : ''
    def cur_arg = params.curation ? "--curation ${file(params.curation)}" : ''
    def gd_arg  = params.gene_diseases ? "--gene-diseases ${file(params.gene_diseases)}" : ''
    """
    make_submission.py --ranked ${ranked} ${phase_arg} ${cur_arg} ${gd_arg} \\
        --min-gene-score ${params.min_gene_score} --vcf ${vcf} --out submission.csv
    """
}

process MULTIQC {
    label 'low'
    publishDir "${params.outdir}/qc", mode: 'copy'

    input:
    path reports

    output:
    path "multiqc_report.html", emit: report

    script:
    """
    multiqc . -n multiqc_report.html
    """
}

/*
 * Tool versions for the methods section.
 */
process VERSIONS {
    label 'low'
    publishDir "${params.outdir}/pipeline_info", mode: 'copy'

    output:
    path "versions.txt"

    script:
    """
    {
      echo "fastp: \$(fastp --version 2>&1 | head -1 || true)"
      echo "bwa: \$(bwa 2>&1 | grep -i version || true)"
      echo "samtools: \$(samtools --version | head -1 || true)"
      echo "bcftools: \$(bcftools --version | head -1 || true)"
      echo "mosdepth: \$(mosdepth --version 2>&1 || true)"
      echo "delly: \$(delly --version 2>&1 | head -1 || true)"
      echo "ExpansionHunter: \$(ExpansionHunter --version 2>&1 | head -1 || true)"
      echo "whatshap: \$(whatshap --version 2>&1 || true)"
      echo "exomiser: ${params.exomiser_jar}"
      echo "deepvariant: ${params.deepvariant_version}"
      echo "vep: ${params.vep_sif} (cache ${params.vep_cache_version})"
      echo "java: \$(java -version 2>&1 | head -1)"
    } > versions.txt
    """
}


// ===========================================================================
// Workflow
// ===========================================================================

workflow {

    if (!params.hpo_ids) {
        error("params.hpo_ids is empty. Set comma-separated HPO IDs in params.yaml.")
    }

    VERSIONS()

    if (params.mode == 'vcf') {

        // ---- Vendor VCF only: bankable submission without alignment ----
        vendor = channel.of(['vendor', file(params.vendor_vcf), file("${params.vendor_vcf}.tbi")])

        NORMALIZE(vendor)
        EXOMISER(NORMALIZE.out.vcf)
        if (params.run_vep) {
            VEP(NORMALIZE.out.vcf)
        }
        MAKE_SUBMISSION(EXOMISER.out.variants, [], NORMALIZE.out.vcf.map { _l, vcf, tbi -> [vcf, tbi] })
    }
    else {

        // ---- Full run from FASTQ ----
        lanes = channel
            .fromFilePairs("${params.fastq_dir}/*_R{1,2}_001.fastq.gz", checkIfExists: true)
            .map { name, reads ->
                def lane = name.tokenize('_').find { tok -> tok.startsWith('L00') }
                [[id: params.sample_id, lane: lane], reads]
            }

        ALIGN_LANE(lanes)
        MARKDUP_CRAM(ALIGN_LANE.out.bam.map { _meta, bam -> bam }.collect())
        cram = MARKDUP_CRAM.out.cram

        ALIGN_QC(cram)
        MOSDEPTH(cram)

        DEEPVARIANT(cram)
        NORMALIZE(DEEPVARIANT.out.vcf.map { vcf, tbi -> ['dv', vcf, tbi] })

        EXOMISER(NORMALIZE.out.vcf)
        if (params.run_vep) {
            VEP(NORMALIZE.out.vcf)
        }

        def phase_ch = channel.value([])
        if (params.panel_bed) {
            panel = file(params.panel_bed)
            if (params.run_mosaic) {
                MUTECT2_PANEL(cram, panel)
            }
            WHATSHAP_PANEL(NORMALIZE.out.vcf, cram, panel)
            // WHATSHAP_PANEL is optional (errorStrategy 'ignore'). An ignored
            // failure leaves its output channel empty, and MAKE_SUBMISSION would
            // then wait forever and silently never run. ifEmpty falls back to
            // "no phasing", which make_submission.py already handles.
            phase_ch = WHATSHAP_PANEL.out.blocks.ifEmpty([])
        }

        if (params.run_sv) {
            DELLY(cram)
        }
        if (params.run_str && params.str_catalog) {
            EXPANSIONHUNTER(cram, file(params.str_catalog))
        }
        if (params.run_cnv) {
            ANEUPLOIDY(MOSDEPTH.out.bins, DEEPVARIANT.out.vcf)
        }

        MAKE_SUBMISSION(EXOMISER.out.variants, phase_ch, NORMALIZE.out.vcf.map { _l, vcf, tbi -> [vcf, tbi] })

        MULTIQC(
            ALIGN_LANE.out.json
                .mix(MARKDUP_CRAM.out.metrics, ALIGN_QC.out.reports.flatten(), MOSDEPTH.out.summary)
                .collect()
        )
    }
}
