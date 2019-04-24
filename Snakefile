configfile: "config.json"

def get_bam_files():
    IDS, = glob_wildcards(config["bam_path"] + "/{id}.bam")
    bamfiles = [base + ".bam" for base in IDS]
    return bamfiles

def get_fastq_files():
    IDS, = glob_wildcards(config["fastq_path"] + "/{id}.fastq.gz")
    fastq_files = [base + ".fastq.gz" for base in IDS]
    return fastq_files

rule bam_links:
    input:
        [config["bam_path"] + "/" + bam for bam in get_bam_files()]
    output:
        ["data/BAM/" + bam for bam in get_bam_files()]
    run:
        for bam in get_bam_files():
            shell("ln -s " + config["bam_path"] + "/" + bam + " data/BAM/" + bam)

rule fastq_links:
    input:
        [config["fastq_path"] + "/" + fastq_file for fastq_file in get_fastq_files()]
    output:
        ["data/FASTQ/" + fastq for fastq in get_fastq_files()]
    run:
        for fastq in get_fastq_files():
            shell("ln -s " + config["fastq_path"] + "/" + fastq + " data/FASTQ/" + fastq)

rule download_genome:
    output:
        "resources/genome.fasta.gz"
    shell:
        'wget ' + config["genome_url"] + ' -O {output}'

rule decompress_genome:
    input:
        "resources/genome.fasta.gz"
    output:
        "resources/genome.fasta"
    shell:
        "gzcat {input} > {output}"

rule bwa_index:
    input:
        "resources/genome.fasta"
    output:
        ["resources/genome.fasta" + ending for ending in ['.amb','.ann','.bwt','.pac','.sa']]
    shell:
        "bwa index -a bwtsw {input}"

rule bwa_single:
    input:
        ref="resources/genome.fasta",
        index=["resources/genome.fasta" + ending for ending in ['.amb','.ann','.bwt','.pac','.sa']],
        forward="data/FASTQ/{sample}_R1.fastq.gz",
        reverse="data/FASTQ/{sample}_R2.fastq.gz"
    threads: 4
    output:
        temp("output/SAM/{sample}.sam")
    shell:
        "bwa mem -t {threads} {input.ref} {input.forward} {input.reverse} | samtools view -h -F 4 > {output}"

rule bam_single:
    input:
        "output/SAM/{sample}.sam"
    threads: 4
    output:
        "output/BAM/{sample}.bam"
    shell:
        "samtools sort -@ 3 -m 5G {input} > {output}"

rule bam_single_slurm_script:
    input:
        "resources/slurm_script.sh"
    output:
        "output/Slurm_scripts/make_BAM/bwa_{sample}.sh"
    shell:
        "cat {input} | sed 's/SAMPLE/{wildcards.sample}/g' > {output}"

rule featurecounts:
    input:
        ["data/BAM/" + bam for bam in get_bam_files()]
    output:
        report="reports/featureCounts.report.txt",
        counts="data/gene_counts.tsv"
    script:
        "scripts/featureCounts.R"

rule dds:
    input:
        counts="data/gene_counts.tsv"
    output:
        rds="data/RData/dds.rds"
    script:
        "scripts/make_dds.R"

rule normalized_counts:
    input:
        rds="data/RData/dds.rds"
    output:
        counts="output/normalized_counts.tsv"
    script:
        "scripts/normalized_counts.R"

rule size_factor_plot:
    input:
        dds="data/RData/dds.rds",
        Rmd="scripts/SizeFactors.Rmd"
    output:
        html="output/SizeFactors.html"
    script:
        "scripts/SizeFactors.Rmd"

rule DE_single:
    input:
        dds="data/RData/dds.rds"
    output:
        dds="data/RData/dds-{comparison}.rds",
        rds="data/RData/DE-{comparison}.rds"
    script:
        "scripts/DifferentialExpression.R"

rule DE_report_single:
    input:
        dds="data/RData/dds-{comparison}.rds",
        rds="data/RData/DE-{comparison}.rds",
        Rmd="scripts/DifferentialExpression.Rmd"
    output:
        html="output/DifferentialExpression/DiffExp-{comparison}.html"
    script:
        "scripts/DifferentialExpression.Rmd"

rule DE_reports:
    input:
        ["output/DifferentialExpression/DiffExp-" + comparison + ".html" for comparison in config['comparisons'].keys()]

rule PCA_vsd:
    input:
        dds="data/RData/dds.rds"
    output:
        vsd="data/RData/vsd-{pca}.rds"
    script:
        "scripts/PCA_vsd.R"

rule PCA_plot:
    input:
        dds="data/RData/dds.rds",
        Rmd="scripts/PCA.Rmd",
        vsd="data/RData/vsd-{pca}.rds"
    output:
        "output/PCA/PCA-{pca}.html"
    script:
        "scripts/PCA.Rmd"

rule PCA_plots:
    input:
        ["output/PCA/PCA-" + pca + ".html" for pca in config['PCA_plots'].keys()]

rule heatmap_vsd:
    input:
        dds="data/RData/dds.rds"
    output:
        vsd="data/RData/heatmap-vsd-{heatmap}.rds"
    script:
        "scripts/heatmap_vsd.R"

rule heatmap:
    input:
        vsd="data/RData/heatmap-vsd-{heatmap}.rds",
        Rmd="scripts/Heatmap.Rmd"
    output:
        "output/Heatmaps/Heatmap-{heatmap}.html"
    script:
        "scripts/Heatmap.Rmd"

rule heatmaps:
    input:
        ["output/Heatmaps/Heatmap-" + heatmap + ".html" for heatmap in config['heatmaps'].keys()]

rule plots:
    input:
        dds="data/RData/dds.rds",
        Rmd="scripts/Plots.Rmd"
    output:
        html="output/Heatmap_PCA.html"
    script:
        "scripts/Plots.Rmd"

rule fgsea_rds_single:
    input:
        rds="data/RData/DE-{comparison}.rds"
    output:
        fgsea="data/RData/fgsea-{comparison}.rds"
    script:
        "scripts/fgsea.R"

rule fgsea_html:
    input:
        rds="data/RData/fgsea-{comparison}.rds",
        Rmd="scripts/fgsea.Rmd"
    output:
        html="output/fgsea/fgsea-{comparison}.html"
    script:
        "scripts/fgsea.Rmd"

rule fgsea:
    input:
        ["output/fgsea/fgsea-" + comparison + ".html" for comparison in config['comparisons'].keys() if config['comparisons'][comparison]['do_fgsea']]
