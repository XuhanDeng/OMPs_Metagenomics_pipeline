configfile: "config/config.yaml"

localrules: get_kofam_chunks

rule all:
    input:
        expand("results/11_bowtie2/4_coverm/{sample}_TPM.tsv", sample=config["samples"]),
        expand("results/11_bowtie2/4_coverm/{sample}_mean.tsv", sample=config["samples"]),
        expand("results/11_bowtie2/4_coverm/{sample}_count.tsv", sample=config["samples"])


# Step 17: Build bowtie2 index for representative sequences
rule build_cluster_index:
    input:
        representatives="results/10_cluster/3_clustered/cluster_representatives.fasta"
    output:
        flag="results/10_cluster/4_bowtie_index/cluster_index.done"
    conda: "envs/bowtie2.yaml"
    threads: config["bowtie2_index"]["threads"]
    log:
        out="log/17_build_cluster_index/build_index.log",
        err="log/17_build_cluster_index/build_index.err"
    resources:
        slurm_partition=config["regular_partition"],
        runtime=config["bowtie2_index"]["runtime"],
        mem_mb_per_cpu=config["bowtie2_index"]["mem_mb_per_cpu"],
        cpus_per_task=config["bowtie2_index"]["threads"],
        slurm_account=config["account"]
    shell:
        """
        mkdir -p results/10_cluster/4_bowtie_index
        mkdir -p log/17_build_cluster_index

        # Build bowtie2 index
        bowtie2-build {input.representatives} results/10_cluster/4_bowtie_index/cluster_index \
            > {log.out} 2> {log.err}

        touch {output.flag}
        """

# Step 18: Bowtie2 alignment and BAM processing
rule bowtie2_alignment:
    input:
        r1="results/1_fastp/{sample}/{sample}_1P.fq.gz",
        r2="results/1_fastp/{sample}/{sample}_2P.fq.gz",
        index_flag="results/10_cluster/4_bowtie_index/cluster_index.done"
    output:
        bam="results/11_bowtie2/3_samtool_bam/{sample}.f2.sorted.bam",
        bam_index="results/11_bowtie2/3_samtool_bam/{sample}.f2.sorted.bam.bai"
    conda: "envs/bowtie2_samtools.yaml"
    threads: config["bowtie2"]["threads"]
    log:
        out="log/18_bowtie2_alignment/{sample}_alignment.log",
        err="log/18_bowtie2_alignment/{sample}_alignment.err"
    resources:
        slurm_partition=config["bowtie2"]["partition"],
        runtime=config["bowtie2"]["runtime"],
        mem_mb_per_cpu=config["bowtie2"]["mem_mb_per_cpu"],
        cpus_per_task=config["bowtie2"]["threads"],
        slurm_account=config["account"]
    shell:
        """
        mkdir -p results/11_bowtie2/{{2_bowtie2_alignment_sam,3_samtool_bam}}
        mkdir -p log/18_bowtie2_alignment

        # Bowtie2 alignment
        bowtie2 --sensitive -t -p {threads} \
            -x results/10_cluster/4_bowtie_index/cluster_index \
            -1 {input.r1} -2 {input.r2} \
            -S results/11_bowtie2/2_bowtie2_alignment_sam/{wildcards.sample}.sam \
            > {log.out} 2> {log.err}

        # Samtools processing
        samtools view -@ {threads} -hbS -f 2 results/11_bowtie2/2_bowtie2_alignment_sam/{wildcards.sample}.sam \
        | samtools sort -@ {threads} -o {output.bam} - \
            >> {log.out} 2>> {log.err}

        samtools index -@ {threads} {output.bam} \
            >> {log.out} 2>> {log.err}

        # Clean up SAM file to save space
        rm -f results/11_bowtie2/2_bowtie2_alignment_sam/{wildcards.sample}.sam
        """

# Step 19: Calculate abundance using CoverM
rule calculate_abundance:
    input:
        bam="results/11_bowtie2/3_samtool_bam/{sample}.f2.sorted.bam",
        bam_index="results/11_bowtie2/3_samtool_bam/{sample}.f2.sorted.bam.bai"
    output:
        tpm="results/11_bowtie2/4_coverm/{sample}_TPM.tsv",
        mean="results/11_bowtie2/4_coverm/{sample}_mean.tsv",
        read_count="results/11_bowtie2/4_coverm/{sample}_count.tsv"
    conda: "envs/coverm.yaml"
    threads: config["coverm_abundance"]["threads"]
    log:
        out="log/19_coverm/{sample}_coverm.log",
        err="log/19_coverm/{sample}_coverm.err"
    resources:
        slurm_partition=config["regular_partition"],
        runtime=config["coverm_abundance"]["runtime"],
        mem_mb_per_cpu=config["coverm_abundance"]["mem_mb_per_cpu"],
        cpus_per_task=config["coverm_abundance"]["threads"],
        slurm_account=config["account"]
    shell:
        """
        mkdir -p results/11_bowtie2/4_coverm
        mkdir -p log/19_coverm

        # CoverM abundance calculations
        coverm contig -m tpm \
            --bam-files {input.bam} \
            --min-read-percent-identity {config[coverm_abundance][min_read_percent_identity]} \
            --min-read-aligned-percent {config[coverm_abundance][min_read_aligned_percent]} \
            --output-file {output.tpm} \
            --contig-end-exclusion {config[coverm_abundance][contig_end_exclusion]} \
            --no-zeros -t {threads} \
            > {log.out} 2> {log.err}

        coverm contig -m mean \
            --bam-files {input.bam} \
            --min-read-percent-identity {config[coverm_abundance][min_read_percent_identity]} \
            --min-read-aligned-percent {config[coverm_abundance][min_read_aligned_percent]} \
            --output-file {output.mean} \
            --contig-end-exclusion {config[coverm_abundance][contig_end_exclusion]} \
            --no-zeros -t {threads} \
            >> {log.out} 2>> {log.err}

        coverm contig -m count \
            --bam-files {input.bam} \
            --min-read-percent-identity {config[coverm_abundance][min_read_percent_identity]} \
            --min-read-aligned-percent {config[coverm_abundance][min_read_aligned_percent]} \
            --output-file {output.read_count} \
            --contig-end-exclusion {config[coverm_abundance][contig_end_exclusion]} \
            --no-zeros -t {threads} \
            >> {log.out} 2>> {log.err}
        """


