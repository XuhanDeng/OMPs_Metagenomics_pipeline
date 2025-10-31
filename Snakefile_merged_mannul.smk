import os

configfile: "config/config.yaml"

rule all:
    input:
        # Taxonomic classification
        "result/97_braken_report",
        # Functional annotation
        config["kofam_dir"] + "/all_samples_kofam_annotations.txt",
        # Per-sample abundance files
        expand("result/10_abundance/{sample}/{sample}.TPM.tsv", sample=config["samples"]),
        expand("result/10_abundance/{sample}/{sample}.mean.tsv", sample=config["samples"]),
        expand("result/10_abundance/{sample}/{sample}.count.tsv", sample=config["samples"]),
        # Merged abundance files
        "result/10_abundance/merged/all_samples.TPM.tsv",
        "result/10_abundance/merged/all_samples.mean.tsv",
        "result/10_abundance/merged/all_samples.count.tsv",
        # Annotated gene-level abundance files
        "result/11_abundance_annotation/all_samples.TPM.annotated.tsv",
        "result/11_abundance_annotation/all_samples.mean.annotated.tsv",
        "result/11_abundance_annotation/all_samples.count.annotated.tsv",
        # KO-aggregated abundance files
        "result/11_abundance_annotation/all_samples.TPM.KO_aggregated.tsv",
        "result/11_abundance_annotation/all_samples.mean.KO_aggregated.tsv",
        "result/11_abundance_annotation/all_samples.count.KO_aggregated.tsv",


#

# ============================================================================
# METAGENOMICS FUNCTIONAL ANNOTATION PIPELINE
# ============================================================================
# This pipeline processes metagenomic sequencing data to generate functional annotations:
#
# 1. Quality Control (fastp): Filters and trims raw paired-end reads
# 2. Assembly (SPAdes): Performs metagenomic assembly to generate contigs
# 3. ORF Prediction (Pyrodigal): Predicts open reading frames (ORFs) from contigs
# 4. Sequence Merging: Combines all sample sequences (nucleotide and protein)
# 5. Clustering (MMseqs2): Clusters proteins to create non-redundant unigene catalog
# 6. Sequence Extraction: Extracts representative sequences (unigenes)
# 7. Functional Annotation (eggNOG-mapper): Annotates proteins with functional information
# 8. Abundance Quantification (CoverM): Calculates gene abundance across samples
# ============================================================================

rule fastp_quality_control:
    input:
        r1=config["input_dir"] + "/{sample}.R1.raw.fastq.gz",
        r2=config["input_dir"] + "/{sample}.R2.raw.fastq.gz"
    output:
        r1_clean=config["fastp_dir"] + "/{sample}/{sample}_1P.fq.gz",
        r2_clean=config["fastp_dir"] + "/{sample}/{sample}_2P.fq.gz",
        html=config["fastp_dir"] + "/{sample}/{sample}.fastp.html",
        json=config["fastp_dir"] + "/{sample}/{sample}.fastp.json"
    params:
        threads=config["fastp"]["threads"],
        qualified_quality_phred=config["fastp"]["qualified_quality_phred"],
        length_required=config["fastp"]["length_required"],
        input_dir=config["input_dir"],
        fastp_dir=config["fastp_dir"]
    threads: config["fastp"]["threads"]
    resources:
        mem_mb_per_cpu=config["regular_memory"],
        runtime=config["fastp"]["runtime"],
        cpus_per_task=config["fastp"]["threads"],
        slurm_partition=config["regular_partition"],
        slurm_account=config["account"]
    conda:
        "envs/fastp.yaml"
    shell:
        """
        mkdir -p {params.fastp_dir}/{wildcards.sample}
        fastp --thread {params.threads} \
              --in1 {input.r1} \
              --in2 {input.r2} \
              --out1 {output.r1_clean} \
              --out2 {output.r2_clean} \
              -h {output.html} \
              -j {output.json} \
              --trim_poly_g \
              --trim_poly_x \
              --qualified_quality_phred {params.qualified_quality_phred} \
              --length_required {params.length_required}
        """

rule spades_assembly_and_reformat:
    input:
        r1=config["fastp_dir"] + "/{sample}/{sample}_1P.fq.gz",
        r2=config["fastp_dir"] + "/{sample}/{sample}_2P.fq.gz"
    output:
        reformated=config["reformated_scaffolds_dir"] + "/{sample}/{sample}_reformated.fa",
        renaming=config["reformated_scaffolds_dir"] + "/{sample}/{sample}_renaming.tsv",
        original=config["reformated_scaffolds_dir"] + "/{sample}/{sample}_original_contigs.fasta"
    params:
        spades_threads=config["spades"]["threads"],
        spades_memory=config["spades"]["memory"],
        kmers=config["spades"]["kmers"],
        spades_outdir=config["spades_dir"] + "/{sample}",
        spades_dir=config["spades_dir"],
        log_dir=config["log_dir"],
        min_length=config["reformat_scaffolds"]["min_length"],
        reformated_scaffolds_dir=config["reformated_scaffolds_dir"],
        input_dir=config["input_dir"]
    threads: config["spades"]["threads"]
    resources:
        mem_mb_per_cpu=config["regular_memory"],
        runtime=config["spades"]["runtime"] + config["reformat_scaffolds"]["runtime"],
        cpus_per_task=config["spades"]["threads"],
        slurm_partition=config["regular_partition"],
        slurm_account=config["account"]
    conda:
        "envs/seqkit-spade.yaml"
    log:
        out="log/2_spades/{sample}_spades.out",
        err="log/2_spades/{sample}_spades.err"
    shell:
        """
        # SPAdes assembly
        mkdir -p {params.log_dir}/2_spades
        mkdir -p {params.spades_dir}/{wildcards.sample}
        spades.py --meta \
                  -o {params.spades_outdir} \
                  -1 {input.r1} \
                  -2 {input.r2} \
                  -t {params.spades_threads} \
                  -m {params.spades_memory} \
                  -k {params.kmers} \
                  --only-assembler > {log.out} 2> {log.err}
        mkdir -p {params.reformated_scaffolds_dir}/{wildcards.sample}
        echo {params.reformated_scaffolds_dir}/{wildcards.sample}
        echo {params.spades_outdir}/contigs.fasta
        # Reformat contig headers and filter by length
        cp {params.spades_outdir}/contigs.fasta {output.original}
        seqkit seq -m {params.min_length} {params.spades_outdir}/contigs.fasta | \
        seqkit replace -p "(.+)" -r "{wildcards.sample}_\$1" > {output.reformated}
        seqkit fx2tab {output.reformated} | cut -f1 | sed 's/>//' | nl -nln | \
        awk 'BEGIN{{OFS="\t"}} {{print $2, "{wildcards.sample}_" $1}}' > {output.renaming}

        # Cleanup intermediate files and raw data for this sample (after all reformatting is done)
        rm -rf {params.spades_dir}/{wildcards.sample}
        """

rule pyrodigal_orfs:
    input:
        reformated=config["reformated_scaffolds_dir"] + "/{sample}/{sample}_reformated.fa"
    output:
        gff=config["pyrodigal_dir"] + "/{sample}/gff/{sample}.gff",
        faa=config["pyrodigal_dir"] + "/{sample}/amino/{sample}.faa",
        ffn=config["pyrodigal_dir"] + "/{sample}/nucl/{sample}.nucl.ffn"
    params:
        mode=config["pyrodigal"]["mode"],
        pyrodigal_dir=config["pyrodigal_dir"]
    threads: config["pyrodigal"]["threads"]
    resources:
        mem_mb_per_cpu=config["regular_memory"],
        runtime=config["pyrodigal"]["runtime"],
        cpus_per_task=config["pyrodigal"]["threads"],
        slurm_partition=config["regular_partition"],
        slurm_account=config["account"]
    conda:
        "envs/pyrodigal.yaml"
    log:
        out="log/4_pyrodigal/{sample}_pyrodigal.out",
        err="log/4_pyrodigal/{sample}_pyrodigal.err"
    shell:
        """
        mkdir -p log/4_pyrodigal
        mkdir -p {params.pyrodigal_dir}/{wildcards.sample}/{{nucl,amino,gff}}
        pyrodigal -i {input.reformated} \
                  -o {output.gff} \
                  -a {output.faa} \
                  -d {output.ffn} \
                  -f gff \
                  -j {threads} \
                  -p {params.mode} > {log.out} 2> {log.err}
        """

rule merge_sequences:
    input:
        nucl_files=expand(config["pyrodigal_dir"] + "/{sample}/nucl/{sample}.nucl.ffn",
                         sample=config["samples"]),
        amino_files=expand(config["pyrodigal_dir"] + "/{sample}/amino/{sample}.faa",
                          sample=config["samples"])
    output:
        nucl_merged=config["mmseqs2_dir"] + "/all_samples_nucl_merged.fasta",
        protein_merged=config["mmseqs2_dir"] + "/all_samples_protein_merged.fasta"
    params:
        mmseqs2_dir=config["mmseqs2_dir"]
    threads: config["merge_sequences"]["threads"]
    resources:
        mem_mb_per_cpu=config["regular_memory"],
        runtime=config["merge_sequences"]["runtime"],
        cpus_per_task=config["merge_sequences"]["threads"],
        slurm_partition=config["regular_partition"],
        slurm_account=config["account"]
    log:
        out="log/5_merge_sequences/merge_sequences.out",
        err="log/5_merge_sequences/merge_sequences.err"
    shell:
        """
        mkdir -p log/5_merge_sequences
        mkdir -p {params.mmseqs2_dir}
        cat {input.nucl_files} > {output.nucl_merged} 2> {log.err}
        cat {input.amino_files} > {output.protein_merged} 2>> {log.err}
        echo "Merge sequences completed" > {log.out}
        """

rule mmseqs_clustering:
    input:
        protein_merged=config["mmseqs2_dir"] + "/all_samples_protein_merged.fasta"
    output:
        cluster_rep=config["mmseqs2_dir"] + "/all_samples_cluster_rep_seq.fasta",
        cluster_all=config["mmseqs2_dir"] + "/all_samples_cluster_all_seqs.fasta",
        cluster_tsv=config["mmseqs2_dir"] + "/all_samples_cluster_cluster.tsv"
    params:
        min_seq_id=config["mmseqs2"]["min_seq_id"],
        coverage=config["mmseqs2"]["coverage"],
        threads=config["mmseqs2"]["threads"],
        prefix=config["mmseqs2_dir"] + "/all_samples_cluster",
        tmp_dir=config["mmseqs2_dir"] + "/tmp",
        memory=config["mmseqs2"]["memory"]
    threads: config["mmseqs2"]["threads"]
    resources:
        mem_mb_per_cpu=config["regular_memory"],
        runtime=config["mmseqs2"]["runtime"],
        cpus_per_task=config["mmseqs2"]["threads"],
        slurm_partition=config["regular_partition"],
        slurm_account=config["account"]
    conda:
        "envs/mmseqs2.yaml"
    log:
        out="log/6_mmseqs_clustering/mmseqs_clustering.out",
        err="log/6_mmseqs_clustering/mmseqs_clustering.err"
    shell:
        """
        mkdir -p log/6_mmseqs_clustering
        mkdir -p {params.tmp_dir}
        mmseqs easy-linclust {input.protein_merged} \
                           {params.prefix} \
                           {params.tmp_dir} \
                           --min-seq-id {params.min_seq_id} \
                           -c {params.coverage} \
                           --split-memory-limit {params.memory} \
                           --threads {params.threads} > {log.out} 2> {log.err}
        """



rule extract_sequences:
    input:
        cluster_rep=config["mmseqs2_dir"] + "/all_samples_cluster_rep_seq.fasta",
        protein_merged=config["mmseqs2_dir"] + "/all_samples_protein_merged.fasta",
        nucl_merged=config["mmseqs2_dir"] + "/all_samples_nucl_merged.fasta"
    output:
        unigene_id=config["mmseqs2_dir"] + "/unigene_id.txt",
        clustered_protein_file="result/6_cluster_seq/clusterd_protein_file.fasta",
        clustered_acid_file="result/6_cluster_seq/clusterd_acid_file.fasta"
    threads: config["extract_sequences"]["threads"]
    resources:
        mem_mb_per_cpu=config["regular_memory"],
        runtime=config["extract_sequences"]["runtime"],
        cpus_per_task=config["extract_sequences"]["threads"],
        slurm_partition=config["regular_partition"],
        slurm_account=config["account"]
    conda:
        "envs/seqkit.yaml"
    log:
        out="log/7_extract_sequences/extract_sequences.out",
        err="log/7_extract_sequences/extract_sequences.err"
    shell:
        """
        mkdir -p log/7_extract_sequences result/6_cluster_seq
        grep '>' {input.cluster_rep} | awk '{{print $1}}' | sed 's/>//' > {output.unigene_id} 2> {log.err}
        seqkit grep -f {output.unigene_id} {input.protein_merged} -o {output.clustered_protein_file} >> {log.out} 2>> {log.err}
        seqkit grep -f {output.unigene_id} {input.nucl_merged} -o {output.clustered_acid_file} >> {log.out} 2>> {log.err}
        """




# # KOfam annotation rules
rule split_protein_for_kofam:
    input:
        protein="result/6_cluster_seq/clusterd_protein_file.fasta"
    output:
        split_done=config["kofam_dir"] + "/split_sequences/split.done"
    params:
        split_dir=config["kofam_dir"] + "/split_sequences",
        chunk_parts=config["kofamscan"]["chunk_parts"]
    threads: 1
    resources:
        mem_mb_per_cpu=config["regular_memory"],
        runtime=60,
        cpus_per_task=1,
        slurm_partition=config["regular_partition"],
        slurm_account=config["account"]
    conda:
        "envs/seqkit.yaml"
    log:
        out="log/10_kofam_annotation/split_sequences.out",
        err="log/10_kofam_annotation/split_sequences.err"
    shell:
        """
        mkdir -p {params.split_dir} log/10_kofam_annotation
        seqkit split2 -p {params.chunk_parts} -O {params.split_dir} {input.protein} > {log.out} 2> {log.err}
        touch {output.split_done}
        """

checkpoint get_kofam_chunks:
    input:
        split_done=config["kofam_dir"] + "/split_sequences/split.done"
    output:
        chunk_list=config["kofam_dir"] + "/split_sequences/chunk_list.txt"
    params:
        split_dir=config["kofam_dir"] + "/split_sequences"
    shell:
        """
        ls {params.split_dir}/clusterd_protein_file.part_*.fasta | \
        sed 's|.*/clusterd_protein_file.part_||; s|.fasta||' | \
        sed 's/^0*//' | \
        sort -n > {output.chunk_list}
        """

rule kofamscan_annotation_chunk:
    input:
        split_file=lambda wildcards: config["kofam_dir"] + f"/split_sequences/clusterd_protein_file.part_{int(wildcards.chunk):03d}.fasta"
    output:
        annotations=config["kofam_dir"] + "/chunks/chunk_{chunk}.kofam.txt"
    params:
        threads=config["kofamscan"]["threads"],
        evalue=config["kofamscan"]["evalue"],
        kofam_profile=config["databases"]["kofam_profiles"],
        kofam_list=config["databases"]["kofam_list"],
        tmp_dir=lambda wildcards: config["kofam_dir"] + f"/tmp/chunk_{wildcards.chunk}"
    threads: config["kofamscan"]["threads"]
    resources:
        mem_mb_per_cpu=config["regular_memory"],
        runtime=config["kofamscan"]["runtime"],
        cpus_per_task=config["kofamscan"]["threads"],
        slurm_partition=config["regular_partition"],
        slurm_account=config["account"]
    conda:
        "envs/kofamscan.yaml"
    log:
        out="log/10_kofam_annotation/chunk_{chunk}.out",
        err="log/10_kofam_annotation/chunk_{chunk}.err"
    shell:
        """
        mkdir -p $(dirname {output.annotations}) log/10_kofam_annotation
        mkdir -p {params.tmp_dir}
        exec_annotation -o {output.annotations} \
                       -p {params.kofam_profile} \
                       -k {params.kofam_list} \
                       --cpu {params.threads} \
                       -E {params.evalue} \
                       -f detail-tsv \
                       --tmp-dir {params.tmp_dir} \
                       {input.split_file} > {log.out} 2> {log.err}
        """

# def get_kofam_chunk_files(wildcards):
#     checkpoint_output = checkpoints.get_kofam_chunks.get(**wildcards).output.chunk_list
#     with open(checkpoint_output) as f:
#         chunks = [line.strip() for line in f]
#     return expand(config["kofam_dir"] + "/chunks/chunk_{chunk}.kofam.txt",
#                   chunk=chunks)

def get_filtered_kofam_chunk_files(wildcards):
    checkpoint_output = checkpoints.get_kofam_chunks.get(**wildcards).output.chunk_list
    with open(checkpoint_output) as f:
        chunks = [line.strip() for line in f]
    return expand(config["kofam_dir"] + "/chunks/chunk_{chunk}.kofam.filtered.txt",
                  chunk=chunks)

rule filter_kofam_best_hit:
    input:
        annotations=config["kofam_dir"] + "/chunks/chunk_{chunk}.kofam.txt"
    output:
        filtered=config["kofam_dir"] + "/chunks/chunk_{chunk}.kofam.filtered.txt"
    threads: 1
    resources:
        mem_mb_per_cpu=config["regular_memory"],
        runtime=600,
        cpus_per_task=2,
        slurm_partition=config["regular_partition"],
        slurm_account=config["account"]
    log:
        out="log/10_kofam_annotation/filter_chunk_{chunk}.out",
        err="log/10_kofam_annotation/filter_chunk_{chunk}.err"
    shell:
        """
        # Extract header
        head -n 1 {input.annotations} > {output.filtered} 2> {log.err}

        # Filter: keep only lines with asterisks (*), E-value <= 1e-5, and best score per gene
        tail -n +3 {input.annotations} | \
        grep '^\*' | \
        awk 'BEGIN{{OFS="\t"}}
        {{
            # Remove leading asterisk and whitespace, re-parse fields
            sub(/^\*[[:space:]]*/, "")
            $1=$1  # Force field re-splitting

            gene=$1
            score=$4
            evalue=$5

            # Only process if E-value <= 1e-5
            if (evalue <= 1e-5) {{
                # Keep best hit: highest score per gene
                if (!(gene in best_score) || score > best_score[gene]) {{
                    best_score[gene] = score
                    best_line[gene] = $0
                }}
            }}
        }}
        END {{
            for (gene in best_line) {{
                print best_line[gene]
            }}
        }}' | sort -k1,1 >> {output.filtered} 2>> {log.err}

        echo "Filtered to best hits with asterisks, highest score per gene, E-value <= 1e-5" > {log.out}
        """

rule merge_kofam_annotations:
    input:
        get_filtered_kofam_chunk_files
    output:
        annotations=config["kofam_dir"] + "/all_samples_kofam_annotations.txt"
    threads: 1
    resources:
        mem_mb_per_cpu=config["regular_memory"],
        runtime=60,
        cpus_per_task=1,
        slurm_partition=config["regular_partition"],
        slurm_account=config["account"]
    log:
        out="log/10_kofam_annotation/merge_annotations.out",
        err="log/10_kofam_annotation/merge_annotations.err"
    shell:
        """
        # Extract header from first chunk
        head -n 1 {input[0]} > {output.annotations} 2> {log.err}

        # Merge all annotations (skip headers)
        for chunk in {input}; do
            tail -n +2 $chunk | grep -v "^#" >> {output.annotations}
        done 2>> {log.err}

        echo "Merged $(echo {input} | wc -w) KOfam annotation chunks" > {log.out}
        """



# Step 16: Clean sequence names and filter by length
rule clean_sequence_names:
    input:
        representatives="result/6_cluster_seq/clusterd_acid_file.fasta"
    output:
        cleaned="result/6_cluster_seq/clusterd_acid_file.cleaned.fasta"
    conda: "envs/seqkit.yaml"
    threads: 1
    log:
        out="log/16_clean_names/clean_names.log",
        err="log/16_clean_names/clean_names.err"
    resources:
        slurm_partition=config["regular_partition"],
        runtime=60,
        mem_mb_per_cpu=config["regular_memory"],
        cpus_per_task=1,
        slurm_account=config["account"]
    shell:
        """
        mkdir -p log/16_clean_names

        # Filter sequences >= 300bp and remove everything after the first space in headers
        seqkit seq -m 300 {input.representatives} | \
        seqkit replace -p ' .+$' -r '' -o {output.cleaned} \
            > {log.out} 2> {log.err}
        """

# Step 17: Build bowtie2 index for representative sequences
rule build_cluster_index:
    input:
        representatives="result/6_cluster_seq/clusterd_acid_file.cleaned.fasta"
    output:
        flag="result/6_cluster_seq/bowtie_index/cluster_index.done"
    conda: "envs/bowtie2.yaml"
    threads: config["bowtie2_index"]["threads"]
    log:
        out="log/17_build_cluster_index/build_index.log",
        err="log/17_build_cluster_index/build_index.err"
    resources:
        slurm_partition=config["bowtie2_index"]["partition"],
        runtime=config["bowtie2_index"]["runtime"],
        mem_mb_per_cpu=config["bowtie2_index"]["mem_mb_per_cpu"],
        cpus_per_task=config["bowtie2_index"]["threads"],
        slurm_account=config["account"]
    shell:
        """
        mkdir -p result/6_cluster_seq/bowtie_index
        mkdir -p log/17_build_cluster_index

        # Build bowtie2 index
        bowtie2-build --threads {threads} {input.representatives} result/6_cluster_seq/bowtie_index/cluster_index \
            > {log.out} 2> {log.err}

        touch {output.flag}
        """

# Step 18: Bowtie2 alignment and BAM processing (combined with pipe)
rule bowtie2_alignment:
    input:
        r1=config["fastp_dir"] + "/{sample}/{sample}_1P.fq.gz",
        r2=config["fastp_dir"] + "/{sample}/{sample}_2P.fq.gz",
        index_flag="result/6_cluster_seq/bowtie_index/cluster_index.done"
    output:
        bam=config["coverm_dir"] + "/{sample}/bowtie2/{sample}.f2.sorted.bam",
        bam_index=config["coverm_dir"] + "/{sample}/bowtie2/{sample}.f2.sorted.bam.bai",
        unsorted_bam=config["coverm_dir"] + "/{sample}/bowtie2/{sample}.unsorted.bam"
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
        # Set sample-specific temporary directory
        export TMPDIR=/scratch/$USER/tmp/{wildcards.sample}
        mkdir -p $TMPDIR
        mkdir -p {config[coverm_dir]}/{wildcards.sample}/bowtie2
        mkdir -p log/18_bowtie2_alignment

        # Bowtie2 alignment piped directly to samtools (generate unsorted BAM)
        bowtie2 --sensitive -t -p {threads} \
            -x result/6_cluster_seq/bowtie_index/cluster_index \
            -1 {input.r1} -2 {input.r2} \
            2> {log.err} \
        | samtools view -@ {threads} -bS \
            -T result/6_cluster_seq/clusterd_acid_file.cleaned.fasta \
            -o {output.unsorted_bam} - \
            > {log.out} 2>> {log.err}

        # Sort BAM with explicit temp directory
        samtools sort -@ {threads} -m 1G \
            -T $TMPDIR/samtools_sort \
            -o {output.bam} {output.unsorted_bam} \
            >> {log.out} 2>> {log.err}

        # Clean up temp directory
        rm -rf $TMPDIR

        # Index BAM
        samtools index -@ {threads} {output.bam} \
            >> {log.out} 2>> {log.err}
        """

# Step 18.5: Extract gene lengths from clusterd_acid_file
rule extract_gene_lengths:
    input:
        reference="result/6_cluster_seq/clusterd_acid_file.fasta"
    output:
        lengths="result/10_abundance/gene_lengths.txt"
    conda: "envs/seqkit.yaml"
    threads: 2
    log:
        out="log/18_gene_lengths/gene_lengths.log",
        err="log/18_gene_lengths/gene_lengths.err"
    resources:
        slurm_partition=config["regular_partition"],
        runtime=160,
        mem_mb_per_cpu=config["regular_memory"],
        cpus_per_task=2,
        slurm_account=config["account"]
    shell:
        """
        mkdir -p result/10_abundance log/18_gene_lengths

        # Use seqkit to get sequence lengths
        seqkit fx2tab -l -n -i {input.reference} | \
        awk '{{print $1"\t"$2}}' > {output.lengths} 2> {log.err}

        echo "Gene lengths extracted" > {log.out}
        """

# Step 19: Run samtools idxstats to get read counts per contig
rule samtools_idxstats:
    input:
        bam=config["coverm_dir"] + "/{sample}/bowtie2/{sample}.f2.sorted.bam",
        bam_index=config["coverm_dir"] + "/{sample}/bowtie2/{sample}.f2.sorted.bam.bai"
    output:
        idxstats="result/10_abundance/{sample}/idxstats.txt"
    conda: "envs/samtools.yaml"
    threads: 8
    log:
        out="log/19_idxstats/{sample}_idxstats.log",
        err="log/19_idxstats/{sample}_idxstats.err"
    resources:
        slurm_partition=config["regular_partition"],
        runtime=60,
        mem_mb_per_cpu=config["regular_memory"],
        cpus_per_task=8,
        slurm_account=config["account"]
    shell:
        """
        mkdir -p result/10_abundance/{wildcards.sample}
        mkdir -p log/19_idxstats

        # Get read counts from idxstats: contig, length, mapped_reads, unmapped_reads
        samtools idxstats {input.bam} > {output.idxstats} 2> {log.err}

        echo "idxstats completed" > {log.out}
        """

# Step 20: Calculate abundance metrics (TPM, mean coverage, read count)
rule calculate_abundance:
    input:
        idxstats="result/10_abundance/{sample}/idxstats.txt",
        gene_lengths="result/10_abundance/gene_lengths.txt"
    output:
        tpm="result/10_abundance/{sample}/{sample}.TPM.tsv",
        mean="result/10_abundance/{sample}/{sample}.mean.tsv",
        read_count="result/10_abundance/{sample}/{sample}.count.tsv"
    threads: 16
    log:
        out="log/20_abundance/{sample}_abundance.log",
        err="log/20_abundance/{sample}_abundance.err"
    resources:
        slurm_partition=config["regular_partition"],
        runtime=60,
        mem_mb_per_cpu=config["regular_memory"],
        cpus_per_task=16,
        slurm_account=config["account"]
    shell:
        """
        mkdir -p result/10_abundance/{wildcards.sample}
        mkdir -p log/20_abundance

        # Calculate sum of (Ni/Li) for all genes - the denominator in TPM formula
        sum_rpk=$(awk 'NR==FNR {{len[$1]=$2; next}}
                       $1 != "*" && $3 > 0 && $1 in len {{
                           sum += $3 / len[$1]
                       }}
                       END {{print sum}}' {input.gene_lengths} \
                           {input.idxstats})

        # Calculate TPM: (Ni/Li * 10^6) / sum(N1/L1 + N2/L2 + ... + Nn/Ln)
        echo -e "contigName\t{wildcards.sample} TPM" > {output.tpm}
        awk -v sum_rpk="$sum_rpk" 'NR==FNR {{len[$1]=$2; next}}
             $1 != "*" && $3 > 0 && $1 in len {{
                 tpm = (($3 / len[$1]) * 1000000) / sum_rpk
                 print $1 "\t" tpm
             }}' {input.gene_lengths} \
                 {input.idxstats} >> {output.tpm} 2>> {log.err}

        # Calculate mean coverage: (mapped_reads * read_length) / gene_length
        echo -e "contigName\t{wildcards.sample} Mean" > {output.mean}
        awk -v read_len={config[read_length]} 'NR==FNR {{len[$1]=$2; next}}
             $1 != "*" && $3 > 0 && $1 in len {{
                 mean = ($3 * read_len) / len[$1]
                 print $1 "\t" mean
             }}' {input.gene_lengths} \
                 {input.idxstats} >> {output.mean} 2>> {log.err}

        # Calculate read count
        echo -e "contigName\t{wildcards.sample} Read Count" > {output.read_count}
        awk '$1 != "*" && $3 > 0 {{
            print $1 "\t" $3
        }}' {input.idxstats} >> {output.read_count} 2>> {log.err}

        echo "Abundance calculations completed" > {log.out}
        """

# Step 21a: Merge TPM files
rule merge_tpm:
    input:
        tpm_files=expand("result/10_abundance/{sample}/{sample}.TPM.tsv", sample=config["samples"])
    output:
        tpm_merged="result/10_abundance/merged/all_samples.TPM.tsv"
    conda: "envs/pandas.yaml"
    threads: 8
    log:
        out="log/21_merge_abundance/merge_tpm.log",
        err="log/21_merge_abundance/merge_tpm.err"
    resources:
        slurm_partition="memory",
        runtime=1800,  # 3 hours
        mem_mb_per_cpu=25000,  # 10GB - reduced with pairwise merging
        cpus_per_task=8,
        slurm_account=config["account"]
    shell:
        """
        mkdir -p result/10_abundance/merged
        mkdir -p log/21_merge_abundance

        python script/merge_abundance.py TPM {output.tpm_merged} {input.tpm_files} > {log.out} 2> {log.err}
        """

# Step 21b: Merge Mean coverage files
rule merge_mean:
    input:
        mean_files=expand("result/10_abundance/{sample}/{sample}.mean.tsv", sample=config["samples"])
    output:
        mean_merged="result/10_abundance/merged/all_samples.mean.tsv"
    conda: "envs/pandas.yaml"
    threads: 8
    log:
        out="log/21_merge_abundance/merge_mean.log",
        err="log/21_merge_abundance/merge_mean.err"
    resources:
        slurm_partition="memory",
        runtime=1800,  # 3 hours
        mem_mb_per_cpu=25000,  # 10GB - reduced with pairwise merging
        cpus_per_task=8,
        slurm_account=config["account"]
    shell:
        """
        mkdir -p result/10_abundance/merged
        mkdir -p log/21_merge_abundance

        python script/merge_abundance.py Mean {output.mean_merged} {input.mean_files} > {log.out} 2> {log.err}
        """

# Step 21c: Merge Read count files
rule merge_count:
    input:
        count_files=expand("result/10_abundance/{sample}/{sample}.count.tsv", sample=config["samples"])
    output:
        count_merged="result/10_abundance/merged/all_samples.count.tsv"
    conda: "envs/pandas.yaml"
    threads: 8
    log:
        out="log/21_merge_abundance/merge_count.log",
        err="log/21_merge_abundance/merge_count.err"
    resources:
        slurm_partition="memory",
        runtime=1800,  # 3 hours
        mem_mb_per_cpu=25000,  # 10GB - reduced with pairwise merging
        cpus_per_task=8,
        slurm_account=config["account"]
    shell:
        """
        mkdir -p result/10_abundance/merged
        mkdir -p log/21_merge_abundance

        python script/merge_abundance.py Count {output.count_merged} {input.count_files} > {log.out} 2> {log.err}
        """

# Step 22: Add KOfam annotations to abundance matrices and aggregate by KO
rule add_annotations_tpm:
    input:
        abundance="result/10_abundance/merged/all_samples.TPM.tsv",
        annotation="result/8_kofam_annotation/all_samples_kofam_annotations.txt"
    output:
        annotated="result/11_abundance_annotation/all_samples.TPM.annotated.tsv",
        ko_aggregated="result/11_abundance_annotation/all_samples.TPM.KO_aggregated.tsv"
    conda: "envs/pandas.yaml"
    threads: 32
    log:
        out="log/22_add_annotations/tpm.log",
        err="log/22_add_annotations/tpm.err"
    resources:
        slurm_partition="compute",
        runtime=600,
        mem_mb_per_cpu=3900,
        cpus_per_task=32,
        slurm_account=config["account"]
    shell:
        """
        mkdir -p result/11_abundance_annotation log/22_add_annotations
        python script/add_annotation_to_abundance.py {input.abundance} {input.annotation} {output.annotated} {output.ko_aggregated} > {log.out} 2> {log.err}
        """

rule add_annotations_mean:
    input:
        abundance="result/10_abundance/merged/all_samples.mean.tsv",
        annotation="result/8_kofam_annotation/all_samples_kofam_annotations.txt"
    output:
        annotated="result/11_abundance_annotation/all_samples.mean.annotated.tsv",
        ko_aggregated="result/11_abundance_annotation/all_samples.mean.KO_aggregated.tsv"
    conda: "envs/pandas.yaml"
    threads: 32
    log:
        out="log/22_add_annotations/mean.log",
        err="log/22_add_annotations/mean.err"
    resources:
        slurm_partition="compute",
        runtime=600,
        mem_mb_per_cpu=3900,
        cpus_per_task=32,
        slurm_account=config["account"]
    shell:
        """
        mkdir -p result/11_abundance_annotation log/22_add_annotations
        python script/add_annotation_to_abundance.py {input.abundance} {input.annotation} {output.annotated} {output.ko_aggregated} > {log.out} 2> {log.err}
        """

rule add_annotations_count:
    input:
        abundance="result/10_abundance/merged/all_samples.count.tsv",
        annotation="result/8_kofam_annotation/all_samples_kofam_annotations.txt"
    output:
        annotated="result/11_abundance_annotation/all_samples.count.annotated.tsv",
        ko_aggregated="result/11_abundance_annotation/all_samples.count.KO_aggregated.tsv"
    conda: "envs/pandas.yaml"
    threads: 32
    log:
        out="log/22_add_annotations/count.log",
        err="log/22_add_annotations/count.err"
    resources:
        slurm_partition="compute",
        runtime=600,
        mem_mb_per_cpu=3900,
        cpus_per_task=32,
        slurm_account=config["account"]
    shell:
        """
        mkdir -p result/11_abundance_annotation log/22_add_annotations
        python script/add_annotation_to_abundance.py {input.abundance} {input.annotation} {output.annotated} {output.ko_aggregated} > {log.out} 2> {log.err}
        """




# ============================================================================


rule kraken2:
    input:
        r1_clean=config["fastp_dir"] + "/{sample}/{sample}_1P.fq.gz",
        r2_clean=config["fastp_dir"] + "/{sample}/{sample}_2P.fq.gz"
    output:
        report=config["kraken2_dir"] + "/{sample}/{sample}_kraken2.report",
        out=config["kraken2_dir"] + "/{sample}/{sample}_kraken2.out"
    log:
        log=config["log_dir"] + "/kraken2/{sample}_kraken2.log",
        err=config["log_dir"] + "/kraken2/{sample}_kraken2.err"
    conda:
        "envs/kraken2.yaml"
    resources:
        mem_mb_per_cpu=config["regular_memory"],
        slurm_partition=config["regular_partition"],
        runtime=config["kraken2"]["runtime"],
        cpus_per_task=config["kraken2"]["threads"],
        slurm_account=config["account"]
    threads: config["kraken2"]["threads"]
    params:
        out_dir=config["kraken2_dir"] + "/{sample}",
        db=config["databases"]["kraken2_db"],
        confidence=config["kraken2"]["confidence"]

    shell:
        """
        mkdir -p {params.out_dir} && \
        mkdir -p $(dirname {log.log}) && \
        kraken2 --db {params.db} \
                --threads {threads} \
                --confidence {params.confidence} \
                --report {output.report} \
                --paired {input.r1_clean} {input.r2_clean} \
                --output {output.out} \
                > {log.log} 2> {log.err}
        """
rule bracken_build:
    input:
        db=config["databases"]["kraken2_db"]
    output:
        marker=config["databases"]["kraken2_db"] + "/bracken_build.done"
    log:
        log=config["log_dir"] + "/bracken_build/bracken_build.log",
        err=config["log_dir"] + "/bracken_build/bracken_build.err"
    conda:
        "envs/kraken2.yaml"
    #resources are for slurm, using the same as kraken2
    resources:
        mem_mb_per_cpu=config["regular_memory"],
        slurm_partition=config["regular_partition"],
        runtime=config["kraken2"]["runtime"],
        cpus_per_task=config["kraken2"]["threads"],
        slurm_account=config["account"]
    params:
        db=config["databases"]["kraken2_db"],
        kmer_length=config["bracken"]["kmer_length"],
        read_length=config["bracken"]["read_length"]
    threads: config["bracken"]["threads"]
    shell:
        """
        mkdir -p $(dirname {log.log}) && \
        bracken-build -d {input.db} -t {threads} -k 35 -l 150 > {log.log} 2> {log.err} && \
        touch {output.marker}
        """
rule bracken:
    input:
        report=config["kraken2_dir"] + "/{sample}/{sample}_kraken2.report",
        db=config["databases"]["kraken2_db"],
        bracken_db=config["databases"]["kraken2_db"] + "/bracken_build.done"
    output:
        bracken_report_P=config["kraken2_dir"] + "/{sample}/{sample}_bracken_phylum.report",
        bracken_report_C=config["kraken2_dir"] + "/{sample}/{sample}_bracken_classes.report",
        bracken_report_O=config["kraken2_dir"] + "/{sample}/{sample}_bracken_orders.report",
        bracken_report_F=config["kraken2_dir"] + "/{sample}/{sample}_bracken_family.report",
        bracken_report_G=config["kraken2_dir"] + "/{sample}/{sample}_bracken_genus.report",
        bracken_report_S=config["kraken2_dir"] + "/{sample}/{sample}_bracken_species.report"
    log:
        log=config["log_dir"] + "/bracken/{sample}_bracken.log",
        err=config["log_dir"] + "/bracken/{sample}_bracken.err"
    conda:
        "envs/kraken2.yaml"
    resources:
        mem_mb_per_cpu=config["regular_memory"],
        slurm_partition=config["regular_partition"],
        runtime=120,  # minutes
        cpus_per_task=2,
        slurm_account=config["account"]
    threads: 2
    params:
        db=config["databases"]["kraken2_db"],
        kmer_length=config["bracken"]["kmer_length"],
        read_length=config["bracken"]["read_length"]
    shell:
        """
        mkdir -p $(dirname {log.log}) && \
        bracken -d /scratch/qou/database/kraken2 -i {input.report} -o {output.bracken_report_P} -r 150 -l P -t 2 && \
        bracken -d /scratch/qou/database/kraken2 -i {input.report} -o {output.bracken_report_C} -r 150 -l C -t 2 && \
        bracken -d /scratch/qou/database/kraken2 -i {input.report} -o {output.bracken_report_O} -r 150 -l O -t 2 && \
        bracken -d /scratch/qou/database/kraken2 -i {input.report} -o {output.bracken_report_F} -r 150 -l F -t 2 && \
        bracken -d /scratch/qou/database/kraken2 -i {input.report} -o {output.bracken_report_G} -r 150 -l G -t 2 && \
        bracken -d /scratch/qou/database/kraken2 -i {input.report} -o {output.bracken_report_S} -r 150 -l S -t 2 \
                > {log.log} 2> {log.err}
        """

rule merge_bracken:
    input:
        P_reports=expand(config["kraken2_dir"] + "/{sample}/{sample}_bracken_phylum.report", sample=config["samples"]),
        C_reports=expand(config["kraken2_dir"] + "/{sample}/{sample}_bracken_classes.report", sample=config["samples"]),
        O_reports=expand(config["kraken2_dir"] + "/{sample}/{sample}_bracken_orders.report", sample=config["samples"]),
        F_reports=expand(config["kraken2_dir"] + "/{sample}/{sample}_bracken_family.report", sample=config["samples"]),
        G_reports=expand(config["kraken2_dir"] + "/{sample}/{sample}_bracken_genus.report", sample=config["samples"]),
        S_reports=expand(config["kraken2_dir"] + "/{sample}/{sample}_bracken_species.report", sample=config["samples"]),
        P_kraken_reports=expand(config["kraken2_dir"] + "/{sample}/{sample}_kraken2_bracken_phylums.report", sample=config["samples"]),
        C_kraken_reports=expand(config["kraken2_dir"] + "/{sample}/{sample}_kraken2_bracken_classes.report", sample=config["samples"]),
        O_kraken_reports=expand(config["kraken2_dir"] + "/{sample}/{sample}_kraken2_bracken_orders.report", sample=config["samples"]),
        F_kraken_reports=expand(config["kraken2_dir"] + "/{sample}/{sample}_kraken2_bracken_families.report", sample=config["samples"]),
        G_kraken_reports=expand(config["kraken2_dir"] + "/{sample}/{sample}_kraken2_bracken_genuses.report", sample=config["samples"]),
        S_kraken_reports=expand(config["kraken2_dir"] + "/{sample}/{sample}_kraken2_bracken_species.report", sample=config["samples"])
    output:
        directory("result/98_bracken_kraken")
    shell:
        """
        mkdir -p {output}/krakenstykle/P {output}/krakenstykle/C {output}/krakenstykle/O {output}/krakenstykle/F {output}/krakenstykle/G {output}/krakenstykle/S && \
        mkdir -p {output}/brackenstykle/P {output}/brackenstykle/C {output}/brackenstykle/O {output}/brackenstykle/F {output}/brackenstykle/G {output}/brackenstykle/S && \
        for file in {input.P_kraken_reports}; do
            cp $file {output}/krakenstykle/P/
        done && \
        for file in {input.C_kraken_reports}; do
            cp $file {output}/krakenstykle/C/
        done && \
        for file in {input.O_kraken_reports}; do
            cp $file {output}/krakenstykle/O/
        done && \
        for file in {input.F_kraken_reports}; do
            cp $file {output}/krakenstykle/F/
        done && \
        for file in {input.G_kraken_reports}; do
            cp $file {output}/krakenstykle/G/
        done && \
        for file in {input.S_kraken_reports}; do
            cp $file {output}/krakenstykle/S/
        done && \
        for file in {input.P_reports}; do
            cp $file {output}/brackenstykle/P/
        done && \
        for file in {input.C_reports}; do
            cp $file {output}/brackenstykle/C/
        done && \
        for file in {input.O_reports}; do
            cp $file {output}/brackenstykle/O/
        done && \
        for file in {input.F_reports}; do
            cp $file {output}/brackenstykle/F/
        done && \
        for file in {input.G_reports}; do
            cp $file {output}/brackenstykle/G/
        done && \
        for file in {input.S_reports}; do
            cp $file {output}/brackenstykle/S/
        done
        """

rule merge_bracken_report:
    input:
        dir="result/98_bracken_kraken"
    output:
        dir=directory("result/97_braken_report")
    conda:
        "envs/pandas.yaml"
    params:
        outdir=lambda wildcards, output: os.path.abspath(output.dir)
    shell:
        """
        mkdir -p {params.outdir} && \
        python script/merge_profiling_reports.py -i {input.dir}/krakenstykle/P -o {params.outdir}/phylum && \
        python script/merge_profiling_reports.py -i {input.dir}/krakenstykle/C -o {params.outdir}/classes && \
        python script/merge_profiling_reports.py -i {input.dir}/krakenstykle/O -o {params.outdir}/orders && \
        python script/merge_profiling_reports.py -i {input.dir}/krakenstykle/F -o {params.outdir}/family && \
        python script/merge_profiling_reports.py -i {input.dir}/krakenstykle/G -o {params.outdir}/genus && \
        python script/merge_profiling_reports.py -i {input.dir}/krakenstykle/S -o {params.outdir}/species
        """

rule pcoa_analysis:
    input:
        csv="result/97_braken_report/{level}_rel_abund.csv"
    output:
        svg="result/99_pcoa_plots/{level}_pcoa.svg"
    conda:
        "envs/R.yaml"
    log:
        "log/pcoa/{level}_pcoa.log"
    shell:
        """
        mkdir -p result/99_pcoa_plots log/pcoa && \
        Rscript script/pcoa_analysis.R {input.csv} {output.svg} > {log} 2>&1
        """

