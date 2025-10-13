configfile: "config/config.yaml"

localrules: get_kofam_chunks

rule all:
    input:
        config["kofam_dir"] + "/all_samples_kofam_annotations.txt",
        expand(config["coverm_dir"] + "/{sample}/coverm_{sample}.TPM.tsv", sample=config["samples"])


#"result/6_cluster_seq/clusterd_protein_file.fasta"
#config["functional_annotation_dir"] + "/all_samples_eggnogresult.emapper.annotations"
#        config["functional_annotation_dir"] + "/all_samples_eggnogresult.emapper.annotations",


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

# rule fastp_quality_control:
#     input:
#         r1=config["input_dir"] + "/{sample}.R1.raw.fastq.gz",
#         r2=config["input_dir"] + "/{sample}.R2.raw.fastq.gz"
#     output:
#         r1_clean=config["fastp_dir"] + "/{sample}/{sample}_1P.fq.gz",
#         r2_clean=config["fastp_dir"] + "/{sample}/{sample}_2P.fq.gz",
#         html=config["fastp_dir"] + "/{sample}/{sample}.fastp.html",
#         json=config["fastp_dir"] + "/{sample}/{sample}.fastp.json"
#     params:
#         threads=config["fastp"]["threads"],
#         qualified_quality_phred=config["fastp"]["qualified_quality_phred"],
#         length_required=config["fastp"]["length_required"],
#         input_dir=config["input_dir"],
#         fastp_dir=config["fastp_dir"]
#     threads: config["fastp"]["threads"]
#     resources:
#         mem_mb_per_cpu=config["regular_memory"],
#         runtime=config["fastp"]["runtime"],
#         cpus_per_task=config["fastp"]["threads"],
#         slurm_partition=config["regular_partition"],
#         slurm_account=config["account"]
#     conda:
#         "envs/fastp.yaml"
#     shell:
#         """
#         mkdir -p {params.fastp_dir}/{wildcards.sample}
#         fastp --thread {params.threads} \
#               --in1 {input.r1} \
#               --in2 {input.r2} \
#               --out1 {output.r1_clean} \
#               --out2 {output.r2_clean} \
#               -h {output.html} \
#               -j {output.json} \
#               --trim_poly_g \
#               --trim_poly_x \
#               --qualified_quality_phred {params.qualified_quality_phred} \
#               --length_required {params.length_required}
#         """

# rule spades_assembly_and_reformat:
#     input:
#         r1=config["fastp_dir"] + "/{sample}/{sample}_1P.fq.gz",
#         r2=config["fastp_dir"] + "/{sample}/{sample}_2P.fq.gz"
#     output:
#         reformated=config["reformated_scaffolds_dir"] + "/{sample}/{sample}_reformated.fa",
#         renaming=config["reformated_scaffolds_dir"] + "/{sample}/{sample}_renaming.tsv",
#         original=config["reformated_scaffolds_dir"] + "/{sample}/{sample}_original_contigs.fasta"
#     params:
#         spades_threads=config["spades"]["threads"],
#         spades_memory=config["spades"]["memory"],
#         kmers=config["spades"]["kmers"],
#         spades_outdir=config["spades_dir"] + "/{sample}",
#         spades_dir=config["spades_dir"],
#         log_dir=config["log_dir"],
#         min_length=config["reformat_scaffolds"]["min_length"],
#         reformated_scaffolds_dir=config["reformated_scaffolds_dir"],
#         input_dir=config["input_dir"]
#     threads: config["spades"]["threads"]
#     resources:
#         mem_mb_per_cpu=config["regular_memory"],
#         runtime=config["spades"]["runtime"] + config["reformat_scaffolds"]["runtime"],
#         cpus_per_task=config["spades"]["threads"],
#         slurm_partition=config["regular_partition"],
#         slurm_account=config["account"]
#     conda:
#         "envs/seqkit-spade.yaml"
#     log:
#         out="log/2_spades/{sample}_spades.out",
#         err="log/2_spades/{sample}_spades.err"
#     shell:
#         """
#         # SPAdes assembly
#         mkdir -p {params.log_dir}/2_spades
#         mkdir -p {params.spades_dir}/{wildcards.sample}
#         spades.py --meta \
#                   -o {params.spades_outdir} \
#                   -1 {input.r1} \
#                   -2 {input.r2} \
#                   -t {params.spades_threads} \
#                   -m {params.spades_memory} \
#                   -k {params.kmers} \
#                   --only-assembler > {log.out} 2> {log.err}
#         mkdir -p {params.reformated_scaffolds_dir}/{wildcards.sample}
#         echo {params.reformated_scaffolds_dir}/{wildcards.sample}
#         echo {params.spades_outdir}/contigs.fasta
#         # Reformat contig headers and filter by length
#         cp {params.spades_outdir}/contigs.fasta {output.original}
#         seqkit seq -m {params.min_length} {params.spades_outdir}/contigs.fasta | \
#         seqkit replace -p "(.+)" -r "{wildcards.sample}_\$1" > {output.reformated}
#         seqkit fx2tab {output.reformated} | cut -f1 | sed 's/>//' | nl -nln | \
#         awk 'BEGIN{{OFS="\t"}} {{print $2, "{wildcards.sample}_" $1}}' > {output.renaming}

#         # Cleanup intermediate files and raw data for this sample (after all reformatting is done)
#         rm -rf {params.spades_dir}/{wildcards.sample}
#         """



# rule pyrodigal_orfs:
#     input:
#         reformated=config["reformated_scaffolds_dir"] + "/{sample}/{sample}_reformated.fa"
#     output:
#         gff=config["pyrodigal_dir"] + "/{sample}/gff/{sample}.gff",
#         faa=config["pyrodigal_dir"] + "/{sample}/amino/{sample}.faa",
#         ffn=config["pyrodigal_dir"] + "/{sample}/nucl/{sample}.nucl.ffn"
#     params:
#         mode=config["pyrodigal"]["mode"],
#         pyrodigal_dir=config["pyrodigal_dir"]
#     threads: config["pyrodigal"]["threads"]
#     resources:
#         mem_mb_per_cpu=config["regular_memory"],
#         runtime=config["pyrodigal"]["runtime"],
#         cpus_per_task=config["pyrodigal"]["threads"],
#         slurm_partition=config["regular_partition"],
#         slurm_account=config["account"]
#     conda:
#         "envs/pyrodigal.yaml"
#     log:
#         out="log/4_pyrodigal/{sample}_pyrodigal.out",
#         err="log/4_pyrodigal/{sample}_pyrodigal.err"
#     shell:
#         """
#         mkdir -p log/4_pyrodigal
#         mkdir -p {params.pyrodigal_dir}/{wildcards.sample}/{{nucl,amino,gff}}
#         pyrodigal -i {input.reformated} \
#                   -o {output.gff} \
#                   -a {output.faa} \
#                   -d {output.ffn} \
#                   -f gff \
#                   -j {threads} \
#                   -p {params.mode} > {log.out} 2> {log.err}
#         """

# rule merge_sequences:
#     input:
#         nucl_files=expand(config["pyrodigal_dir"] + "/{sample}/nucl/{sample}.nucl.ffn",
#                          sample=config["samples"]),
#         amino_files=expand(config["pyrodigal_dir"] + "/{sample}/amino/{sample}.faa",
#                           sample=config["samples"])
#     output:
#         nucl_merged=config["mmseqs2_dir"] + "/all_samples_nucl_merged.fasta",
#         protein_merged=config["mmseqs2_dir"] + "/all_samples_protein_merged.fasta"
#     params:
#         mmseqs2_dir=config["mmseqs2_dir"]
#     threads: config["merge_sequences"]["threads"]
#     resources:
#         mem_mb_per_cpu=config["regular_memory"],
#         runtime=config["merge_sequences"]["runtime"],
#         cpus_per_task=config["merge_sequences"]["threads"],
#         slurm_partition=config["regular_partition"],
#         slurm_account=config["account"]
#     log:
#         out="log/5_merge_sequences/merge_sequences.out",
#         err="log/5_merge_sequences/merge_sequences.err"
#     shell:
#         """
#         mkdir -p log/5_merge_sequences
#         mkdir -p {params.mmseqs2_dir}
#         cat {input.nucl_files} > {output.nucl_merged} 2> {log.err}
#         cat {input.amino_files} > {output.protein_merged} 2>> {log.err}
#         echo "Merge sequences completed" > {log.out}
#         """

# rule mmseqs_clustering:
#     input:
#         protein_merged=config["mmseqs2_dir"] + "/all_samples_protein_merged.fasta"
#     output:
#         cluster_rep=config["mmseqs2_dir"] + "/all_samples_cluster_rep_seq.fasta",
#         cluster_all=config["mmseqs2_dir"] + "/all_samples_cluster_all_seqs.fasta",
#         cluster_tsv=config["mmseqs2_dir"] + "/all_samples_cluster_cluster.tsv"
#     params:
#         min_seq_id=config["mmseqs2"]["min_seq_id"],
#         coverage=config["mmseqs2"]["coverage"],
#         threads=config["mmseqs2"]["threads"],
#         prefix=config["mmseqs2_dir"] + "/all_samples_cluster",
#         tmp_dir=config["mmseqs2_dir"] + "/tmp",
#         memory=config["mmseqs2"]["memory"]
#     threads: config["mmseqs2"]["threads"]
#     resources:
#         mem_mb_per_cpu=config["regular_memory"],
#         runtime=config["mmseqs2"]["runtime"],
#         cpus_per_task=config["mmseqs2"]["threads"],
#         slurm_partition=config["regular_partition"],
#         slurm_account=config["account"]
#     conda:
#         "envs/mmseqs2.yaml"
#     log:
#         out="log/6_mmseqs_clustering/mmseqs_clustering.out",
#         err="log/6_mmseqs_clustering/mmseqs_clustering.err"
#     shell:
#         """
#         mkdir -p log/6_mmseqs_clustering
#         mkdir -p {params.tmp_dir}
#         mmseqs easy-linclust {input.protein_merged} \
#                            {params.prefix} \
#                            {params.tmp_dir} \
#                            --min-seq-id {params.min_seq_id} \
#                            -c {params.coverage} \
#                            --split-memory-limit {params.memory} \
#                            --threads {params.threads} > {log.out} 2> {log.err}
#         """



# rule extract_sequences:
#     input:
#         cluster_rep=config["mmseqs2_dir"] + "/all_samples_cluster_rep_seq.fasta",
#         protein_merged=config["mmseqs2_dir"] + "/all_samples_protein_merged.fasta",
#         nucl_merged=config["mmseqs2_dir"] + "/all_samples_nucl_merged.fasta"
#     output:
#         unigene_id=config["mmseqs2_dir"] + "/unigene_id.txt",
#         clustered_protein_file="result/6_cluster_seq/clusterd_protein_file.fasta",
#         clustered_acid_file="result/6_cluster_seq/clusterd_acid_file.fasta"
#     threads: config["extract_sequences"]["threads"]
#     resources:
#         mem_mb_per_cpu=config["regular_memory"],
#         runtime=config["extract_sequences"]["runtime"],
#         cpus_per_task=config["extract_sequences"]["threads"],
#         slurm_partition=config["regular_partition"],
#         slurm_account=config["account"]
#     conda:
#         "envs/seqkit.yaml"
#     log:
#         out="log/7_extract_sequences/extract_sequences.out",
#         err="log/7_extract_sequences/extract_sequences.err"
#     shell:
#         """
#         mkdir -p log/7_extract_sequences result/6_cluster_seq
#         grep '>' {input.cluster_rep} | awk '{{print $1}}' | sed 's/>//' > {output.unigene_id} 2> {log.err}
#         seqkit grep -f {output.unigene_id} {input.protein_merged} -o {output.clustered_protein_file} >> {log.out} 2>> {log.err}
#         seqkit grep -f {output.unigene_id} {input.nucl_merged} -o {output.clustered_acid_file} >> {log.out} 2>> {log.err}
#         """






#----------------------------------------------------------------------------
# eggNOG annotation rules

# rule split_protein_for_eggnog:
#     input:
#         protein="result/6_cluster_seq/clusterd_protein_file.fasta"
#     output:
#         split_done=config["functional_annotation_dir"] + "/split_sequences/split.done"
#     params:
#         split_dir=config["functional_annotation_dir"] + "/split_sequences",
#         chunk_size=config["eggnog"]["chunk_size"]
#     threads: config["extract_sequences"]["threads"]
#     resources:
#         mem_mb_per_cpu=config["regular_memory"],
#         runtime=config["extract_sequences"]["runtime"],
#         cpus_per_task=config["extract_sequences"]["threads"],
#         slurm_partition=config["regular_partition"],
#         slurm_account=config["account"]
#     conda:
#         "envs/seqkit.yaml"
#     log:
#         out="log/9_eggnog_annotation/split_sequences.out",
#         err="log/9_eggnog_annotation/split_sequences.err"
#     shell:
#         """
#         mkdir -p {params.split_dir} log/9_eggnog_annotation
#         seqkit split2 -s {params.chunk_size} -O {params.split_dir} {input.protein} > {log.out} 2> {log.err}
#         touch {output.split_done}
#         """

# checkpoint get_eggnog_chunks:
#     input:
#         split_done=config["functional_annotation_dir"] + "/split_sequences/split.done"
#     output:
#         chunk_list=config["functional_annotation_dir"] + "/split_sequences/chunk_list.txt"
#     params:
#         split_dir=config["functional_annotation_dir"] + "/split_sequences"
#     shell:
#         """
#         ls {params.split_dir}/clusterd_protein_file.part_*.fasta | \
#         sed 's|.*/clusterd_protein_file.part_||; s|.fasta||' | \
#         sort -n > {output.chunk_list}
#         """

# rule eggnog_annotation_chunk:
#     input:
#         split_file=config["functional_annotation_dir"] + "/split_sequences/clusterd_protein_file.part_{chunk}.fasta"
#     output:
#         annotations=config["functional_annotation_dir"] + "/chunks/chunk_{chunk}.emapper.annotations"
#     params:
#         threads_per_job=config["eggnog"]["threads_per_job"],
#         evalue=config["eggnog"]["evalue"],
#         data_dir=config["databases"]["eggnog_data_dir"],
#         output_prefix=config["functional_annotation_dir"] + "/chunks/chunk_{chunk}"
#     threads: config["eggnog"]["threads_per_job"]
#     resources:
#         mem_mb_per_cpu=config["regular_memory"],
#         runtime=config["eggnog"]["runtime"],
#         cpus_per_task=config["eggnog"]["threads_per_job"],
#         slurm_partition=config["regular_partition"],
#         slurm_account=config["account"]
#     conda:
#         "envs/eggnog.yaml"
#     log:
#         out="log/9_eggnog_annotation/chunk_{chunk}.out",
#         err="log/9_eggnog_annotation/chunk_{chunk}.err"
#     shell:
#         """
#         mkdir -p $(dirname {params.output_prefix}) log/9_eggnog_annotation
#         emapper.py -m diamond \
#                   --cpu {params.threads_per_job} \
#                   --seed_ortholog_evalue {params.evalue} \
#                   --data_dir {params.data_dir} \
#                   -i {input.split_file} \
#                   -o {params.output_prefix} \
#                   --resume > {log.out} 2> {log.err}
#         """

# def get_eggnog_chunk_files(wildcards):
#     checkpoint_output = checkpoints.get_eggnog_chunks.get(**wildcards).output.chunk_list
#     with open(checkpoint_output) as f:
#         chunks = [line.strip() for line in f]
#     return expand(config["functional_annotation_dir"] + "/chunks/chunk_{chunk}.emapper.annotations",
#                   chunk=chunks)

# rule merge_eggnog_annotations:
#     input:
#         get_eggnog_chunk_files
#     output:
#         annotations=config["functional_annotation_dir"] + "/all_samples_eggnogresult.emapper.annotations"
#     threads: 1
#     resources:
#         mem_mb_per_cpu=config["regular_memory"],
#         runtime=60,
#         cpus_per_task=1,
#         slurm_partition=config["regular_partition"],
#         slurm_account=config["account"]
#     log:
#         out="log/9_eggnog_annotation/merge_annotations.out",
#         err="log/9_eggnog_annotation/merge_annotations.err"
#     shell:
#         """
#         # Extract header from first chunk
#         head -n 4 {input[0]} > {output.annotations} 2> {log.err}

#         # Merge all annotations (skip headers)
#         for chunk in {input}; do
#             tail -n +5 $chunk | grep -v "^#" >> {output.annotations}
#         done 2>> {log.err}

#         echo "Merged $(echo {input} | wc -w) annotation chunks" > {log.out}
#         """

#----------------------------------------------------------------------------


# rule coverm_abundance:
#     input:
#         r1_clean=config["fastp_dir"] + "/{sample}/{sample}_1P.fq.gz",
#         r2_clean=config["fastp_dir"] + "/{sample}/{sample}_2P.fq.gz",
#         reference="result/6_cluster_seq/clusterd_acid_file.fasta"
#     output:
#         tpm=config["coverm_dir"] + "/{sample}/coverm_{sample}.TPM.tsv"
#     params:
#         threads=config["coverm"]["threads"],
#         min_read_percent_identity=config["coverm"]["min_read_percent_identity"],
#         min_read_aligned_percent=config["coverm"]["min_read_aligned_percent"],
#         min_covered_fraction=config["coverm"]["min_covered_fraction"],
#         contig_end_exclusion=config["coverm"]["contig_end_exclusion"],
#         coverm_dir=config["coverm_dir"]
#     threads: config["coverm"]["threads"]
#     resources:
#         mem_mb_per_cpu=config["regular_memory"],
#         runtime=config["coverm"]["runtime"],
#         cpus_per_task=config["coverm"]["threads"],
#         slurm_partition=config["regular_partition"],
#         slurm_account=config["account"]
#     conda:
#         "envs/coverm.yaml"
#     log:
#         out="log/10_coverm_abundance/{sample}_coverm.out",
#         err="log/10_coverm_abundance/{sample}_coverm.err"
#     shell:
#         """
#         mkdir -p log/10_coverm_abundance
#         mkdir -p {params.coverm_dir}/{wildcards.sample}
#         coverm contig -1 {input.r1_clean} \
#                      -2 {input.r2_clean} \
#                      -t {params.threads} \
#                      --reference {input.reference} \
#                      -m tpm \
#                      -o {output.tpm} \
#                      --exclude-supplementary \
#                      --min-read-percent-identity {params.min_read_percent_identity} \
#                      --min-read-aligned-percent {params.min_read_aligned_percent} \
#                      --min-covered-fraction {params.min_covered_fraction} \
#                      --contig-end-exclusion {params.contig_end_exclusion} > {log.out} 2> {log.err}
#         """



# KOfam annotation rules
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

def get_kofam_chunk_files(wildcards):
    checkpoint_output = checkpoints.get_kofam_chunks.get(**wildcards).output.chunk_list
    with open(checkpoint_output) as f:
        chunks = [line.strip() for line in f]
    return expand(config["kofam_dir"] + "/chunks/chunk_{chunk}.kofam.txt",
                  chunk=chunks)

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

        # Filter for best hit per gene where score > threshold and E-value <= 1e-5
        # Process all lines (with or without *), apply filters
        tail -n +3 {input.annotations} | \
        awk 'BEGIN{{OFS="\t"}}
        {{
            # Remove leading * or tab if present
            sub(/^[\*\t]+/, "")
            gene=$1
            threshold=$3
            score=$4
            evalue=$5
            diff=score-threshold
            # Only keep hits where score > threshold AND E-value <= 1e-5
            if (diff > 0 && evalue <= 1e-5) {{
                # Keep best hit: highest (score - threshold)
                if (!(gene in best_diff) || diff > best_diff[gene]) {{
                    best_diff[gene] = diff
                    best_line[gene] = $0
                }}
            }}
        }}
        END {{
            for (gene in best_line) {{
                print best_line[gene]
            }}
        }}' | sort -k1,1 >> {output.filtered} 2>> {log.err}

        echo "Filtered to best hits for chunk (score > threshold)" > {log.out}
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

# '''
# Untill here, all the rules are finished for abundance calculation and functional annotation. the following scripts for result analyze.m
# '''