configfile: "config/config.yaml"

rule all:
    input:
        expand(config["gtdb_dir"] + "/{group}/annotation.txt", 
               group=config["sample_groups"].keys()),
        expand(config["functional_annotation_dir"] + "/{group}_eggnogresult.emapper.annotations", 
               group=config["sample_groups"].keys()),
        expand(config["coverm_dir"] + "/{group}/coverm_{sample}.TPM.tsv", 
               group=config["sample_groups"].keys(),
               sample=config["samples"])

rule fastp_quality_control:
    input:
        r1=config["input_dir"] + "/{sample}/{sample}.R1.fq.gz",
        r2=config["input_dir"] + "/{sample}/{sample}.R2.fq.gz"
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

rule spades_assembly:
    input:
        r1=rules.fastp_quality_control.output.r1_clean,
        r2=rules.fastp_quality_control.output.r2_clean
    output:
        contigs=config["spades_dir"] + "/{sample}/contigs.fasta"
    params:
        threads=config["spades"]["threads"],
        memory=config["spades"]["memory"],
        kmers=config["spades"]["kmers"],
        outdir=config["spades_dir"] + "/{sample}",
        spades_dir=config["spades_dir"],
        log_dir=config["log_dir"]
    threads: config["spades"]["threads"]
    resources:
        mem_mb_per_cpu=config["spades"]["memory"] * 1000 // config["spades"]["threads"],
        runtime=config["spades"]["runtime"],
        cpus_per_task=config["spades"]["threads"],
        slurm_partition=config["regular_partition"],
        slurm_account=config["account"]
    conda:
        "envs/spades.yaml"
    log:
        "log/2_spades/{sample}_spades.txt"
    shell:
        """
        mkdir -p {params.log_dir}/2_spades
        mkdir -p {params.spades_dir}/{wildcards.sample}
        spades.py --meta \
                  -o {params.outdir} \
                  -1 {input.r1} \
                  -2 {input.r2} \
                  -t {params.threads} \
                  -m {params.memory} \
                  -k {params.kmers} \
                  --only-assembler > {log} 2>&1
        """

rule reformat_scaffolds:
    input:
        contigs=rules.spades_assembly.output.contigs
    output:
        reformated=config["reformated_scaffolds_dir"] + "/{sample}/{sample}_reformated.fa",
        renaming=config["reformated_scaffolds_dir"] + "/{sample}/{sample}_renaming.tsv"
    params:
        min_length=config["reformat_scaffolds"]["min_length"],
        reformated_scaffolds_dir=config["reformated_scaffolds_dir"]
    threads: config["reformat_scaffolds"]["threads"]
    resources:
        mem_mb_per_cpu=config["regular_memory"],
        runtime=config["reformat_scaffolds"]["runtime"],
        cpus_per_task=config["reformat_scaffolds"]["threads"],
        slurm_partition=config["regular_partition"],
        slurm_account=config["account"]
    conda:
        "envs/anvio.yaml"
    shell:
        """
        mkdir -p {params.reformated_scaffolds_dir}/{wildcards.sample}
        anvi-script-reformat-fasta {input.contigs} \
                                   -l {params.min_length} \
                                   -o {output.reformated} \
                                   --simplify-names \
                                   -r {output.renaming} \
                                   --prefix {wildcards.sample}
        """

rule prodigal_orfs:
    input:
        reformated=rules.reformat_scaffolds.output.reformated
    output:
        gff=config["prodigal_dir"] + "/{sample}/gff/{sample}.gff",
        faa=config["prodigal_dir"] + "/{sample}/amino/{sample}.faa",
        ffn=config["prodigal_dir"] + "/{sample}/nucl/{sample}.nucl.ffn"
    params:
        mode=config["prodigal"]["mode"],
        prodigal_dir=config["prodigal_dir"]
    threads: config["prodigal"]["threads"]
    resources:
        mem_mb_per_cpu=config["regular_memory"],
        runtime=config["prodigal"]["runtime"],
        cpus_per_task=config["prodigal"]["threads"],
        slurm_partition=config["regular_partition"],
        slurm_account=config["account"]
    conda:
        "envs/prodigal.yaml"
    shell:
        """
        mkdir -p {params.prodigal_dir}/{wildcards.sample}/{{nucl,amino,gff}}
        prodigal -i {input.reformated} \
                 -o {output.gff} \
                 -a {output.faa} \
                 -d {output.ffn} \
                 -f gff \
                 -p {params.mode}
        """

def get_samples_for_group(group):
    return config["sample_groups"][group]

rule merge_sequences:
    input:
        nucl_files=lambda wildcards: expand("{prodigal_dir}/{sample}/nucl/{sample}.nucl.ffn",
                                           prodigal_dir=config["prodigal_dir"],
                                           sample=get_samples_for_group(wildcards.group)),
        amino_files=lambda wildcards: expand("{prodigal_dir}/{sample}/amino/{sample}.faa",
                                            prodigal_dir=config["prodigal_dir"],
                                            sample=get_samples_for_group(wildcards.group))
    output:
        nucl_merged=config["cdhit_dir"] + "/{group}/{group}_nucl_merged.fasta",
        protein_merged=config["cdhit_dir"] + "/{group}/{group}_protein_merge.fasta"
    params:
        cdhit_dir=config["cdhit_dir"]
    threads: config["merge_sequences"]["threads"]
    resources:
        mem_mb_per_cpu=config["regular_memory"],
        runtime=config["merge_sequences"]["runtime"],
        cpus_per_task=config["merge_sequences"]["threads"],
        slurm_partition=config["regular_partition"],
        slurm_account=config["account"]
    shell:
        """
        mkdir -p {params.cdhit_dir}/{wildcards.group}
        cat {input.nucl_files} > {output.nucl_merged}
        cat {input.amino_files} > {output.protein_merged}
        """

rule mmseqs_clustering:
    input:
        nucl_merged=rules.merge_sequences.output.nucl_merged
    output:
        cluster_rep=config["cdhit_dir"] + "/{group}/{group}_cluster_rep_seq.fasta",
        cluster_all=config["cdhit_dir"] + "/{group}/{group}_cluster_all_seqs.fasta",
        cluster_tsv=config["cdhit_dir"] + "/{group}/{group}_cluster.tsv"
    params:
        min_seq_id=config["mmseqs2"]["min_seq_id"],
        coverage=config["mmseqs2"]["coverage"],
        threads=config["mmseqs2"]["threads"],
        prefix=config["cdhit_dir"] + "/{group}/{group}_cluster",
        tmp_dir=config["cdhit_dir"] + "/{group}/tmp"
    threads: config["mmseqs2"]["threads"]
    resources:
        mem_mb_per_cpu=config["regular_memory"],
        runtime=config["mmseqs2"]["runtime"],
        cpus_per_task=config["mmseqs2"]["threads"],
        slurm_partition=config["regular_partition"],
        slurm_account=config["account"]
    conda:
        "envs/mmseqs2.yaml"
    shell:
        """
        mmseqs easy-cluster {input.nucl_merged} \
                           {params.prefix} \
                           {params.tmp_dir} \
                           --min-seq-id {params.min_seq_id} \
                           -c {params.coverage} \
                           --threads {params.threads}
        """

rule extract_sequences:
    input:
        cluster_rep=rules.mmseqs_clustering.output.cluster_rep,
        protein_merged=rules.merge_sequences.output.protein_merged
    output:
        unigene_id=config["cdhit_dir"] + "/{group}/unigene_id.txt",
        unigene_protein=config["cdhit_dir"] + "/{group}/unigene_protein.fasta"
    threads: config["extract_sequences"]["threads"]
    resources:
        mem_mb_per_cpu=config["regular_memory"],
        runtime=config["extract_sequences"]["runtime"],
        cpus_per_task=config["extract_sequences"]["threads"],
        slurm_partition=config["regular_partition"],
        slurm_account=config["account"]
    conda:
        "envs/seqkit.yaml"
    shell:
        """
        grep '>' {input.cluster_rep} | awk '{{print $1}}' | sed 's/>//' > {output.unigene_id}
        seqkit grep -f {output.unigene_id} {input.protein_merged} -o {output.unigene_protein}
        """

rule gtdb_annotation:
    input:
        protein=rules.extract_sequences.output.unigene_protein
    output:
        annotation=config["gtdb_dir"] + "/{group}/annotation.txt"
    params:
        gtdb_db=config["databases"]["gtdb_diamond"],
        threads=config["diamond"]["threads"],
        evalue=config["diamond"]["evalue"],
        max_hits=config["diamond"]["max_hits"],
        gtdb_dir=config["gtdb_dir"]
    threads: config["diamond"]["threads"]
    resources:
        mem_mb_per_cpu=config["regular_memory"],
        runtime=config["diamond"]["runtime"],
        cpus_per_task=config["diamond"]["threads"],
        slurm_partition=config["regular_partition"],
        slurm_account=config["account"]
    conda:
        "envs/diamond.yaml"
    shell:
        """
        mkdir -p {params.gtdb_dir}/{wildcards.group}
        diamond blastp -d {params.gtdb_db} \
                      -q {input.protein} \
                      -o {output.annotation} \
                      --threads {params.threads} \
                      --evalue {params.evalue} \
                      --sensitive \
                      -k {params.max_hits} \
                      --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen qcovhsp
        """

rule eggnog_annotation:
    input:
        protein=rules.extract_sequences.output.unigene_protein
    output:
        annotations=config["functional_annotation_dir"] + "/{group}_eggnogresult.emapper.annotations"
    params:
        threads=config["eggnog"]["threads"],
        evalue=config["eggnog"]["evalue"],
        data_dir=config["databases"]["eggnog_data_dir"],
        output_prefix=config["functional_annotation_dir"] + "/{group}_eggnogresult",
        functional_annotation_dir=config["functional_annotation_dir"]
    threads: config["eggnog"]["threads"]
    resources:
        mem_mb_per_cpu=config["regular_memory"],
        runtime=config["eggnog"]["runtime"],
        cpus_per_task=config["eggnog"]["threads"],
        slurm_partition=config["regular_partition"],
        slurm_account=config["account"]
    conda:
        "envs/eggnog.yaml"
    shell:
        """
        mkdir -p {params.functional_annotation_dir}
        emapper.py -m diamond \
                  --cpu {params.threads} \
                  --seed_ortholog_evalue {params.evalue} \
                  --data_dir {params.data_dir} \
                  -i {input.protein} \
                  -o {params.output_prefix} \
                  --override
        """

rule coverm_abundance:
    input:
        r1_clean=config["fastp_dir"] + "/{sample}/{sample}_1P.fq.gz",
        r2_clean=config["fastp_dir"] + "/{sample}/{sample}_2P.fq.gz",
        reference=lambda wildcards: config["cdhit_dir"] + "/{group}/{group}_cluster_rep_seq.fasta".format(
            group=[g for g, samples in config["sample_groups"].items() if wildcards.sample in samples][0]
        )
    output:
        tpm=config["coverm_dir"] + "/{group}/coverm_{sample}.TPM.tsv"
    params:
        threads=config["coverm"]["threads"],
        min_read_percent_identity=config["coverm"]["min_read_percent_identity"],
        min_read_aligned_percent=config["coverm"]["min_read_aligned_percent"],
        min_covered_fraction=config["coverm"]["min_covered_fraction"],
        contig_end_exclusion=config["coverm"]["contig_end_exclusion"],
        coverm_dir=config["coverm_dir"]
    threads: config["coverm"]["threads"]
    resources:
        mem_mb_per_cpu=config["regular_memory"],
        runtime=config["coverm"]["runtime"],
        cpus_per_task=config["coverm"]["threads"],
        slurm_partition=config["regular_partition"],
        slurm_account=config["account"]
    conda:
        "envs/coverm.yaml"
    shell:
        """
        mkdir -p {params.coverm_dir}/{wildcards.group}
        coverm contig -1 {input.r1_clean} \
                     -2 {input.r2_clean} \
                     -t {params.threads} \
                     --reference {input.reference} \
                     -m tpm \
                     -o {output.tpm} \
                     --exclude-supplementary \
                     --min-read-percent-identity {params.min_read_percent_identity} \
                     --min-read-aligned-percent {params.min_read_aligned_percent} \
                     --min-covered-fraction {params.min_covered_fraction} \
                     --contig-end-exclusion {params.contig_end_exclusion}
        """