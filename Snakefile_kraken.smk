import os

configfile: "config/config.yaml"

rule all:
    input:
        "result/96_taxpasta/report_bracken_with_lineage.tsv"

'''
rule all:
    input:
        config["databases"]["kraken2_db"] + "/bracken_build.done",
        expand(config["kraken2_dir"] + "/{sample}/{sample}_bracken_{level}.report",
               sample=config["samples"], level=['phylum','classes','orders','family','genus','species']),
        "result/98_bracken_kraken",
        "result/97_braken_report",
        "result/96_taxpasta/report_bracken_with_lineage.tsv"
'''




# rule kraken2:
#     input:
#         r1_clean=config["fastp_dir"] + "/{sample}/{sample}_1P.fq.gz",
#         r2_clean=config["fastp_dir"] + "/{sample}/{sample}_2P.fq.gz"
#     output:
#         report=config["kraken2_dir"] + "/{sample}/{sample}_kraken2.report",
#         out=config["kraken2_dir"] + "/{sample}/{sample}_kraken2.out"
#     log:
#         log=config["log_dir"] + "/kraken2/{sample}_kraken2.log",
#         err=config["log_dir"] + "/kraken2/{sample}_kraken2.err"
#     conda:
#         "envs/kraken2.yaml"
#     resources:
#         mem_mb_per_cpu=config["regular_memory"],
#         slurm_partition=config["regular_partition"],
#         runtime=config["kraken2"]["runtime"],
#         cpus_per_task=config["kraken2"]["threads"],
#         slurm_account=config["account"]
#     threads: config["kraken2"]["threads"]
#     params:
#         out_dir=config["kraken2_dir"] + "/{sample}",
#         db=config["databases"]["kraken2_db"],
#         confidence=config["kraken2"]["confidence"]

#     shell:
#         """
#         mkdir -p {params.out_dir} && \
#         mkdir -p $(dirname {log.log}) && \
#         kraken2 --db {params.db} \
#                 --threads {threads} \
#                 --confidence {params.confidence} \
#                 --report {output.report} \
#                 --paired {input.r1_clean} {input.r2_clean} \
#                 --output {output.out} \
#                 > {log.log} 2> {log.err}
#         """
# rule bracken_build:
#     input:
#         db=config["databases"]["kraken2_db"]
#     output:
#         marker=config["databases"]["kraken2_db"] + "/bracken_build.done"
#     log:
#         log=config["log_dir"] + "/bracken_build/bracken_build.log",
#         err=config["log_dir"] + "/bracken_build/bracken_build.err"
#     conda:
#         "envs/kraken2.yaml"
#     #resources are for slurm, using the same as kraken2
#     resources:
#         mem_mb_per_cpu=config["regular_memory"],
#         slurm_partition=config["regular_partition"],
#         runtime=config["kraken2"]["runtime"],
#         cpus_per_task=config["kraken2"]["threads"],
#         slurm_account=config["account"]
#     params:
#         db=config["databases"]["kraken2_db"],
#         kmer_length=config["bracken"]["kmer_length"],
#         read_length=config["bracken"]["read_length"]
#     threads: config["bracken"]["threads"]
#     shell:
#         """
#         mkdir -p $(dirname {log.log}) && \
#         bracken-build -d {input.db} -t {threads} -k 35 -l 150 > {log.log} 2> {log.err} && \
#         touch {output.marker}
#         """
# rule bracken:
#     input:
#         report=config["kraken2_dir"] + "/{sample}/{sample}_kraken2.report",
#         db=config["databases"]["kraken2_db"],
#         bracken_db=config["databases"]["kraken2_db"] + "/bracken_build.done"
#     output:
#         bracken_report_P=config["kraken2_dir"] + "/{sample}/{sample}_bracken_phylum.report",
#         bracken_report_C=config["kraken2_dir"] + "/{sample}/{sample}_bracken_classes.report",
#         bracken_report_O=config["kraken2_dir"] + "/{sample}/{sample}_bracken_orders.report",
#         bracken_report_F=config["kraken2_dir"] + "/{sample}/{sample}_bracken_family.report",
#         bracken_report_G=config["kraken2_dir"] + "/{sample}/{sample}_bracken_genus.report",
#         bracken_report_S=config["kraken2_dir"] + "/{sample}/{sample}_bracken_species.report"
#     log:
#         log=config["log_dir"] + "/bracken/{sample}_bracken.log",
#         err=config["log_dir"] + "/bracken/{sample}_bracken.err"
#     conda:
#         "envs/kraken2.yaml"
#     resources:
#         mem_mb_per_cpu=config["regular_memory"],
#         slurm_partition=config["regular_partition"],
#         runtime=120,  # minutes
#         cpus_per_task=2,
#         slurm_account=config["account"]
#     threads: 2
#     params:
#         db=config["databases"]["kraken2_db"],
#         kmer_length=config["bracken"]["kmer_length"],
#         read_length=config["bracken"]["read_length"]
#     shell:
#         """
#         mkdir -p $(dirname {log.log}) && \
#         bracken -d /scratch/qou/database/kraken2 -i {input.report} -o {output.bracken_report_P} -r 150 -l P -t 2 && \
#         bracken -d /scratch/qou/database/kraken2 -i {input.report} -o {output.bracken_report_C} -r 150 -l C -t 2 && \
#         bracken -d /scratch/qou/database/kraken2 -i {input.report} -o {output.bracken_report_O} -r 150 -l O -t 2 && \
#         bracken -d /scratch/qou/database/kraken2 -i {input.report} -o {output.bracken_report_F} -r 150 -l F -t 2 && \
#         bracken -d /scratch/qou/database/kraken2 -i {input.report} -o {output.bracken_report_G} -r 150 -l G -t 2 && \
#         bracken -d /scratch/qou/database/kraken2 -i {input.report} -o {output.bracken_report_S} -r 150 -l S -t 2 \
#                 > {log.log} 2> {log.err}
#         """

# rule merge_bracken:
#     input:
#         P_reports=expand(config["kraken2_dir"] + "/{sample}/{sample}_bracken_phylum.report", sample=config["samples"]),
#         C_reports=expand(config["kraken2_dir"] + "/{sample}/{sample}_bracken_classes.report", sample=config["samples"]),
#         O_reports=expand(config["kraken2_dir"] + "/{sample}/{sample}_bracken_orders.report", sample=config["samples"]),
#         F_reports=expand(config["kraken2_dir"] + "/{sample}/{sample}_bracken_family.report", sample=config["samples"]),
#         G_reports=expand(config["kraken2_dir"] + "/{sample}/{sample}_bracken_genus.report", sample=config["samples"]),
#         S_reports=expand(config["kraken2_dir"] + "/{sample}/{sample}_bracken_species.report", sample=config["samples"]),
#         P_kraken_reports=expand(config["kraken2_dir"] + "/{sample}/{sample}_kraken2_bracken_phylums.report", sample=config["samples"]),
#         C_kraken_reports=expand(config["kraken2_dir"] + "/{sample}/{sample}_kraken2_bracken_classes.report", sample=config["samples"]),
#         O_kraken_reports=expand(config["kraken2_dir"] + "/{sample}/{sample}_kraken2_bracken_orders.report", sample=config["samples"]),
#         F_kraken_reports=expand(config["kraken2_dir"] + "/{sample}/{sample}_kraken2_bracken_families.report", sample=config["samples"]),
#         G_kraken_reports=expand(config["kraken2_dir"] + "/{sample}/{sample}_kraken2_bracken_genuses.report", sample=config["samples"]),
#         S_kraken_reports=expand(config["kraken2_dir"] + "/{sample}/{sample}_kraken2_bracken_species.report", sample=config["samples"])
#     output:
#         directory("result/98_bracken_kraken")
#     shell:
#         """
#         mkdir -p {output}/krakenstykle/P {output}/krakenstykle/C {output}/krakenstykle/O {output}/krakenstykle/F {output}/krakenstykle/G {output}/krakenstykle/S && \
#         mkdir -p {output}/brackenstykle/P {output}/brackenstykle/C {output}/brackenstykle/O {output}/brackenstykle/F {output}/brackenstykle/G {output}/brackenstykle/S && \
#         for file in {input.P_kraken_reports}; do
#             cp $file {output}/krakenstykle/P/
#         done && \
#         for file in {input.C_kraken_reports}; do
#             cp $file {output}/krakenstykle/C/
#         done && \
#         for file in {input.O_kraken_reports}; do
#             cp $file {output}/krakenstykle/O/
#         done && \
#         for file in {input.F_kraken_reports}; do
#             cp $file {output}/krakenstykle/F/
#         done && \
#         for file in {input.G_kraken_reports}; do
#             cp $file {output}/krakenstykle/G/
#         done && \
#         for file in {input.S_kraken_reports}; do
#             cp $file {output}/krakenstykle/S/
#         done && \
#         for file in {input.P_reports}; do
#             cp $file {output}/brackenstykle/P/
#         done && \
#         for file in {input.C_reports}; do
#             cp $file {output}/brackenstykle/C/
#         done && \
#         for file in {input.O_reports}; do
#             cp $file {output}/brackenstykle/O/
#         done && \
#         for file in {input.F_reports}; do
#             cp $file {output}/brackenstykle/F/
#         done && \
#         for file in {input.G_reports}; do
#             cp $file {output}/brackenstykle/G/
#         done && \
#         for file in {input.S_reports}; do
#             cp $file {output}/brackenstykle/S/
#         done
#         """

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

# abortion rule taxpasta_merge_bracken:
#     input:
#         reports=expand(config["kraken2_dir"] + "/{sample}/{sample}_bracken_species.report", sample=config["samples"])
#     output:
#         merged="result/96_taxpasta/report_bracken_with_lineage.tsv"
#     log:
#         log=config["log_dir"] + "/taxpasta/taxpasta_merge.log",
#         err=config["log_dir"] + "/taxpasta/taxpasta_merge.err"
#     conda:
#         "envs/taxpasta.yaml"
#     params:
#         taxonomy=config["databases"]["taxa_db"]
#     shell:
#         """
#         mkdir -p $(dirname {output.merged}) && \
#         mkdir -p $(dirname {log.log}) && \
#         taxpasta merge -p bracken -o {output.merged} --taxonomy {params.taxonomy} --add-lineage {input.reports} > {log.log} 2> {log.err}
#         """

