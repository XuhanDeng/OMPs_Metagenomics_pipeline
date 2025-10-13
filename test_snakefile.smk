configfile: "test_config.yaml"

# Test Snakefile for debugging spades_assembly_and_reformat rule
# This focuses on testing just the spades assembly and reformatting steps

rule all:
    input:
        expand(config["reformated_scaffolds_dir"] + "/{sample}/{sample}_reformated.fa",
               sample=config["samples"])

# Dummy rule to create test input files
rule create_test_fastq:
    output:
        r1=config["fastp_dir"] + "/{sample}/{sample}_1P.fq.gz",
        r2=config["fastp_dir"] + "/{sample}/{sample}_2P.fq.gz"
    params:
        fastp_dir=config["fastp_dir"]
    shell:
        """
        mkdir -p {params.fastp_dir}/{wildcards.sample}

        # Create small test FASTQ files with synthetic sequences
        # R1 file
        cat > {params.fastp_dir}/{wildcards.sample}/temp_r1.fq << 'EOF'
@read1/1
ATGCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATC
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
@read2/1
GCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGA
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
@read3/1
TCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGT
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
@read4/1
CGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGAT
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
@read5/1
GATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATC
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
EOF

        # R2 file
        cat > {params.fastp_dir}/{wildcards.sample}/temp_r2.fq << 'EOF'
@read1/2
GATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGAT
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
@read2/2
TCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGA
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
@read3/2
CGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGAT
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
@read4/2
GATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATC
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
@read5/2
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
EOF

        # Compress the files
        gzip {params.fastp_dir}/{wildcards.sample}/temp_r1.fq
        gzip {params.fastp_dir}/{wildcards.sample}/temp_r2.fq

        # Rename to expected names
        mv {params.fastp_dir}/{wildcards.sample}/temp_r1.fq.gz {output.r1}
        mv {params.fastp_dir}/{wildcards.sample}/temp_r2.fq.gz {output.r2}
        """

# This is the rule we want to test - copied exactly from your main Snakefile
rule spades_assembly_and_reformat:
    input:
        r1=rules.create_test_fastq.output.r1,
        r2=rules.create_test_fastq.output.r2
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
        mem_mb_per_cpu=config["spades"]["memory"] * 1000 // config["spades"]["threads"],
        runtime=config["spades"]["runtime"] + config["reformat_scaffolds"]["runtime"],
        cpus_per_task=config["spades"]["threads"],
        slurm_partition=config["spades"]["parition"],
        slurm_account=config["account"]
    conda:
        "envs/seqkit-spade.yaml"
    log:
        out="test_log/2_spades/{sample}_spades.out",
        err="test_log/2_spades/{sample}_spades.err"
    shell:
        """
        # SPAdes assembly
        mkdir -p {params.log_dir}/2_spades
        mkdir -p {params.spades_dir}/{wildcards.sample}

        echo "Starting SPAdes assembly for {wildcards.sample}..."
        echo "Command: spades.py --meta -o {params.spades_outdir} -1 {input.r1} -2 {input.r2} -t {params.spades_threads} -m {params.spades_memory} -k {params.kmers} --only-assembler"

        spades.py --meta \
                  -o {params.spades_outdir} \
                  -1 {input.r1} \
                  -2 {input.r2} \
                  -t {params.spades_threads} \
                  -m {params.spades_memory} \
                  -k {params.kmers} \
                  --only-assembler > {log.out} 2> {log.err}

        echo "SPAdes completed. Checking output..."

        # Check if contigs.fasta was created
        if [ ! -f {params.spades_outdir}/contigs.fasta ]; then
            echo "ERROR: SPAdes did not produce contigs.fasta"
            echo "SPAdes output directory contents:"
            ls -la {params.spades_outdir}/
            echo "SPAdes stdout log:"
            cat {log.out}
            echo "SPAdes stderr log:"
            cat {log.err}
            exit 1
        fi

        echo "✓ contigs.fasta found"

        # Reformat scaffolds
        mkdir -p {params.reformated_scaffolds_dir}/{wildcards.sample}
        echo "Reformating scaffolds directory: {params.reformated_scaffolds_dir}/{wildcards.sample}"
        echo "Source contigs file: {params.spades_outdir}/contigs.fasta"

        # Copy original contigs
        cp {params.spades_outdir}/contigs.fasta {output.original}
        echo "✓ Original contigs copied"

        # Check if seqkit is available
        if ! command -v seqkit &> /dev/null; then
            echo "ERROR: seqkit command not found"
            exit 1
        fi

        # Reformat contig headers and filter by length
        echo "Filtering sequences ≥ {params.min_length}bp and renaming..."
        seqkit seq -m {params.min_length} {params.spades_outdir}/contigs.fasta | \
        seqkit replace -p "(.+)" -r "{wildcards.sample}_${{1}}" > {output.reformated}

        echo "✓ Sequences reformatted"

        # Create renaming table
        seqkit fx2tab {output.reformated} | cut -f1 | sed 's/>//' | nl -nln | \
        awk 'BEGIN{{OFS="\t"}} {{print $2, "{wildcards.sample}_" $1}}' > {output.renaming}

        echo "✓ Renaming table created"

        # Show statistics
        original_count=$(grep -c "^>" {output.original} || echo "0")
        reformated_count=$(grep -c "^>" {output.reformated} || echo "0")

        echo "Assembly statistics:"
        echo "  Original contigs: $original_count"
        echo "  Filtered contigs (≥{params.min_length}bp): $reformated_count"

        # Cleanup intermediate files (commented out for debugging)
        # rm -rf {params.spades_dir}/{wildcards.sample}

        echo "✓ Rule completed successfully!"
        """