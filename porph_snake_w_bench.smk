#!/usr/bin/python3
# modify path to folder if necessary
base = "/mnt/disk1/RUNS/sonec_iv/porhyromonas/for_benchmark/"
suffix_1 = "_R1.fastq.gz"
suffix_2 = "_R2.fastq.gz"

SAMPLES, = glob_wildcards(base +"{sample}"+suffix_1)

want_all = []
want_all.append(expand("results/kraken2/{sample}.report", sample=SAMPLES))
want_all.append(expand("results/megahit_final/{sample}.final_contigs.fa",sample=SAMPLES))
want_all.append(expand("results/aln/bwa2/{sample}_contigs_sorted.bam",sample=SAMPLES))
want_all.append(expand("results/aln/strobealign/{sample}_contigs_sorted.bam",sample=SAMPLES))
want_all.append(expand("results/binning/bwa2_semibin2/{sample}",sample=SAMPLES))
want_all.append(expand("results/binning/strobealign_semibin2/{sample}",sample=SAMPLES))
want_all.append(expand("results/binning_check/bwa2_semibin2/{sample}/CheckM",sample=SAMPLES))
want_all.append(expand("results/binning_check/strobealign_semibin2/{sample}/CheckM",sample=SAMPLES))

rule all:
    input: want_all

rule kraken2:
    input: read1=base+"{sample}"+suffix_1,
           read2=base+"{sample}"+suffix_2
    output: kraken2_report="results/kraken2/{sample}.report",
            kraken2_out="results/kraken2/{sample}.out"
    params: kraken2_db="/mnt/disk1/DATABASES/kraken2/pluspf"
    log:
        "results/logs/kraken2_{sample}.log"
    benchmark:
        "results/benchmarks/kraken2_{sample}.benchmark.txt"
    threads: 16
    shell: """
        kraken2 --threads {threads} --confidence 0.5 --db {params.kraken2_db} {input.read1} {input.read2} --use-names --report {output.kraken2_report} --output {output.kraken2_out}
    """
# megahit rule may be unstable, proceed with caution, see issue here: https://stackoverflow.com/questions/69929526/snakemake-megahit-output-issue
# nowadays megahit results folder created inside rule with mkdir
rule megahit:
    input: read1=base+"{sample}"+suffix_1,
           read2=base+"{sample}"+suffix_2
    output: megahit_asm_folder = directory("results/megahit_asm_{sample}/"),
            megahit_asm = "results/megahit_final/{sample}.final_contigs.fa",
    threads: 16
    log:
        "results/logs/megahit_{sample}.log"
    benchmark:
        "results/benchmarks/megahit_{sample}.benchmark.txt"
    shell: """
        megahit -t {threads} -1 {input.read1} -2 {input.read2} -o {output.megahit_asm_folder}
        mkdir -p results/megahit_final/
        cp {output.megahit_asm_folder}/final.contigs.fa {output.megahit_asm}
    """
rule bwa2_align_to_asm:
    input: megahit_asm="results/megahit_final/{sample}.final_contigs.fa",
           read1=base+"{sample}"+suffix_1,
           read2=base+"{sample}"+suffix_2,
    output: bwa2_aln_sort="results/aln/bwa2/{sample}_contigs_sorted.bam",
    params: prefix=directory("results/aln/bwa2/idx/"),
    threads: 16
    log:
        "results/logs/bwa2_{sample}.log"
    benchmark:
        "results/benchmarks/bwa2_{sample}.benchmark.txt"
    shell: """
          mkdir -p {params.prefix}
          bwa-mem2 index -p {params.prefix}/{wildcards.sample} {input.megahit_asm}
          bwa-mem2 mem -t {threads} {params.prefix}/{wildcards.sample} {input.read1} {input.read2} \
          | samtools view -b -@ {threads} - | samtools sort -@ {threads} - > {output.bwa2_aln_sort}
    """
rule strobealign_to_asm:
    input: megahit_asm="results/megahit_final/{sample}.final_contigs.fa",
           read1=base+"{sample}"+suffix_1,
           read2=base+"{sample}"+suffix_2,
    output: strobealign_sorted="results/aln/strobealign/{sample}_contigs_sorted.bam",
            strobealign_aemb_res="results/aln/strobealign/{sample}_abundance.tsv",
    threads: 16
    log:
        "results/logs/strobealign_{sample}.log"
    benchmark:
        "results/benchmarks/strobealign_{sample}.benchmark.txt"
    shell: """
          strobealign -t {threads} {input.megahit_asm} \
          {input.read1} {input.read2} | samtools sort -@ {threads} - > {output.strobealign_sorted}
          strobealign -t {threads} --aemb {input.megahit_asm} \
          {input.read1} {input.read2} > {output.strobealign_aemb_res}
    """
rule semibin2_bwa2:
    input: bwa2_aln_sort="results/aln/bwa2/{sample}_contigs_sorted.bam",
           megahit_asm="results/megahit_final/{sample}.final_contigs.fa",
    output: semibin2_bwa2_res = directory("results/binning/bwa2_semibin2/{sample}"),
    threads: 16
    log:
        "results/logs/semibin2_bwa2_{sample}.log"
    benchmark:
        "results/benchmarks/semibin2_bwa2_{sample}.benchmark.txt"

    shell: """
           SemiBin2 single_easy_bin --threads {threads} --input-fasta {input.megahit_asm} \
           --input-bam {input.bwa2_aln_sort} --output {output.semibin2_bwa2_res} --environment human_oral
    """

rule semibin2_strobealign:
    input: strobealign_sorted="results/aln/strobealign/{sample}_contigs_sorted.bam",
           strobealign_aemb_res="results/aln/strobealign/{sample}_abundance.tsv",
           megahit_asm="results/megahit_final/{sample}.final_contigs.fa",
    output: semibin2_strobealign_res = directory("results/binning/strobealign_semibin2/{sample}")
    threads: 16
    log:
        "results/logs/semibin2_strobealign_{sample}.log"
    benchmark:
        "results/benchmarks/semibin2_strobealign_{sample}.benchmark.txt"

    shell: """
           SemiBin2  single_easy_bin --threads {threads} --input-fasta {input.megahit_asm} \
           --input-bam {input.strobealign_sorted} --abundance {input.strobealign_aemb_res} \
           --output {output.semibin2_strobealign_res} --environment human_oral
    """
rule checkM_bwa2:
    input: semibin2_bwa2_res="results/binning/bwa2_semibin2/{sample}",
    output: checkM_bwa2_txt="results/binning_check/bwa2_semibin2/{sample}/CheckM/checkm.txt",
            checkM_bwa2_res=directory("results/binning_check/bwa2_semibin2/{sample}/CheckM"),
    threads:16
    params: bwa2_lineage_ms="results/binning_check/bwa2_semibin2/{sample}/CheckM/lineage.ms",
    log:
        "results/logs/checkM_bwa2_{sample}.log"
    benchmark:
        "results/benchmarks/checkM_bwa2_{sample}.benchmark.txt"

    shell:"""
          checkm lineage_wf {input.semibin2_bwa2_res}/sorted_bins -t {threads} -x fa {output.checkM_bwa2_res}
          checkm qa {params.bwa2_lineage_ms} {output.checkM_bwa2_res} -f {output.checkM_bwa2_txt} -o 2 --tab_table -t {threads}
    """

rule checkM_strobealign:
    input: semibin2_strobealign_res="results/binning/strobealign_semibin2/{sample}",
    output: checkM_strobealign_res=directory("results/binning_check/strobealign_semibin2/{sample}/CheckM"),
            checkM_strobealign_txt="results/binning_check/strobealign_semibin2/{sample}/CheckM/checkm.txt",
    threads:16
    params: strobealign_lineage_ms="results/binning_check/strobealign_semibin2/{sample}/CheckM/lineage.ms",
    log:
        "results/logs/checkM_strobealign_{sample}.log"
    benchmark:
        "results/benchmarks/checkM_strobealign_{sample}.benchmark.txt"

    shell:"""
          checkm lineage_wf {input.semibin2_strobealign_res}/sorted_bins -t {threads} -x fa {output.checkM_strobealign_res}
          checkm qa {params.strobealign_lineage_ms} {output.checkM_strobealign_res} -f {output.checkM_strobealign_txt} -o 2 --tab_table -t {threads}
    """
