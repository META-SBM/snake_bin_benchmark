#!/usr/bin/python3
#base = "/mnt/disk1/RUNS/sonec_iv/porphyromonas/for_benchmark/"
base = "/mnt/disk1/PROJECTS/METAGENOMICS/Porphyromonas/external_data/PRJNA934596/reads/raw"
suffix_1 = "_R1.fastq.gz"
suffix_2 = "_R2.fastq.gz"

SAMPLES, = glob_wildcards(base +"{sample}"+suffix_1)

want_all = []
want_all.append(expand("results/kraken2/{sample}.report", sample=SAMPLES))
want_all.append(expand("results/megahit_final/{sample}.final_contigs.fa",sample=SAMPLES))
want_all.append(expand("results/metaspades/{sample}", sample=SAMPLES))
want_all.append(expand("results/metaspades_final/{sample}.assembly.fasta", sample=SAMPLES))
want_all.append(expand("results/aln/bwa2/{sample}_contigs_sorted.bam",sample=SAMPLES))
want_all.append(expand("results/aln/bwa2_spades/{sample}_contigs_sorted.bam",sample=SAMPLES))
want_all.append(expand("results/aln/strobealign/{sample}_contigs_sorted.bam",sample=SAMPLES))
want_all.append(expand("results/aln/strobealign_spades/{sample}_contigs_sorted.bam",sample=SAMPLES))
want_all.append(expand("results/binning/bwa2_semibin2/{sample}",sample=SAMPLES))
want_all.append(expand("results/binning/strobealign_semibin2/{sample}",sample=SAMPLES))
want_all.append(expand("results/binning/bwa2_semibin2_spades/{sample}",sample=SAMPLES))
want_all.append(expand("results/binning/strobealign_semibin2_spades/{sample}",sample=SAMPLES))
want_all.append(expand("results/aln/bwa2/{sample}_contigs_sorted.depth.txt",sample=SAMPLES))
want_all.append(expand("results/binning/metabat2_bwa2/{sample}/bin",sample=SAMPLES))
want_all.append(expand("results/aln/strobealign/{sample}_contigs_sorted.depth.txt",sample=SAMPLES))
want_all.append(expand("results/binning/metabat2_strobealign/{sample}/bin",sample=SAMPLES))
want_all.append(expand("results/binning_check/bwa2_semibin2/{sample}/CheckM",sample=SAMPLES))
want_all.append(expand("results/binning_check/strobealign_semibin2/{sample}/CheckM",sample=SAMPLES))

rule all:
    input: want_all

rule kraken2:
    input: read1=base+"{sample}"+suffix_1,
           read2=base+"{sample}"+suffix_2
    output: kraken2_report="results/kraken2/{sample}.report",
            kraken2_out="results/kraken2/{sample}.out"
    params: kraken2_db="/mnt/disk1/DATABASES/kraken2/k2_pluspf_20240112/"
    log:
        "results/logs/kraken2_{sample}.log"
    benchmark:
        "results/benchmarks/kraken2_{sample}.benchmark.txt"
    threads: 32
    shell: """
       (kraken2 --threads {threads} --confidence 0.5 --db {params.kraken2_db} {input.read1} {input.read2} --use-names --report {output.kraken2_report} --output {output.kraken2_out}) > {log}
    """
rule megahit:
    input: read1=base+"{sample}"+suffix_1,
           read2=base+"{sample}"+suffix_2
    output: megahit_asm_folder=directory("results/megahit_asm_{sample}/"),
            megahit_asm="results/megahit_final/{sample}.final_contigs.fa",
    threads: 32
    log:
        "results/logs/megahit_{sample}.log"
    benchmark:
        "results/benchmarks/megahit_{sample}.benchmark.txt"
    shell: """
        (megahit -t {threads} -1 {input.read1} -2 {input.read2} -o {output.megahit_asm_folder}) > {log}
        mkdir -p results/megahit_final/
        cp {output.megahit_asm_folder}/final.contigs.fa {output.megahit_asm}
    """
rule spades:
    input: read1=base+"{sample}"+suffix_1,
           read2=base+"{sample}"+suffix_2
    output: spades_res=directory("results/metaspades/{sample}"),
            spades_asm="results/metaspades_final/{sample}.assembly.fasta",
    log:
        "results/logs/spades_{sample}.log"
    benchmark:
        "results/benchmarks/spades_{sample}.benchmark.txt"
    threads: 32
    shell: """
            (spades.py -t {threads} --meta -1 {input.read1} -2 {input.read2} -o {output.spades_res}) > {log}
            cp {output.spades_res}/contigs.fasta {output.spades_asm}
    """

rule bwa2_align_to_megahit_asm:
    input: megahit_asm="results/megahit_final/{sample}.final_contigs.fa",
           read1=base+"{sample}"+suffix_1,
           read2=base+"{sample}"+suffix_2,
    output: bwa2_aln_sort="results/aln/bwa2/{sample}_contigs_sorted.bam",
    params: prefix=directory("results/aln/bwa2/idx/"),
    threads: 32
    log:
        "results/logs/bwa2_{sample}.log"
    benchmark:
        "results/benchmarks/bwa2_{sample}.benchmark.txt"
    shell: """
          mkdir -p {params.prefix}
          bwa-mem2 index -p {params.prefix}/{wildcards.sample} {input.megahit_asm}
          (bwa-mem2 mem -t {threads} {params.prefix}/{wildcards.sample} {input.read1} {input.read2} | samtools view -b -@ {threads} - | samtools sort -@ {threads} - > {output.bwa2_aln_sort}) > {log}
    """
rule strobealign_to_megahit_asm:
    input: megahit_asm="results/megahit_final/{sample}.final_contigs.fa",
           read1=base+"{sample}"+suffix_1,
           read2=base+"{sample}"+suffix_2,
    output: strobealign_sorted="results/aln/strobealign/{sample}_contigs_sorted.bam",
    threads: 32
    log:
        "results/logs/strobealign_{sample}.log"
    benchmark:
        "results/benchmarks/strobealign_{sample}.benchmark.txt"
    shell: """
          (strobealign -t {threads} {input.megahit_asm} {input.read1} {input.read2} | samtools sort -@ {threads} - > {output.strobealign_sorted}) > {log}
    """
rule strobealign_aemb_megahit:
    input: megahit_asm="results/megahit_final/{sample}.final_contigs.fa",
           read1=base+"{sample}"+suffix_1,
           read2=base+"{sample}"+suffix_2,
    output: strobealign_aemb_megahit="results/aln/strobealign/{sample}_abundance.tsv",
    threads: 32
    log:
        "results/logs/strobealign_aemb_{sample}.log"
    benchmark:
        "results/benchmarks/strobealign_aemb_{sample}.benchmark.txt"
    shell: """
          (strobealign -t {threads} --aemb {input.megahit_asm} {input.read1} {input.read2} > {output.strobealign_aemb_megahit}) > {log}
    """
rule bwa2_align_to_spades_asm:
    input: spades_asm="results/metaspades_final/{sample}.assembly.fasta",
           read1=base+"{sample}"+suffix_1,
           read2=base+"{sample}"+suffix_2,
    output: bwa2_aln_sort="results/aln/bwa2_spades/{sample}_contigs_sorted.bam",
    params: prefix=directory("results/aln/bwa2_spades/idx/"),
    threads: 32
    log:
        "results/logs/bwa2_spades_{sample}.log"
    benchmark:
        "results/benchmarks/bwa2_spades_{sample}.benchmark.txt"
    shell: """
          mkdir -p {params.prefix}
          bwa-mem2 index -p {params.prefix}/{wildcards.sample} {input.spades_asm}
          (bwa-mem2 mem -t {threads} {params.prefix}/{wildcards.sample} {input.read1} {input.read2} | samtools view -b -@ {threads} - | samtools sort -@ {threads} - > {output.bwa2_aln_sort}) > {log}
    """
rule strobealign_to_spades_asm:
    input: spades_asm="results/metaspades_final/{sample}.assembly.fasta",
           read1=base+"{sample}"+suffix_1,
           read2=base+"{sample}"+suffix_2,
    output: strobealign_sorted="results/aln/strobealign_spades/{sample}_contigs_sorted.bam",
    threads: 32
    log:
        "results/logs/strobealign_spades{sample}.log"
    benchmark:
        "results/benchmarks/strobealign_spades_{sample}.benchmark.txt"
    shell: """
          (strobealign -t {threads} {input.spades_asm} {input.read1} {input.read2} | samtools sort -@ {threads} - > {output.strobealign_sorted}) > {log}
    """
rule strobealign_to_spades_asm_aemb:
    input: spades_asm="results/metaspades_final/{sample}.assembly.fasta",
           read1=base+"{sample}"+suffix_1,
           read2=base+"{sample}"+suffix_2,
    output:strobealign_aemb_spades="results/aln/strobealign_spades/{sample}_abundance.tsv",
    threads: 32
    log:
        "results/logs/strobealign_spades_aemb_{sample}.log"
    benchmark:
        "results/benchmarks/strobealign_spades_aemb_{sample}.benchmark.txt"
    shell: """
          (strobealign -t {threads} --aemb {input.spades_asm} {input.read1} {input.read2} > {output.strobealign_aemb_spades}) > {log}
    """

rule semibin2_bwa2:
    input: bwa2_aln_sort="results/aln/bwa2/{sample}_contigs_sorted.bam",
           megahit_asm="results/megahit_final/{sample}.final_contigs.fa",
    output: semibin2_bwa2_res = directory("results/binning/bwa2_semibin2/{sample}"),
    threads: 32
    log:
        "results/logs/semibin2_bwa2_{sample}.log"
    benchmark:
        "results/benchmarks/semibin2_bwa2_{sample}.benchmark.txt"

    shell: """
           (SemiBin2 single_easy_bin --threads {threads} --input-fasta {input.megahit_asm} --input-bam {input.bwa2_aln_sort} --output {output.semibin2_bwa2_res} --environment human_oral) > {log}
    """

rule semibin2_strobealign:
    input: strobealign_sorted="results/aln/strobealign/{sample}_contigs_sorted.bam",
           megahit_asm="results/megahit_final/{sample}.final_contigs.fa",
    output: semibin2_strobealign_res = directory("results/binning/strobealign_semibin2/{sample}")
    threads: 32
    log:
        "results/logs/semibin2_strobealign_{sample}.log"
    benchmark:
        "results/benchmarks/semibin2_strobealign_{sample}.benchmark.txt"

    shell: """
          (SemiBin2 single_easy_bin --threads {threads} --input-fasta {input.megahit_asm} --input-bam {input.strobealign_sorted} --output {output.semibin2_strobealign_res} --environment human_oral) > {log}
    """
rule semibin2_bwa2_spades:
    input: bwa2_aln_sort="results/aln/bwa2_spades/{sample}_contigs_sorted.bam",
           spades_asm="results/metaspades_final/{sample}.assembly.fasta",
    output: semibin2_bwa2_res = directory("results/binning/bwa2_semibin2_spades/{sample}"),
    threads: 32
    log:
        "results/logs/semibin2_bwa2_spades_{sample}.log"
    benchmark:
        "results/benchmarks/semibin2_bwa2_spades_{sample}.benchmark.txt"

    shell: """
           (SemiBin2 single_easy_bin --threads {threads} --input-fasta {input.spades_asm} --input-bam {input.bwa2_aln_sort} --output {output.semibin2_bwa2_res} --environment human_oral) > {log}
    """

rule semibin2_strobealign_spades:
    input: strobealign_sorted="results/aln/strobealign_spades/{sample}_contigs_sorted.bam",
           spades_asm="results/metaspades_final/{sample}.assembly.fasta",
    output: semibin2_strobealign_res = directory("results/binning/strobealign_semibin2_spades/{sample}")
    threads: 32
    log:
        "results/logs/semibin2_strobealign_spades_{sample}.log"
    benchmark:
        "results/benchmarks/semibin2_strobealign_spades_{sample}.benchmark.txt"

    shell: """
          (SemiBin2 single_easy_bin --threads {threads} --input-fasta {input.spades_asm} --input-bam {input.strobealign_sorted} --output {output.semibin2_strobealign_res} --environment human_oral) > {log}
    """
rule metabat2_bwa2:
    input: bwa2_aln_sort="results/aln/bwa2/{sample}_contigs_sorted.bam",
           megahit_asm="results/megahit_final/{sample}.final_contigs.fa",
    output: bwa2_aln_depth = "results/aln/bwa2/{sample}_contigs_sorted.depth.txt",
            metabat2_bwa2_res = directory("results/binning/metabat2_bwa2/{sample}/bin")
    threads: 32
    params: path="/mnt/disk1/RUNS/sonec_iv/porphyromonas:/data"
    log:
        "results/logs/metabat2_bwa2_{sample}.log"
    benchmark:
        "results/benchmarks/metabat2_bwa2_{sample}.benchmark.txt"
    shell: """
          jgi_summarize_bam_contig_depths --outputDepth {output.bwa2_aln_depth} {input.bwa2_aln_sort}
          (metabat2 -t {threads} -i {input.megahit_asm} -a {output.bwa2_aln_depth} -o {output.metabat2_bwa2_res}) > {log}
    """

rule metabat2_strobealign:
    input: strobealign_sorted="results/aln/strobealign/{sample}_contigs_sorted.bam",
           megahit_asm="results/megahit_final/{sample}.final_contigs.fa",
    output: strobealign_depth = "results/aln/strobealign/{sample}_contigs_sorted.depth.txt",
            metabat2_strobealign_res = directory("results/binning/metabat2_strobealign/{sample}/bin")
    threads:32
    params: path="/mnt/disk1/RUNS/sonec_iv/porphyromonas:/data"
    log:
        "results/logs/metabat2_strobealign_{sample}.log"
    benchmark:
        "results/benchmarks/metabat2_strobealign_{sample}.benchmark.txt"
    shell: """
           jgi_summarize_bam_contig_depths --outputDepth {output.strobealign_depth} {input.strobealign_sorted}
           (metabat2 -t {threads} -i {input.megahit_asm} -a {output.strobealign_depth} -o {output.metabat2_strobealign_res}) > {log}
    """

rule checkM_bwa2:
    input: semibin2_bwa2_res="results/binning/bwa2_semibin2/{sample}",
    output: checkM_bwa2_txt="results/binning_check/bwa2_semibin2/{sample}/CheckM/checkm.txt",
            checkM_bwa2_res=directory("results/binning_check/bwa2_semibin2/{sample}/CheckM"),
    threads: 32
    params: bwa2_lineage_ms="results/binning_check/bwa2_semibin2/{sample}/CheckM/lineage.ms",
    log:
        "results/logs/checkM_bwa2_{sample}.log"
    benchmark:
        "results/benchmarks/checkM_bwa2_{sample}.benchmark.txt"

    shell:"""
          (checkm lineage_wf {input.semibin2_bwa2_res}/output_bins -t {threads} -x fa.gz {output.checkM_bwa2_res}) > {log}
          (checkm qa {params.bwa2_lineage_ms} {output.checkM_bwa2_res} -f {output.checkM_bwa2_txt} -o 2 --tab_table -t {threads}) >> {log}
    """

rule checkM_strobealign:
    input: semibin2_strobealign_res="results/binning/strobealign_semibin2/{sample}",
    output: checkM_strobealign_res=directory("results/binning_check/strobealign_semibin2/{sample}/CheckM"),
            checkM_strobealign_txt="results/binning_check/strobealign_semibin2/{sample}/CheckM/checkm.txt",
    threads: 32
    params: strobealign_lineage_ms="results/binning_check/strobealign_semibin2/{sample}/CheckM/lineage.ms",
    log:
        "results/logs/checkM_strobealign_{sample}.log"
    benchmark:
        "results/benchmarks/checkM_strobealign_{sample}.benchmark.txt"

    shell:"""
          (checkm lineage_wf {input.semibin2_strobealign_res}/output_bins -t {threads} -x fa.gz {output.checkM_strobealign_res}) > {log}
          (checkm qa {params.strobealign_lineage_ms} {output.checkM_strobealign_res} -f {output.checkM_strobealign_txt} -o 2 --tab_table -t {threads}) >> {log}
    """
