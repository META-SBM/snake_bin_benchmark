# snake_bin_benchmark
Simple Snakemake pipeline for WGS metagenomics assembly, binning and benchmarking

**Prerequisites:**  
 - bwa-mem2;
 - strobealign;
 - semibin2;
 - checkM;
 - kraken2 and its DB;
 - megahit.

All tools should be installed in conda fresh env.

**Ready:**

 - assembly with MEGAHIT;  
 - kraken2 with pluspf base;  
 - binning with SemiBin2 in mode ```human_oral```;  
 - realignment of reads to assembly with bwa-mem2 and strobealign;  
 - storbealign with AEMB mode;  
 - saving all logs and benchmarks files.

**Not ready:**  
- automatic creation of benchmark plots(WIP).
