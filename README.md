# snake_bin_benchmark
Simple Snakemake pipeline for WGS metagenomics assembly, binning, quality assessing along with benchmarking all programs used in this pipeline.

**Prerequisites:**  
 - bwa-mem2;
 - strobealign;
 - semibin2;
 - checkM;
 - kraken2 and its DB;
 - megahit;
 - SPAdes (added 09/07/24);
 - MetaBat2 (added 09/07/24).

All tools(except MetaBat2) should be installed in conda fresh env. Use MetaBat2 as a Docker image. 

**Ready:**
 - assembly with MEGAHIT;  
 - kraken2 with pluspf base;  
 - binning with SemiBin2 in mode ```human_oral```;  
 - realignment of reads to assembly with bwa-mem2 and strobealign;  
 - storbealign with AEMB mode;  
 - saving all logs and benchmarks files.

**UPD 09/07/2024**  
  - added MetaBat2;
  - added SPAdes;
  - SemiBin fails for some samples; possibly because of megahit bad assemblies;
  - checkM not working for SPAdes assemblies and bins(WIP);
  - prereqisites updated.

**Notes**  

MetaBat2 usage as a Docker image brings its own difficulties. DO NOT forget mount the folder so Docker can access all files.

**TBA and ideas:**  
- automatic creation of benchmark plots(WIP);
- .yaml file for conda environment;
- maybe deploy this pipeline as a Docker image or as a package?
