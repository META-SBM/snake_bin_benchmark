# snake_bin_benchmark
Simple Snakemake pipeline for WGS metagenomics assembly, binning, quality assessing along with benchmarking all programs used in this pipeline.

**Prerequisites:**  
 - bwa-mem2;
 - strobealign;
 - semibin2;
 - checkM;
 - kraken2 and its pluspf DB;
 - megahit;
 - SPAdes (added 09/07/24);
 - MetaBat2 (added 09/07/24).

~~All tools (except MetaBat2) should be installed in conda fresh env. Use MetaBat2 as a Docker image.~~
**UPD 17/09/24**: Metabat2 now works through conda

**Ready:**
 - assembly with MEGAHIT;  
 - kraken2 with pluspf base;  
 - binning with SemiBin2 in mode ```human_oral```;  
 - realignment of reads to assembly with bwa-mem2 and strobealign;  
 - storbealign with AEMB mode;  
 - ~~saving all logs and benchmarks files~~ Fails a lot, needs major update.

**UPD 09/07/2024**  
  - added MetaBat2;
  - added SPAdes;
  - SemiBin fails for some samples; possibly because of megahit bad assemblies;
  - checkM not working for SPAdes assemblies and bins(WIP);
  - prereqisites updated.

**UPD 17/09/24**
  - added logging of each step (now works!);
  - made small script to create a table with sample name and read count (later it will be used for making benchmark plots);
  - minor fixes across all steps;
  - checkM works as intended (except some files with 0 bins found);

**WIP**
  - making script to automatically create and save benchmarking plots;
  - some files contains 0 bins, is it some error or just bad quality reads? Check and/or fix.

**Notes**  

~~MetaBat2 usage as a Docker image brings its own difficulties. DO NOT forget mount the folder so Docker can access all files.~~ Metabat2 now works through conda.

**TBA and ideas:**  
- automatic creation of benchmark plots (WIP);
- .yaml file for conda environment (almost ready; commit ASAP);
- ~~- maybe deploy this pipeline as a Docker image or as a package?~~ pip package would be more effective; but it's *on hold*
