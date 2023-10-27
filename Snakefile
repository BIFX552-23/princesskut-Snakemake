# count number of reads
rule countreads: 
  output: "{indir}.{myfile}.fq.count"
  input: "{indir}/{myfile}.fq"
  shell: 
    "echo $(( $(wc -l <{input}) / 4)) > {output}"

#trim reads
rule trimreads: 
  output:"trimmed/{myfile}.fq"
  input:"reads/{myfile}.fq"
  shell:
    "fastq_quality_trimmer -t 20 -l 100 -o {output} < {input}"

# align reads
rule kallisto_quant:
  output:
    h5   = "kallisto.{sample}/abundance.h5",
    tsv  = "kallisto.{sample}/abundance.tsv",
    json = "kallisto.{sample}/run_info.json"
  input:
    idx = "transcriptome/Saccharomyces_cerevisiae.R64-1-1.kallisto_index", 
    fq1 = "trimmed/{sample}_1.fq",
    fq2 = "trimmed/{sample}_2.fq"
  shell:
    "kallisto quant -i {index.idx} -o kallisto.{wildcards.sample} {input.fq1} {input.fq2}"

# kallisto index 
rule kallisto_index: 
  output:
    idx = "transcriptome/{strain}.kallisto_index",
    log = "transcriptome/{strain}.kallisto_log"
  input: "transcriptome/{strain}.cdna.all.fa.gz"
  shell:
    "kallisto index -i {output.idx} {input} >& {output.log}"
  

  