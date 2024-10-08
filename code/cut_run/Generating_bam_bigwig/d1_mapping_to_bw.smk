# specify config file

configfile: "config.yaml"

# Define a global variable for the sleep duration (in seconds)
sleep_duration = 10  # Adjust as needed

rule all:
    input:
        expand("bam_shigeki/sort_{sample}.bam", sample=config["samples"]), 
        expand("bam_shigeki/sort_{sample}.bam.bai", sample=config["samples"]),
        expand("bam_shigeki/bt2_{sample}.log", sample=config["samples"]),
        expand("bigwig_shigeki/{sample}.bw", sample=config["samples"])

# Read mapping and converting SAM to BAM
rule bowtie2:
    input:
        fq1=lambda wildcards: config["samples"][wildcards.sample]["fq1"],
        fq2=lambda wildcards: config["samples"][wildcards.sample]["fq2"]
    output:
        bam1 = "bam_shigeki/ori_{sample}.bam"
    params:
        index = "bt2_index/mouse38-index"  # change
    log:
        log = "bam_shigeki/bt2_{sample}.log"
    resources:
        mem_mb = 40000  # 40 GB
    shell:
        """
        module purge
        module load Bioinformatics
        module load python/3.9.12
        module load bowtie2/2.4.2-3pufpzz
        module load samtools
        bowtie2 -p 16 -q --local \
        -x {params.index} \
        -1 {input.fq1} \
        -2 {input.fq2} | samtools view -h -S -b -q 10 --threads 16 > {output.bam1} 
        samtools stats {output.bam1} > bam/{wildcards.sample}_raw_stats.txt
        sleep {sleep_duration}  # Introduce a delay to avoid overloading the server
        """

# Define the rule for filtering BAM to remove reads < 140 nt
rule filter_bam:
    input:
        bam1 = "bam_shigeki/ori_{sample}.bam"
    output:
        bam_filtered = "bam_shigeki/filtered_{sample}.bam"
    shell:
        """
        module load Bioinformatics
        module load python/3.9.12
        module load samtools
        samtools view -h {input.bam1} | 
        awk 'BEGIN {{OFS="\t"}} {{if ($1 ~ /^@/) {{print $0}} else if ($9 >= 140 || $9 <= -140) {{print $0}} }}' | 
        samtools view -Sb - > {output.bam_filtered}
        sleep {sleep_duration}  # Introduce a delay to avoid overloading the server
        """

# Define the rule for sorting BAM
rule bam_sort:
    input:
        bam_filtered = "bam_shigeki/filtered_{sample}.bam"
    output:
        bam_sorted = "bam_shigeki/sort_{sample}.bam"
    resources:
        mem_mb = 40000  # 40 GB
    shell:
        """
        module load Bioinformatics
        module load python/3.9.12
        module load samtools
        samtools sort -o {output.bam_sorted} -T {wildcards.sample}_temp --threads 16 {input.bam_filtered} -m 1G
        sleep {sleep_duration}  # Introduce a delay to avoid overloading the server
        """

# Define the rule for indexing BAM
rule bam_index:
    input:
        bam_sorted = "bam_shigeki/sort_{sample}.bam"
    output:
        bai1 = "bam_shigeki/sort_{sample}.bam.bai"
    shell:
        """
        module load Bioinformatics
        module load python/3.9.12
        module load samtools
        samtools index {input.bam_sorted}
        sleep {sleep_duration}  # Introduce a delay to avoid overloading the server
        """

# Define the rule for generating bigwig files
rule bam_to_bigwig:
    input:
        bam_sorted = "bam_shigeki/sort_{sample}.bam",
        bai1 = "bam_shigekim/sort_{sample}.bam.bai"
    output:
        bigwig = "bigwig_shigeki/{sample}.bw"
    shell:
        """
        source ~/miniconda3/etc/profile.d/conda.sh  # change if needed
        conda activate deeptools
        bamCoverage -b {input.bam_sorted} -of bigwig -o {output.bigwig} --normalizeUsing BPM -p 16
        sleep {sleep_duration}  # Introduce a delay to avoid overloading the server
        """

