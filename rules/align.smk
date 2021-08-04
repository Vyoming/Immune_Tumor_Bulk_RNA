rule align:
    input:
        unpack(get_fq)
    output:
        "star/{sample}-{unit}/Aligned.out.bam",
        "star/{sample}-{unit}/ReadsPerGene.out.tab"
    log:
        "logs/star/{sample}-{unit}.log"
    params:
        readcmd=get_fileend,
        # path to STAR reference genome index
        index=config["ref"]["index"],
        # optional parameters
        extra="--quantMode GeneCounts --sjdbGTFfile {} {}".format(
              config["ref"]["annotation"], config["params"]["star"])
    conda:
        "../envs/align.yaml" 
    resources:
        mem_mb = 120000
    threads: 4
    shell:
        """
        STAR \
        {params.extra} \
            --runThreadN {threads} \
            --runMode alignReads \
            --genomeDir {params.index} \
            --readFilesIn {input.fq1} {input.fq2} {params.readcmd} \
            --outReadsUnmapped Fastq \
            --outSAMtype BAM Unsorted \
            --outFileNamePrefix star/{wildcards.sample}-{wildcards.unit}/
        """

rule sort:
    input:
        "star/{sample}-{unit}/Aligned.out.bam"
    output:
        "star/{sample}-{unit}/Aligned.sortedByCoord.out.bam"
    log:
        "logs/star/{sample}-{unit}-sam-sort.log"
    conda:
        "../envs/align.yaml"
    resources:
        mem_mb = 120000
    threads: 4
    shell:
        """
        samtools sort {input} -o {output} -@ {threads} -m 20G
        """

rule index:
    input:
        "star/{sample}-{unit}/Aligned.sortedByCoord.out.bam"
    output:
        "star/{sample}-{unit}/Aligned.sortedByCoord.out.bam.bai"
    log:
        "logs/star/{sample}-{unit}-sam-index.log"
    conda:
        "../envs/align.yaml"
    resources:
        mem_mb = 10000
    threads: 4
    shell:
        """
        samtools index {input} -@ {threads}
        """         
