rule cutadapt_pe:
    input:
        get_fastq
    output:
        fastq1="trimmed/{sample}-{unit}.1.fastq.gz",
        fastq2="trimmed/{sample}-{unit}.2.fastq.gz",
        qc="trimmed/{sample}-{unit}.qc.txt"
    params:
        adapters="-a {}".format(config["trimming"]["adapter"]),
        others=config["params"]["cutadapt-pe"]
    conda:
        "../envs/trim.yaml"
    resources:
        mem_mb = 10000
    threads: 1
    log:
        "logs/cutadapt/{sample}-{unit}.log"
    shell:
        """
        cutadapt \
            {params.adapters} \
            {params.others} \
            -o {output.fastq1} \
            -p {output.fastq2} \
            -j {threads} \
            {input} \
        > {output.qc}
        """

rule cutadapt:
    input:
        get_fastq
    output:
        fastq="trimmed/{sample}-{unit}.fastq.gz",
        qc="trimmed/{sample}-{unit}.qc.txt"
    params:
        adapters="-a {}".format(config["trimming"]["adapter"]),
        others=config["params"]["cutadapt-se"]
    conda:
        "../envs/trim.yaml"
    resources:
        mem_mb = 10000
    threads: 1
    log:
        "logs/cutadapt/{sample}-{unit}.log"
    shell:
        """
        cutadapt \
            {params.adapters} \
            {params.others} \
            -o {output.fastq} \
            -j {threads} \
            {input} \
        > {output.qc}
        """
