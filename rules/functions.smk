def get_fastq(wildcards):
    return samples.loc[(wildcards.sample, wildcards.unit), ["fq1", "fq2"]].dropna()

def is_single_end(sample, unit):
    return pd.isnull(samples.loc[(sample, unit), "fq2"])

def get_fq(wildcards):
    if config["trimming"]["skip"]:
        # no trimming, use raw reads
        return {"fq1":samples.loc[(wildcards.sample, wildcards.unit), ["fq1"]].dropna(),
                "fq2":samples.loc[(wildcards.sample, wildcards.unit), ["fq2"]].dropna()}
    else:
        # yes trimming, use trimmed data
        if not is_single_end(**wildcards):
            # paired-end sample
            return {"fq1": expand("trimmed/{sample}-{unit}.{group}.fastq.gz",
                          group=1, **wildcards),
                    "fq2": expand("trimmed/{sample}-{unit}.{group}.fastq.gz",
                          group=2, **wildcards)}
        # single end sample
        return {"fq1":"trimmed/{sample}-{unit}.fastq.gz".format(**wildcards)}


def get_fileend(wildcards):
    if config["trimming"]["skip"]:
        fq1=samples.loc[(wildcards.sample, wildcards.unit), ["fq1"]].dropna()
    else:
        if not is_single_end(**wildcards):
            fq1=expand("trimmed/{sample}-{unit}.{group}.fastq.gz",
                          group=1, **wildcards)
        else:
            fq1=expand("trimmed/{sample}-{unit}.fastq.gz", **wildcards),
    if fq1[0].endswith(".gz"):
        readcmd = "--readFilesCommand zcat"
    else:
        readcmd = ""
    return readcmd

def get_strandness(samples):
    if "strandedness" in samples.columns:
        return samples["strandedness"].tolist()
    else:
        strand_list=["none"]
        return strand_list*samples.shape[0]

def get_deseq2_threads(wildcards=None):
    # https://twitter.com/mikelove/status/918770188568363008
    few_coeffs = False if wildcards is None else len(get_contrast(wildcards)) < 10
    return 1 if len(samples) < 100 or few_coeffs else 6

def get_contrast(wildcards):
    return config["diffexp"]["contrasts"][wildcards.contrast]


