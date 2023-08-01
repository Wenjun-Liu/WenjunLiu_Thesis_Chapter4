rule count:
    input:
        bams = expand(["data/aligned/bam/{sample}{tag}/Aligned.sortedByCoord.out.bam"],
                      tag = [tag], sample = samples['sample']),
        gtf = rules.get_annotation.output
    output:
        temp(counts_file)
    conda:
        "../envs/subread.yml"
    threads: 4
    params:
        minOverlap = config['featureCounts']['minOverlap'],
        fracOverlap = config['featureCounts']['fracOverlap'],
        q = config['featureCounts']['minQual'],
        s = config['featureCounts']['strandedness'],
        extra = config['featureCounts']['extra']
    shell:
       """
       featureCounts \
         {params.extra} \
         -Q {params.q} \
         -s {params.s} \
         --minOverlap {params.minOverlap} \
         --fracOverlap {params.fracOverlap} \
         -T {threads} \
         -a {input.gtf} \
         -o {output} \
         {input.bams}
       """

rule zip_counts:
    input: counts_file
    output: counts_file + ".gz"
    threads: 1
    shell: "gzip -c {input} > {output}"
