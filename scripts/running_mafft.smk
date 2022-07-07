IDS, = glob_wildcards("top_clusters/{id}.fasta")

rule all:
    input:  expand("top_clusters_msa/{id}.fasta", id=IDS)

rule:
    input: "top_clusters/{id}.fasta"
    output: "top_clusters_msa/{id}.fasta"
    threads: 10
    shell:
        "/usr/bin/./mafft --anysymbol --thread {threads} --auto {input} > {output}"

