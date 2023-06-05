configfile: "config/config.yaml"
configfile: "config/samples.yaml"

rule all:
    input: 
        expand("results/bsmooth_fit/{samples}/{samples}.fitted.rds", samples=config["samples"])

rule bsmooth_fit:
    input:
        filePath=lambda wildcards: config["samples"][wildcards.samples]
    output:
        bsfit="results/{rule}/{samples}/{samples}.fitted.rds"
    params:
        rscript=config["BSmooth_fit"]
    log:
        "logs/{rule}/{samples}/{samples}.log"
    shell:
        "Rscript {params.rscript} --sample {wildcards.samples} --input_file {input.filePath} --output_file {output.bsfit} 2> {log}"

#this runs ok
#snakemake -np #test run with
#snakemake -s bsmooth_Snakefile.smk --cores 12 --forcerun -np #deploy code
#nohup snakemake -s bsmooth_Snakefile.smk --latency-wait 60 --restart-times 2 --keep-going --forceall --cluster "bsub -J {rule} -R "rusage[mem=32]" -W 5:00 -n 12 -o logs/cluster/{rule}.%J.out -e logs/cluster/{rule}.%J.err" -j 3 &
