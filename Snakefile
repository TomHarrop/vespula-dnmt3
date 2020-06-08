#!/usr/bin/env python3

from pathlib import Path
import itertools


def pairwise_combinations(x, y):
    '''
    Has to be output that snakemake likes. Each iteration of combo should
    yield something like: 
    (('ref', 'ma_FR_norm_k71_diplo1'), ('query', 'ma_MA_norm_k71_diplo0'))
    '''
    for combo in itertools.combinations(x, 2):
        yield(((combo[0]), (y[0][0], combo[1][1])))


def resolve_input(wildcards):
    my_ref = spec_to_fna[wildcards.ref]
    my_query = spec_to_fna[wildcards.query]
    return ({'ref': my_ref, 'query': my_query})


def resolve_path(x):
    return Path(x).resolve().as_posix()


bbmap = 'shub://TomHarrop/seq-utils:bbmap_38.76'
bioconductor = 'shub://TomHarrop/r-containers:bioconductor_3.11'
minimap2 = 'shub://TomHarrop/singularity-containers:minimap2_2.17r941'
samtools = 'shub://TomHarrop/align-utils:samtools_1.10'

spec_to_fna = {
    'Pcad': 'data/GCF_001313835.1_ASM131383v1_genomic.fna',
    'Pdomi': 'data/Pdomi.assembly.fna',
    'Pdor': 'data/GCA_010416905.1_CU_Pdor_10_genomic.fna',
    'Pfus': 'data/GCA_010416935.1_CU_Pfus_HIC_genomic.fna',
    'Pmet': 'data/GCA_010416925.1_CU_Pmet_PB_genomic.fna',
    'Vgerm': 'data/Vgerm.assembly.fna',
    'Vpens': 'data/Vpens.assembly.fna',
    'Vvulg': 'data/Vvulg.assembly.fna',
    # 'Vvulg': 'output/00_chrs/Chr11.fa',
    }

all_spec = sorted(set([
    x for x in spec_to_fna.keys()
    if x != "Vvulg"]))

# And Dnmt3 is Vvulg11g01820.t1
# >Pcad|NW_014569547.1:1110528-1117352
# >Pdomi|NW_015149055.1:96797-102031
# >Pdor|QUOG01000006.1:5707517-5716197
# >Pfus|QUOI01000010.1:1610914-1619584
# >Pmet|QUOH01000011.1:3318883-3323970

rule target:
    input:
        'output/03_gviz/gviz_bamfiles.pdf'


rule gviz_bamfiles:
    input:
        bamfiles = expand('output/02_minimap/{ref}.{query}/align.bam',
                          ref=['Vvulg'],
                          query=all_spec),
        fa = 'data/Vvulg.assembly.fna',
        gff = 'data/Vvulg_final_sorted.gff3'
    output:
        txdb = 'output/03_gviz/txdb.sqlite',
        plot = 'output/03_gviz/gviz_bamfiles.pdf'
    params:
        width = 483,        # point
        height = 483 // 2   # point
    log:
        'output/logs/gviz_bamfiles.log'
    singularity:
        bioconductor
    script:
        'src/gviz_bamfiles.R'

rule minimap2:
    input:
        unpack(resolve_input)
    output:
        temp('output/02_minimap/{ref}.{query}/align.sam')
    log:
        'output/logs/minimap2.{ref}.{query}.log'
    threads:
        workflow.cores
    singularity:
        minimap2
    shell:
        'minimap2 '
        '-a '
        # '-x asm20 '
        '-t {threads} '
        '{input.ref} '
        # '<( cat {input.vesp} {input.pol} ) '
        '{input.query} '
        '> {output} '
        '2> {log}'

rule extract_chr:
    input:
        'data/Vvulg.assembly.fna'
    output:
        'output/00_chrs/{chr}.fa'
    singularity:
        bbmap
    shell:
        'filterbyname.sh '
        'in={input} '
        'names={wildcards.chr} '
        'substring=name '
        'include=t '
        'out={output} '

rule sam_to_bam:
    input:
        Path('{folder}', '{samfile}.sam')
    output:
        bam = Path('{folder}', '{samfile}.bam'),
        bai = Path('{folder}', '{samfile}.bam.bai')
    singularity:
        samtools
    shell:
        'samtools view -S -b {input} '
        '| '
        'samtools sort - '
        '> {output.bam} ; '
        'samtools index {output.bam}'
