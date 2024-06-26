
# Snakemake rules below
# suitable for snm3C-seq

# example 1: run on local HPC:

# snakemake -d /gale/ddn/bican/Wubin/run_pipeline/mapping/AMB_220510_8wk_12D_13B_2_P3-5-A11 \
# --snakefile /gale/ddn/bican/Wubin/run_pipeline/mapping/AMB_220510_8wk_12D_13B_2_P3-5-A11/Snakefile \
# -j 10 --default-resources mem_mb=100 --resources mem_mb=50000

# example 2: run on GCP
# sync Snakefile to GCP:~/sky_workdir

# conda activate yap
# prefix="DATASET/mapping/Pool1/AMB_220510_8wk_12D_13B_2_P3-5-A11/"
# snakemake --snakefile ~/sky_workdir/Snakefile -j 10 --default-resources mem_mb=100 --resources mem_mb=50000 
# --config gcp=True --default-remote-prefix ${prefix} \
# --default-remote-provider GS --google-lifesciences-region us-west1 -np

if "gcp" in config:
    gcp=config["gcp"] # if the fastq files stored in GCP cloud, set gcp=True in snakemake: --config gcp=True
else:
    gcp=False

if gcp:
    from snakemake.remote.GS import RemoteProvider as GSRemoteProvider
    GS = GSRemoteProvider()
    os.environ['GOOGLE_APPLICATION_CREDENTIALS'] =os.path.expanduser('~/.config/gcloud/application_default_credentials.json')
    bam_dir=workflow.default_remote_prefix+"/bam"
    allc_dir=workflow.default_remote_prefix+"/allc"
    hic_dir=workflow.default_remote_prefix+"/hic"
else:
    bam_dir="bam"
    allc_dir="allc"
    hic_dir="hic"

fastq_dir=os.path.abspath("fastq")

# the summary rule is the final target
rule summary:
    input:
        expand("allc/{cell_id}.allc.tsv.gz", cell_id=CELL_IDS),
        # also add all the stats path here, so they won't be deleted until summary is generated
        expand("allc/{cell_id}.allc.tsv.gz.count.csv", cell_id=CELL_IDS),
        local(expand("fastq/{cell_id}-{read_type}.trimmed.stats.tsv", 
                        cell_id=CELL_IDS,read_type=['R1','R2'])),
        local(expand(bam_dir+"/{cell_id}-{read_type}.two_mapping.deduped.matrix.txt", 
                        cell_id=CELL_IDS,read_type=['R1','R2'])),
        local(expand(bam_dir+"/{cell_id}-{read_type}.two_mapping.filter.bam", 
                        cell_id=CELL_IDS,read_type=['R1','R2'])),
        local(expand(bam_dir+"/{cell_id}-{read_type}.two_mapping.deduped.bam", 
                        cell_id=CELL_IDS,read_type=['R1','R2'])),
        expand("hic/{cell_id}.3C.contact.tsv.gz", cell_id=CELL_IDS),
        expand("hic/{cell_id}.3C.contact.tsv.counts.txt", cell_id=CELL_IDS)

        #expand("allc/{cell_id}.allc.tsv.gz", cell_id=CELL_IDS),
        # also add all the stats path here, so they won't be deleted until summary is generated
        #expand("allc/{cell_id}.allc.tsv.gz.count.csv", cell_id=CELL_IDS),
        #local(expand("fastq/{cell_id}-{read_type}.trimmed.stats.tsv", 
        #                cell_id=CELL_IDS,read_type=['R1','R2'])),
        #local(expand(bam_dir+"/{cell_id}-{read_type}.m3C.deduped.bam", 
        #                cell_id=CELL_IDS,read_type=['R1','R2'])),
        #expand(bam_dir+"/{cell_id}-{read_type}.m3C.sorted.bam", 
        #                cell_id=CELL_IDS,read_type=['R1','R2']),
        #expand("hic/{cell_id}.3C.contact.tsv.gz", cell_id=CELL_IDS),
        #expand("hic/{cell_id}.3C.contact.tsv.counts.txt", cell_id=CELL_IDS)
    output:
        "MappingSummary.csv.gz"
    params:
        outdir=os.path.abspath("./") if not gcp else workflow.default_remote_prefix,
    shell:
        """
        yap-internal summary --output_dir {params.outdir} --fastq_dir {fastq_dir} --mode {mode} --barcode_version {barcode_version} \
                    --mc_stat_feature "{mc_stat_feature}" --mc_stat_alias "{mc_stat_alias}" \
                    --num_upstr_bases {num_upstr_bases}
        """

# Trim reads
rule trim_r1:
    input:
        fq=local("fastq/{cell_id}-R1.fq.gz") #if not gcp else GS.remote("gs://"+workflow.default_remote_prefix+"/fastq/{cell_id}-R1.fq.gz")
    output:
        fq=local(temp("fastq/{cell_id}-R1.trimmed.fq.gz")),
        stats=local(temp("fastq/{cell_id}-R1.trimmed.stats.tsv"))
    threads:
        2
    shell:
        "cutadapt --report=minimal -a {r1_adapter} {input.fq} 2> {output.stats} | "
        "cutadapt --report=minimal -O 6 -q 20 -u {r1_left_cut} -u -{r1_right_cut} -m 30 "
        "-o {output.fq} - >> {output.stats}"

rule trim_r2:
    input:
        fq=local("fastq/{cell_id}-R2.fq.gz") #if not gcp else GS.remote("gs://"+workflow.default_remote_prefix+"/fastq/{cell_id}-R2.fq.gz")
    output:
        fq=local(temp("fastq/{cell_id}-R2.trimmed.fq.gz")),
        stats=local(temp("fastq/{cell_id}-R2.trimmed.stats.tsv"))
    threads:
        2
    shell:
        "cutadapt --report=minimal -a {r2_adapter} {input.fq} 2> {output.stats} | "
        "cutadapt --report=minimal -O 6 -q 20 -u {r2_left_cut} -u -{r2_right_cut} -m 30 "
        "-o {output.fq} - >> {output.stats}"

# bismark mapping
rule bismark:
    input:
        local("fastq/{cell_id}-{read_type}.trimmed.fq.gz")
    output:
        bam=local(temp(bam_dir+"/{cell_id}-{read_type}.trimmed_bismark_bt2.bam")),
        um=local(temp(bam_dir+"/{cell_id}-{read_type}.trimmed.fq.gz_unmapped_reads.fq.gz")),
        stats=local(temp(bam_dir+"/{cell_id}-{read_type}.trimmed_bismark_bt2_SE_report.txt"))
    params:
        mode=lambda wildcards: "--pbat" if wildcards.read_type=="R1" else ""
    threads:
        3
    resources:
        mem_mb=14000
    shell:
        # map R1 with --pbat mode; map R2 with normal SE mode
        """
        mkdir -p {bam_dir}
        bismark {bismark_reference} -un --bowtie2 {input} {params.mode} -o {bam_dir} --temp_dir {bam_dir}
        """

# split unmapped fastq
rule split_um_fastq:
    input:
        local(bam_dir+"/{cell_id}-{read_type}.trimmed.fq.gz_unmapped_reads.fq.gz")
    output:
        local(temp(bam_dir+"/{cell_id}-{read_type}.trimmed.fq.gz_unmapped_reads.split.fq.gz"))
    threads:
        1
    shell:
        "yap-internal m3c-split-reads --fastq_path {input} --output_path {output} "
        "--size_l {split_left_size} --size_r {split_right_size} "
        "--size_m {split_middle_min_size} --trim_b {trim_on_both_end}"

# map split fastq again
rule bismark_split:
    input:
        local(bam_dir+"/{cell_id}-{read_type}.trimmed.fq.gz_unmapped_reads.split.fq.gz")
    output:
        bam=local(temp(bam_dir+"/{cell_id}-{read_type}.trimmed.fq.gz_unmapped_reads.split_bismark_bt2.bam")),
        stats=local(temp(bam_dir+"/{cell_id}-{read_type}.trimmed.fq.gz_unmapped_reads.split_bismark_bt2_SE_report.txt"))
    params:
        mode=lambda wildcards: "--pbat" if wildcards.read_type=="R1" else ""
    threads:
        3
    resources:
        mem_mb=14000
    shell:
        """
        bismark {bismark_reference} --bowtie2 {input} {params.mode} -o {bam_dir} --temp_dir {bam_dir}
        """

# merge two bam files
rule merge_raw_bam:
    input:
        local(bam_dir+"/{cell_id}-{read_type}.trimmed_bismark_bt2.bam"),
        local(bam_dir+"/{cell_id}-{read_type}.trimmed.fq.gz_unmapped_reads.split_bismark_bt2.bam")
    output:
        local(temp("{indir}/{cell_id}-{read_type}.two_mapping.bam"))
    shell:
        "samtools merge -f {output} {input}"

# filter bam
rule filter_bam:
    input:
        local("{sname}.bam")
    output:
        bam=local(temp("{sname}.filter.bam")),
    shell:
        """
        samtools view -q 10 -b -h -o {output.bam} {input}
        """

# merge R1 and R2, get final bam for mC calling
rule merge_mc_bam:
    input:
        local(bam_dir+"/{cell_id}-R1.two_mapping.bam"),
        local(bam_dir+"/{cell_id}-R2.two_mapping.bam")
    output:
        bam=local(temp(bam_dir+"/{cell_id}.m3C.bam")),
    shell:
        "samtools merge -f {output.bam} {input} && samtools sort -n -o {output.bam} {output.bam}"


# remove PCR duplicates
rule dedup_bam:
    input:
        local(bam_dir+"/{cell_id}.m3C.bam")
    output:
        bam=local(temp(bam_dir+"/{cell_id}.m3C.dedup.bam")),
    params:
        tmp_dir=os.path.abspath("bam/temp") if not gcp else workflow.default_remote_prefix+"/bam/temp"
    resources:
        mem_mb=1000
    shell:
        """
         yap-internal mark-duplicates --bam_path {input} --output_path {output}
        """


rule sort_bam:
    input:
        local(bam_dir+"/{cell_id}.m3C.dedup.bam")
    output:
        bam=bam_dir+"/{cell_id}.m3C.sorted.bam",
        bai=bam_dir+"/{cell_id}.m3C.sorted.bam.bai",
    shell:
        """
        samtools sort -@ 20 -o {ouptut.bam} {input} && samtools index -@ 20 {output.bam}
        """


# generate ALLC
rule allc:
    input:
        bam=bam_dir+"/{cell_id}.m3C.sorted.bam",
        index=bam_dir+"/{cell_id}.m3C.sorted.bam.bai"
    output:
        allc="{indir}/{cell_id}.allc.tsv.gz",
        tbi="{indir}/{cell_id}.allc.tsv.gz.tbi",
        stats="{indir}/{cell_id}.allc.tsv.gz.count.csv"
    threads:
        2
    resources:
        mem_mb=500
    shell:
        """
        mkdir -p {allc_dir}
        allcools bam-to-allc \
                --bam_path {input.bam} \
                --reference_fasta {reference_fasta} \
                --output_path {output.allc} \
                --cpu 1 \
                --num_upstr_bases {num_upstr_bases} \
                --num_downstr_bases {num_downstr_bases} \
                --compress_level {compress_level} \
                --save_count_df
        """

rule generate_contact:
    input:
       local(bam_dir+"/{cell_id}.m3C.dedup.bam"),
    output:
        contact="{indir}/{cell_id}.3C.contact.tsv.gz",
        stats="{indir}/{cell_id}.3C.contact.tsv.counts.txt"
    resources:
        mem_mb=300
    shell:
        """
        mkdir -p {hic_dir}
        yap-internal generate-contacts --bam_path {input} --output_path {output.contact} \
                    --chrom_size_path {chrom_size_path} --min_gap {min_gap}
        """
