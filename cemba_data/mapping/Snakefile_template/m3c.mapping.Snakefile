
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

# the sort_bam rule is the final target
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


# merge and sort (by read name) bam before dedup for generating contact
# contact dedup happen within generate contact
rule merge_3c_bam_for_contact:
    input:
        local(bam_dir+"/{cell_id}-R1.two_mapping.sorted.bam"),
        local(bam_dir+"/{cell_id}-R2.two_mapping.sorted.bam")
    output:
        local(temp(bam_dir+"/{cell_id}.3C.bam"))
    shell:
        "samtools merge -f {output} {input}"

