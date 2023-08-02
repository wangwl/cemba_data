
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

# the summary rule is the final target
rule summary:
    input:
        expand("allc/{cell_id}.allc.tsv.gz", cell_id=CELL_IDS),
        # also add all the stats path here, so they won't be deleted until summary is generated
        expand("allc/{cell_id}.allc.tsv.gz.count.csv", cell_id=CELL_IDS),
        expand("fastq/{cell_id}-R1.trimmed.stats.tsv", cell_id=CELL_IDS),
        expand("fastq/{cell_id}-R2.trimmed.stats.tsv", cell_id=CELL_IDS),
        expand("bam/{cell_id}-R1.two_mapping.deduped.matrix.txt", cell_id=CELL_IDS),
        expand("bam/{cell_id}-R2.two_mapping.deduped.matrix.txt", cell_id=CELL_IDS),
        expand("bam/{cell_id}-R1.two_mapping.filter.bam", cell_id=CELL_IDS),
        expand("bam/{cell_id}-R2.two_mapping.filter.bam", cell_id=CELL_IDS),
        expand("bam/{cell_id}-R1.two_mapping.deduped.bam", cell_id=CELL_IDS),
        expand("bam/{cell_id}-R2.two_mapping.deduped.bam", cell_id=CELL_IDS),
        expand("hic/{cell_id}.3C.contact.tsv.gz", cell_id=CELL_IDS),
        expand("hic/{cell_id}.3C.contact.tsv.counts.txt", cell_id=CELL_IDS)
    output:
        "MappingSummary.csv.gz"
    params:
        outdir=os.path.abspath("./") if not gcp else workflow.default_remote_prefix
    shell:
        """
        yap-internal summary --output_dir {params.outdir} --mode {mode} --barcode_version {barcode_version} \
                    --mc_stat_feature {mc_stat_feature} --mc_stat_alias {mc_stat_alias} \
                    --num_upstr_bases {num_upstr_bases}
        """

# Trim reads
rule trim_r1:
    input:
        fq=local("fastq/{cell_id}-R1.fq.gz") #if not gcp else GS.remote("gs://"+workflow.default_remote_prefix+"/fastq/{cell_id}-R1.fq.gz")
    output:
        fq=temp("fastq/{cell_id}-R1.trimmed.fq.gz"),
        stats=temp("fastq/{cell_id}-R1.trimmed.stats.tsv")
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
        fq=temp("fastq/{cell_id}-R2.trimmed.fq.gz"),
        stats=temp("fastq/{cell_id}-R2.trimmed.stats.tsv")
    threads:
        2
    shell:
        "cutadapt --report=minimal -a {r2_adapter} {input.fq} 2> {output.stats} | "
        "cutadapt --report=minimal -O 6 -q 20 -u {r2_left_cut} -u -{r2_right_cut} -m 30 "
        "-o {output.fq} - >> {output.stats}"

# bismark mapping, R1 and R2 separately
rule bismark_r1:
    input:
        "fastq/{cell_id}-R1.trimmed.fq.gz"
    output:
        bam=temp("bam/{cell_id}-R1.trimmed_bismark.bam"),
        um=temp("bam/{cell_id}-R1.trimmed.fq.gz_unmapped_reads.fq.gz"),
        stats=temp("bam/{cell_id}-R1.trimmed_bismark_SE_report.txt")
    params:
        bam_dir=os.path.abspath("bam") if not gcp else workflow.default_remote_prefix+"/bam"
    threads:
        3
    resources:
        mem_mb=14000
    shell:
        # map R1 with --pbat mode
        "bismark {bismark_reference} -un --bowtie1 {input} "
        "--pbat -o {params.bam_dir} --temp_dir {params.bam_dir}"

rule bismark_r2:
    input:
        "fastq/{cell_id}-R2.trimmed.fq.gz"
    output:
        bam=temp("bam/{cell_id}-R2.trimmed_bismark.bam"),
        um=temp("bam/{cell_id}-R2.trimmed.fq.gz_unmapped_reads.fq.gz"),
        stats=temp("bam/{cell_id}-R2.trimmed_bismark_SE_report.txt")
    params:
        bam_dir=os.path.abspath("bam") if not gcp else workflow.default_remote_prefix+"/bam"
    threads:
        3
    resources:
        mem_mb=14000
    shell:
        # map R2 with normal SE mode
        "bismark {bismark_reference} -un --bowtie1 {input} "
        "-o {params.bam_dir} --temp_dir {params.bam_dir}"


# split unmapped fastq
rule split_um_fastq_r1:
    input:
        "bam/{cell_id}-R1.trimmed.fq.gz_unmapped_reads.fq.gz"
    output:
        temp("bam/{cell_id}-R1.trimmed.fq.gz_unmapped_reads.split.fq.gz")
    threads:
        1
    shell:
        "yap-internal m3c-split-reads --fastq_path {input} --output_path {output} "
        "--size_l {split_left_size} --size_r {split_right_size} "
        "--size_m {split_middle_min_size} --trim_b {trim_on_both_end}"

rule split_um_fastq_r2:
    input:
        "bam/{cell_id}-R2.trimmed.fq.gz_unmapped_reads.fq.gz"
    output:
        temp("bam/{cell_id}-R2.trimmed.fq.gz_unmapped_reads.split.fq.gz")
    threads:
        1
    shell:
        "yap-internal m3c-split-reads --fastq_path {input} --output_path {output} "
        "--size_l {split_left_size} --size_r {split_right_size} "
        "--size_m {split_middle_min_size} --trim_b {trim_on_both_end}"

# map split fastq again
rule bismark_split_r1:
    input:
        "bam/{cell_id}-R1.trimmed.fq.gz_unmapped_reads.split.fq.gz"
    output:
        bam=temp("bam/{cell_id}-R1.trimmed.fq.gz_unmapped_reads.split_bismark.bam"),
        stats=temp("bam/{cell_id}-R1.trimmed.fq.gz_unmapped_reads.split_bismark_SE_report.txt")
    params:
        bam_dir=os.path.abspath("bam") if not gcp else workflow.default_remote_prefix+"/bam"
    threads:
        3
    resources:
        mem_mb=14000
    shell:
        # map R1 with --pbat mode
        "bismark {bismark_reference} --bowtie1 {input} "
        "--pbat -o {params.bam_dir} --temp_dir {params.bam_dir}"

rule bismark_split_r2:
    input:
        "bam/{cell_id}-R2.trimmed.fq.gz_unmapped_reads.split.fq.gz"
    output:
        bam=temp("bam/{cell_id}-R2.trimmed.fq.gz_unmapped_reads.split_bismark.bam"),
        stats=temp("bam/{cell_id}-R2.trimmed.fq.gz_unmapped_reads.split_bismark_SE_report.txt")
    params:
        bam_dir=os.path.abspath("bam") if not gcp else workflow.default_remote_prefix+"/bam"
    threads:
        3
    resources:
        mem_mb=14000
    shell:
        # map R2 with normal SE mode
        "bismark {bismark_reference} --bowtie1 {input} "
        "-o {params.bam_dir} --temp_dir {params.bam_dir}"

# merge two bam files
rule merge_r1_raw_bam:
    input:
        "bam/{cell_id}-R1.trimmed_bismark.bam",
        "bam/{cell_id}-R1.trimmed.fq.gz_unmapped_reads.split_bismark.bam"
    output:
        temp("bam/{cell_id}-R1.two_mapping.bam")
    shell:
        "samtools merge -f {output} {input}"

rule merge_r2_raw_bam:
    input:
        "bam/{cell_id}-R2.trimmed_bismark.bam",
        "bam/{cell_id}-R2.trimmed.fq.gz_unmapped_reads.split_bismark.bam"
    output:
        temp("bam/{cell_id}-R2.two_mapping.bam")
    shell:
        "samtools merge -f {output} {input}"


# filter bam
rule filter_r1_bam:
    input:
        "bam/{cell_id}-R1.two_mapping.bam"
    output:
        temp("bam/{cell_id}-R1.two_mapping.filter.bam")
    shell:
        "samtools view -b -h -q 10 -o {output} {input}"

rule filter_r2_bam:
    input:
        "bam/{cell_id}-R2.two_mapping.bam"
    output:
        temp("bam/{cell_id}-R2.two_mapping.filter.bam")
    shell:
        "samtools view -b -h -q 10 -o {output} {input}"

# sort bam by coords
rule sort_r1_bam:
    input:
        "bam/{cell_id}-R1.two_mapping.filter.bam"
    output:
        temp("bam/{cell_id}-R1.two_mapping.sorted.bam")
    resources:
        mem_mb=1000
    shell:
        "samtools sort -o {output} {input}"

rule sort_r2_bam:
    input:
        "bam/{cell_id}-R2.two_mapping.filter.bam"
    output:
        temp("bam/{cell_id}-R2.two_mapping.sorted.bam")
    resources:
        mem_mb=1000
    shell:
        "samtools sort -o {output} {input}"

# remove PCR duplicates
rule dedup_r1_bam:
    input:
        "bam/{cell_id}-R1.two_mapping.sorted.bam"
    output:
        bam=temp("bam/{cell_id}-R1.two_mapping.deduped.bam"),
        stats=temp("bam/{cell_id}-R1.two_mapping.deduped.matrix.txt")
    params:
        tmp_dir=os.path.abspath("bam/temp") if not gcp else workflow.default_remote_prefix+"/bam/temp"
    resources:
        mem_mb=1000
    shell:
        "picard MarkDuplicates I={input} O={output.bam} M={output.stats} "
        "REMOVE_DUPLICATES=true TMP_DIR={params.tmp_dir}"

rule dedup_r2_bam:
    input:
        "bam/{cell_id}-R2.two_mapping.sorted.bam"
    output:
        bam=temp("bam/{cell_id}-R2.two_mapping.deduped.bam"),
        stats=temp("bam/{cell_id}-R2.two_mapping.deduped.matrix.txt")
    params:
        tmp_dir=os.path.abspath("bam/temp") if not gcp else workflow.default_remote_prefix+"/bam/temp"
    resources:
        mem_mb=1000
    shell:
        "picard MarkDuplicates I={input} O={output.bam} M={output.stats} "
        "REMOVE_DUPLICATES=true TMP_DIR={params.tmp_dir}"

# merge R1 and R2, get final bam for mC calling
rule merge_mc_bam:
    input:
        "bam/{cell_id}-R1.two_mapping.deduped.bam",
        "bam/{cell_id}-R2.two_mapping.deduped.bam"
    output:
        bam=temp("bam/{cell_id}.mC.bam"),
        bai=temp("bam/{cell_id}.mC.bam.bai")
    shell:
        "samtools merge -f {output.bam} {input} && samtools index {output.bam}"

# generate ALLC
rule allc:
    input:
        bam="bam/{cell_id}.mC.bam",
        index="bam/{cell_id}.mC.bam.bai"
    output:
        allc="allc/{cell_id}.allc.tsv.gz",
        stats=temp("allc/{cell_id}.allc.tsv.gz.count.csv")
    threads:
        2
    resources:
        mem_mb=500
    shell:
        'allcools bam-to-allc '
        '--bam_path {input.bam} '
        '--reference_fasta {reference_fasta} '
        '--output_path {output.allc} '
        '--cpu 1 '
        '--num_upstr_bases {num_upstr_bases} '
        '--num_downstr_bases {num_downstr_bases} '
        '--compress_level {compress_level} '
        '--save_count_df'


# merge and sort (by read name) bam before dedup for generating contact
# contact dedup happen within generate contact
rule merge_3c_bam_for_contact:
    input:
        "bam/{cell_id}-R1.two_mapping.sorted.bam",
        "bam/{cell_id}-R2.two_mapping.sorted.bam"
    output:
        temp("bam/{cell_id}.3C.bam")
    shell:
        "samtools merge -f {output} {input}"

rule sort_bam_for_contact:
    input:
        "bam/{cell_id}.3C.bam"
    output:
        "bam/{cell_id}.3C.sorted.bam"
    resources:
        mem_mb=1000
    shell:
        "samtools sort -n -o {output} {input}"

rule generate_contact:
    input:
        "bam/{cell_id}.3C.sorted.bam"
    output:
        contact="hic/{cell_id}.3C.contact.tsv.gz",
        stats=temp("hic/{cell_id}.3C.contact.tsv.counts.txt")
    resources:
        mem_mb=300
    shell:
        "yap-internal generate-contacts --bam_path {input} --output_path {output.contact} "
        "--chrom_size_path {chrom_size_path} --min_gap {min_gap}"
