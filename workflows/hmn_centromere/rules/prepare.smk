# ### Compute fai idx
# rule index_fasta:
#     input:
#         fa = os.path.join(DATADIR, "{sample}.fa.gz")
#     output:
#         fai_idx = os.path.join(WORKDIR, "data/asm/fa_idx/{sample}.fa.gz.fai"),
#         gzi_idx = os.path.join(WORKDIR, "data/asm/fa_idx/{sample}.fa.gz.gzi")
#     log:
#         "logs/annotations/asm/index_{sample}_fa.log"
#     conda:
#         "../envs/standalone/sambedtools.yaml"
#     shell:
#         r"""
#         samtools faidx {input.fa} --fai-idx {output.fai_idx} --gzi-idx {output.gzi_idx} 2> {log}
#         """

# ### Download mapping to CHM13
# rule download_mashmaps:
#     output:
#         mashmap = os.path.join(WORKDIR, "data/annotations/mashmap/{sample}.pi95.paf")
#     params:
#         hprc_url = config["hprc_url"]
#     log:
#         "logs/download/mashmap/download_{sample}_mashmap.log"
#     shell:
#         r"""
#         exec > {log} 2>&1
#         S3="{params.hprc_url}"

#         BASE={wildcards.sample}
#         IFS='_' read -r F1 F2 _ <<< "$BASE"
#         IFS='.' read -r S1 S2 _ <<< "$F2"
#         SAMPLE="$S1"
#         SAMPLE_HAP="${{S1}}_${{S2}}"
#         SAMPLE_FULL="${{F1}}_${{S1}}.${{S2}}"

#         PFX="submissions/DC27718F-5F38-43B0-9A78-270F395F13E8--INT_ASM_PRODUCTION/${{SAMPLE}}/assemblies/freeze_2/annotation/chrom_assignment/mashmap/asm_to_chm13v2/"
#         KEY=$(curl -fsS "$S3/?list-type=2&prefix=$PFX" \
#             | grep -o '<Key>[^<]*</Key>' | sed 's|</*Key>||g' \
#             | grep "/mashmap_${{SAMPLE_HAP}}_.*\_pi95\.paf$" | head -1)

#         curl -fsSL -o {output.mashmap} "$S3/$KEY"
#         """

### Download and merge chromosome aliases
rule download_chrom_alias:
    output:
        chromalias = os.path.join(WORKDIR, "data/annotations/chromalias/{sample}.chromAlias.txt")
    params:
        hprc_url = config["hprc_url"]
    log:
        "logs/download/chromalias/download_{sample}_chromalias.log"
    shell:
        r"""
        exec > {log} 2>&1
        S3="{params.hprc_url}"

        BASE={wildcards.sample}
        IFS='_' read -r F1 F2 _ <<< "$BASE"
        IFS='.' read -r S1 S2 _ <<< "$F2"
        SAMPLE="$S1"
        SAMPLE_HAP="${{S1}}_${{S2}}"
        SAMPLE_FULL="${{F1}}_${{S1}}.${{S2}}"

        PFX="submissions/DC27718F-5F38-43B0-9A78-270F395F13E8--INT_ASM_PRODUCTION/${{SAMPLE}}/assemblies/freeze_2/annotation/chrom_assignment/"
        KEY=$(curl -fsS "$S3/?list-type=2&prefix=$PFX" \
            | grep -o '<Key>[^<]*</Key>' | sed 's|</*Key>||g' \
            | grep "/${{SAMPLE_HAP}}_.*\.chromAlias\.txt$" | head -1)

        curl -fsSL -o {output.chromalias} "$S3/$KEY"
        """

rule create_chromalias_table:
    input:
        chromalias = expand(
            os.path.join(WORKDIR, "data/annotations/chromalias/{sample}.chromAlias.txt"),
            sample = sample_list
        )
    output:
        chrom_table = os.path.join(WORKDIR, "results/annotations/chromalias_table.txt")
    log:
        "logs/annotations/chromalias/create_chromalias_table.txt"
    shell:
        """
        printf 'assembly\\tucsc\\tgenbank\\n' > {output.chrom_table} 2> {log}
        awk 'FNR>1' {input.chromalias} >> {output.chrom_table} 2>> {log}
        """

### Download centromeric and pericentromeric coordinates
# rule download_cenSat_bed:
#     output:
#         bed = os.path.join(WORKDIR, "data/annotations/cenSat/{sample}.active.centromeres.bed")
#     params:
#         hprc_url = config["hprc_url"]
#     log:
#         "logs/download/cenSat/download_{sample}_cenSat_anno.log"
#     shell:
#         r"""
#         exec > {log} 2>&1
#         S3="{params.hprc_url}"

#         BASE={wildcards.sample}
#         IFS='_' read -r F1 F2 _ <<< "$BASE"
#         IFS='.' read -r S1 S2 _ <<< "$F2"
#         SAMPLE="$S1"
#         SAMPLE_HAP="${{S1}}_${{S2}}"
#         SAMPLE_FULL="${{F1}}_${{S1}}.${{S2}}"

#         PFX="submissions/DC27718F-5F38-43B0-9A78-270F395F13E8--INT_ASM_PRODUCTION/${{SAMPLE}}/assemblies/freeze_2/annotation/censat/"
#         KEY=$(curl -fsS "$S3/?list-type=2&prefix=$PFX" \
#             | grep -o '<Key>[^<]*</Key>' | sed 's|</*Key>||g' \
#             | grep "/${{SAMPLE_HAP}}_.*\.active.centromeres\.bed$" | head -1)

#         curl -fsSL -o {output.bed} "$S3/$KEY"
#         """

rule download_cenSat_bed:
    output:
        bed = os.path.join(WORKDIR, "data/annotations/cenSat/{sample}.cenSat.bed")
    params:
        hprc_url = config["hprc_url"]
    log:
        "logs/download/cenSat/download_{sample}_cenSat_anno.log"
    shell:
        r"""
        exec > {log} 2>&1
        S3="{params.hprc_url}"

        BASE={wildcards.sample}
        IFS='_' read -r F1 F2 _ <<< "$BASE"
        IFS='.' read -r S1 S2 _ <<< "$F2"
        SAMPLE="$S1"
        SAMPLE_HAP="${{S1}}_${{S2}}"
        SAMPLE_FULL="${{F1}}_${{S1}}.${{S2}}"

        PFX="submissions/DC27718F-5F38-43B0-9A78-270F395F13E8--INT_ASM_PRODUCTION/${{SAMPLE}}/assemblies/freeze_2/annotation/censat/"
        KEY=$(curl -fsS "$S3/?list-type=2&prefix=$PFX" \
            | grep -o '<Key>[^<]*</Key>' | sed 's|</*Key>||g' \
            | grep "/${{SAMPLE_HAP}}_.*\.cenSat\.bed$" | head -1)

        curl -fsSL -o {output.bed} "$S3/$KEY"
        """

rule process_cenSat_bed:
    input:
        bed = os.path.join(WORKDIR, "data/annotations/cenSat/{sample}.cenSat.bed")
    output:
        cm_bed = os.path.join(WORKDIR, "results/cenSat/chrALL/{sample}.centromeric.bed"),
        pcm_bed = os.path.join(WORKDIR, "results/cenSat/chrALL/{sample}.pericentromeric.bed")
    log:
        "logs/annotations/cenSat/process_{sample}_cenSat_bed.log"
    shell:
        r"""
        awk -F'\t' '!/^track|^#|^browser/ {{
            f4 = tolower($4)
            if (f4 ~ /^active_hor\(/) print
        }}' {input.bed} > {output.cm_bed} 2> {log}

        awk -F'\t' '!/^track|^#|^browser/ {{
            f4 = tolower($4)
            if (f4 !~ /^active_hor\(/) print
        }}' {input.bed} > {output.pcm_bed} 2> {log}
        """

rule fix_bed_sf_names:
    input:
        bed = os.path.join(WORKDIR, "results/cenSat/chrALL/{sample}.{region}.bed"),
        sf = os.path.join(WORKDIR, "data/sfalias/sf_alias_table.txt")
    output:
        bed = os.path.join(WORKDIR, "results/cenSat/chrALL/{sample}.{region}.sf.bed")
    log:
        "logs/annotations/centromeres/chrALL/add_{sample}_{region}_sf_names_bed.txt"
    script:
        "../scripts/py/add_bed_sf.py"

rule split_cenSat_bed_chr_names:
    input:
        bed = os.path.join(WORKDIR, "results/cenSat/chrALL/{sample}.{region}.sf.bed"),
        chrom_table = os.path.join(WORKDIR, "results/annotations/chromalias_table.txt"),
        # mashmap = os.path.join(WORKDIR, "data/annotations/mashmap/{sample}.pi95.paf"),
        # fai_idx = os.path.join(WORKDIR, "data/asm/fa_idx/{sample}.fa.gz.fai")
    output:
        bed = os.path.join(WORKDIR, "results/cenSat/{chr}/{sample}.{region}.sf.bed")
    params:
        chrom = "{chr}"
    log:
        "logs/annotations/centromeres/{chr}/split_{sample}_{region}_chr_names_bed.txt"
    script:
        "../scripts/py/split_bed_chr.py"

rule create_bed_list_cm:
    input:
        bed_paths = expand(
            os.path.join(WORKDIR, "results/cenSat/{chr}/{sample}.{region}.sf.bed"),
            chr="{chr}",
            sample=sample_list,
            region="{region}"
        )
    output:
        bed_list = os.path.join(WORKDIR, "results/annotations/bed_lists/{chr}/bed_list.{region}.txt")
    log:
        "logs/annotations/centromeres/{chr}/create_{region}_bed_list.txt"
    shell:
        """
        printf '%s\\n' {input.bed_paths} > {output.bed_list} 2> {log}
        """

# rule extend_centromeres_bed:
#     input:
#         bed = os.path.join(WORKDIR, "results/cenSat/{chr}/{sample}.centromeric.bed"),
#         fai_idx = os.path.join(WORKDIR, "data/asm/fa_idx/{sample}.fa.gz.fai")
#     output:
#         bed = os.path.join(WORKDIR, "results/cenSat/{chr}/{sample}.pericentromeric.{span}.bed")
#     params:
#         span = lambda wildcards: config["ranges"][wildcards.span]
#     log:
#         "logs/annotations/centromeres/{chr}/{sample}.extend_{span}_centromeric_bed.txt"
#     script:
#         "../scripts/py/extend_centromeres_boundaries.py"


# rule create_bed_list_pcm:
#     input:
#         bed_paths = expand(
#             os.path.join(WORKDIR, "results/cenSat/{chr}/{sample}.pericentromeric.{span}.bed"),
#             chr="{chr}",
#             sample=sample_list,
#             span="{span}"
#         )
#     output:
#         bed_list = os.path.join(WORKDIR, "results/annotations/bed_lists/{chr}/bed_list.pericentromeric.{span}.txt")
#     log:
#         "logs/annotations/centromeres/{chr}/create_{span}_pericentromeric_bed_list.txt"
#     shell:
#         """
#         printf '%s\\n' {input.bed_paths} > {output.bed_list} 2> {log}
#         """

# rule split_fa_chr:
#     input:
#         fa_file = os.path.join(DATADIR, "{sample}.fa.gz"),
#         chromalias = os.path.join(WORKDIR, "data/annotations/chromalias/{sample}.chromAlias.txt")
#     output:
#         fa_file = os.path.join(WORKDIR, "results/asm/{chr}/{sample}.fa.gz")
#     params:
#         chrom = "{chr}"
#     log:
#         "logs/asm/{chr}/split_{sample}_fasta_by_chr.log"
#     conda:
#         "../envs/py/biopython.yaml"
#     script:
#         "../scripts/py/split_fa_chr.py"






# rule download_pcm_regions_cenhap_paper:
#     output:
#         tds = os.path.join(WORKDIR, "data/reference/chrALL/cenhap_anno.tds")
#     params:
#         url = config["cenhap_regions_url"]
#     log:
#         "logs/ref/chrALL/download_cenhap_regions.txt"
#     shell:
#         """
#         mkdir -p $(dirname {output.tds}) 2> {log};
#         wget -O {output.tds} {params.url} 2> {log};
#         """

# rule process_pcm_regions_cenhap_paper:
#     input:
#         tds = os.path.join(WORKDIR, "data/reference/chrALL/cenhap_anno.tds")
#     output:
#         bed = os.path.join(WORKDIR, "data/reference/chrALL/cenhap_anno.bed")
#     log:
#         "logs/ref/chrALL/process_cenhap_regions.txt"
#     conda:
#         "../envs/py/ml.yaml"
#     script:
#         "../scripts/py/format_cenhap_regions.py"

# rule split_pcm_regions_cenhap_paper:
#     input:
#         bed = os.path.join(WORKDIR, "data/reference/chrALL/cenhap_anno.bed"),
#         chrom_table = os.path.join(WORKDIR, "results/annotations/chromalias_table.txt")
#     output:
#         bed = os.path.join(WORKDIR, "data/reference/{chr}/cenhap_anno.bed")
#     params:
#         chrom = "{chr}"
#     log:
#         "logs/ref/{chr}/split_cenhap_regions.txt"
#     conda:
#         "../envs/py/ml.yaml"
#     script:
#         "../scripts/py/split_bed_chr.py"

# rule download_hs19:
#     output:
#         fa = os.path.join(WORKDIR, "data/reference/GRCh37.fa.gz")
#     params:
#         url = config["hs19_url"]
#     log:
#         "logs/ref/download_hs19.txt"
#     shell:
#         """
#         mkdir -p $(dirname {output.fa}) 2> {log};
#         wget -O {output.fa} {params.url} 2> {log};
#         """