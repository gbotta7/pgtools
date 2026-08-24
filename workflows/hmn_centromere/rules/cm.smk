# Count SNP-mers and rare k-mers in centromeres
rule detect_snpmers_cm:
    input:
        fa_files = expand(
            os.path.join(DATADIR, "{sample}.fa.gz"),
            sample=sample_list
        ),
        bed_list = os.path.join(WORKDIR, "results/annotations/bed_lists/{chr}/bed_list.centromeric.txt")
    output:
        snpmers = os.path.join(WORKDIR, "results/snps/centromeres/{chr}/hmn462_unfilt.centromeric.snpmers.txt")
    params:
        k = config["kmer_length"],
        msf = 0,
        f = 0
    threads:
        config["pgtools_threads"]
    log:
        "logs/snps/centromeres/{chr}/detect_centromeric_snpmers_hmn462_unfilt.log"
    shell:
        """
        pgtools detect -v \
            --snp \
            -k {params.k} \
            --msf {params.msf} \
            -t {threads} \
            -f {params.f} \
            -b {input.bed_list} \
            -o {output.snpmers} \
            {input.fa_files} \
            2> {log}
        """

rule count_snpmers_cm:
    input:
        fa_files = expand(
            os.path.join(DATADIR, "{sample}.fa.gz"),
            sample=sample_list
        ),
        bed_list = os.path.join(WORKDIR, "results/annotations/bed_lists/{chr}/bed_list.centromeric.txt"),
        snpmers = os.path.join(WORKDIR, "results/snps/centromeres/{chr}/hmn462_unfilt.centromeric.snpmers.txt")
    output:
        cts = os.path.join(WORKDIR, "results/snps/centromeres/{chr}/hmn462_unfilt.centromeric.snpmers.tsv")
    params:
        k = config["kmer_length"]
    threads:
        config["pgtools_threads"] * 4
    log:
        "logs/snps/centromeres/{chr}/count_centromeric_snpmers_hmn462_unfilt.log"
    shell:
        """
        pgtools count -v \
            --snp \
            -k {params.k} \
            -t {threads} \
            -b {input.bed_list} \
            --kmers {input.snpmers} \
            -o {output.cts} \
            {input.fa_files} \
            2> {log}
        """

rule detect_rare_kmers_cm:
    input:
        fa_files = expand(
            os.path.join(DATADIR, "{sample}.fa.gz"),
            sample=sample_list
        ),
        bed_list = os.path.join(WORKDIR, "results/annotations/bed_lists/{chr}/bed_list.centromeric.txt")
    output:
        rare_kmers = os.path.join(WORKDIR, "results/snps/centromeres/{chr}/hmn462_unfilt.centromeric.rare_kmers.txt")
    params:
        k = config["kmer_length"],
        msf = 0,
        f = 0,
        mko = 2
    threads:
        config["pgtools_threads"]
    log:
        "logs/snps/centromeres/{chr}/detect_centromeric_rare_kmers_hmn462_unfilt.log"
    shell:
        """
        pgtools detect -v \
            -k {params.k} \
            --msf {params.msf} \
            --mko {params.mko} \
            -t {threads} \
            -f {params.f} \
            -b {input.bed_list} \
            -o {output.rare_kmers} \
            {input.fa_files} \
            2> {log}
        """

rule count_rare_kmers_cm:
    input:
        fa_files = expand(
            os.path.join(DATADIR, "{sample}.fa.gz"),
            sample=sample_list
        ),
        bed_list = os.path.join(WORKDIR, "results/annotations/bed_lists/{chr}/bed_list.centromeric.txt"),
        rare_kmers = os.path.join(WORKDIR, "results/snps/centromeres/{chr}/hmn462_unfilt.centromeric.rare_kmers.txt")
    output:
        cts = os.path.join(WORKDIR, "results/snps/centromeres/{chr}/hmn462_unfilt.centromeric.rare_kmers.tsv")
    params:
        k = config["kmer_length"]
    threads:
        config["pgtools_threads"] * 4
    log:
        "logs/snps/centromeres/{chr}/count_centromeric_rare_kmers_hmn462_unfilt.log"
    shell:
        """
        pgtools count -v \
            --snp \
            -k {params.k} \
            -t {threads} \
            -b {input.bed_list} \
            --kmers {input.rare_kmers} \
            -o {output.cts} \
            {input.fa_files} \
            2> {log}
        """