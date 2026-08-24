# Call SNPs in centromeres
rule call_SNPs_pcm:
    input:
        fa_files = expand(
            os.path.join(DATADIR, "{sample}.fa.gz"),
            sample=sample_list
        ),
        bed_list = os.path.join(WORKDIR, "results/annotations/bed_lists/{chr}/bed_list.pericentromeric.{span}.txt")
    output:
        vcf = os.path.join(WORKDIR, "results/snps/centromeres/{chr}/hmn462_unfilt.pericentromeric.{span}.vcf")
    params:
        k = config["kmer_length"],
        m = 0,
        f = 0
    threads:
        config["pgtools_threads"]
    log:
        "logs/snps/centromeres/{chr}/call_{span}_pericentromeric_snps_hmn462_unfilt.log"
    shell:
        """
        pgtools count -v \
            -k {params.k} \
            -m {params.m} \
            -t {threads} \
            -f {params.f} \
            -b {input.bed_list} \
            -o {output.vcf} \
            {input.fa_files} \
            2> {log}
        """