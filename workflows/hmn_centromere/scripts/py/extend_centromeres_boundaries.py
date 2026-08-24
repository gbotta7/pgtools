import sys
import csv

sys.stdout = sys.stderr = open(snakemake.log[0], "w")

ACROCENTRIC = {"chr13", "chr14", "chr15", "chr21", "chr22"}


def load_contig_length(fai_idx, contig):
    with open(fai_idx) as fh:
        for line in fh:
            fields = line.rstrip("\n").split("\t")
            if fields[0] == contig:
                return int(fields[1])
    return None


def extend_centromere_bed(input_bed, fai_idx, output_bed, span):
    span = int(float(span))

    with open(input_bed) as fh:
        reader = csv.reader(fh, delimiter="\t")
        rows = [row for row in reader if row and not row[0].startswith("#")]

    if not rows:
        open(output_bed, "w").close()
        return

    first_row = rows[0]
    last_row = rows[-1]

    first_contig = first_row[0]
    last_contig = last_row[0]
    first_start = int(first_row[1])
    last_end = int(last_row[2])
    chrom = last_row[3]
    last_contig_len = load_contig_length(fai_idx, last_contig)

    preceding_contig = first_row[4] if len(first_row) > 4 else ""
    following_contig = last_row[5] if len(last_row) > 5 else ""

    out_lines = []

    if chrom in ACROCENTRIC:
        if following_contig:
            following_len = load_contig_length(fai_idx, following_contig)
            new_end = min(span, following_len)
            out_lines.append([following_contig, "0", str(new_end), chrom])
        else:
            new_end = min(last_end + span, last_contig_len)
            out_lines.append([last_contig, str(last_end), str(new_end), chrom])
    else:
        if preceding_contig:
            preceding_len = load_contig_length(fai_idx, preceding_contig)
            new_start = max(0, preceding_len - span)
            out_lines.append([preceding_contig, str(new_start), str(preceding_len), chrom])
        else:
            new_start = max(0, first_start - span)
            out_lines.append([first_contig, str(new_start), str(first_start), chrom])

        if following_contig:
            following_len = load_contig_length(fai_idx, following_contig)
            new_end = min(span, following_len)
            out_lines.append([following_contig, "0", str(new_end), chrom])
        else:
            new_end = min(last_end + span, last_contig_len)
            out_lines.append([last_contig, str(last_end), str(new_end), chrom])

    with open(output_bed, "w") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerows(out_lines)


extend_centromere_bed(
    input_bed=snakemake.input["bed"],
    fai_idx=snakemake.input["fai_idx"],
    output_bed=snakemake.output["bed"],
    span=snakemake.params["span"],
)