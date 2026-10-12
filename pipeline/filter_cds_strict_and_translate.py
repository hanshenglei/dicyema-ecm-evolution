#!/usr/bin/env python3
import argparse
import re

STOP_CODONS = {"TAA", "TAG", "TGA"}

CODON_TABLE = {
    "TTT":"F","TTC":"F","TTA":"L","TTG":"L",
    "TCT":"S","TCC":"S","TCA":"S","TCG":"S",
    "TAT":"Y","TAC":"Y","TAA":"*","TAG":"*",
    "TGT":"C","TGC":"C","TGA":"*","TGG":"W",
    "CTT":"L","CTC":"L","CTA":"L","CTG":"L",
    "CCT":"P","CCC":"P","CCA":"P","CCG":"P",
    "CAT":"H","CAC":"H","CAA":"Q","CAG":"Q",
    "CGT":"R","CGC":"R","CGA":"R","CGG":"R",
    "ATT":"I","ATC":"I","ATA":"I","ATG":"M",
    "ACT":"T","ACC":"T","ACA":"T","ACG":"T",
    "AAT":"N","AAC":"N","AAA":"K","AAA":"K",
    "AGT":"S","AGC":"S","AGA":"R","AGG":"R",
    "GTT":"V","GTC":"V","GTA":"V","GTG":"V",
    "GCT":"A","GCC":"A","GCA":"A","GCG":"A",
    "GAT":"D","GAC":"D","GAA":"E","GAG":"E",
    "GGT":"G","GGC":"G","GGA":"G","GGG":"G"
}

# 修正 AAG
CODON_TABLE["AAG"] = "K"

def read_fasta(path):
    seqs = {}
    order = []
    name = None
    chunks = []

    with open(path) as f:
        for line in f:
            line = line.rstrip()
            if not line:
                continue

            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(chunks)
                    order.append(name)

                name = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line.strip())

    if name is not None:
        seqs[name] = "".join(chunks)
        order.append(name)

    return seqs, order

def write_fasta(records, path, width=60):
    with open(path, "w") as out:
        for seq_id, seq in records:
            out.write(f">{seq_id}\n")
            for i in range(0, len(seq), width):
                out.write(seq[i:i+width] + "\n")

def translate(cds):
    aa = []
    cds = cds.upper().replace("U", "T")

    for i in range(0, len(cds), 3):
        codon = cds[i:i+3]
        if len(codon) != 3:
            aa.append("X")
        elif re.search("[^ACGT]", codon):
            aa.append("X")
        else:
            aa.append(CODON_TABLE.get(codon, "X"))

    return "".join(aa)

def main():
    parser = argparse.ArgumentParser(
        description="Strictly filter complete CDS and translate to high-confidence protein sequences."
    )

    parser.add_argument("--cds", required=True, help="Input CDS fasta")
    parser.add_argument("--out_prefix", required=True, help="Output prefix")
    parser.add_argument("--min_aa", type=int, default=50, help="Minimum protein length, default 50")
    args = parser.parse_args()

    seqs, order = read_fasta(args.cds)

    pass_cds = []
    pass_faa = []

    report_path = args.out_prefix + ".strict_cds_filter.report.tsv"

    with open(report_path, "w") as report:
        report.write(
            "id\tstatus\treasons\tcds_len\taa_len\tstart_codon\tstop_codon\n"
        )

        for seq_id in order:
            cds = seqs[seq_id].upper().replace("U", "T")
            reasons = []

            cds_len = len(cds)

            if cds_len == 0:
                reasons.append("empty_cds")

            # 最严格：CDS 只能包含 A/C/G/T
            if re.search("[^ACGT]", cds):
                reasons.append("non_ACGT_base")

            # 必须是 3 的倍数
            if cds_len % 3 != 0:
                reasons.append("cds_not_multiple_of_3")

            start_codon = cds[:3] if cds_len >= 3 else "NA"
            stop_codon = cds[-3:] if cds_len >= 3 else "NA"

            # 必须 ATG 起始
            if start_codon != "ATG":
                reasons.append("no_ATG_start")

            # 必须标准终止密码子结尾
            if stop_codon not in STOP_CODONS:
                reasons.append("no_terminal_stop")

            aa_with_stop = translate(cds)

            # 理论上完整 CDS 最后一位应该是 *
            if not aa_with_stop.endswith("*"):
                reasons.append("translated_protein_no_terminal_stop")

            # 去掉末尾终止符号
            aa = aa_with_stop[:-1] if aa_with_stop.endswith("*") else aa_with_stop

            # 内部不能有 *
            if "*" in aa:
                reasons.append("premature_stop_codon")

            # 蛋白中不能有 X
            if "X" in aa:
                reasons.append("unknown_amino_acid_X")

            aa_len = len(aa)

            # 最小蛋白长度
            if aa_len < args.min_aa:
                reasons.append("too_short")

            status = "PASS" if len(reasons) == 0 else "FAIL"

            report.write(
                f"{seq_id}\t{status}\t{';'.join(reasons) if reasons else 'NA'}\t"
                f"{cds_len}\t{aa_len}\t{start_codon}\t{stop_codon}\n"
            )

            if status == "PASS":
                pass_cds.append((seq_id, cds))
                pass_faa.append((seq_id, aa))

    write_fasta(pass_cds, args.out_prefix + ".strict.complete.cds.fa")
    write_fasta(pass_faa, args.out_prefix + ".strict.complete.faa")

    print("Done.")
    print(f"Input CDS: {len(order)}")
    print(f"PASS CDS: {len(pass_cds)}")
    print(f"Output CDS: {args.out_prefix}.strict.complete.cds.fa")
    print(f"Output FAA: {args.out_prefix}.strict.complete.faa")
    print(f"Report: {report_path}")

if __name__ == "__main__":
    main()
