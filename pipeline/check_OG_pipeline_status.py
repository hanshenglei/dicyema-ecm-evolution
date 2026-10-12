#!/usr/bin/env python3
import os
import argparse
import re

def nonempty(path):
    return os.path.exists(path) and os.path.getsize(path) > 0

def grep_error(path):
    if not os.path.exists(path):
        return "no_log"

    err_patterns = re.compile(
        r"(error|failed|killed|cannot|segmentation|exception|aborted|not found|no such|exceeded|maxsize|fatal)",
        re.IGNORECASE
    )

    hits = []
    with open(path, errors="ignore") as f:
        for line in f:
            if err_patterns.search(line):
                hits.append(line.strip())

    if hits:
        return " | ".join(hits[:3])
    return "NA"

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--root_dir", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    with open(args.out, "w") as out:
        out.write(
            "Orthogroup\tmember_count\tmafft\ttrimal\ttreefile\tmeme_xml\t"
            "pipeline_log\tstatus\terror_hint\n"
        )

        for og in sorted(os.listdir(args.root_dir)):
            og_dir = os.path.join(args.root_dir, og)
            if not os.path.isdir(og_dir):
                continue

            faa = os.path.join(og_dir, f"{og}.faa")
            mafft = os.path.join(og_dir, f"{og}.mafft.fa")
            trim = os.path.join(og_dir, f"{og}.mafft.trim.fa")
            tree = os.path.join(og_dir, f"{og}.mafft.trim.fa.treefile")
            meme = os.path.join(og_dir, "meme", "meme.xml")
            log = os.path.join(og_dir, f"{og}.pipeline.log")

            member_count = 0
            if os.path.exists(faa):
                with open(faa) as f:
                    member_count = sum(1 for line in f if line.startswith(">"))

            has_mafft = nonempty(mafft)
            has_trim = nonempty(trim)
            has_tree = nonempty(tree)
            has_meme = nonempty(meme)
            has_log = nonempty(log)

            if has_tree and has_meme:
                status = "complete_tree_meme"
            elif has_tree and not has_meme:
                status = "tree_only_or_meme_failed"
            elif has_trim and not has_tree:
                status = "iqtree_failed_or_running"
            elif has_mafft and not has_trim:
                status = "trimal_failed_or_running"
            elif not has_mafft:
                status = "mafft_failed_or_not_started"
            else:
                status = "unknown"

            err = grep_error(log)

            out.write(
                f"{og}\t{member_count}\t{has_mafft}\t{has_trim}\t{has_tree}\t"
                f"{has_meme}\t{has_log}\t{status}\t{err}\n"
            )

if __name__ == "__main__":
    main()
