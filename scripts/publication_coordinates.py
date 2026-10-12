"""Validate model identities before adding genomic coordinate evidence."""
import csv
import json
import os
from collections import defaultdict
from itertools import combinations
from pathlib import Path
from urllib.parse import unquote
from Bio import SeqIO

ROOT = Path(os.environ.get('ECM_PROJECT_ROOT', Path(__file__).resolve().parents[1] / 'workspace')).resolve()
SOURCE = ROOT / 'publication_package' / 'remote_sources'
OUT = ROOT / 'publication_package' / 'coordinate_validation'

def write_tsv(path, rows):
    if not rows:
        path.write_text('No_records\n', encoding='utf-8')
        return
    with path.open('w', encoding='utf-8', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter='\t')
        writer.writeheader(); writer.writerows(rows)

def merged(intervals):
    result=[]
    for start,end in sorted(intervals):
        if result and start<=result[-1][1]+1:
            result[-1]=(result[-1][0],max(result[-1][1],end))
        else:
            result.append((start,end))
    return result

def intersection(a,b):
    return sum(max(0,min(ae,be)-max(ast,bst)+1) for ast,ae in merged(a) for bst,be in merged(b))

def main():
    OUT.mkdir(parents=True, exist_ok=True)
    proteins = {r.id: str(r.seq).rstrip('*') for r in SeqIO.parse(
        SOURCE / '05.OrthoFinder_input/D_japonicum.fa', 'fasta')}
    gff_proteins = {r.id: str(r.seq).rstrip('*') for r in SeqIO.parse(
        SOURCE / 'Dicyema_japonicum/Dicyema_japonicum.strict.complete.renamed.faa', 'fasta')}
    by_sequence = defaultdict(list)
    for key, seq in gff_proteins.items():
        by_sequence[seq].append(key)
    coords, cds = {}, defaultdict(list)
    for line in (SOURCE / 'Dicyema_japonicum/Dicyema_japonicum.strict.complete.renamed.gff3').open(encoding='utf-8'):
        if line.startswith('#'):
            continue
        fields = line.rstrip().split('\t')
        if len(fields) != 9:
            continue
        attrs = dict((k,unquote(v)) for k,v in (a.split('=',1) for a in fields[8].split(';') if '=' in a))
        if fields[2] in ('mRNA','transcript'):
            coords[attrs['ID']] = dict(contig=fields[0], start=int(fields[3]), end=int(fields[4]), strand=fields[6])
        elif fields[2] == 'CDS':
            for parent in attrs.get('Parent','').split(','):
                cds[parent].append((int(fields[3]), int(fields[4])))
    with (ROOT / 'remote_results_cache/07.ECM_Muscle_candidates/D_japonicum.ECM.related.og_map.tsv').open(encoding='utf-8-sig') as f:
        mapping = [line.strip().split('\t')[:2] for line in f if line.strip() and not line.startswith('Orthogroup')]
    results = []
    for og, gene in mapping:
        seq = proteins.get(gene)
        same_id_valid = seq is not None and gff_proteins.get(gene) == seq
        matches = by_sequence.get(seq, [])
        match = gene if same_id_valid else (matches[0] if len(matches)==1 else '')
        coord = coords.get(match, {})
        results.append(dict(Orthogroup=og, analysis_gene_id=gene, coordinate_gene_id=match,
            identity_validation='same_id_exact_sequence' if same_id_valid else ('unique_sequence_remapping' if match else 'unresolved'),
            contig=coord.get('contig',''), start=coord.get('start',''), end=coord.get('end',''),
            strand=coord.get('strand',''), CDS_segments=len(cds.get(match,[])),
            protein_length=len(seq) if seq else 0))
    write_tsv(OUT/'ECM_gene_coordinates.tsv',results)
    groups = defaultdict(list)
    for row in results:
        if row['contig']:
            groups[row['Orthogroup']].append(row)
    pairs=[]
    for og, rows in groups.items():
        for a,b in combinations(rows,2):
            if a['contig'] != b['contig']:
                continue
            overlap=max(0,min(a['end'],b['end'])-max(a['start'],b['start'])+1)
            gap=0 if overlap else max(a['start'],b['start'])-min(a['end'],b['end'])-1
            if gap<=50000:
                ac=cds[a['coordinate_gene_id']]; bc=cds[b['coordinate_gene_id']]
                shared=intersection(ac,bc)
                alen=sum(e-s+1 for s,e in merged(ac)); blen=sum(e-s+1 for s,e in merged(bc))
                pairs.append(dict(Orthogroup=og,gene1=a['analysis_gene_id'],gene2=b['analysis_gene_id'],
                    contig=a['contig'],gap_bp=gap,overlap_bp=overlap,
                    same_strand=a['strand']==b['strand'],
                    shared_CDS_bp=shared, CDS_overlap_fraction_shorter=shared/min(alen,blen) if min(alen,blen) else 0,
                    evidence='overlapping_models_require_review' if overlap else 'distinct_neighboring_models'))
    write_tsv(OUT/'ECM_neighboring_model_pairs.tsv',pairs)
    sensitivity=[]
    locus_groups=[]
    for threshold in [0.25,0.5,0.8]:
        for og,rows in groups.items():
            parents={x['analysis_gene_id']:x['analysis_gene_id'] for x in rows}
            def find(k):
                while parents[k]!=k:
                    parents[k]=parents[parents[k]]; k=parents[k]
                return k
            for pair in pairs:
                if pair['Orthogroup']==og and pair['same_strand'] and pair['CDS_overlap_fraction_shorter']>=threshold:
                    parents[find(pair['gene2'])]=find(pair['gene1'])
            components=defaultdict(list)
            for row in rows:
                components[find(row['analysis_gene_id'])].append(row['analysis_gene_id'])
            sensitivity.append(dict(Orthogroup=og,threshold=threshold,raw_gene_models=len(rows),
                                    overlap_collapsed_components=len(components)))
            if threshold==0.5:
                for idx,genes in enumerate(components.values(),1):
                    representative=max(genes,key=lambda k:len(proteins[k]))
                    locus_groups.append(dict(Orthogroup=og,component=f'{og}_C{idx}',
                        representative=representative,member_models=','.join(genes),n_models=len(genes)))
    write_tsv(OUT/'ECM_copy_sensitivity.tsv',sensitivity)
    write_tsv(OUT/'ECM_overlap_components.tsv',locus_groups)
    with (OUT/'ECM_representatives.faa').open('w',encoding='utf-8') as f:
        for row in locus_groups:
            f.write('>'+row['representative']+'\n'+proteins[row['representative']]+'\n')
    summary=dict(candidate_genes=len(results),coordinate_resolved=sum(bool(x['contig']) for x in results),
                 same_id_exact_sequence=sum(x['identity_validation']=='same_id_exact_sequence' for x in results),
                 unique_remapping=sum(x['identity_validation']=='unique_sequence_remapping' for x in results),
                 neighboring_pairs_50kb=len(pairs),overlapping_pairs=sum(x['overlap_bp']>0 for x in pairs),
                 neighboring_nonoverlap_pairs_10kb=sum(x['gap_bp']<=10000 and x['overlap_bp']==0 for x in pairs),
                 same_strand_shared_CDS_pairs=sum(x['same_strand'] and x['shared_CDS_bp']>0 for x in pairs),
                 overlap_collapsed_components_by_threshold={str(t):sum(x['overlap_collapsed_components'] for x in sensitivity if x['threshold']==t) for t in [0.25,0.5,0.8]},
                 multi_model_OGs_after_50pct_collapse=sum(x['overlap_collapsed_components']>1 for x in sensitivity if x['threshold']==0.5),
                 caution='Coordinates support distinct gene-model locations, not independently validated duplication or expression.')
    (OUT/'summary.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')
    print(json.dumps(summary,indent=2))

if __name__=='__main__':
    main()
