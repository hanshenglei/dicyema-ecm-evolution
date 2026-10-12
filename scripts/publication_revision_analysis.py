"""Summarize matched annotation QC and map MEME occurrences to Pfam intervals."""
import argparse
import collections
import csv
import hashlib
import json
import re
from pathlib import Path


def read(path):
    with path.open(encoding='utf-8-sig', newline='') as f:
        return list(csv.DictReader(f, delimiter='\t'))


def write(path, rows):
    with path.open('w', encoding='utf-8', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter='\t')
        w.writeheader(); w.writerows(rows)


def intervals(text):
    return [(m.group(1), int(m.group(2)), int(m.group(3))) for m in re.finditer(r'([^;]+):(\d+)-(\d+)', text)]


def marker_status(path):
    status = {}
    for line in path.read_text().splitlines():
        if line.startswith('#') or not line.strip():
            continue
        t = line.split('\t')
        if t[0] in status and status[t[0]] != t[1]:
            raise ValueError('Inconsistent marker status '+t[0])
        status[t[0]] = t[1]
    return status


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--package-root', required=True, type=Path)
    parser.add_argument('--revision-root', required=True, type=Path)
    args = parser.parse_args()
    pkg, rev = args.package_root, args.revision_root
    out = rev/'analysis'
    out.mkdir(parents=True, exist_ok=True)
    genes = read(pkg/'analysis/gene_evidence.tsv')
    mapped = []
    for r in genes:
        domains, motifs = intervals(r['domain_architecture']), intervals(r['motif_architecture'])
        dom_res = set()
        for _, lo, hi in domains:
            dom_res.update(range(lo, hi+1))
        for motif, lo, hi in motifs:
            residues = set(range(lo, hi+1))
            overlap = sorted({name for name, a, b in domains if max(lo, a) <= min(hi, b)})
            mapped.append({'Orthogroup': r['Orthogroup'], 'gene_id': r['gene_id'], 'motif': motif,
                           'start': lo, 'end': hi, 'overlap_residues': len(residues & dom_res),
                           'motif_length': len(residues), 'overlap_fraction': len(residues & dom_res)/len(residues),
                           'overlapping_domains': ','.join(overlap)})
    write(out/'motif_domain_overlap.tsv', mapped)
    per_gene=[]
    for r in genes:
        selected=[m for m in mapped if m['Orthogroup']==r['Orthogroup'] and m['gene_id']==r['gene_id']]
        per_gene.append({'Orthogroup':r['Orthogroup'],'gene_id':r['gene_id'],'motif_occurrences':len(selected),
                         'domain_overlapping_occurrences':sum(m['overlap_residues']>0 for m in selected),
                         'majority_domain_overlapping_occurrences':sum(m['overlap_fraction']>=.5 for m in selected),
                         'total_motif_residues':sum(m['motif_length'] for m in selected),
                         'domain_overlapping_motif_residues':sum(m['overlap_residues'] for m in selected)})
    write(out/'motif_domain_per_model.tsv',per_gene)
    fasta = pkg/'remote_sources/05.OrthoFinder_input/D_japonicum.fa'
    sequences = {}
    for line in fasta.open():
        if line.startswith('>'):
            ident = line[1:].split()[0]
            sequences[ident] = ''
        else:
            sequences[ident] += line.strip()
    checks=[]
    for r in {r['gene_id']:r for r in genes}.values():
        seq=sequences[r['gene_id']]
        assert len(seq)==int(r['protein_length'])
        for label,pattern in [('metzincin_HEXXHXXGXXH','HE..H..G..H'),('astacin_HEXXHXXGXXHE','HE..H..G..HE'),('canonical_NPxY','NP.Y')]:
            found=list(re.finditer(pattern,seq))
            checks.append({'gene_id':r['gene_id'],'sequence_sha256':hashlib.sha256(seq.encode()).hexdigest(),
                           'sequence_length':len(seq),'pattern_name':label,'regular_expression':pattern,
                           'matches':len(found),'one_based_intervals':';'.join(str(m.start()+1)+'-'+str(m.end()) for m in found),
                           'matched_sequences':';'.join(m.group() for m in found),
                           'interpretation':'sequence-pattern check; distinct from de novo MEME motifs; no enzyme or receptor activity assay'})
    write(out/'sequence_motif_checks.tsv',checks)
    import scipy
    from scipy.stats import fisher_exact
    previous=json.loads((pkg/'analysis/summary.json').read_text())['fisher_disjoint_candidates']
    recomputed=fisher_exact(previous['table'],alternative='two-sided')
    assert abs(recomputed.statistic-previous['odds_ratio'])<1e-14
    assert abs(recomputed.pvalue-previous['pvalue'])<1e-18
    (out/'statistical_reproduction.json').write_text(json.dumps({'scipy_version':scipy.__version__,'table':previous['table'],'alternative':'two-sided','odds_ratio':recomputed.statistic,'pvalue':recomputed.pvalue,'matches_archived_results':True},indent=2))
    (out/'motif_source_provenance.json').write_text(json.dumps({'protein_file':str(fasta),'protein_sha256':hashlib.sha256(fasta.read_bytes()).hexdigest(),'domain_motif_source_sha256':hashlib.sha256((pkg/'analysis/gene_evidence.tsv').read_bytes()).hexdigest(),'models':len(genes),'distinct_proteins':len({r['gene_id'] for r in genes}),'intervals':'one-based inclusive; minimum-overlap sensitivity defined at 50 percent of each motif occurrence'},indent=2))
    hmms = read(pkg/'domain_validation/models.tsv')
    group_map={r['NAME']:r['group'] for r in hmms}
    published=rev/'published_profiles/hits.tsv'
    if published.exists():
        hits=read(published)
        rows=[]
        for mode in ['GA','E1e-5']:
            for group in sorted(set(group_map.values())):
                sub=[h for h in hits if h['mode']==mode and group_map[h['model']]==group]
                rows.append({'dataset':'Published_5012','mode':mode,'group':group,
                             'hit_proteins':len({h['protein_id'] for h in sub}),
                             'hit_domains':len({(h['protein_id'],h['model'],h['ali_from'],h['ali_to']) for h in sub})})
        write(out/'published_domain_counts.tsv',rows)
        comparisons=[]
        allcounts=read(pkg/'domain_validation/counts_long.tsv')
        for r in rows:
            for name in ['Dicyema_unfiltered','Dicyema_strict_original']:
                s=next(z for z in allcounts if z['dataset']==name and z['group']==r['group'] and z['mode']==r['mode'])
                comparisons.append({'domain_group':r['group'],'mode':r['mode'],'published_models':r['hit_proteins'],
                                    'current_dataset':name,'current_models':int(s['hit_proteins']),
                                    'difference_in_hit_bearing_predictions':int(s['hit_proteins'])-r['hit_proteins']})
        write(out/'published_current_domain_comparison.tsv',comparisons)
    busco=rev/'matched_busco'
    if busco.exists():
        rows=[];status={}
        for label in ['Published_5012','Integrated_17260','Strict_12402']:
            directory=busco/label
            path=next(directory.glob('short_summary.specific.*.json'))
            data=json.loads(path.read_text())
            r=data['results']
            assert data['versions']['busco']=='5.8.2'
            assert data['lineage_dataset']['creation_date']=='2024-01-08'
            assert r['n_markers']==255
            rows.append({'dataset':label, 'complete':r['Complete BUSCOs'],'single_copy':r['Single copy BUSCOs'],
                         'duplicated':r['Multi copy BUSCOs'],'fragmented':r['Fragmented BUSCOs'],'missing':r['Missing BUSCOs'],
                         'C_percent':r['Complete percentage'],'S_percent':r['Single copy percentage'],
                         'D_percent':r['Multi copy percentage'],'F_percent':r['Fragmented percentage'],'M_percent':r['Missing percentage']})
            status[label]=marker_status(directory/'run_eukaryota_odb10/full_table.tsv')
        assert all(len(v)==255 for v in status.values())
        write(out/'matched_BUSCO_summary.tsv',rows)
        pairs=[]
        for marker in sorted(status['Published_5012']):
            pairs.append({'marker':marker, **{k:v[marker] for k,v in status.items()}})
        write(out/'matched_BUSCO_marker_transitions.tsv',pairs)
        counts={}
        for name in ['Integrated_17260','Strict_12402']:
            c=collections.Counter((r['Published_5012'],r[name]) for r in pairs)
            counts[name]={'transitions':{a+' -> '+b:n for (a,b),n in c.items()},
                          'new_complete_markers':sum(n for (a,b),n in c.items() if a not in {'Complete','Duplicated'} and b in {'Complete','Duplicated'}),
                          'lost_complete_markers':sum(n for (a,b),n in c.items() if a in {'Complete','Duplicated'} and b not in {'Complete','Duplicated'})}
        (out/'matched_BUSCO_transitions.json').write_text(json.dumps(counts,indent=2))
    records=[{'file':p.name,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()} for p in out.glob('*') if p.is_file() and p.name!='checksums.tsv']
    write(out/'checksums.tsv',records)
    print(json.dumps({'output':str(out),'models':len(genes),'motif_occurrences':len(mapped),'written':len(records)},indent=2))


if __name__ == '__main__':
    main()
