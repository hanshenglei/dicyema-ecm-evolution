"""Attach original domain evidence and global OG mapping without new searches."""
import collections
import csv
import json
from publication_domain_validation import BASE, OUT, sha, table

def main():
    if not (OUT / 'SUCCESS.json').exists():
        raise RuntimeError('Primary computation is not complete')
    source = BASE / '05.OrthoFinder_input/OrthoFinder/Results_Jun11/Orthogroups/Orthogroups.tsv'
    ogs = list(csv.DictReader(source.open(), delimiter='\t'))
    membership = {}
    for row in ogs:
        for ident in row['D_japonicum'].split(','):
            if ident.strip():
                membership[ident.strip()] = row['Orthogroup']
    (OUT / 'source_Orthogroups.tsv').write_bytes(source.read_bytes())
    legacy = BASE / '08.ECM_Muscle_annotation/pfam/ECM_Muscle.with_Dj.pfam.domtblout'
    (OUT / 'source_legacy.hmmscan.domtblout').write_bytes(legacy.read_bytes())
    hits = list(csv.DictReader((OUT / 'hits.tsv').open(), delimiter='\t'))
    families = {r['Orthogroup']: r for r in csv.DictReader((OUT / 'source_OG_family_summary.tsv').open(), delimiter='\t')}
    support = []
    index = collections.defaultdict(list)
    for hit in hits:
        if hit['dataset'] == 'D_japonicum':
            index[(hit['protein_id'], hit['model'], hit['mode'])].append(hit)
    for (ident, model, mode), hs in sorted(index.items()):
        og = membership.get(ident, 'unassigned')
        support.append({'protein_id': ident, 'Orthogroup': og, 'mode': mode, 'model': model,
                        'group': hs[0]['group'], 'domain_iE_min': min(float(h['domain_independent_evalue']) for h in hs),
                        'protein_length': hs[0]['protein_length'], 'legacy_assigned_family': families.get(og, {}).get('assigned_family', 'not_in_legacy_family_summary'),
                        'legacy_quality_class': families.get(og, {}).get('quality_class', '')})
    table('Dicyema_all_OG_domain_support.tsv', list(support[0]), support)
    counts = []
    for mode in ['GA', 'E1e-5']:
        for model in sorted({h['model'] for h in support}):
            subset = [h for h in support if h['mode'] == mode and h['model'] == model]
            counts.append({'mode': mode, 'model': model, 'proteins': len(subset),
                           'Orthogroups': ','.join(sorted({h['Orthogroup'] for h in subset})),
                           'protein_ids': ','.join(h['protein_id'] for h in subset)})
    table('Dicyema_all_OG_domain_counts.tsv', list(counts[0]), counts)
    paths = [OUT / 'source_Orthogroups.tsv', OUT / 'source_legacy.hmmscan.domtblout',
             OUT / 'Dicyema_all_OG_domain_support.tsv', OUT / 'Dicyema_all_OG_domain_counts.tsv',
             OUT / 'publication_domain_genome_check.py', OUT / 'publication_domain_supplement.py']
    paths += [p for p in sorted((OUT / 'genome_check').glob('*')) if p.is_file()]
    table('supplement_checksums.tsv', ['file', 'bytes', 'sha256'],
          [{'file': str(p.relative_to(OUT)), 'bytes': p.stat().st_size, 'sha256': sha(p)} for p in paths])
    print(json.dumps(counts, indent=2))

if __name__ == '__main__':
    main()
