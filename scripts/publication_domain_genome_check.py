"""Bounded, optional six-query TBLASTN against the original genome (read-only)."""
import csv
import json
import os
import shutil
import subprocess
from publication_domain_validation import BASE, OUT, fasta, sha, table

BLAST = shutil.which(os.environ.get('ECM_TBLASTN', 'tblastn')) or os.environ.get('ECM_TBLASTN', 'tblastn')

def main():
    if not (OUT / 'SUCCESS.json').exists():
        raise RuntimeError('Primary whole-proteome HMMER must finish first')
    dest = OUT / 'genome_check'
    dest.mkdir(exist_ok=False)
    hits = list(csv.DictReader((OUT / 'hits.tsv').open(), delimiter='\t'))
    queries = []
    for ds in ['C_teleta', 'C_gigas']:
        proteins = dict(fasta(OUT / 'inputs' / (ds + '.faa')))
        for model in ['C4', 'Laminin_N']:
            candidates = [h for h in hits if h['dataset'] == ds and h['mode'] == 'GA' and h['model'] == model
                          and float(h['model_coverage']) >= 0.7]
            best = min(candidates, key=lambda h: float(h['domain_independent_evalue']))
            lo, hi = max(1, int(best['ali_from']) - 20), min(len(proteins[best['protein_id']]), int(best['ali_to']) + 20)
            queries.append({'query_id': ds + '_' + model, 'role': 'cross_species_basement_domain', 'dataset': ds,
                            'model': model, 'protein_id': best['protein_id'], 'from_aa': lo, 'to_aa': hi,
                            'profile_iE': best['domain_independent_evalue'], 'sequence': proteins[best['protein_id']][lo-1:hi]})
    proteins = dict(fasta(OUT / 'inputs/D_japonicum.faa'))
    for model in ['Actin', 'Integrin_beta']:
        candidates = [h for h in hits if h['dataset'] == 'D_japonicum' and h['mode'] == 'GA' and h['model'] == model]
        best = min(candidates, key=lambda h: float(h['domain_independent_evalue']))
        lo, hi = max(1, int(best['ali_from']) - 20), min(len(proteins[best['protein_id']]), int(best['ali_to']) + 20)
        queries.append({'query_id': 'D_japonicum_' + model, 'role': 'within_species_positive_control',
                        'dataset': 'D_japonicum', 'model': model, 'protein_id': best['protein_id'],
                        'from_aa': lo, 'to_aa': hi, 'profile_iE': best['domain_independent_evalue'],
                        'sequence': proteins[best['protein_id']][lo-1:hi]})
    query_path = dest / 'queries.faa'
    query_path.write_text(''.join('>' + q['query_id'] + '\n' + q['sequence'] + '\n' for q in queries))
    table('genome_check/query_provenance.tsv', list(queries[0]), queries)
    genome = BASE / 'Dicyema_japonicum/Dicyema_japonicum.genome.fa'
    version = subprocess.check_output([BLAST, '-version'], text=True)
    (dest / 'tblastn_version.txt').write_text(version)
    output = dest / 'hits.outfmt6.tsv'
    fields = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend',
              'sstart', 'send', 'evalue', 'bitscore', 'qlen', 'qcovhsp', 'sframe']
    cmd = [BLAST, '-query', str(query_path), '-subject', str(genome), '-num_threads', '4',
           '-seg', 'yes', '-evalue', '1e-5', '-max_target_seqs', '100', '-max_hsps', '20',
           '-outfmt', '6 ' + ' '.join(fields), '-out', str(output)]
    provenance = {'argv': cmd, 'genome': str(genome), 'genome_sha256': sha(genome),
                  'query_sha256': sha(query_path), 'tblastn_sha256': sha(__import__('pathlib').Path(BLAST)),
                  'timeout_seconds': 180, 'status': 'started',
                  'scope': 'Four domains from two positive-control spiralians; two within-Dicyema positive controls. No absence/loss proof and no general genome-wide annotation.'}
    try:
        with (dest / 'run.log').open('w') as log:
            subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT, check=True, timeout=180)
        rows = list(csv.DictReader(output.open(), fieldnames=fields, delimiter='\t'))
        table('genome_check/hits.tsv', fields, rows)
        summary = []
        for q in queries:
            hs = [r for r in rows if r['qseqid'] == q['query_id']]
            summary.append({'query_id': q['query_id'], 'role': q['role'], 'model': q['model'],
                            'HSPs_E1e5': len(hs), 'contigs': len({h['sseqid'] for h in hs}),
                            'best_E': min((float(h['evalue']) for h in hs), default=''),
                            'best_bitscore': max((float(h['bitscore']) for h in hs), default='')})
        table('genome_check/query_summary.tsv', list(summary[0]), summary)
        for row in summary:
            if row['role'] == 'within_species_positive_control' and row['HSPs_E1e5'] == 0:
                raise RuntimeError('Within-species genomic positive control failed: ' + row['query_id'])
        provenance['status'] = 'completed_positive_controls_passed'
        print(json.dumps(summary, indent=2), flush=True)
    except subprocess.TimeoutExpired:
        provenance['status'] = 'bounded_timeout_no_negative_inference'
    finally:
        (dest / 'provenance.json').write_text(json.dumps(provenance, indent=2) + '\n')

if __name__ == '__main__':
    main()
