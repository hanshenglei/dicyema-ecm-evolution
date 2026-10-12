"""Independently verify downloaded domain validation, then write compact reports."""
import collections
import csv
import hashlib
import json
import os
from pathlib import Path

ROOT = Path(os.environ.get('ECM_PROJECT_ROOT', Path(__file__).resolve().parents[1] / 'workspace')).resolve()
OUT = ROOT / 'publication_package/domain_validation'
KEYS = ['Collagen', 'C4', 'Laminin_N', 'Laminin_G_1', 'Laminin_G_2', 'Laminin_EGF',
        'Integrin_beta', 'FG-GAP', 'Cadherin', 'Actin', 'Myosin_head', 'Tropomyosin', 'Troponin']

def read(name):
    with (OUT / name).open() as handle:
        return list(csv.DictReader(handle, delimiter='\t'))

def write(name, rows, fields=None):
    with (OUT / name).open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields or list(rows[0]), delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)

def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()

def fasta(path):
    result, ident = {}, None
    with path.open() as handle:
        for line in handle:
            if line.startswith('>'):
                ident = line[1:].split()[0]
                assert ident not in result
                result[ident] = ''
            else:
                result[ident] += line.strip()
    return result

def main():
    success = json.loads((OUT / 'SUCCESS.json').read_text())
    assert success['completed_searches'] == 28
    for item in read('checksums.tsv'):
        path = OUT / item['file']
        assert path.stat().st_size == int(item['bytes']), path
        assert sha(path) == item['sha256'], path
    models = {r['NAME']: r for r in read('models.tsv')}
    inputs = read('inputs.tsv')
    fastas = {r['dataset']: fasta(OUT / 'inputs' / (r['dataset'] + '.faa')) for r in inputs}
    independent = []
    commands = [json.loads(line) for line in (OUT / 'commands.jsonl').read_text().splitlines()]
    assert len(commands) == 28
    for cmd in commands:
        raw_name = Path(cmd['argv'][cmd['argv'].index('--domtblout') + 1]).name
        path = OUT / 'raw' / raw_name
        lines = path.read_text().splitlines()
        assert lines[-1] == '# [ok]', path
        for line in lines:
            if not line or line.startswith('#'):
                continue
            p = line.split()
            ds, mode = cmd['dataset'], cmd['mode']
            name = p[3]
            assert p[0] in fastas[ds]
            assert int(p[2]) == len(fastas[ds][p[0]])
            assert p[4] == models[name]['ACC']
            if mode == 'GA':
                ga1, ga2 = [float(x) for x in models[name]['GA'].strip(';').split()]
                assert float(p[7]) >= ga1 - 0.11
                assert float(p[13]) >= ga2 - 0.11
            elif float(p[6]) > 1e-5 or float(p[12]) > 1e-5:
                continue
            independent.append((ds, mode, p[0], name, int(p[17]), int(p[18]), float(p[6]), float(p[12])))
    hits = read('hits.tsv')
    parsed = [(r['dataset'], r['mode'], r['protein_id'], r['model'], int(r['ali_from']),
               int(r['ali_to']), float(r['sequence_evalue']), float(r['domain_independent_evalue'])) for r in hits]
    assert collections.Counter(independent) == collections.Counter(parsed), 'Independent raw parsing mismatch'
    indexed = collections.defaultdict(set)
    group_indexed = collections.defaultdict(set)
    domain_indexed = collections.defaultdict(set)
    for r in hits:
        indexed[(r['dataset'], r['mode'], r['model'])].add(r['protein_id'])
        key = r['dataset'], r['mode'], r['group']
        group_indexed[key].add(r['protein_id'])
        domain_indexed[key].add((r['protein_id'], r['model'], r['ali_from'], r['ali_to']))
    for row in read('counts_long.tsv'):
        key = row['dataset'], row['mode'], row['group']
        assert len(group_indexed[key]) == int(row['hit_proteins'])
        assert len(domain_indexed[key]) == int(row['hit_domains'])
    models_long = []
    core_wide = []
    for ds in inputs:
        assert len(fastas[ds['dataset']]) == int(ds['proteins'])
        for mode in ['GA', 'E1e-5']:
            row = {'dataset': ds['dataset'], 'species': ds['species'], 'proteins': ds['proteins'], 'mode': mode}
            for name in KEYS:
                row[name] = len(indexed[(ds['dataset'], mode, name)])
            core_wide.append(row)
            for name in sorted(models):
                models_long.append({'dataset': ds['dataset'], 'species': ds['species'], 'mode': mode,
                                    'model': name, 'accession': models[name]['ACC'],
                                    'hit_proteins': len(indexed[(ds['dataset'], mode, name)])})
    write('per_species_core_counts.tsv', core_wide)
    write('per_species_model_counts.tsv', models_long)
    full, strict, renamed = (fastas[k] for k in ['Dicyema_unfiltered', 'Dicyema_strict_original', 'D_japonicum'])
    strict_sequence_set, full_sequence_set = set(strict.values()), set(full.values())
    missing_strict = [ident for ident, seq in strict.items() if seq not in full_sequence_set]
    assert set(strict.values()) == set(renamed.values()), 'Strict-original and OrthoFinder proteins differ'
    report = {r['id']: r for r in read('source_Dicyema_filter_report.tsv')}
    for r in read('filtered_rescue_ids.tsv'):
        ident = r['unfiltered_id']
        assert full[ident] not in strict_sequence_set
        assert ident in report
        assert report[ident]['status'] == 'FAIL', (ident, report[ident])
    controls = []
    for name in ['C4', 'Laminin_N']:
        for mode in ['GA', 'E1e-5']:
            positive = [r['dataset'] for r in inputs if r['role'] == 'OrthoFinder_strict'
                        and r['dataset'] != 'D_japonicum' and indexed[(r['dataset'], mode, name)]]
            assert len(positive) >= 2, (name, mode, positive)
            controls.append({'model': name, 'mode': mode, 'positive_species': positive})
    legacy = read('legacy_Dicyema_reparsed.tsv')
    legacy_evalue_discordance = sum(float(r['sequence_evalue']) <= 1e-5 < float(r['domain_independent_evalue']) for r in legacy)
    ga_cad = [r for r in read('adhesion_OG_audit.tsv') if r['legacy_family'] == 'Cadherin' and r['mode'] == 'GA']
    cad_core_ogs = [r['Orthogroup'] for r in ga_cad if r['core_supported_ids']]
    supplementary_only = [r['Orthogroup'] for r in ga_cad if not r['core_supported_ids'] and r['new_models']]
    no_support = [r['Orthogroup'] for r in ga_cad if not r['new_models']]
    verification = {'status': 'PASS', 'checksums_verified': len(read('checksums.tsv')),
                    'independent_raw_domain_rows': len(independent), 'searches_completed_verified': len(commands),
                    'all_432_group_count_rows_verified': len(read('counts_long.tsv')),
                    'full_proteins': len(full), 'strict_proteins': len(strict),
                    'strict_original_and_OrthoFinder_sequence_sets_identical': True,
                    'strict_ids_without_exact_full_sequence': missing_strict,
                    'filtered_rescue_rows_verified_against_FAIL_report': len(read('filtered_rescue_ids.tsv')),
                    'positive_controls': controls, 'legacy_sequence_significant_but_domain_iE_not_1e5_rows': legacy_evalue_discordance,
                    'cadherin_GA_core_supported_OGs': cad_core_ogs,
                    'cadherin_GA_auxiliary_or_shared_only_OGs': supplementary_only,
                    'cadherin_GA_no_tested_profile_OGs': no_support}
    (OUT / 'verification.json').write_text(json.dumps(verification, indent=2) + '\n')
    selected = [r for r in core_wide if r['mode'] == 'GA' and r['dataset'] != 'Dicyema_strict_original']
    columns = ['dataset', 'proteins'] + KEYS
    md = ['# New Dicyema Domain-Profile Validation', '',
          'Completed and independently checked 28 HMMER searches: 12 inputs (unfiltered Dicyema, original-ID strict Dicyema, ten OrthoFinder proteomes), GA and E<=1e-5, plus four exhaustive core-profile searches in unfiltered and strict Dicyema.', '',
          '## GA Protein Counts', '', '| ' + ' | '.join(columns) + ' |',
          '| ' + ' | '.join(['---'] * len(columns)) + ' |']
    md += ['| ' + ' | '.join(str(r[c]) for c in columns) + ' |' for r in selected]
    md += ['', 'Counts are unique protein IDs per model, not genes or independent domains. Default and --max results are unioned for focal Dicyema; alternative/clan-related models can overlap. E-value counts are in per_species_core_counts.tsv.', '',
           '## Evidence Boundaries', '',
           '- GA calls use curated per-sequence and per-domain bit scores. Sensitivity calls require BOTH whole-sequence E<=1e-5 and independent domain iE<=1e-5 (domtblout zero-based field 12); conditional domain E is not the final criterion.',
           '- Non-detection is limited to the supplied predicted proteins and tested models; it does not demonstrate genomic deletion, evolutionary loss, or causal body-plan simplification.',
           '- Collagen-repeat hits alone are not collagen-IV identity. C4 is the type-IV NC1-associated model; Laminin_N is more diagnostic than shared EGF/Laminin_G domains. Full protein architecture, signal peptides, transmembrane domains, phylogeny and contamination checks remain separate validation steps.',
           '- The Troponin and Troponin-I_N models do not comprehensively test the entire troponin complex. Actin/Myosin_head/Tropomyosin domains do not demonstrate an organized contractile muscle system.',
           '- Strict filtering is controlled by exact normalized sequence matching and original FAIL reasons; rescue IDs can be fragments, isoforms or alternative predictions rather than additional genes.',
           '- Legacy Integrin-family labels include unrelated or integrin-binding families. Do not aggregate TSP_1, Reprolysin, PP2C, Disintegrin, EF-hand or ITFG2 as integrin receptor retention.', '',
           '## Cadherin OG Audit', '',
           'GA core-domain-supported OGs: ' + ', '.join(cad_core_ogs) + '.',
           'GA auxiliary/shared-only OGs: ' + (', '.join(supplementary_only) or 'none') + '.',
           'GA no tested profile support: ' + (', '.join(no_support) or 'none') + '.',
           'See adhesion_OG_audit.tsv for exact genes, new model names, GA/E calls and legacy independent domain E-values. Motifs and OG membership alone are not domain support.', '',
           '## Provenance', '',
           'Existing HMMER 3.4 executable and local Pfam-A.hmm were read without modification. Exact database SHA-256, all model accession versions/GA scores, input hashes, sequence normalization and argv are in provenance.json, models.tsv, inputs.tsv, input_normalization.tsv and commands.jsonl. The local database release was not explicitly labeled and is intentionally not guessed.', '',
           'Re-run remotely in a fresh authorized directory using publication_domain_validation.py; the worker refuses to overwrite an existing command log. Transport uses publication_domain_remote.py and reads ECM_SSH_PASSWORD only from the process environment. publication_remote_audit.py is untouched.', '',
           'Official threshold reference: [HMMER User Guide](https://eddylab.org/software/hmmer/CURRENT/Userguide.pdf).', '']
    (OUT / 'RESULTS.md').write_text('\n'.join(md))
    print(json.dumps(verification, indent=2))
    print('\nGA CORE COUNTS')
    for row in core_wide:
        if row['mode'] == 'GA':
            print(json.dumps(row))

if __name__ == '__main__':
    main()
