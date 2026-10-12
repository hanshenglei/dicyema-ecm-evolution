"""Targeted Pfam validation with explicit, configurable input/output roots."""
import collections
import csv
import datetime
import hashlib
import json
import os
import platform
import shutil
import shlex
import subprocess
import time
from pathlib import Path

PROJECT = Path(os.environ.get('ECM_PROJECT_ROOT', Path(__file__).resolve().parents[1] / 'workspace')).resolve()
BASE = Path(os.environ.get('ECM_INPUT_ROOT', PROJECT / 'inputs' / 'ecm')).resolve()
OUT = Path(os.environ.get('ECM_DOMAIN_OUTPUT', PROJECT / 'publication_package' / 'domain_validation')).resolve()
DB = Path(os.environ.get('ECM_PFAM_HMM', PROJECT / 'inputs' / 'Pfam-A.hmm')).resolve()
BIN = Path(shutil.which(os.environ.get('ECM_HMMSEARCH', 'hmmsearch')) or os.environ.get('ECM_HMMSEARCH', 'hmmsearch'))
GROUPS = {
    'collagen_repeat': ['Collagen'],
    'collagenIV_NC1': ['C4'],
    'fibrillar_collagen_C': ['COLFI'],
    'laminin_N': ['Laminin_N'],
    'laminin_G': ['Laminin_G_1', 'Laminin_G_2'],
    'laminin_EGF': ['Laminin_EGF'],
    'laminin_other': ['Laminin_B', 'Laminin_I', 'Laminin_II', 'LAMININ_IV_B'],
    'shared_EGF': ['EGF', 'EGF_CA', 'hEGF'],
    'integrin_beta_core': ['Integrin_beta'],
    'integrin_alpha_FG_GAP': ['FG-GAP', 'FG-GAP_2'],
    'integrin_alpha_aux': ['Integrin_alpha', 'Integrin_A_Ig_1', 'Integrin_A_Ig_2', 'Integrin_A_Ig_3'],
    'integrin_beta_aux': ['Integrin_b_cyt', 'Integrin_B_tail', 'I-EGF_1', 'PSI_integrin'],
    'cadherin_core': ['Cadherin', 'Cadherin_2', 'Cadherin_3', 'Cadherin_5', 'Cadherin-like', 'ECD'],
    'cadherin_aux': ['Cadherin_pro', 'Cadherin_tail', 'Cadherin_C_2', 'Protocadherin'],
    'actin': ['Actin'],
    'myosin_head': ['Myosin_head'],
    'troponin': ['Troponin', 'Troponin-I_N'],
    'tropomyosin': ['Tropomyosin'],
}
CORE = {'Collagen', 'C4', 'COLFI', 'Laminin_N', 'Laminin_G_1', 'Laminin_G_2',
        'Laminin_EGF', 'Integrin_beta', 'FG-GAP', 'FG-GAP_2', 'Cadherin',
        'Cadherin_2', 'Troponin', 'Tropomyosin'}
MODEL_GROUP = {name: group for group, names in GROUPS.items() for name in names}
SPECIES = {'C_gigas': 'Crassostrea_gigas', 'C_teleta': 'Capitella_teleta',
           'D_japonicum': 'Dicyema_japonicum', 'H_robusta': 'Helobdella_robusta',
           'I_linei': 'Intoshia_linei', 'L_gigantea': 'Lottia_gigantea',
           'O_sinensis': 'Octopus_sinensis', 'O_vulgaris': 'Octopus_vulgaris',
           'P_canaliculata': 'Pomacea_canaliculata', 'S_mansoni': 'Schistosoma_mansoni'}

def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as handle:
        for chunk in iter(lambda: handle.read(4 * 1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()

def save_json(name, value):
    (OUT / name).write_text(json.dumps(value, indent=2) + '\n')

def table(name, fields, rows):
    with (OUT / name).open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t', extrasaction='ignore')
        writer.writeheader()
        writer.writerows(rows)

def fasta(path):
    ident, pieces = None, []
    with path.open() as handle:
        for line in handle:
            if line.startswith('>'):
                if ident is not None:
                    yield ident, ''.join(pieces)
                ident, pieces = line[1:].split()[0], []
            else:
                pieces.append(line.strip())
    if ident is not None:
        yield ident, ''.join(pieces)

def models():
    metadata, blocks, selected = [], [], {}
    digest = hashlib.sha256()
    count = 0
    with DB.open('rb') as handle:
        for line in handle:
            digest.update(line)
            blocks.append(line)
            if line.strip() != b'//':
                continue
            count += 1
            header = {}
            for raw in blocks:
                fields = raw.decode().strip().split(maxsplit=1)
                if fields and fields[0] in {'NAME', 'ACC', 'DESC', 'LENG', 'GA', 'DATE'}:
                    header[fields[0]] = fields[1]
                if raw.startswith(b'HMM '):
                    break
            if header['NAME'] in MODEL_GROUP:
                name = header['NAME']
                header['group'] = MODEL_GROUP[name]
                block = b''.join(blocks)
                header['sha256'] = hashlib.sha256(block).hexdigest()
                selected[name] = block
                metadata.append(header)
            blocks = []
    missing = set(MODEL_GROUP) - set(selected)
    if missing:
        raise RuntimeError('Requested model NAMES absent in database: ' + repr(missing))
    for filename, names in [('target_profiles.hmm', set(selected)), ('core_profiles.hmm', CORE)]:
        (OUT / filename).write_bytes(b''.join(selected[name] for name in sorted(names)))
    table('models.tsv', ['group', 'NAME', 'ACC', 'DESC', 'LENG', 'GA', 'DATE', 'sha256'], metadata)
    return metadata, {'path': str(DB), 'sha256': digest.hexdigest(), 'bytes': DB.stat().st_size,
                      'model_count': count, 'release': 'Not explicitly documented locally; do not infer release from count',
                      'model_metadata': 'Exact accession versions and GA scores recorded in models.tsv'}

def prepare_inputs():
    datasets, normalization = [], []
    paths = [('Dicyema_unfiltered', 'Dicyema_japonicum', 'unfiltered',
              BASE / 'Dicyema_japonicum/Dicyema_japonicum.faa'),
             ('Dicyema_strict_original', 'Dicyema_japonicum', 'strict_original_ids',
              BASE / 'Dicyema_japonicum/Dicyema_japonicum.strict.complete.faa')]
    paths += [(path.stem, SPECIES[path.stem], 'OrthoFinder_strict', path)
              for path in sorted((BASE / '05.OrthoFinder_input').glob('*.fa'))]
    if len(paths) != 12:
        raise RuntimeError('Expected full + original strict Dicyema and ten OrthoFinder datasets')
    (OUT / 'inputs').mkdir(exist_ok=True)
    for name, species, role, source in paths:
        dest = OUT / 'inputs' / (name + '.faa')
        ids, residues = set(), 0
        with dest.open('w') as handle:
            for ident, seq in fasta(source):
                if ident in ids:
                    raise RuntimeError('Duplicate FASTA identifier: ' + ident)
                ids.add(ident)
                stops = seq.count('*')
                clean = seq.rstrip('*').replace('*', 'X').replace('-', 'X').replace('.', 'X').upper()
                if not clean:
                    raise RuntimeError('Empty protein ' + ident)
                if stops or clean != seq:
                    normalization.append({'dataset': name, 'id': ident, 'raw_length': len(seq),
                                          'scan_length': len(clean), 'stop_count': stops,
                                          'rule': 'strip terminal stops; internal stops/gaps -> X; uppercase'})
                handle.write('>' + ident + '\n' + clean + '\n')
                residues += len(clean)
        datasets.append({'dataset': name, 'species': species, 'role': role, 'source': str(source),
                         'source_sha256': sha(source), 'scan_file': str(dest), 'scan_sha256': sha(dest),
                         'proteins': len(ids), 'residues': residues})
        print('prepared', name, len(ids), flush=True)
    table('input_normalization.tsv', ['dataset', 'id', 'raw_length', 'scan_length', 'stop_count', 'rule'], normalization)
    table('inputs.tsv', list(datasets[0]), datasets)
    return datasets

def scan(dataset, mode, profile='target_profiles.hmm', exhaustive=False):
    name = dataset['dataset'] + '.' + mode
    if exhaustive:
        name += '.max'
    prefix = OUT / 'raw' / name
    command = [str(BIN), '--cpu', '4', '--noali', '--seed', '42', '--domtblout', str(prefix) + '.domtblout',
               '--tblout', str(prefix) + '.tblout']
    if mode == 'GA':
        command += ['--cut_ga']
    else:
        command += ['-E', '1e-5', '--domE', '1e-5', '--incE', '1e-5', '--incdomE', '1e-5']
    if exhaustive:
        command += ['--max']
    command += [str(OUT / profile), dataset['scan_file']]
    with (OUT / 'commands.jsonl').open('a') as handle:
        handle.write(json.dumps({'dataset': dataset['dataset'], 'mode': mode, 'exhaustive': exhaustive,
                                 'argv': command, 'shell': shlex.join(command), 'started_utc': datetime.datetime.utcnow().isoformat()}) + '\n')
    print('START', name, flush=True)
    start = time.time()
    with Path(str(prefix) + '.log').open('w') as log:
        subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, check=True)
    if not Path(str(prefix) + '.domtblout').read_text().rstrip().endswith('# [ok]'):
        raise RuntimeError('Incomplete HMMER domtblout ' + name)
    rows = parse_dom(Path(str(prefix) + '.domtblout'), dataset['dataset'], mode, 'max' if exhaustive else 'default')
    print('DONE', name, 'domains', len(rows), 'seconds', round(time.time() - start, 2), flush=True)
    return rows

def parse_dom(path, dataset, mode, heuristic, hmmscan=False):
    rows = []
    for line in path.read_text().splitlines():
        if not line or line.startswith('#'):
            continue
        p = line.split(maxsplit=22)
        protein, name, acc, plen, mlen = (p[3], p[0], p[1], p[5], p[2]) if hmmscan else (p[0], p[3], p[4], p[2], p[5])
        if name not in MODEL_GROUP:
            continue
        if mode == 'E1e-5' and (float(p[6]) > 1e-5 or float(p[12]) > 1e-5):
            continue
        rows.append({'dataset': dataset, 'mode': mode, 'heuristic': heuristic, 'protein_id': protein,
                     'model': name, 'accession': acc, 'group': MODEL_GROUP[name], 'protein_length': int(plen),
                     'model_length': int(mlen), 'sequence_evalue': float(p[6]), 'sequence_score': float(p[7]),
                     'domain_conditional_evalue': float(p[11]), 'domain_independent_evalue': float(p[12]),
                     'domain_score': float(p[13]), 'domain_bias': float(p[14]), 'hmm_from': int(p[15]),
                     'hmm_to': int(p[16]), 'ali_from': int(p[17]), 'ali_to': int(p[18]),
                     'model_coverage': (int(p[16]) - int(p[15]) + 1) / int(mlen), 'source_domtblout': str(path)})
    return rows

def summarize(datasets, hits):
    fields = list(hits[0])
    table('hits.tsv', fields, hits)
    counters = []
    for ds in datasets:
        for mode in ['GA', 'E1e-5']:
            for group in GROUPS:
                subset = [h for h in hits if h['dataset'] == ds['dataset'] and h['mode'] == mode and h['group'] == group]
                # Union default and --max evidence, deduplicating domain coordinates.
                domains = {(h['protein_id'], h['model'], h['ali_from'], h['ali_to']) for h in subset}
                counters.append({'dataset': ds['dataset'], 'species': ds['species'], 'role': ds['role'],
                                 'proteins_scanned': ds['proteins'], 'mode': mode, 'group': group,
                                 'hit_proteins': len({h['protein_id'] for h in subset}), 'hit_domains': len(domains)})
    table('counts_long.tsv', list(counters[0]), counters)
    wide = []
    for ds in datasets:
        row = {'dataset': ds['dataset'], 'species': ds['species'], 'proteins_scanned': ds['proteins']}
        for mode in ['GA', 'E1e-5']:
            for group in GROUPS:
                row[group + '_' + mode] = next(c['hit_proteins'] for c in counters if c['dataset'] == ds['dataset'] and c['mode'] == mode and c['group'] == group)
        wide.append(row)
    table('species_counts.tsv', list(wide[0]), wide)
    strict = dict(fasta(OUT / 'inputs/Dicyema_strict_original.faa'))
    full = dict(fasta(OUT / 'inputs/Dicyema_unfiltered.faa'))
    seq_to_strict = collections.defaultdict(list)
    for ident, seq in strict.items():
        seq_to_strict[hashlib.sha256(seq.encode()).hexdigest()].append(ident)
    filter_report = BASE / 'Dicyema_japonicum/Dicyema_japonicum.strict_cds_filter.report.tsv'
    reports = {r['id']: r for r in csv.DictReader(filter_report.open(), delimiter='\t')}
    rescues, concordance = [], []
    for mode in ['GA', 'E1e-5']:
        for group in GROUPS:
            raw_hits = {h['protein_id'] for h in hits if h['dataset'] == 'Dicyema_unfiltered' and h['mode'] == mode and h['group'] == group}
            strict_hits = {h['protein_id'] for h in hits if h['dataset'] == 'Dicyema_strict_original' and h['mode'] == mode and h['group'] == group}
            matched = set()
            for ident in sorted(raw_hits):
                mapping = seq_to_strict.get(hashlib.sha256(full[ident].encode()).hexdigest(), [])
                matched.update(mapping)
                if not mapping:
                    report = reports.get(ident, {})
                    rescues.append({'mode': mode, 'group': group, 'unfiltered_id': ident,
                                    'filter_status': report.get('status', 'ID_not_in_report'),
                                    'filter_reasons': report.get('reasons', ''), 'scan_length': len(full[ident]),
                                    'evidence': 'protein sequence absent from strict set; not necessarily a new gene'})
            concordance.append({'mode': mode, 'group': group, 'full_hits': len(raw_hits), 'strict_hits': len(strict_hits),
                                'full_hits_matching_strict_sequence': len(matched),
                                'strict_hits_not_in_full_exact_sequence': ','.join(sorted(strict_hits - matched))})
    table('filtered_rescue_ids.tsv', ['mode', 'group', 'unfiltered_id', 'filter_status', 'filter_reasons', 'scan_length', 'evidence'], rescues)
    table('strict_full_concordance.tsv', list(concordance[0]), concordance)
    # Rename mapping is checked by exact normalized protein sequence, not guessed ID equivalence.
    renamed = dict(fasta(OUT / 'inputs/D_japonicum.faa'))
    renamed_hashes = collections.defaultdict(list)
    for ident, seq in renamed.items():
        renamed_hashes[hashlib.sha256(seq.encode()).hexdigest()].append(ident)
    mapping = [{'original_strict_id': ident, 'renamed_ids': ','.join(renamed_hashes.get(hashlib.sha256(seq.encode()).hexdigest(), []))}
               for ident, seq in strict.items()]
    table('Dicyema_strict_sequence_id_mapping.tsv', list(mapping[0]), mapping)
    legacy_path = BASE / '08.ECM_Muscle_annotation/pfam/ECM_Muscle.with_Dj.pfam.domtblout'
    legacy = parse_dom(legacy_path, 'legacy_selected_OG', 'GA', 'legacy_hmmscan', hmmscan=True)
    legacy = [h for h in legacy if h['protein_id'].startswith('D_japonicum')]
    table('legacy_Dicyema_reparsed.tsv', fields, legacy)
    families = list(csv.DictReader((BASE / '12.family_groups/OG_family_summary.tsv').open(), delimiter='\t'))
    og_rows = list(csv.DictReader((BASE / '10.family_analysis/selected_orthogroups.info.tsv').open(), delimiter='\t'))
    og_genes = {r['Orthogroup']: r['D_japonicum_genes'].split(',') for r in og_rows}
    (OUT / 'source_OG_family_summary.tsv').write_bytes((BASE / '12.family_groups/OG_family_summary.tsv').read_bytes())
    (OUT / 'source_selected_orthogroups.info.tsv').write_bytes((BASE / '10.family_analysis/selected_orthogroups.info.tsv').read_bytes())
    (OUT / 'source_Dicyema_filter_report.tsv').write_bytes(filter_report.read_bytes())
    audits = []
    # Family column name is obtained from the actual source header.
    print('family source columns', list(families[0]), flush=True)
    for fam in families:
        if not any(v in {'Cadherin', 'Integrin'} for v in fam.values()):
            continue
        og = fam.get('Orthogroup', next(iter(fam.values())))
        genes = set(og_genes.get(og, []))
        dh = [h for h in hits if h['dataset'] == 'D_japonicum' and h['protein_id'] in genes]
        lh = [h for h in legacy if h['protein_id'] in genes]
        for mode in ['GA', 'E1e-5']:
            hs = [h for h in dh if h['mode'] == mode]
            core = [h for h in hs if h['group'] in {'cadherin_core', 'integrin_beta_core', 'integrin_alpha_FG_GAP'}]
            audits.append({'Orthogroup': og, 'legacy_family': 'Cadherin' if 'Cadherin' in fam.values() else 'Integrin',
                           'mode': mode, 'Dicyema_genes': ','.join(sorted(genes)), 'new_models': ','.join(sorted({h['model'] for h in hs})),
                           'core_supported_ids': ','.join(sorted({h['protein_id'] for h in core})),
                           'legacy_models': ','.join(sorted({h['model'] for h in lh})),
                           'legacy_min_domain_iE': min((h['domain_independent_evalue'] for h in lh), default=''),
                           'new_min_domain_iE': min((h['domain_independent_evalue'] for h in hs), default=''),
                           'interpretation': 'core profile supported' if core else ('auxiliary/shared profiles only' if hs else 'no tested profile support; motifs/OG membership alone are not domain validation')})
    table('adhesion_OG_audit.tsv', list(audits[0]), audits)
    return counters

def main():
    if OUT == BASE or BASE in OUT.parents:
        raise RuntimeError('Keep validation outputs outside the input tree')
    OUT.mkdir(parents=True, exist_ok=True)
    if not DB.is_file() or not BIN.is_file():
        raise RuntimeError('Set ECM_PFAM_HMM and install hmmsearch or set ECM_HMMSEARCH')
    (OUT / 'raw').mkdir(exist_ok=True)
    if (OUT / 'commands.jsonl').exists():
        raise RuntimeError('Existing run detected; do not overwrite provenance')
    version = subprocess.check_output([str(BIN), '-h'], text=True)
    (OUT / 'hmmsearch_version.txt').write_text(version)
    metadata, db_info = models()
    datasets = prepare_inputs()
    provenance = {'started_utc': datetime.datetime.utcnow().isoformat(), 'host': platform.node(),
                  'python': platform.python_version(), 'hmmsearch': str(BIN), 'hmmsearch_sha256': sha(BIN),
                  'database': db_info, 'datasets': datasets, 'models': metadata,
                  'thresholds': {'GA': 'HMMER --cut_ga (per sequence and domain model-specific bit scores)',
                                 'E1e-5': 'sequence E <= 1e-5 AND independent domain E <= 1e-5; conditional E is not used as final evidence'},
                  'cpu': 4, 'sensitivity': 'Default all datasets; --max GA and E1e-5 for 14 core models in full and original strict Dicyema',
                  'normalization': 'Terminal stop symbols stripped; internal stop/gap symbols replaced with X; no length or quality filtering',
                  'scope': 'Profile non-detection in supplied predicted proteins is not definitive genomic loss; shared domains do not establish full family identity; no clan overlap competition'}
    save_json('provenance.json', provenance)
    hits = []
    for ds in datasets:
        for mode in ['GA', 'E1e-5']:
            hits += scan(ds, mode)
    for ds in datasets[:2]:
        for mode in ['GA', 'E1e-5']:
            hits += scan(ds, mode, 'core_profiles.hmm', True)
    counts = summarize(datasets, hits)
    provenance['finished_utc'] = datetime.datetime.utcnow().isoformat()
    provenance['completed_searches'] = 28
    save_json('provenance.json', provenance)
    manifest = [{'file': str(p.relative_to(OUT)), 'bytes': p.stat().st_size, 'sha256': sha(p)}
                for p in sorted(OUT.rglob('*')) if p.is_file() and p.name not in {'run.log', 'checksums.tsv', 'SUCCESS.json'}]
    table('checksums.tsv', ['file', 'bytes', 'sha256'], manifest)
    save_json('SUCCESS.json', {'completed_searches': 28, 'datasets': len(datasets), 'models': len(metadata),
                              'domain_rows': len(hits), 'count_rows': len(counts), 'files_hashed': len(manifest)})
    print('SUCCESS', json.dumps(json.loads((OUT / 'SUCCESS.json').read_text())), flush=True)

if __name__ == '__main__':
    main()
