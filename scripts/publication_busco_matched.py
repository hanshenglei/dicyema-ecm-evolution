"""Run matched protein-mode BUSCO benchmarks without altering original inputs."""
import argparse
import hashlib
import json
import os
import subprocess
import sys
import tarfile
import time
import urllib.request
from pathlib import Path


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(2**20), b''):
            h.update(block)
    return h.hexdigest()


def retrieve(url, dest):
    if not dest.exists():
        print('Downloading', url, flush=True)
        request = urllib.request.Request(url, headers={'User-Agent': 'Dicyema-annotation-benchmark/1.0'})
        with urllib.request.urlopen(request, timeout=90) as response, dest.open('wb') as f:
            while block := response.read(2**20):
                f.write(block)
    return {'url': url, 'sha256': sha(dest), 'bytes': dest.stat().st_size}


def extract(archive, destination):
    root = destination.resolve()
    with tarfile.open(archive) as tf:
        for m in tf.getmembers():
            target = (root/m.name).resolve()
            if not target.is_relative_to(root) or m.issym() or m.islnk():
                raise ValueError('Unsafe archive member '+m.name)
        tf.extractall(root)


def run(command, cwd, env, log):
    print('Running', ' '.join(command[:5]), flush=True)
    with log.open('w') as out:
        p = subprocess.run(command, cwd=cwd, env=env, stdout=out, stderr=subprocess.STDOUT, timeout=7200)
    if p.returncode:
        print(log.read_text()[-6000:], flush=True)
        raise RuntimeError('Command failed with status '+str(p.returncode))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--project-root', required=True, type=Path)
    parser.add_argument('--output', required=True, type=Path)
    parser.add_argument('--hmmer-bin', required=True, type=Path)
    parser.add_argument('--cpu', type=int, default=8)
    args = parser.parse_args()
    out = args.output.resolve()
    out.mkdir(parents=True, exist_ok=True)
    software = out/'software'
    software.mkdir(exist_ok=True)
    source = software/'busco-5.8.2.tar.gz'
    dataset = software/'eukaryota_odb10.2024-01-08.tar.gz'
    provenance = {'started_utc': time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()), 'python': sys.version, 'software': {}, 'inputs': [], 'commands': []}
    provenance['software']['BUSCO'] = retrieve('https://gitlab.com/ezlab/busco/-/archive/5.8.2/busco-5.8.2.tar.gz', source)
    provenance['software']['lineage'] = retrieve('https://busco-data.ezlab.org/v5/data/lineages/eukaryota_odb10.2024-01-08.tar.gz', dataset)
    if not (software/'busco-5.8.2').exists():
        extract(source, software)
    if not (software/'eukaryota_odb10').exists():
        extract(dataset, software)
    busco_source = software/'busco-5.8.2'
    env = os.environ.copy()
    env['PATH'] = str(args.hmmer_bin)+os.pathsep+env.get('PATH', '')
    env['PYTHONPATH'] = str(busco_source/'src')+os.pathsep+env.get('PYTHONPATH', '')
    env['OMP_NUM_THREADS'] = '1'
    version = subprocess.run([sys.executable, str(busco_source/'bin/busco'), '--version'], env=env, capture_output=True, text=True)
    provenance['busco_version'] = version.stdout.strip()
    if version.returncode:
        raise RuntimeError(version.stderr)
    provenance['hmmer'] = {'binary': str(args.hmmer_bin/'hmmsearch'), 'sha256': sha(args.hmmer_bin/'hmmsearch')}
    lineage = software/'eukaryota_odb10'
    provenance['dataset_cfg'] = (lineage/'dataset.cfg').read_text()
    provenance['dataset_cfg_sha256'] = sha(lineage/'dataset.cfg')
    input_paths = {
        'Published_5012': '00.genome/dicyema_japonicum_genemodel_prot.fa',
        'Integrated_17260': '14.ecm/Dicyema_japonicum/Dicyema_japonicum.faa',
        'Strict_12402': '14.ecm/Dicyema_japonicum/Dicyema_japonicum.strict.complete.faa',
    }
    if (out/'provenance.json').exists():
        previous=json.loads((out/'provenance.json').read_text())
        recorded={r['name']:r['sha256'] for r in previous.get('inputs',[])}
        for name, relative in input_paths.items():
            if name in recorded and sha(args.project_root/relative)!=recorded[name]:
                raise ValueError('Input changed for existing benchmark; use a new output directory')
        if previous.get('software')!=provenance['software']:
            raise ValueError('Software or lineage changed; use a new output directory')
    for name, relative in input_paths.items():
        path = args.project_root/relative
        models = sum(line.startswith(b'>') for line in path.open('rb'))
        provenance['inputs'].append({'name': name, 'path': str(path), 'sha256': sha(path), 'models': models, 'normalization': 'none; original input bytes passed directly to BUSCO'})
        command = [sys.executable, str(busco_source/'bin/busco'), '-i', str(path), '-o', name, '-m', 'proteins', '-l', str(lineage), '-c', str(args.cpu), '--offline', '--opt-out-run-stats', '--download_path', str(software/'downloads')]
        provenance['commands'].append(command)
        (out/'provenance.json').write_text(json.dumps(provenance, indent=2))
        complete = list((out/name).glob('short_summary.specific.*.json'))
        if not complete:
            run(command, out, env, out/(name+'.console.log'))
        summary = json.loads(next((out/name).glob('short_summary.specific.*.json')).read_text())
        print(name, json.dumps(summary['results']), flush=True)
    provenance['finished_utc'] = time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())
    (out/'provenance.json').write_text(json.dumps(provenance, indent=2))
    import zipfile
    with zipfile.ZipFile(out/'matched_busco_evidence.zip', 'w', compression=zipfile.ZIP_DEFLATED, compresslevel=6) as z:
        for path in out.rglob('*'):
            if path.is_file() and ('software' not in path.relative_to(out).parts) and path.name != z.filename.split('/')[-1]:
                z.write(path, path.relative_to(out).as_posix())
    print('Evidence archive', out/'matched_busco_evidence.zip', flush=True)


if __name__ == '__main__':
    main()
