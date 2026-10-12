"""Add the published protein inventory to the existing controlled Pfam searches."""
import argparse
import csv
import hashlib
import json
import subprocess
import zipfile
from pathlib import Path


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--input', type=Path, required=True)
    p.add_argument('--profiles', type=Path, required=True)
    p.add_argument('--core-profiles', type=Path, required=True)
    p.add_argument('--hmmer', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    args = p.parse_args()
    out = args.output
    out.mkdir(exist_ok=True, parents=True)
    if (out/'provenance.json').exists():
        previous=json.loads((out/'provenance.json').read_text())
        for key,path in [('input_sha256',args.input),('profiles_sha256',args.profiles),('core_profiles_sha256',args.core_profiles),('hmmer_sha256',args.hmmer)]:
            if previous[key]!=sha(path):
                raise ValueError('Inputs or program changed; use a new output directory')
    seqs = {}
    ident = None
    for line in args.input.open():
        if line.startswith('>'):
            ident = line[1:].split()[0]
            if ident in seqs:
                raise ValueError('Duplicate FASTA ID')
            seqs[ident] = ''
        else:
            seqs[ident] += line.strip()
    normalized = out/'Published_5012.normalized.faa'
    with normalized.open('w') as f:
        for ident, seq in seqs.items():
            clean = seq.rstrip('*').replace('*', 'X').replace('-', 'X').replace('.', 'X').upper()
            f.write('>'+ident+'\n'+clean+'\n')
    prov = {'input':str(args.input),'input_sha256':sha(args.input),'normalized_sha256':sha(normalized),'models':len(seqs),'normalization':'strip terminal stops; internal stops/gaps to X; uppercase','hmmer_sha256':sha(args.hmmer),'profiles_sha256':sha(args.profiles),'core_profiles_sha256':sha(args.core_profiles),'commands':[]}
    hits = []
    for exhaustive in [False, True]:
        for mode in ['GA', 'E1e-5']:
            name = mode+('.max' if exhaustive else '')
            dom = out/(name+'.domtblout')
            command = [str(args.hmmer),'--cpu','4','--noali','--seed','42','--domtblout',str(dom),'--tblout',str(out/(name+'.tblout'))]
            command += ['--cut_ga'] if mode == 'GA' else ['-E','1e-5','--domE','1e-5','--incE','1e-5','--incdomE','1e-5']
            if exhaustive:
                command += ['--max']
            command += [str(args.core_profiles if exhaustive else args.profiles),str(normalized)]
            prov['commands'].append(command)
            if not dom.exists() or not dom.read_text().rstrip().endswith('# [ok]'):
                with (out/(name+'.log')).open('w') as log:
                    subprocess.run(command,stdout=log,stderr=subprocess.STDOUT,check=True,timeout=3600)
            for line in dom.read_text().splitlines():
                if not line or line.startswith('#'):
                    continue
                t = line.split(maxsplit=22)
                if mode == 'E1e-5' and (float(t[6]) > 1e-5 or float(t[12]) > 1e-5):
                    continue
                hits.append({'mode':mode,'search':'max' if exhaustive else 'default','protein_id':t[0],'model':t[3],'accession':t[4],'sequence_evalue':t[6],'independent_evalue':t[12],'ali_from':t[17],'ali_to':t[18]})
            print('Completed published protein scan',name,flush=True)
    (out/'provenance.json').write_text(json.dumps(prov,indent=2))
    with (out/'hits.tsv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(hits[0]),delimiter='\t')
        w.writeheader();w.writerows(hits)
    with zipfile.ZipFile(out/'published_profile_evidence.zip','w',compression=zipfile.ZIP_DEFLATED) as z:
        for path in out.iterdir():
            if path.is_file() and path.suffix != '.zip':
                z.write(path,path.name)
    print('Finished published protein-domain validation',flush=True)


if __name__ == '__main__':
    main()
