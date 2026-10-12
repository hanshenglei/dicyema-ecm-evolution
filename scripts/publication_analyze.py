"""Recompute comparative evidence without transferring seed family names to Dicyema."""
from __future__ import annotations
import csv
import json
import os
from pathlib import Path
from collections import Counter, defaultdict
from itertools import product
import numpy as np
from scipy.stats import fisher_exact
from Bio import SeqIO, Phylo
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(os.environ.get('ECM_PROJECT_ROOT', Path(__file__).resolve().parents[1] / 'workspace')).resolve()
CACHE=ROOT/'remote_results_cache'
SOURCE=ROOT/'publication_package/remote_sources'
OUT=ROOT/'publication_package/analysis'
ORTH='05.OrthoFinder_input/OrthoFinder/Results_Jun11'

def read(path):
    with Path(path).open(encoding='utf-8-sig',newline='') as f:
        return list(csv.DictReader(f,delimiter='\t'))

def write(name,rows):
    with (OUT/name).open('w',encoding='utf-8',newline='') as f:
        if rows:
            w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
        else:
            f.write('No_records\n')

def classification(names):
    n=set(names)
    if 'Integrin_beta' in n: return 'Integrin beta related','receptor_domain'
    if 'FG-GAP' in n: return 'FG-GAP containing alpha integrin candidate','shared_or_partial'
    if n & {'Cadherin','Cadherin_2'}: return 'Cadherin domain containing','adhesion_domain'
    if n & {'Reprolysin','Reprolysin_2','Reprolysin_3','Reprolysin_4','Reprolysin_5'}:
        return 'Reprolysin disintegrin thrombospondin modules','proteolysis_domain'
    if n & {'TSP_1','TSP1_ADAMTS'}: return 'Thrombospondin repeat containing','extracellular_associated'
    if n & {'CUB','Astacin'}: return 'CUB astacin containing','extracellular_associated'
    if 'VWD' in n: return 'VWD EGF TIL modular proteins','extracellular_associated'
    if 'Lectin_C' in n: return 'C type lectin domain containing','extracellular_associated'
    if 'Sushi' in n: return 'Sushi repeat containing','extracellular_associated'
    if 'VWC' in n: return 'VWC domain containing','extracellular_associated'
    if 'VIT' in n: return 'VIT domain containing','extracellular_associated'
    if n & {'TIG','G8'}: return 'TIG G8 modular proteins','shared_or_partial'
    if 'VWA_3' in n: return 'VWA related modules','shared_or_partial'
    if n & {'EGF_CA','EGF','Laminin_EGF'}: return 'EGF repeat containing','shared_or_partial'
    if 'Kazal_2' in n: return 'Kazal domain containing','shared_or_partial'
    if n & {'FERM_N','FERM_M','FERM_C'}: return 'FERM cytoplasmic adhesion regulators','intracellular_associated'
    if n & {'GPS','7tm_2'}: return 'GPS seven transmembrane receptor modules','shared_or_partial'
    if n & {'SMC_N','SMC_hinge','MCM_bind','Homeodomain','GCN5L1','Paf1','RRS1','Mt_ATP-synt_B','PP2C','AAA'}:
        return 'Intracellular domain supported proteins','intracellular_or_nonECM'
    if n: return 'Shared domains without specific ECM assignment','shared_or_partial'
    return 'Sequence candidate without Dicyema Pfam support','unresolved'

def main():
    OUT.mkdir(parents=True, exist_ok=True); (OUT/'figures').mkdir(exist_ok=True)
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42})
    ecm=read(CACHE/'09.presence_absence/ECM.orthogroup.presence_absence.tsv')
    muscle=read(CACHE/'09.presence_absence/Muscle.orthogroup.presence_absence.tsv')
    es={r['Orthogroup'] for r in ecm};ms={r['Orthogroup'] for r in muscle}
    counts=read(CACHE/ORTH/'Orthogroups/Orthogroups.GeneCount.tsv')
    cm={r['Orthogroup']:r for r in counts};species=[x for x in counts[0] if x not in ['Orthogroup','Total']]
    families={r['Orthogroup']:r for r in read(CACHE/'12.family_groups/OG_family_summary.tsv')}
    domains=defaultdict(list)
    for r in read(CACHE/'08.ECM_Muscle_annotation/pfam/ECM_Muscle.with_Dj.pfam.domain.tsv'):
        domains[r['gene_id']].append(r)
    independent_evalues={}
    raw_domtbl=SOURCE/'08.ECM_Muscle_annotation/pfam/ECM_Muscle.with_Dj.pfam.domtblout'
    for line in raw_domtbl.open(encoding='utf-8'):
        if line.startswith('#'):continue
        fields=line.split(maxsplit=22)
        if len(fields)>=22:
            independent_evalues[(fields[3],fields[0],fields[17],fields[18])]=fields[12]
    proteins={r.id:str(r.seq) for r in SeqIO.parse(SOURCE/'05.OrthoFinder_input/D_japonicum.fa','fasta')}
    pairs=read(CACHE/ORTH/'Orthologues/D_japonicum.tsv')
    partner=defaultdict(list)
    for r in pairs:
        if r['Orthogroup'] in es:
            partner[r['Orthogroup']].append(r)
    dup=defaultdict(list)
    for r in read(CACHE/ORTH/'Gene_Duplication_Events/Duplications.tsv'):
        if r['Orthogroup'] in es:
            dup[r['Orthogroup']].append(r)
    coord=read(ROOT/'publication_package/coordinate_validation/ECM_gene_coordinates.tsv')
    co={r['analysis_gene_id']:r for r in coord}
    comp=read(ROOT/'publication_package/coordinate_validation/ECM_overlap_components.tsv')
    components=defaultdict(list)
    component_by_gene={}
    for r in comp:
        components[r['Orthogroup']].append(r)
        for gene in r['member_models'].split(','):component_by_gene[gene]=r['component']
    maps=defaultdict(list)
    for line in (CACHE/'07.ECM_Muscle_candidates/D_japonicum.ECM.related.og_map.tsv').read_text(encoding='utf-8-sig').splitlines():
        if line and not line.startswith('Orthogroup'):
            og,gene=line.split('\t')[:2];maps[og].append(gene)
    gene_rows=[];og_rows=[];pair_rows=[];event_rows=[]
    for entry in ecm:
        if entry['has_D_japonicum']!='yes': continue
        og=entry['Orthogroup'];genes=maps[og];family=families.get(og,{})
        selected=SOURCE/'10.family_analysis/selected_OGs'/og
        motiffile=CACHE/f'10.family_analysis/selected_OGs/{og}/{og}.motif.tsv'
        motifs=defaultdict(list)
        if motiffile.exists():
            for r in read(motiffile): motifs[r['gene_id']].append(r)
        allnames=[];domain_genes=[]
        for gene in genes:
            doms=sorted(domains[gene],key=lambda x:(int(x['ali_start']),int(x['ali_end'])))
            names=[r['domain_name'] for r in doms];allnames+=names
            if names: domain_genes.append(gene)
            module,status=classification(names)
            gene_rows.append(dict(Orthogroup=og,gene_id=gene,protein_length=len(proteins.get(gene,'')),
                original_family=family.get('assigned_family','unplotted'),original_quality=family.get('quality_class','unplotted'),
                domain_architecture=';'.join(f"{r['domain_name']}:{r['ali_start']}-{r['ali_end']}" for r in doms),
                domain_names=','.join(sorted(set(names))),domain_hits=len(doms),
                domain_independent_evalues=';'.join(independent_evalues.get((gene,r['domain_name'],r['ali_start'],r['ali_end']),'NA') for r in doms),
                domain_acceptance='original_hmmscan_model_gathering_thresholds',
                motif_unique_types=len({r['motif'] for r in motifs[gene]}),motif_hits=len(motifs[gene]),
                motif_table_available=motiffile.exists(),
                motif_architecture=';'.join(f"{r['motif']}:{r['start']}-{r['end']}" for r in motifs[gene]),
                evidence_module=module,evidence_status=status,contig=co[gene]['contig'],
                genomic_start=co[gene]['start'],genomic_end=co[gene]['end'],strand=co[gene]['strand']))
        module,status=classification(allnames)
        terminal=[r for r in dup[og] if r['Species Tree Node']=='D_japonicum']
        internal=[r for r in dup[og] if r['Species Tree Node'].startswith('N')]
        disjoint_terminal=0;overlap_sensitive_terminal=0
        for r in dup[og]:
            component_sets=[]
            for side in ['Genes 1','Genes 2']:
                ids=[x.strip().replace('D_japonicum_D_japonicum','D_japonicum') for x in r[side].split(',') if 'D_japonicum' in x]
                component_sets.append({component_by_gene[k] for k in ids if k in component_by_gene})
            overlap_sensitive=bool(component_sets[0]&component_sets[1])
            disjoint=bool(component_sets[0] and component_sets[1] and not overlap_sensitive)
            if r['Species Tree Node']=='D_japonicum':
                disjoint_terminal+=disjoint;overlap_sensitive_terminal+=overlap_sensitive
            event_rows.append({**r,'event_scope':'Dicyema_terminal' if r['Species Tree Node']=='D_japonicum' else ('internal_species_node' if r['Species Tree Node'].startswith('N') else 'other_species_terminal'),
                'same_overlap_component_in_both_descendants':overlap_sensitive,
                'Dicyema_descendants_in_disjoint_overlap_components':disjoint})
        ortho_species=set();ortho_genes=set();actualpairs=set()
        for r in partner[og]:
            ortho_species.add(r['Species'])
            dj=[x.strip() for x in r['D_japonicum'].split(',') if x.strip()]
            others=[x.strip() for x in r['Orthologs'].split(',') if x.strip()]
            ortho_genes.update(dj)
            for a,b in product(dj,others):
                actualpairs.add((r['Species'],a,b))
        for sp,a,b in sorted(actualpairs):
            pair_rows.append(dict(Orthogroup=og,Dicyema_gene=a,partner_species=sp,partner_gene=b,
                                 relationship='OrthoFinder_inferred_orthology',source='Orthologues/D_japonicum.tsv'))
        treepath=selected/f'{og}.mafft.trim.fa.treefile'
        tree=Phylo.read(treepath,'newick') if treepath.exists() else None
        tree_tips=len(tree.get_terminals()) if tree else 0
        support=[]
        if tree:
            for cl in tree.get_nonterminals():
                if cl.name and '/' in cl.name:
                    try: support.append(tuple(float(x) for x in cl.name.split('/')))
                    except ValueError: pass
        og_rows.append(dict(Orthogroup=og,original_family=family.get('assigned_family','unplotted'),
            original_quality=family.get('quality_class','unplotted'),corrected_module=module,evidence_status=status,
            Dicyema_gene_models=len(genes),overlap_components_50pct=len(components[og]),
            species_present=entry['total_species_present'],total_models=entry['total_gene_number'],
            Dicyema_domain_supported_models=len(domain_genes),Dicyema_domain_names=','.join(sorted(set(allnames))),
            Dicyema_genes=','.join(genes),orthology_partner_species=','.join(sorted(ortho_species)),
            orthology_Dicyema_models=len(ortho_genes),orthology_pair_count=len(actualpairs),
            Dicyema_terminal_duplications=len(terminal),internal_node_duplications=len(internal),
            Dicyema_terminal_overlap_sensitive_events=overlap_sensitive_terminal,
            Dicyema_terminal_disjoint_component_events=disjoint_terminal,
            other_terminal_duplications=len(dup[og])-len(terminal)-len(internal),
            IQTREE_available=bool(tree),IQTREE_tips=tree_tips,
            IQTREE_supported_branches_80_95=sum(a>=80 and b>=95 for a,b in support)))
    write('gene_evidence.tsv',gene_rows);write('OG_evidence.tsv',og_rows)
    write('orthology_pairs.tsv',pair_rows);write('duplication_events.tsv',event_rows)
    retention=[]
    for label,ogs in [('ECM_seed',es),('Muscle_seed',ms),('ECM_only',es-ms),('Muscle_only',ms-es),('Shared_seed',es&ms),('All_orthogroups',set(cm))]:
        for sp in species:
            n=sum(int(cm[og][sp])>0 for og in ogs)
            retention.append(dict(category=label,species=sp,denominator=len(ogs),present_OGs=n,representation_percent=n/len(ogs)*100 if ogs else 0))
    write('retention.tsv',retention)
    er=[x for x in retention if x['category']=='ECM_only' and x['species']=='D_japonicum'][0]
    mr=[x for x in retention if x['category']=='Muscle_only' and x['species']=='D_japonicum'][0]
    a,b=er['present_OGs'],er['denominator']-er['present_OGs'];c,d=mr['present_OGs'],mr['denominator']-mr['present_OGs']
    odds,pvalue=fisher_exact([[a,b],[c,d]])
    logse=np.sqrt(1/a+1/b+1/c+1/d);ci=np.exp(np.log(odds)+np.array([-1,1])*1.96*logse)
    supported_cad=[r for r in gene_rows if set(r['domain_names'].split(','))&{'Cadherin','Cadherin_2'}]
    summary=dict(ECM_OGs=len(es),Muscle_OGs=len(ms),shared_candidate_OGs=len(es&ms),
        ECM_retained_OGs=len(og_rows),ECM_gene_models=len(gene_rows),
        ECM_overlap_components=sum(int(r['overlap_components_50pct']) for r in og_rows),
        ECM_quality_counts=dict(Counter(r['original_quality'] for r in og_rows)),
        ECM_domain_status_counts=dict(Counter(r['evidence_status'] for r in og_rows)),
        ECM_domain_supported_gene_models=sum(r['domain_hits']>0 for r in gene_rows),
        cadherin_domain_OGs=sorted({r['Orthogroup'] for r in supported_cad}),cadherin_domain_models=len(supported_cad),
        ECM_with_orthology_OGs=sum(int(r['orthology_pair_count'])>0 for r in og_rows),
        ECM_orthology_pairs=len(pair_rows),ECM_Dicyema_terminal_duplication_events=sum(int(r['Dicyema_terminal_duplications']) for r in og_rows),
        ECM_Dicyema_overlap_sensitive_terminal_events=sum(int(r['Dicyema_terminal_overlap_sensitive_events']) for r in og_rows),
        ECM_Dicyema_disjoint_component_terminal_events=sum(int(r['Dicyema_terminal_disjoint_component_events']) for r in og_rows),
        ECM_IQTREE_count=sum(r['IQTREE_available'] for r in og_rows),
        fisher_disjoint_candidates=dict(table=[[a,b],[c,d]],odds_ratio=float(odds),pvalue=float(pvalue),
            approximate_Wald_95pct_CI=ci.tolist(),scope='conditional_on_keyword_seed_selection_not_phylogenetically_independent'),
        primary_interpretation='Residual domain-defined modules; raw candidate representation is not a direct ECM function or gene-loss rate.')
    (OUT/'summary.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')
    figures(og_rows,gene_rows,retention)
    print(json.dumps(summary,indent=2))

def save(fig,name):
    fig.savefig(OUT/'figures'/f'{name}.png',dpi=240,bbox_inches='tight')
    fig.savefig(OUT/'figures'/f'{name}.pdf',bbox_inches='tight');plt.close(fig)

def figures(ogs,genes,retention):
    species=[r['species'] for r in retention if r['category']=='ECM_seed']
    fig,axs=plt.subplots(1,2,figsize=(10,4),gridspec_kw={'width_ratios':[1.3,1]})
    y=np.arange(len(species));labels=[s.replace('_',' ') for s in species]
    for offset,cat,col in [(-.17,'ECM_seed','#007A87'),(.17,'Muscle_seed','#A3486C')]:
        v=[next(r['representation_percent'] for r in retention if r['category']==cat and r['species']==s) for s in species]
        axs[0].barh(y+offset,v,height=.32,color=col,label=cat.replace('_',' '))
    axs[0].set_yticks(y,labels,fontsize=8);axs[0].invert_yaxis();axs[0].set_xlabel('Candidate orthogroup representation (%)');axs[0].legend(fontsize=8)
    cats=Counter(r['evidence_status'] for r in ogs)
    items=list(cats.items());axs[1].barh(range(len(items)),[v for k,v in items],color='#578557')
    axs[1].set_yticks(range(len(items)),[k.replace('_',' ') for k,v in items],fontsize=8);axs[1].set_xlabel('Retained ECM seed orthogroups')
    fig.tight_layout();save(fig,'Figure1_representation_and_domain_evidence')
    fig,axs=plt.subplots(1,2,figsize=(9,3.6))
    labels=['Raw models','25% overlap','50% overlap','80% overlap'];values=[137,100,101,105]
    axs[0].bar(range(4),values,color=['#777777','#007A87','#007A87','#007A87']);axs[0].set_xticks(range(4),labels,rotation=25,ha='right');axs[0].set_ylabel('ECM models or overlap components')
    for i,v in enumerate(values):axs[0].text(i,v+2,str(v),ha='center')
    raw=[int(r['Dicyema_gene_models']) for r in ogs];collapsed=[int(r['overlap_components_50pct']) for r in ogs]
    axs[1].scatter(raw,collapsed,color='#A3486C',s=26);axs[1].plot([0,18],[0,18],color='#999999',linestyle='--');axs[1].set_xlabel('Gene models per orthogroup');axs[1].set_ylabel('50% overlap components');axs[1].set_xlim(0,18);axs[1].set_ylim(0,18)
    fig.tight_layout();save(fig,'Figure2_coordinate_sensitivity')
    domain_names=sorted({name for g in genes for name in g['domain_names'].split(',') if name})
    palette={name:plt.get_cmap('tab20')(i%20) for i,name in enumerate(domain_names)}
    for og in ogs:
        selected=[g for g in genes if g['Orthogroup']==og['Orthogroup']]
        height=max(2.8,.6*len(selected)+1.2)
        fig,ax=plt.subplots(figsize=(9,height));maxlen=max(g['protein_length'] for g in selected)
        used=set()
        for i,g in enumerate(selected):
            ax.hlines(i,0,g['protein_length'],color='#888888',lw=1)
            for item in g['domain_architecture'].split(';'):
                if not item:continue
                name,interval=item.split(':');start,end=map(int,interval.split('-'));used.add(name)
                ax.broken_barh([(start,end-start+1)],(i-.17,.34),facecolors=palette[name])
            for item in g['motif_architecture'].split(';'):
                if not item:continue
                name,interval=item.split(':');start,end=map(int,interval.split('-'))
                idx=int(name.split('-')[-1]);ax.broken_barh([(start,end-start+1)],(i+.21,.08),facecolors=plt.get_cmap('tab10')((idx-1)%10))
        ax.set_yticks(range(len(selected)),[g['gene_id'] for g in selected],fontsize=8)
        ax.invert_yaxis();ax.set_xlim(0,maxlen*1.05);ax.set_xlabel('Protein position (amino acids)')
        ax.set_title(og['Orthogroup']+'  '+og['corrected_module'],fontsize=11,loc='left')
        handles=[plt.Line2D([0],[0],color=palette[n],lw=5,label=n) for n in sorted(used)]
        if handles:ax.legend(handles=handles,ncol=3,loc='upper center',bbox_to_anchor=(.5,-.2),fontsize=7,frameon=False)
        fig.tight_layout();fig.savefig(OUT/'figures'/f"ecm_{og['Orthogroup']}.png",dpi=200,bbox_inches='tight');plt.close(fig)
    fig,ax=plt.subplots(figsize=(9,4));ordered=sorted(ogs,key=lambda r:int(r['orthology_pair_count']),reverse=True)[:15]
    ax.barh(range(len(ordered)),[int(r['orthology_pair_count']) for r in ordered],color='#007A87');ax.set_yticks(range(len(ordered)),[r['Orthogroup'] for r in ordered]);ax.invert_yaxis();ax.set_xlabel('OrthoFinder inferred orthology pairs (model level)');fig.tight_layout();save(fig,'Figure3_orthology_evidence')

if __name__=='__main__':main()
