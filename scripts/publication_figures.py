"""Publication figures from tabular searches and original maximum likelihood trees."""
import csv
import copy
import re
import os
from pathlib import Path
from collections import defaultdict
import numpy as np
from Bio import SeqIO, Phylo
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(os.environ.get('ECM_PROJECT_ROOT', Path(__file__).resolve().parents[1] / 'workspace')).resolve()
PKG=ROOT/'publication_package'
OUT=PKG/'analysis/figures'
SOURCE=PKG/'remote_sources'

def read(path):
    with path.open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f,delimiter='\t'))

def save(fig,name):
    for ext in ['png','pdf','tiff']:
        kw={'dpi':300} if ext!='pdf' else {}
        if ext=='tiff':kw['pil_kwargs']={'compression':'tiff_lzw'}
        fig.savefig(OUT/f'{name}.{ext}',bbox_inches='tight',**kw)
    plt.close(fig)

def domain_overview():
    rows=read(PKG/'domain_validation/species_counts.tsv')
    lookup={r['dataset']:r for r in rows}
    sets=['D_japonicum','Dicyema_unfiltered','I_linei','C_teleta','H_robusta','L_gigantea','P_canaliculata','C_gigas','O_sinensis','O_vulgaris','S_mansoni']
    cols=['collagen_repeat_GA','collagenIV_NC1_GA','fibrillar_collagen_C_GA','laminin_N_GA','laminin_G_GA','laminin_EGF_GA','integrin_beta_core_GA','integrin_alpha_FG_GAP_GA','cadherin_core_GA','actin_GA','myosin_head_GA','troponin_GA','tropomyosin_GA']
    labels=['Collagen','C4 NC1','COLFI','Laminin N','Laminin G','Laminin EGF','Integrin beta','FG GAP','Cadherin core','Actin','Myosin head','Troponin','Tropomyosin']
    counts=np.array([[int(lookup[s][c]) for c in cols] for s in sets])
    fig,ax=plt.subplots(figsize=(10,5.5))
    im=ax.imshow(np.log2(counts+1),cmap='YlGnBu',vmin=0,vmax=np.log2(counts.max()+1),aspect='auto')
    ax.set_xticks(range(len(cols)),labels,rotation=55,ha='right',fontsize=9)
    names=['D japonicum strict','D japonicum unfiltered']+[s.replace('_',' ') for s in sets[2:]]
    ax.set_yticks(range(len(sets)),names,fontsize=9)
    for i in range(len(sets)):
        for j in range(len(cols)):
            ax.text(j,i,str(counts[i,j]),ha='center',va='center',fontsize=9,color='white' if counts[i,j]>10 else 'black')
    ax.set_title('Profile supported protein models at Pfam gathering thresholds',loc='left',fontsize=12)
    cb=fig.colorbar(im,ax=ax,fraction=.025,pad=.025);cb.set_label('log2 model count plus 1',fontsize=9)
    fig.tight_layout();save(fig,'Figure2_whole_proteome_domains')

def inventory_overview():
    rows=read(PKG/'analysis/retention.tsv')
    lookup={(r['category'],r['species']):r for r in rows}
    species=['D_japonicum','I_linei','C_teleta','H_robusta','L_gigantea','P_canaliculata','C_gigas','O_sinensis','O_vulgaris','S_mansoni']
    fig,axes=plt.subplots(1,3,figsize=(12,5.2),gridspec_kw={'width_ratios':[1.35,1,.8]})
    y=np.arange(len(species))
    for offset,category,color,label in [(-.18,'ECM_seed','#007A87','ECM seeds'),(.18,'Muscle_seed','#A34770','Muscle seeds')]:
        values=[float(lookup[category,s]['representation_percent']) for s in species]
        axes[0].barh(y+offset,values,height=.35,color=color,label=label)
    axes[0].set_yticks(y,[s.replace('_',' ') for s in species],fontsize=8)
    axes[0].invert_yaxis();axes[0].set_xlabel('Seeded OG representation (%)',fontsize=9)
    axes[0].legend(fontsize=8,frameon=False,loc='lower right')
    ogs=read(PKG/'analysis/OG_evidence.tsv')
    categories=['adhesion_domain','receptor_domain','extracellular_associated','proteolysis_domain','shared_or_partial','intracellular_associated','intracellular_or_nonECM','unresolved']
    labels=['Cadherin core','Integrin beta','Extracellular associated','Proteolysis associated','Shared or partial','Intracellular adhesion','Intracellular or other','Unresolved']
    values=[sum(r['evidence_status']==c for r in ogs) for c in categories]
    axes[1].barh(range(len(categories)),values,color='#578559')
    axes[1].set_yticks(range(len(categories)),labels,fontsize=8);axes[1].invert_yaxis()
    axes[1].set_xlabel('Retained ECM seeded OGs',fontsize=9)
    axes[2].bar(range(4),[137,100,101,105],color=['#A34770','#007A87','#007A87','#007A87'])
    axes[2].set_xticks(range(4),['Raw\nmodels','0.25','0.50','0.80'],fontsize=8)
    axes[2].set_ylabel('Models or overlap components',fontsize=9)
    axes[2].set_xlabel('Shared CDS threshold',fontsize=9)
    for i,v in enumerate([137,100,101,105]):axes[2].text(i,v+2,str(v),ha='center',fontsize=9)
    axes[2].set_ylim(0,155)
    for ax,label in zip(axes,'ABC'):ax.set_title(label,loc='left',fontweight='bold')
    fig.tight_layout();save(fig,'Figure1_candidate_inventory')

def architecture_examples():
    rows=read(PKG/'analysis/gene_evidence.tsv')
    selected=['OG0000187','OG0000510','OG0000136','OG0000188','OG0000240','OG0000034']
    fig,axes=plt.subplots(3,2,figsize=(12,9))
    palette=['#007A87','#A34770','#578559','#D28031','#566FB1','#8C5CA6','#558FAD','#767676','#C83F41','#A6A135','#41663B','#B762A0']
    for ax,og,label in zip(axes.flat,selected,'ABCDEF'):
        genes=[r for r in rows if r['Orthogroup']==og]
        if len(genes)>4:genes=sorted(genes,key=lambda r:-int(r['protein_length']))[:4]
        names=sorted({part.split(':')[0] for r in genes for part in r['domain_architecture'].split(';') if part})
        colors={n:palette[i%len(palette)] for i,n in enumerate(names)}
        for y,r in enumerate(genes):
            ax.hlines(y,0,int(r['protein_length']),color='#AAAAAA',lw=1)
            for part in r['domain_architecture'].split(';'):
                if not part:continue
                n,coords=part.split(':');start,end=map(int,coords.split('-'))
                ax.broken_barh([(start,end-start+1)],(y-.18,.36),facecolor=colors[n])
            for part in r['motif_architecture'].split(';'):
                if not part:continue
                n,coords=part.split(':');start,end=map(int,coords.split('-'));idx=int(n.split('-')[-1])
                ax.broken_barh([(start,end-start+1)],(y+.22,.09),facecolor=plt.get_cmap('tab10')((idx-1)%10))
        ax.set_yticks(range(len(genes)),[r['gene_id'].replace('D_japonicum','Dj ') for r in genes],fontsize=8)
        ax.set_ylim(len(genes)-.4,-.7);ax.set_xlabel('Protein position (aa)',fontsize=8)
        ax.set_title(label+'  '+og,loc='left',fontsize=10,fontweight='bold')
        handles=[plt.Line2D([0],[0],color=colors[n],lw=5,label=n) for n in names]
        ax.legend(handles=handles,ncol=3,loc='upper left',bbox_to_anchor=(0,-.35),frameon=False,fontsize=7)
    fig.subplots_adjust(hspace=1.1,wspace=.38,top=.94,bottom=.10)
    save(fig,'Figure3_selected_architectures')

def selected_tree(og,focus,ax_tree,ax_domain,prune=False):
    d=SOURCE/'10.family_analysis/selected_OGs'/og
    original=Phylo.read(d/f'{og}.mafft.trim.fa.treefile','newick')
    tree=copy.deepcopy(original)
    proteins={r.id:str(r.seq) for r in SeqIO.parse(d/f'{og}.faa','fasta')}
    if prune:
        chosen={focus}
        by_species=defaultdict(list)
        for tip in tree.get_terminals():by_species[re.sub(r'\d+$','',tip.name)].append(tip.name)
        for candidates in by_species.values():chosen.add(min(candidates,key=lambda k:tree.distance(focus,k)))
        for tip in list(tree.get_terminals()):
            if tip.name not in chosen:tree.prune(tip)
    tips=tree.get_terminals();ys={tip:i for i,tip in enumerate(tips)}
    depths=tree.depths()
    if not max(depths.values()):depths=tree.depths(unit_branch_lengths=True)
    def draw(cl):
        if cl in ys:return ys[cl]
        yvalues=[draw(child) for child in cl.clades]
        mid=sum(yvalues)/len(yvalues)
        ax_tree.vlines(depths[cl],min(yvalues),max(yvalues),color='#444444',lw=.8)
        for child,y in zip(cl.clades,yvalues):
            ax_tree.hlines(y,depths[cl],depths[child],color='#444444',lw=.8)
            if child.name and '/' in child.name:
                ax_tree.text((depths[cl]+depths[child])/2,y-.15,child.name,fontsize=6,ha='center')
        return mid
    draw(tree.root)
    for tip in tips:
        label=tip.name
        ax_tree.text(depths[tip]+.04,ys[tip],label,fontsize=8,va='center',color='#007A87' if label.startswith('D_japonicum') else '#333333')
    ax_tree.set_ylim(len(tips)-.5,-.7);ax_tree.set_xlim(0,max(depths.values())*1.65);ax_tree.set_yticks([])
    ax_tree.set_xlabel('Substitutions per site',fontsize=8)
    ax_tree.set_title(og+' original tree'+(' with selected representatives' if prune else ''),loc='left',fontsize=10)
    doms=defaultdict(list)
    for row in read(d/f'{og}.pfam.domain.tsv'):doms[row['gene_id']].append(row)
    motif=defaultdict(list)
    for row in read(d/f'{og}.motif.tsv'):motif[row['gene_id']].append(row)
    names=sorted({r['domain_name'] for t in tips for r in doms[t.name]})
    colors={n:plt.get_cmap('tab20')(i%20) for i,n in enumerate(names)}
    for tip in tips:
        gene=tip.name;y=ys[tip]
        ax_domain.hlines(y,0,len(proteins[gene]),color='#999999',lw=.8)
        for row in doms[gene]:
            start,end=int(row['ali_start']),int(row['ali_end'])
            ax_domain.broken_barh([(start,end-start+1)],(y-.18,.36),facecolors=colors[row['domain_name']])
        for row in motif[gene]:
            start,end=int(row['start']),int(row['end']);idx=int(row['motif'].split('-')[-1])
            ax_domain.broken_barh([(start,end-start+1)],(y+.2,.08),facecolors=plt.get_cmap('tab10')((idx-1)%10))
    ax_domain.set_ylim(len(tips)-.5,-.7);ax_domain.set_yticks([]);ax_domain.set_xlabel('Protein position in amino acids',fontsize=8)
    handles=[plt.Line2D([0],[0],color=colors[n],lw=5,label=n) for n in names]
    ax_domain.legend(handles=handles,ncol=3,loc='lower left',bbox_to_anchor=(-.05,1.01),fontsize=7,frameon=False)

def trees():
    fig,axes=plt.subplots(2,2,figsize=(11,8),gridspec_kw={'height_ratios':[2,1],'width_ratios':[1.5,1]})
    selected_tree('OG0000767','D_japonicum05019',axes[0,0],axes[0,1],prune=True)
    selected_tree('OG0008174','D_japonicum03307',axes[1,0],axes[1,1])
    fig.subplots_adjust(hspace=.7,wspace=.06,top=.9,bottom=.08)
    save(fig,'Figure4_trees_domains_motifs')

if __name__=='__main__':
    OUT.mkdir(exist_ok=True)
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42})
    inventory_overview();domain_overview();architecture_examples();trees()
