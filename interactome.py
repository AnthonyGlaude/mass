# ──────────────────────────────────────────────────────────────────────────────
# 0) Imports et chemins
# ──────────────────────────────────────────────────────────────────────────────
import requests
import pandas as pd
from collections import defaultdict

import networkx as nx
import matplotlib.pyplot as plt

from goatools.obo_parser import GODag
from goatools.go_enrichment import GOEnrichmentStudy
from goatools.associations import read_gaf

# Chemins vers tes fichiers
INTERACT = r"C:\Users\Antho\OneDrive - USherbrooke\Documents\0\YAhoop\interaction.txt"
FASTA    = r"D:\Anthony\Mass_spect\human_proteoforme_ref\uniprotkb_proteome_UP000005640_2025_04_29.fasta"
OBO = r"D:\Anthony\Mass_spect\Final_output\interpretation\data\go-basic.obo"
GAF      = r"D:\Anthony\Mass_spect\Final_output\interpretation\data\goa_human.gaf"
# Score minimal STRING à 400 (0.4)
STRING_SCORE = 400  
# Seuil Fold_MS pour définir tes "hits" à enrichir
FOLD_THRESHOLD = 5

# ──────────────────────────────────────────────────────────────────────────────
# 1) Charger et préparer les données SAINT
# ──────────────────────────────────────────────────────────────────────────────
# interaction.txt : IP_name, bait, prey, count
int_df = pd.read_csv(
    INTERACT, sep='\t', header=None,
    names=['IP_name','bait','prey','count']
)

# Calcul des SpecSum (somme des counts) par bait/prey
spec_sum = (
    int_df
    .groupby(['bait','prey'])['count']
    .sum()
    .reset_index(name='SpecSum')
)

# Extraction du SpecSum dans tes contrôles IgG
ctrl_sum = (
    spec_sum
    .query("bait=='ctrl_IgG'")
    .loc[:, ['prey','SpecSum']]
    .rename(columns={'SpecSum':'CtrlSum'})
)

# Merge bait vs ctrl + calcul Fold_MS
ratio_df = (
    spec_sum
    .merge(ctrl_sum, on='prey', how='left')
    .fillna({'CtrlSum':0})
)
ratio_df['Fold_MS'] = (ratio_df['SpecSum'] + 1) / (ratio_df['CtrlSum'] + 1)

# Sélectionne tes interacteurs d’intérêt (Fold_MS ≥ seuil)
study_proteins = (
    ratio_df
    .query("Fold_MS >= @FOLD_THRESHOLD")
    ['prey']
    .unique()
    .tolist()
)

# Background = toutes les protéines identifiées
background = int_df['prey'].unique().tolist()

# Nettoyage : retire les suffixes d'isoformes (le -1, -2 etc.)
def clean_uniprot_ids(protein_list):
    return [prot.split('-')[0] for prot in protein_list]

study_proteins = clean_uniprot_ids(study_proteins)
background = clean_uniprot_ids(background)
# ➔ Map UniProt -> gene symbols pour les listes d'IDs
def uniprot_to_gene(protein_list, gaf_file):
    conv = {}
    with open(gaf_file) as f:
        for l in f:
            if l.startswith('!'): continue
            parts = l.split('\t')
            conv[parts[1]] = parts[2]
    return [conv.get(prot, prot) for prot in protein_list]

study_proteins = uniprot_to_gene(study_proteins, GAF)
background = uniprot_to_gene(background, GAF)

print(f"Gene symbols pour enrichissement (extrait) : {study_proteins[:5]}")
print(f"Gene symbols background (extrait) : {background[:5]}")



print(f"{len(study_proteins)} interacteurs retenus pour l’enrichissement (Fold_MS ≥ {FOLD_THRESHOLD})")
print(f"{len(background)} protéines background")


# ──────────────────────────────────────────────────────────────────────────────
# 2) Enrichissement GO avec goatools
# ──────────────────────────────────────────────────────────────────────────────
def ids_to_gene2gos(id2gos, gaf_file):
    # map UniProt -> gene symbol
    conv = {}
    with open(gaf_file) as f:
        for l in f:
            if l.startswith('!'): continue
            parts = l.split('\t')
            conv[parts[1]] = parts[2]
    # swap keys to gene names
    g2gos = {}
    for uid, gos in id2gos.items():
        gene = conv.get(uid)
        if gene:
            g2gos[gene] = gos
    return g2gos

def gene_enrichment(study, background, namespace):
    dag    = GODag(OBO)
    id2gos = read_gaf(GAF, godag=dag, namespace=namespace)
    g2gos  = ids_to_gene2gos(id2gos, GAF)
    gee = GOEnrichmentStudy(
        background, g2gos, dag,
        propagate_counts=False,
        alpha=0.05,
        methods=['fdr_bh']
    )
    results = gee.run_study(study)
    # filtre p < 0.05
    sig = []
    for r in results:
        p = getattr(r, 'p_fdr_bh')
        if p < 0.05:
            sig.append({
                'GO': r.GO,
                'Name': r.name,
                'NS': r.NS,
                'p_fdr_bh': p,
                'Study_Count': r.study_count,
                'Pop_Count': r.pop_count
            })
    df = pd.DataFrame(sig)
    if not df.empty:
        return df.sort_values('p_fdr_bh')
    else:
        print(f"@@ Aucun terme significatif trouvé pour le namespace {namespace}.")
        return pd.DataFrame()


# BP et MF
go_bp = gene_enrichment(study_proteins, background, namespace='BP')
go_mf = gene_enrichment(study_proteins, background, namespace='MF')

print("\n➤ Top 10 GO Biological Process:")
print(go_bp.head(10))
print("\n➤ Top 10 GO Molecular Function:")
print(go_mf.head(10))


# sauvegarde éventuelle
go_bp.to_csv("GO_BP_enrichment.csv", index=False)
go_mf.to_csv("GO_MF_enrichment.csv", index=False)

# ──────────────────────────────────────────────────────────────────────────────
# 3) Interactions connues via STRING et graphe NetworkX
# ──────────────────────────────────────────────────────────────────────────────
from io import StringIO  # ajoute cet import en haut

def string_api(method, **params):
    url    = "https://version-12-0.string-db.org/api/tsv/" + method
    resp   = requests.post(url, data=params)
    df     = pd.read_csv(StringIO(resp.text), sep='\t')
    return df


# récupère les edges pour tous tes study_proteins
ids = "%0d".join(study_proteins)
sf = string_api("network",
                identifiers=ids,
                species=9606,
                required_score=STRING_SCORE)

# on ne garde que les colonnes utiles
edges = sf[['preferredName_A','preferredName_B','score']]
edges['score'] = edges['score'].astype(float)

# construit le graphe
G = nx.from_pandas_edgelist(edges, 'preferredName_A','preferredName_B',
                            edge_attr='score')

# positionnement
pos = nx.spring_layout(G, k=0.5, weight='score', seed=42)

# trace
plt.figure(figsize=(10,10))
nx.draw_networkx_nodes(G, pos,
                       node_color='skyblue',
                       node_size=300,
                       alpha=0.8)
# edges colorées selon le score
nx.draw_networkx_edges(G, pos,
                       edge_color=['black' if w>=0.7 else 'grey'
                                   for _,_,w in G.edges(data='score')],
                       width=1.5, alpha=0.7)
nx.draw_networkx_labels(G, pos, font_size=8)
plt.title("Interactome de ton Fold_MS ≥ 5 via STRING (score≥0.4)")
plt.axis('off')
plt.tight_layout()
plt.show()
