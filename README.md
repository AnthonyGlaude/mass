# mass# Pipeline MaxQuant + SAINTq

Ce document sert de mémo personnel pour reproduire l'analyse des données de spectrométrie de masse avec **MaxQuant** et l'analyse des interactions protéiques avec **SAINTq**.

---

## Données brutes

- Charger le dossier contenant les fichiers bruts (`.raw` ou `.d`).
- **Pas de fraction**.
- Définir un **parameter group** par condition (ex. : lignée cellulaire).

### Paramètres spécifiques au groupe

- **Type** : `TIMS-DDA`
- **Label-free quantification (LFQ)** : ✅ Activé

### Paramètres globaux

- **Identification** : cocher `Match between runs`

>  **Mémoire requise** : 1 analyse = 1 thread ≈ **6–8 Go de RAM**

---

## Fichier de sortie important

- `combined/txt/proteinGroups.txt`

---

##  Traitement de `proteinGroups.txt`

1. Ouvrir dans Excel.
2. Filtrer pour supprimer les lignes contenant :
   - `Only identified by site`
   - `Reverse`
   - `Potential contaminant`
3. Sauvegarder sous :  
   ➤ `proteinGroups_filtered.tsv`

---

##  Fichiers requis pour SAINTq

Ces fichiers sont créés par un script de préparation :

- `proteinGroups_filtered.tsv`  
- `bait.txt`
- `prey.txt`
- `interaction.txt`

---

##  Lancement de SAINTq

Utiliser WSL pour exécuter la commande suivante dans le répertoire contenant SAINTq :

```bash
./SAINTexpress-spc ../interaction.txt ../prey.txt ../bait.txt
```
---
---
## Analyse de l’interactome et enrichissement fonctionnel

Ce script Python permet d’analyser les résultats de l’interactome obtenus avec SAINTq en combinant :

1. Le calcul du **Fold_MS** pour identifier les interacteurs spécifiques (vs. contrôle IgG),
2. L’identification des **gènes enrichis en fonctions GO** via le package `goatools` (fichiers GAF et OBO),
3. L’interrogation de la **base STRING** pour obtenir les interactions connues entre les protéines sélectionnées,
4. La **visualisation de l’interactome** sous forme de graphe à l’aide de `networkx` et `matplotlib`.

### Entrées requises

- `interaction.txt` : par SAINTq (colonnes : IP_name, bait, prey, count)
- Fichier GAF (`goa_human.gaf`) : annotations gène/GO
- Fichier OBO (`go-basic.obo`) : hiérarchie des termes GO
- Fichier FASTA de référence


### Sorties

- `GO_BP_enrichment.csv` : enrichissement des processus biologiques (BP)
- `GO_MF_enrichment.csv` : enrichissement des fonctions moléculaires (MF)
- Une figure matplotlib représentant l’interactome entre les protéines enrichies (score STRING ≥ 0.4)

### Seuils utilisés

- **Fold_MS ≥ 5** pour définir les interacteurs d’intérêt (ajustable)
- **Score STRING ≥ 400** (correspond à 0.4)

---

### Exemple d’exécution

Ce script peut être exécuté tel quel après modification des chemins des fichiers :

```bash
python analyse_interactome.py

