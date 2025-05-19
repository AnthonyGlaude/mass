# mass# Pipeline MaxQuant + SAINTq

Ce document sert de mémo personnel pour reproduire l'analyse des données de spectrométrie de masse avec **MaxQuant** et l'analyse des interactions protéiques avec **SAINTq**.

---

## 🧪 Données brutes

- Charger le dossier contenant les fichiers bruts (`.raw` ou `.d`).
- **Pas de fraction**.
- Définir un **parameter group** par condition (ex. : lignée cellulaire).

### Paramètres spécifiques au groupe

- **Type** : `TIMS-DDA`
- **Label-free quantification (LFQ)** : ✅ Activé

### Paramètres globaux

- **Identification** : cocher `Match between runs`

> 💡 **Mémoire requise** : 1 analyse = 1 thread ≈ **6–8 Go de RAM**

---

## 📂 Fichier de sortie important

- `combined/txt/proteinGroups.txt`

---

## 🧼 Traitement de `proteinGroups.txt`

1. Ouvrir dans Excel.
2. Filtrer pour supprimer les lignes contenant :
   - `Only identified by site`
   - `Reverse`
   - `Potential contaminant`
3. Sauvegarder sous :  
   ➤ `proteinGroups_filtered.tsv`

---

## 🧰 Fichiers requis pour SAINTq

Ces fichiers sont créés par un script de préparation :

- `proteinGroups_filtered.tsv`  
- `bait.txt`
- `prey.txt`
- `interaction.txt`

---

## 🚀 Lancement de SAINTq

Utiliser WSL pour exécuter la commande suivante dans le répertoire contenant SAINTq :

```bash
./SAINTexpress-spc ../interaction.txt ../prey.txt ../bait.txt
