import pandas as pd
from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator
from rdkit.SimDivFilters.rdSimDivPickers import MaxMinPicker

# Configuración
INPUT = "tus_moleculas.csv"
OUTPUT = "moleculas_para_docking.csv"
SMILES_COL = "SMILES"
CLASS_COL = "Clase"
TARGET_BINDERS = 2000
SEED = 42

# 1. Cargar y validar columnas
df = pd.read_csv(INPUT)

required = {SMILES_COL, CLASS_COL}
missing = required - set(df.columns)
if missing:
    raise ValueError(f"Faltan columnas requeridas: {sorted(missing)}")

labels = df[CLASS_COL].astype("string").str.strip().str.casefold()
binder_mask = labels.eq("binder").fillna(False)

df_binders = df.loc[binder_mask]
df_other = df.loc[~binder_mask]

if df_binders.empty:
    raise ValueError("No se encontraron filas con Clase = 'binder'.")

# 2. Parsear SMILES, eliminar duplicados estructurales exactos
# y mantener el índice original de cada molécula válida.
fpgen = rdFingerprintGenerator.GetMorganGenerator(
    radius=2,
    fpSize=2048
)

valid_mols = []
valid_indices = []
seen_smiles = set()
invalid_indices = []
duplicate_indices = []

for idx, smi in df_binders[SMILES_COL].items():
    mol = Chem.MolFromSmiles(str(smi)) if pd.notna(smi) else None

    if mol is None:
        invalid_indices.append(idx)
        continue

    canonical = Chem.MolToSmiles(
        mol,
        canonical=True,
        isomericSmiles=True
    )

    if canonical in seen_smiles:
        duplicate_indices.append(idx)
        continue

    seen_smiles.add(canonical)
    valid_mols.append(mol)
    valid_indices.append(idx)

if not valid_mols:
    raise ValueError("No hay binders válidos y únicos para seleccionar.")

# 3. Generar fingerprints
fps = [fpgen.GetFingerprint(mol) for mol in valid_mols]

# 4. Seleccionar diversidad química con MaxMin
n_pick = min(TARGET_BINDERS, len(fps))

picker = MaxMinPicker()
picks = list(
    picker.LazyBitVectorPick(
        fps,
        len(fps),
        n_pick,
        seed=SEED
    )
)

# Correspondencia correcta fingerprint -> fila original
selected_indices = [valid_indices[i] for i in picks]

# 5. Conservar todos los no-binders y agregar los binders seleccionados
df_selected_binders = df.loc[selected_indices]

df_final = pd.concat(
    [df_other, df_selected_binders]
).sort_index().reset_index(drop=True)

df_final.to_csv(OUTPUT, index=False)

# 6. Auditoría: guardar filas descartadas sin alterar el original
df.loc[invalid_indices].to_csv(
    "binders_invalidos.csv", index=False
)
df.loc[duplicate_indices].to_csv(
    "binders_duplicados_exactos.csv", index=False
)

print(f"Filas iniciales: {len(df)}")
print(f"Binders iniciales: {len(df_binders)}")
print(f"Binders inválidos: {len(invalid_indices)}")
print(f"Duplicados exactos excluidos: {len(duplicate_indices)}")
print(f"Binders elegibles: {len(valid_mols)}")
print(f"Binders seleccionados: {len(picks)}")
print(f"Filas no-binder conservadas: {len(df_other)}")
print(f"Total final: {len(df_final)}")
print(f"Archivo generado: {OUTPUT}")
