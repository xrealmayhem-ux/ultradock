# Selección y control de calidad de ligandos para docking 5-HT2A

Generado: 2026-10-28T20:56:42+00:00

## Método

Selección calculada a partir de `ligands_master.csv`, `experimental_evidence.csv`, `functional_labels.csv` y `label_conflicts.csv`. Se reutilizaron los CSV locales; no se descargaron fuentes ni se modificaron los datos raw. La unidad de selección es la entidad molecular con InChIKey completo, que conserva estereoisómeros distintos.

Se incluyen anotaciones funcionales explícitas solo cuando proceden de evidencia funcional curada, y binders solo con evidencia pertinente de clase binding. Las moléculas sin clase funcional explícita se admiten como `binder` únicamente si tienen evidencia de unión; una anotación funcional conflictiva se representa como `binder` solo cuando además existe evidencia de unión, manteniendo `conflict_flag=true`. Las clases permitidas son exclusivamente `agonist`, `antagonist` y `binder`.

## Recuentos

- Moléculas del inventario maestro revisadas: **10157**.
- Moléculas seleccionadas: **6565**.
- Estructuras válidas entre seleccionadas: **6565**.
- Elegibles para pasar a preparación estructural: **6565**.
- Conflictos entre clases funcionales incluidos como binder con binding confirmado: **7**.
- Conflictos entre clases funcionales sin binding y excluidos: **3**.
- Discrepancias de subtipo incluidas manteniendo la clase analítica: **32**.
- Moléculas marcadas en `label_conflicts.csv` revisadas: **43**.
- Moléculas antes sin clase funcional, reclasificadas como binder por evidencia de unión: **629**.
- Moléculas antes sin clase funcional, excluidas por falta de unión: **1443**.
- Moléculas excluidas por falta de clase admitida y binding: **1619**.
- Registros de evidencia sin estructura enlazable, excluidos de docking: **98**.
- Moléculas con etiqueta/evidencia suficiente, pero diferidas por revisión química: **2149**.
- Binders antes sin clase concluyente que pasan ahora a revisión química: **247**.
- Moléculas diferidas con múltiples fragmentos: **739**.
- Moléculas diferidas con estereoquímica no especificada: **1612**.

| Etiqueta | Moléculas |
|---|---:|
| agonist | 186 |
| antagonist | 191 |
| binder | 6188 |

### Cobertura de evidencia integrada

| Fuente | Registros | Moléculas con enlace químico |
|---|---:|---:|
| chembl | 21219 | 10078 |
| guide_to_pharmacology | 139 | 130 |
| pubchem | 2754 | 1230 |

### Solapamientos y duplicación de evidencia

| Fuentes | Moléculas compartidas por InChIKey | Evidencias espejo exactas marcadas |
|---|---:|---:|
| chembl / guide_to_pharmacology | 81 | 0 |
| chembl / pubchem | 1200 | 231 |
| guide_to_pharmacology / pubchem | 36 | 0 |

Los espejos se conservan en `experimental_evidence.csv` y se señalan mediante `duplicate_of`; no cuentan como experimentos independientes. Las moléculas compartidas se mantienen una sola vez en `docking_ligands.csv`.

## Calidad química y limitaciones

- Identificadores duplicados por InChIKey entre ligand_id: **0**; la exportación se aborta si aparece una colisión.
- Estructuras inválidas entre moléculas seleccionadas: **0**.
- Estructuras inválidas de moléculas del inventario maestro excluidas: **0**.
- SMILES original ausente entre seleccionadas: **0**.
- SMILES original se recupera de los CSV raw cuando se encuentra una representación que RDKit mapea al mismo ligand_id; si no está disponible queda vacío y el estado lo indica.
- Se conservan sales/múltiples fragmentos y se marcan para revisión; no se eliminan fragmentos, no se neutralizan estructuras, y no se fusionan tautómeros ni estados de protonación.
- `eligible_for_structure_preparation` significa estructura parseable y evidencia admisible; no significa PDBQT preparado ni docking listo. `requires_chemical_review` sigue en el CSV seleccionado, pero se separa de la cuenta apta hasta resolver fragmentos/estereoquímica. Las tablas no contienen scores de docking.
- Las 2.319 moléculas antes `unclassified` se revisan individualmente mediante sus filas de evidencia: solo pasan como binder si hay evidencia binding; el resto se exporta en exclusiones.
- El inventario de evidencia disponible para PubChem se limita a los tres datatables locales iniciales; los metadatos de AIDs adicionales no aportan tablas de resultados.
- Un ChEMBL assay con organismo experimental vacío se retiene solo si el target anotado es humano; la incertidumbre queda visible en la evidencia curada.

### Exclusiones

| Motivo | Filas |
|---|---:|
| evidence_record_missing_or_invalid_structure | 98 |
| no_supported_functional_or_binding_label | 1443 |
| requires_chemical_review:multiple_fragments | 537 |
| requires_chemical_review:multiple_fragments|stereochemistry_unspecified | 202 |
| requires_chemical_review:stereochemistry_unspecified | 1410 |

## Archivos

- `data/processed/docking_ligands.csv`
- `results/docking_ligands_qc.csv`
- `results/docking_ligands_exclusions.csv`
- `results/dataset_manifest.json` actualizado aditivamente.
