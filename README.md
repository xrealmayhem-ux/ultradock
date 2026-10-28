# UltraDock — paquete de entrada 5-HT2A

Paquete autocontenido para transferir al clúster los ligandos operativos, el generador y validador SDF, los informes de selección y las estructuras receptoras disponibles.

## Contenido

- `data/processed/docking_ligands.csv`: 8.714 ligandos con estructura parseable y etiqueta respaldada, con SMILES, etiquetas y trazabilidad. 6.565 no requieren revisión química; 2.149 están incluidos con `chemical_review_required=true` y deben revisarse antes de preparar estructuras.
- `scripts/prepare_docking_sdf.py`, `scripts/rawstore.py`: generación de un conformero 3D por ligando y utilidades de checksum.
- `receptors/7WC8.cif`, `receptors/9LL8.cif`, `receptors/9LLB.cif`: copias de los CIF ya presentes en `../docking/proteins/`.
- `results/docking_ligands_summary.md`, `docking_ligands_qc.csv`, `docking_ligands_exclusions.csv` y `dataset_manifest.json`: selección, controles por molécula, exclusiones y manifiesto de la curación original.
- `bundle_manifest.json`: checksums de los archivos incluidos en este paquete.

## Generar y validar el SDF en el clúster

Requiere Python, pandas y RDKit. Desde la raíz de `ultradock/`:

```bash
python scripts/prepare_docking_sdf.py
```

La ejecución crea `structures/docking_ligands_3d.sdf` y `structures/docking_ligands_sdf_qc.csv`. Preserva el estado representado por el SMILES; usa RDKit ETKDGv3 con semilla fija `20261028`, añade hidrógenos y minimiza con MMFF94s, usando UFF solo como fallback. El script valida identidad/InChIKey, estereoquímica, hidrógenos, coordenadas y convergencia antes de promover el SDF final. Si algún ligando falla, registra el QC y no presenta el SDF parcial como válido.

El SDF es un conformero inicial para preparación posterior: **no es PDBQT y no está listo por sí solo para docking**.

## Receptores

Los tres CIF son estructuras de HTR2A (7WC8: complejo con lumateperona; 9LL8: complejo con psilocina-Gi; 9LLB: complejo con DOI-NBOMe-Gq). Se incluyen como opciones, sin modificar ni convertir. Aún hay que seleccionar el receptor y definir preparación, estado/cofactores y caja de docking antes de acoplar ligandos. Este paquete no ejecuta docking.

## Integridad

Verifica `bundle_manifest.json` después de copiar el directorio al clúster. `dataset_manifest.json` conserva el manifiesto original de la curación y no se sustituye por el manifiesto del paquete.
