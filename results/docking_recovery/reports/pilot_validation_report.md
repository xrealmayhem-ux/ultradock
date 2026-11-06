# Validación piloto de docking y ProLIF 5-HT2A

Ligandos elegidos: 12 (6 agonistas, 6 antagonistas). No se ejecutó el conjunto de 423 ligandos.
Dockings piloto completados: 36/36. ProLIF: 259 modelos de pose. Reproducibilidad de un complejo repetido: True.

Nota de vigencia: después del piloto se cambió la plantilla de LYS146 de 9LL8 para resolver su contacto estrecho con ASP218. Por ello, los 12 dockings y 84 fingerprints ProLIF previos de 9LL8 son históricos y no validan el receptor actual. No se repetirán como piloto aparte; se recalcularán una sola vez dentro del lote final. Tampoco se reutilizarán poses piloto sin hashes de entrada por ejecución.

| Receptor | Vina runs OK | Poses ProLIF | Scores Vina mediana (kcal/mol) |
|---|---:|---:|---:|
| 7WC8 | 12/12 | 83 | -10.415 |
| 9LL8 | 12/12 | 84 | -9.733 |
| 9LLB | 12/12 | 92 | -9.8775 |

## Frecuencias ProLIF (eventos)

| Receptor | Hydrophobic | HBDonor | HBAcceptor | PiStacking | Anionic | Cationic | CationPi | PiCation | VdWContact |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 7WC8 | 512 | 34 | 47 | 40 | 0 | 24 | 1 | 5 | 1229 |
| 9LL8 | 408 | 43 | 31 | 28 | 0 | 28 | 0 | 5 | 1075 |
| 9LLB | 500 | 39 | 36 | 43 | 0 | 30 | 0 | 0 | 1224 |

Se conservaron las nueve categorías, incluso cuando la frecuencia es cero. ASP155 se inspeccionó como diagnóstico y no se impuso como filtro. Los scores/poses nuevos no se comparan como equivalentes con los antiguos porque cambiaron el receptor y la protonación de los ligandos; los resultados históricos se conservaron.

Parámetros homogéneos del piloto: AutoDock Vina 1.2.7 CPU, scoring vina, exhaustiveness=8, num_modes=9, energy_range=3, CPU=4, seed explícito por combinación. Centro y dimensiones de caja iguales a las configuraciones originales.

Controles químicos y del pipeline se registran en `pilot_prolif/pilot_prolif_summary.json`. Errores de preparación/docking: 0.
