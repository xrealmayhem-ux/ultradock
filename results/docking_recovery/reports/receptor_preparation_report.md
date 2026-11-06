# Preparación localizada de receptores 5-HT2A

Se preservaron los archivos originales. Solo se transfirieron átomos de cadenas laterales faltantes o de las tres mutaciones de 7WC8 desde receptores experimentales locales, tras ajuste rígido por backbone/CB; no se modeló el resto de la proteína.

| Receptor | Átomos pesados antes/después | Reparaciones aplicadas | Problemas no reconstruidos | Dictamen |
|---|---:|---:|---|---|
| 7WC8 | 2798/2787 | 3 | ninguno | preparado para el piloto |
| 9LL8 | 1993/2019 | 7 | ninguno | preparado para el piloto |
| 9LLB | 2024/2059 | 11 | ninguno | preparado para el piloto |

## Residuos injertados y evidencia

| Receptor | Posición UniProt | Cambio | Plantilla | RMSD backbone/CB (Å) | Clearance mínimo (Å) | Evidencia de selección | Estado |
|---|---:|---|---|---:|---:|---|---|
| 7WC8 | 162 | LYS→SER (OG) | 9LLB | 0.027 | 2.825 | lowest local fit with acceptable steric clearance OG:LEU116:CD1 | grafted_from_local_experimental_template |
| 7WC8 | 164 | TRP→MET (CE;CG;SD) | 9LLB | 0.048 | 3.376 | lowest local fit with acceptable steric clearance SD:PRO246:CB | grafted_from_local_experimental_template |
| 7WC8 | 372 | ASN→SER (OG) | 9LL8 | 0.065 | 3.837 | lowest local fit with acceptable steric clearance OG:LEU331:O | grafted_from_local_experimental_template |
| 9LL8 | 146 | LYS→LYS (CD;CE;CG;NZ) | 9LLB | 0.17 | 4.11 | joint sidechain fit with ASP218 minimizes a close contact CG:LEU215:CD2 | grafted_from_local_experimental_template |
| 9LL8 | 216 | GLN→GLN (CD;CG;NE2;OE1) | 7WC8 | 0.186 | 3.088 | lowest local fit with acceptable steric clearance NE2:VAL212:O | grafted_from_local_experimental_template |
| 9LL8 | 217 | ASP→ASP (CG;OD1;OD2) | 7WC8 | 0.599 | 2.874 | lowest local fit with acceptable steric clearance OD2:LYS220:CE | grafted_from_local_experimental_template |
| 9LL8 | 218 | ASP→ASP (CG;OD1;OD2) | 7WC8 | 0.068 | 2.673 | lowest local fit with acceptable steric clearance OD1:GLY214:O | grafted_from_local_experimental_template |
| 9LL8 | 223 | LYS→LYS (CD;CE;CG;NZ) | 7WC8 | 0.02 | 3.115 | lowest local fit with acceptable steric clearance CE:LEU228:CD1 | grafted_from_local_experimental_template |
| 9LL8 | 224 | GLU→GLU (CD;CG;OE1;OE2) | 7WC8 | 0.616 | 3.62 | lowest local fit with acceptable steric clearance CG:GLY225:N | grafted_from_local_experimental_template |
| 9LL8 | 355 | GLU→GLU (CD;CG;OE1;OE2) | 7WC8 | 0.097 | 3.166 | lowest local fit with acceptable steric clearance CG:ASP356:N | grafted_from_local_experimental_template |
| 9LLB | 217 | ASP→ASP (CG;OD1;OD2) | 7WC8 | 0.284 | 2.347 | lowest local fit with acceptable steric clearance OD1:ASP218:N | grafted_from_local_experimental_template |
| 9LLB | 219 | SER→SER (OG) | 7WC8 | 0.101 | 2.945 | lowest local fit with acceptable steric clearance OG:LYS220:N | grafted_from_local_experimental_template |
| 9LLB | 223 | LYS→LYS (CD;CE;CG;NZ) | 7WC8 | 0.023 | 3.888 | lowest local fit with acceptable steric clearance CD:LEU228:CD1 | grafted_from_local_experimental_template |
| 9LLB | 224 | GLU→GLU (CD;CG;OE1;OE2) | 7WC8 | 0.658 | 3.309 | lowest local fit with acceptable steric clearance CG:GLY225:N | grafted_from_local_experimental_template |
| 9LLB | 231 | ASP→ASP (CG;OD1;OD2) | 7WC8 | 0.084 | 2.657 | lowest local fit with acceptable steric clearance OD1:ASP232:N | grafted_from_local_experimental_template |
| 9LLB | 232 | ASP→ASP (CG;OD1;OD2) | 9LL8 | 0.13 | 3.83 | lowest local fit with acceptable steric clearance OD1:VAL347:CG1 | grafted_from_local_experimental_template |
| 9LLB | 330 | PHE→PHE (CD1;CD2;CE1;CE2;CG;CZ) | 7WC8 | 0.041 | 3.421 | lowest local fit with acceptable steric clearance CZ:VAL334:CG2 | grafted_from_local_experimental_template |
| 9LLB | 345 | MET→MET (CE;CG;SD) | 7WC8 | 0.043 | 3.115 | lowest local fit with acceptable steric clearance CG:ILE341:O | grafted_from_local_experimental_template |
| 9LLB | 349 | CYS→CYS (SG) | 7WC8 | 0.071 | 2.419 | lowest local fit with acceptable steric clearance SG:MET345:O | grafted_from_local_experimental_template |
| 9LLB | 355 | GLU→GLU (CD;CG;OE1;OE2) | 7WC8 | 0.039 | 3.181 | lowest local fit with acceptable steric clearance CG:ASP356:N | grafted_from_local_experimental_template |
| 9LLB | 356 | ASP→ASP (CG;OD1;OD2) | 7WC8 | 0.029 | 4.008 | lowest local fit with acceptable steric clearance OD1:GLU355:CG | grafted_from_local_experimental_template |

Las sustituciones S162K, M164W y S372N de 7WC8 se trataron como mutaciones ingenierizadas, no como identidad WT. Se generó un modelo de trabajo con SER162, MET164 y SER372 transferidos desde los receptores WT 9LL8/9LLB; se conserva la conformación de backbone de 7WC8. Por ello, este receptor es un modelo de constructo WT-like, no una estructura experimental WT.

GLU351 de 7WC8 no se reconstruyó: no hay plantilla experimental local completa, su residuo pesado más cercano está a 12.26 Å de las poses piloto y a 4.80 Å fuera del volumen de caja. Se documenta como una limitación periférica; no se inventó su rotámero.

Para 9LL8 se seleccionó LYS146 de la plantilla experimental 9LLB y ASP218 de 7WC8. La búsqueda combinatoria local mostró que injertar ambos desde 7WC8 deja NZ(LYS146)-OD2(ASP218) a 1.835 Å; la combinación seleccionada aumenta la separación entre esas cadenas laterales a 5.765 Å y conserva clearance mínimo de 2.673 Å frente al resto del receptor. Los ajustes locales backbone/CB fueron 0.170 Å para LYS146 y 0.068 Å para ASP218. Esta corrección cambia el receptor respecto al docking piloto previo de 9LL8; sus 12 poses no se reutilizarán.

Los campos de RMSD son ajustes locales de backbone/CB y no validan por sí solos la rotámera injertada. Se revisó clearance estérico; el modelo permanece sujeto a incertidumbre conformacional. Los residuos distantes fuera de la envolvente de 6 Å no se tocaron.
