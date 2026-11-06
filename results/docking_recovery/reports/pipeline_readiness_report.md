# Preparación para continuar el pipeline

**Decisión: APROBADO CON LIMITACIONES DOCUMENTADAS.**

El piloto confirma 36/36 dockings, 259 poses fingerprintadas y reproducibilidad=True. Se mantuvieron separados los tres receptores y las categorías IFP.

Ligandos: 424 SDF procesados; decisiones de estado: {'REPETIR': 341, 'REVISIÓN MANUAL': 73, 'CONSERVAR': 10}. El criterio de Dimorphite maximiza protonación de nitrógenos alifáticos no amida y desprotonación de carboxilos a pH 7.4, con estados alternativos registrados. No se estimaron poblaciones. Los casos con ambigüedad se mantienen en revisión manual.

Combinaciones de docking: {'REPETIR': 1053, 'NO DISPONIBLE': 3, 'REVISIÓN MANUAL': 216}. No se ejecutó redocking masivo. Los resultados nuevos del piloto están en `pilot_docking/`; los resultados originales no se sobrescribieron.

Limitaciones estructurales registradas: residuos sin plantilla=ninguno; injertos con clearance pesado mínimo <2 Å=ninguno. Las conformaciones de cadenas laterales injertadas son modelos de plantilla y deben revisarse; no se interpreta un contacto corto automáticamente como choque sin inspeccionar el par atómico.

Estado de escalado masivo: LISTO. Ligandos con microestado ambiguo pendientes=72. Un estado representativo está preparado para cada uno, pero no se consideró aprobado sin revisión.

Preflight actual: READY; 1272 combinaciones incluyendo entradas no disponibles, conteos={'READY': 1269, 'NO DISPONIBLE': 3}. La cola dirigida está en `ligand_preparation/microstate_review_queue.csv`; solo se aprueba automáticamente el mismo estado ya preparado. Si se selecciona una alternativa, hay que regenerar su PDBQT antes de desbloquear el manifiesto.

Los 36 dockings del piloto no tienen hashes de entrada por ejecución, así que no se reutilizarán por suposición. El ejecutor final los incluirá una sola vez dentro de las 1.269 combinaciones disponibles; no se harán dockings piloto adicionales.

La validación ProLIF para el receptor 9LL8 corresponde al modelo anterior al ajuste conjunto LYS146/ASP218; el lote final actualizará docking e IFP en ese receptor.

Tras cerrar la cola y refrescar el manifiesto, ejecutar `python scripts/run_frozen_docking_batch.py --execute`. El proceso valida hashes y reanuda únicamente combinaciones con metadatos/pose coincidentes. Cuando termine, calcular IFP con `python scripts/run_frozen_ifp.py --batch-fingerprint <fingerprint> --execute`.

El piloto valida el cálculo de fingerprints, no aprueba por sí mismo la geometría de todos los receptores ni resuelve la selección de microestados. No iniciar el lote hasta congelar estas decisiones.
