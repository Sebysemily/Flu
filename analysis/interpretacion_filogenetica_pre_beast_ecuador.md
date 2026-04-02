# Interpretacion filogenetica pre-BEAST (H5N1 Ecuador)

## Alcance
Este documento resume la interpretacion de los arboles ML ya generados, enfocada en las muestras de Ecuador (Flu-xxxx), y propone un marco de trabajo pre-BEAST sin ejecutar BEAST en esta etapa.

## Hallazgos principales
1. Existe un clado principal de Ecuador 2024 bien soportado en el arbol concatenado, consistente con una expansion local dominante.
2. La muestra Flu-0406 (2023) aparece como una rama separada respecto al nucleo 2024, compatible con al menos una introduccion previa o una trayectoria evolutiva distinta.
3. La posicion de Flu-0407 y Flu-0023 cambia entre segmentos frente al concatenado, lo que sugiere incongruencia inter-segmento compatible con reassortment o con historia incompleta por segmento.

## Lectura epidemiologica conservadora
1. Escenario minimo compatible: dos introducciones hacia Ecuador.
2. Introduccion 1: linaje asociado al clado mayoritario detectado en 2024.
3. Introduccion 2: linaje asociado a Flu-0406 (senal 2023) que no colapsa claramente dentro del nucleo 2024.
4. Reassortment: hay indicios, pero no debe declararse como conclusion unica solo con topologia; debe integrarse con soporte de nodos, consistencia entre genes y contexto temporal.

## Estrategia de subconjuntos para reloj molecular (pre-BEAST)
Se adopta la estrategia acordada: contexto regional cercano por proximidad filogenetica y sin anclas distantes en el conjunto temporal principal.

### Panel A (principal)
1. Incluye nucleo de Ecuador completo de interes epidemiologico.
2. Agrega vecinos regionales mas cercanos al clado principal y a la rama 2023 de interes.
3. Excluye anclas eurasia/usa lejanas para no forzar pendiente temporal por outgroups distantes.

### Panel B (sensibilidad)
1. Mantiene el mismo nucleo de Ecuador.
2. Reemplaza vecinos por segunda capa de cercania regional para evaluar estabilidad de senales de reloj y posicion relativa.

## Interpretacion de root-to-tip esperada
1. Si panel A y panel B mantienen pendiente positiva y outliers similares, la senal temporal es robusta al muestreo regional cercano.
2. Si cambia drásticamente la pendiente o emergen outliers distintos entre paneles, la calibracion temporal es sensible a composicion y debe ajustarse antes de BEAST.
3. El objetivo aqui no es maximizar R2 de forma aislada, sino verificar coherencia biologica y estabilidad de senal para justificar el dataset final de BEAST.

## Relacion con hipotesis de introduccion y reassortment
1. Introducciones: el patron concatenado favorece al menos dos entradas/linajes en el periodo observado.
2. Reassortment: se mantiene como hipotesis de trabajo para Flu-0407 y Flu-0023 hasta contrastar con evidencia temporal y soporte por segmento.
3. Decision pre-BEAST: avanzar con dataset regional cercano y registrar en auditoria cualquier taxon con fecha incompleta o posicion inestable.

## Que queda para validar cuando ejecutes
1. Tabla de seleccion de paneles y auditoria de taxones elegidos.
2. Subarboles y alineamientos podados por panel.
3. Tabla de fechas por panel (sin faltantes).
4. Resultados root-to-tip por panel para comparar estabilidad.
5. Resumen de incongruencia inter-segmento para muestras objetivo (Flu-0406, Flu-0407, Flu-0023).

## Conclusion ejecutiva
Con la evidencia topologica actual, la lectura mas solida es: un clado principal de Ecuador 2024, una senal separada de 2023 (Flu-0406) compatible con introduccion adicional, e incongruencias por segmento que justifican prueba formal pre-BEAST con paneles regionales cercanos. Esta configuracion es metodologicamente consistente con tu objetivo (foco Ecuador) y evita sesgo por anclas filogeneticamente distantes en el modelo temporal principal.
