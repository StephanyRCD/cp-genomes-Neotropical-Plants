# cp-ROIs-Cordillera Escalera
Pipeline reproducible para seleccionar ROIs (regiones de interés) del cloroplasto mediante mapas de diversidad (π), distancias (p-distance), barcoding gap y filogenia comparativa (árbol cp vs. árbol de ROIs). Enfoque aplicado a familias de la Cordillera Escalera.
## Objetivo 
Comparar la conservación y divergencia de 5–15 regiones cloroplásticas (ROIs) entre familias relevantes de la Cordillera Escalera y priorizar las ROIs con mejor poder discriminante en función de su variación (π, p-distance, barcoding gap) y coherencia con la filogenia del cloroplasto. 
## Flujo
1. **Construcción del set de referencia cp** → descarga, estandarización, control de calidad.
2. **Alineamientos y ventanas** → MSA (alineamiento múltiple) por genes/regiones; π y entropía (Shannon) con ventana 600 bp / paso 200 bp + chequeo 400/150 y 800/250;    shortlist 10–15 ROIs (≥200 bp).
3. **Métricas por ROI** → π intra-familia, p-distance inter-género/inter-familia, barcoding gap; tabla ranking (longitud, % sitios informativos, π_intra mediana/P90, p-dist_inter mediana/P10, Δgap, bandera “cerca de IR/repeats”).
4. **Filogenia comparativa** → árbol cp (subconjunto de CDS concatenados) vs. árbol de ROIs concatenadas (particionado); comparar topologías y soporte. 

