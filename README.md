# Patrones filogenómicos y diversidad genética en un conjunto representativo de familias de plantas neotropicales con mayor diversidad y representación cloroplastidial
Este repositorio contiene el flujo completo para analizar patrones filogenómicos y diversidad genética en un conjunto representativo de familias de plantas neotropicales, seleccionadas en función de su diversidad taxonómica y disponibilidad de genomas cloroplastidiales completos.

El análisis se basa en la integración de información taxonómica y geográfica del World Checklist of Vascular Plants (WCVP, Kew Gardens) con genomas cloroplastidiales completos provenientes de NCBI RefSeq, seguido de análisis filogenéticos y estimaciones de divergencia genética.
## Hipótesis
La diversidad nucleotídica intra-familiar es menor que la diversidad inter-familiar en plantas neotropicales, lo que permite recuperar agrupamientos filogenéticos coherentes a nivel de familia utilizando genomas cloroplastidiales completos.
## Objetivo general
Analizar los patrones filogenómicos y la diversidad genética entre familias de plantas neotropicales mediante la inferencia filogenética basada en genomas cloroplastidiales completos y la estimación de divergencia nucleotídica intra- e inter-familiar.
## Objetivos específicos: 
1.	Construir un conjunto de datos genómicos curado y balanceado de plastomas de plantas neotropicales.
2.	Inferir la relación filogenética entre familias mediante métodos de máxima verosimilitud.
3.	Estimar y comparar la diversidad nucleotídica intra-familiar (π) y la distancia genética inter-familiar (p-distance).
## Metodología:
### 1. Construcción y curación de la base de datos
#### 1.1 Filtrado taxonómico neotropical (WCVP)
Se utilizan las tablas wcvp_names y wcvp_distribution del World Checklist of Vascular Plants (WCVP).
El filtrado se realiza en R mediante:
-selección de taxones aceptados y con rango Species,
detección de regiones neotropicales a partir de los campos continent, region y area,
-uso de listas definidas de regiones y países neotropicales (e.g., Amazon, South America, Peru, Ecuador, Brazil).
#### 1.2 Integración con plastomas RefSeq (NCBI)
Se descarga el release completo de plastomas RefSeq.
Desde los encabezados FASTA se extraen accession, género, especie y se construye una llave taxonómica compatible con la del WCVP.
La unión entre ambas bases se realiza usando herramientas de Unix.
#### 1.3 Selección final de plastomas (50 total)
A partir del conjunto integrado se selecciona un subconjunto representativo siguiendo criterios taxonómicos y computacionales:

- Top 10 familias neotropicales con mayor diversidad de géneros (según WCVP).

- Disponibilidad genómica mínima: ≥5 plastomas completos por familia en RefSeq.

- Diseño balanceado: 1 plastoma por género, 5 plastomas por familia.

- Total final: 50 genomas cloroplastidiales completos.

Se priorizarán accesiones completas, curadas (NC_) y familias con alta diversidad de géneros disponibles.
### 2. Alineamiento y filogenia
#### 2.1 Alineamiento
Los plastomas seleccionados se alinean utilizando MAFFT (modo automático), seguido de un trimado automático del alineamiento con trimAl.
Se evaluó el contenido de gaps posterior al trimado para asegurar la calidad del alineamiento final.
#### 2.2 Árbol filogenético
La inferencia filogenética se realiza mediante IQ-TREE 2, utilizando:

- Selección automática del modelo de sustitución (ModelFinder).

- Inferencia por máxima verosimilitud.

- Soporte nodal mediante Ultrafast Bootstrap (UFBoot).
### 3. Diversidad genética
Se estiman métricas de divergencia genética a partir del alineamiento final utilizando R:

- Diversidad nucleotídica (π) a nivel intra-familiar.

- Distancia genética promedio (p-distance) entre familias.

- Evaluación del efecto del tratamiento de gaps mediante:

  pairwise deletion

  complete deletion

# Herramientas computacionales 
## Lenguajes y entornos
-Bash — procesamiento de archivos, unión de tablas, extracción de secuencias FASTA.
-R — análisis de diversidad genética (p-distance, π), manipulación de datos, visualización.
## Software de alineamiento y filogenia
-MAFFT — alineamiento múltiple de genomas cloroplastidiales completos.
-IQ-TREE 2 — inferencia filogenética por máxima verosimilitud, ModelFinder, UFBoot.


<img width="3349" height="2797" alt="Pipeline_gv" src="https://github.com/user-attachments/assets/c6f454bb-23f9-4903-9edf-cd85181d8889" />

