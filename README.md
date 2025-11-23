# Patrones filogenéticos y diversidad genética entre familias de plantas neotropicales basadas en genomas cloroplastidiales
Este repositorio contendrá el flujo completo para analizar la estructura filogenómica y la diversidad genética entre familias de plantas neotropicales, utilizando genomas cloroplastidiales completos (RefSeq) y un filtrado taxonómico basado en WCVP (Kew Gardens).
## Hipótesis
La diversidad nucleotídica dentro de las familias neotropicales es menor que la diversidad entre familias, lo que permite recuperar agrupamientos filogenéticos coherentes mediante genomas cloroplastidiales.
## Objetivo general
Analizar los patrones filogenómicos y la diversidad genética entre familias de plantas neotropicales mediante la construcción de un árbol filogenético de genomas cloroplastidiales completos y la estimación de divergencia nucleotídica intra- e inter-familiar
## Objetivos específicos: 
1.	Construir una base de datos curada de plastomas de plantas neotropicales
2.	Generar alineamientos múltiples de los plastomas seleccionados y realizar la inferencia filogenética mediante máxima verosimilitud
3.	Estimar la diversidad nucleotídica dentro y entre familias, y comparar estos patrones con la topología filogenética obtenida.
## Metodología:
### 1. Construcción y curación de la base de datos
#### 1.1 Filtrado taxonómico neotropical (WCVP)
Se utilizarán las tablas wcvp_names y wcvp_distribution del World Checklist of Vascular Plants (WCVP).
El filtrado se realizará en R mediante:
-selección de taxones aceptados y con rango Species,
detección de regiones neotropicales a partir de los campos continent, region y area,
-uso de listas definidas de regiones y países neotropicales (e.g., Amazon, South America, Peru, Ecuador, Brazil).
#### 1.2 Integración con plastomas RefSeq (NCBI)
Se descargará el release completo de plastomas RefSeq.
Desde los encabezados FASTA se extraerán accession, género, especie y se construirá una llave taxonómica compatible con la del WCVP.
La unión entre ambas bases se realizará usando herramientas de Unix.
#### 1.3 Selección final de plastomas (150 total)
Se seleccionará un conjunto representativo compuesto por:
- 30 familias neotropicales,
-≥ 5 géneros por familia,
-1 plastoma completo RefSeq por género,
-150 plastomas en total.
Se priorizarán accesiones completas, curadas (NC_) y familias con alta diversidad de géneros disponibles.
### 2. Alineamiento y filogenia
#### 2.1 Alineamiento
Los 150 plastomas se alinearán con MAFFT (--auto, --adjustdirection, --reorder).
#### 2.2 Árbol filogenético
El árbol se inferirá mediante máxima verosimilitud usando:
-IQ-TREE 2 o RAxML-NG
-modelos de sustitución seleccionados automáticamente (ModelFinder)
-soporte de ramas mediante 1000 replicaciones de ultrafast bootstrap (UFBoot)
### 3. Diversidad genética
Se calcularán p-distance y π (diversidad nucleotídica) con R o Python.
Los valores intra- e inter-familia se compararán con la estructura del árbol.
# Herramientas computacionales 
## Lenguajes y entornos
-Bash — procesamiento de archivos, filtrado taxonómico (WCVP), unión de tablas, extracción de secuencias FASTA.
-R — análisis de diversidad genética (p-distance, π), manipulación de datos, visualización.
-Python (opcional) — parsing avanzado de secuencias y automatización adicional.
## Software de alineamiento y filogenia
-MAFFT — alineamiento múltiple de genomas cloroplastidiales completos.
-IQ-TREE 2 / RAxML-NG — inferencia filogenética por máxima verosimilitud, ModelFinder, UFBoot, SH-aLRT.

