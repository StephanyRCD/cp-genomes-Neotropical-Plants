# Descripción general:
Este cuaderno documenta el flujo seguido para construir la base de datos de plastomas neotropicales y el alineamiento de 150 genomas cloroplastidiales, a partir de:
- La lista taxonómica del World Checklist of Vascular Plants (WCVP).

- La base de genomas plastidiales de NCBI RefSeq.

- Scripts de filtrado, curación y selección almacenados en la carpeta scripts/ del repositorio.

Los pasos computacionales principales están implementados en:

- scripts/01_filtrado_wcvp.R

- scripts/02_curacion_y_seleccion_plastomas.sh

- scripts/03_alineamiento_filogenia.sh

En este notebook no se reejecuta todo el pipeline, sino que se documentan los resultados clave leyendo los archivos ya generados.
# Filtrado neotropical a partir de WCVP
Este paso está implementado en el script:

- scripts/01_filtrado_wcvp.R

De manera resumida, se:

1. Leen las tablas wcvp_names.csv y wcvp_distribution.csv (formato |-separado).

2. Se filtran únicamente especies aceptadas (taxon_status == "Accepted", taxon_rank == "Species").

3. Se seleccionan taxones con distribución en el Neotrópico mediante:

continentes: South America, Central America, Caribbean;

regiones y áreas que contienen palabras clave como Neotropics, Amazon, South America, etc.

4. Se genera el archivo lista_neotropical_filtrada.csv, que resume las especies neotropicales aceptadas.
### Carga y chequeos básicos de la lista filtrada
lista_neotrop <- read_delim(
"Dataset/lista_neotropical_filtrada.csv",
delim = "|"
)

nrow(lista_neotrop)            # número de filas (accesiones WCVP)
length(unique(lista_neotrop$plant_name_id))  # debe coincidir con nrow (1 fila por especie)
length(unique(lista_neotrop$family))        # número de familias
head(lista_neotrop %>%
select(plant_name_id, family, genus, species, geographic_area))
# Intersección WCVP–NCBI RefSeq

Esta etapa está implementada principalmente en el script:

- scripts/02_curacion_y_seleccion_plastomas.sh

Los pasos fueron:

1. Descarga de todos los plastidios RefSeq (plastid.*.genomic.fna.gz) y concatenación en plastid_refseq.fna.

2. Construcción de llaves taxonómicas:

- WCVP: genus_species o genus (cuando sp./sin epíteto).

- NCBI: idem, a partir de los encabezados FASTA (>NC_... Genus species ...).

3. Unión (join) de ambas tablas para obtener wcvp_plastid_matches.tsv, que contiene:

- información taxonómica de WCVP;

- accesión plastidial RefSeq (NC_*).
  
### Cargar la tabla WCVP–NCBI
matches <- read_tsv(
"Dataset/wcvp_plastid_matches.tsv",
col_names = c(
"key",
"plant_name_id",
"family",
"genus_wcvp",
"species_wcvp",
"taxlevel_wcvp",
"accession",
"genus_ncbi",
"species_ncbi",
"taxlevel_ncbi"
)
)

nrow(matches)      # filas = combinaciones llave–accesión
head(matches)
### Resumen rápido de la cobertura

#### Cuántas accesiones únicas (plastomas RefSeq) se vinculan a la lista neotropical

length(unique(matches$accession))

#### Cuántas especies/taxones WCVP tienen al menos un plastoma

length(unique(matches$plant_name_id))

#### Número de familias representadas

length(unique(matches$family))

# Selección de 150 plastomas (30 familias × 5 géneros)
La selección final se realizó también en:

- scripts/02_curacion_y_seleccion_plastomas.sh

Criterios aplicados:

1. Familias con al menos 5 géneros con plastoma completo.

2. Solo accesiones NC_ (RefSeq curadas).

3. Una accesión por combinación (familia, género).

4. Selección de 5 géneros por familia para construir un conjunto balanceado:

- 30 familias × 5 géneros = 150 plastomas.

El resultado se guarda en plastomas_5porfamilia.tsv.

### Carga de la tabla de 150 plastomas

plastomas150 <- read_tsv(
"Dataset/plastomas_5porfamilia.tsv",
col_names = c("family", "genus", "accession")
)

nrow(plastomas150)                  # debería ser 150
length(unique(plastomas150$family)) # debería ser 30
length(unique(plastomas150$genus))  # 150 (un género por plastoma)

head(plastomas150)

### Distribución por familia

plastomas150 %>%
count(family, name = "n_plastomes") %>%
arrange(desc(n_plastomes))

# Alineamiento de los 150 plastomas

El alineamiento se realizó en:

scripts/03_alineamiento.sh


