# Información de los datasets utilizados en este proyecto
Este repositorio utiliza datasets externos que, por su tamaño, no pueden ser incluidos directamente.
A continuación se describe cada dataset, su contenido, origen y cómo acceder a él.
## 1. World Checklist of Vascular Plants (WCVP)
Esta base de datos onstituye una de las fuentes taxonómicas más autorizadas y exhaustivas para la flora mundial, desarrollada y mantenida por el Royal Botanic Gardens, Kew. Su principal fortaleza radica en ofrecer un marco taxonómico globalmente estandarizado, con información actualizada sobre nombres aceptados, sinonimias y distribución geográfica de las plantas vasculares. 

**Archivos:**
- `wcvp_names.csv`
- `wcvp_distribution.csv`
- `README_WCVP`
  
**Uso en el proyecto:**  
Filtrado de plantas neotropicales y construcción de la lista taxonómica base.

**Contenido:**
- *wcvp_names.csv*: Información taxonómica completa (género, especie, familia, estatus, sinónimos, etc.).
- *wcvp_distribution.csv*: Distribución geográfica detallada (continente, región, área, países).
- *README_WCVP*: Documentación oficial del dataset.
  
**Acceso:**
- Página oficial: https://powo.science.kew.org/
  
## 2. NCBI RefSeq Plastid Genomes

**Archivo utilizado:**
- `plastid.*.genomic.fna.gz` (colección completa)
  
**Uso en el proyecto:**  
Construcción de la base de datos de genomas cloroplastidiales completos.

**Contenido:**
- Secuencias FASTA curadas del conjunto RefSeq para plastidios.
- Incluye cloroplastos de plantas, pero también otros plastidios (apicoplastos, cromatóforos, etc.), por lo que se realizó un filtrado manual.

**Acceso:**
- FTP oficial:  
  ftp://ftp.ncbi.nlm.nih.gov/refseq/release/plastid/
