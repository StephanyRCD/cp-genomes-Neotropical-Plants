#Descarga de DB Kew GardenS (WCVP)
#Generación de una nueva lista filtrada por términos como Neotrópico o Sudamérica o países pertenecientes a estas regiones (en R)
# En bash:
cut -d '|' -f1 wcvp_names.csv | sort -u | wc -l #-1
cut -d '|' -f2 wcvp_distribution.csv | sort -u | wc -l   #-1
cut -d '|' -f1 lista_neotropical_filtrada.csv | sort -u | wc -l   #-1
#names (1440076) distribution (443347) filtrado (127743)
#Descarga de la DB del NCBI
wget ftp://ftp.ncbi.nlm.nih.gov/refseq/release/plastid/plastid.*.genomic.fna.gz #ojo que NO contiene solamente plastidios de plantas, sino todo lo que NCBI considera "plastid" (apicoplastos (protistas como Plasmodium, Babesia, Eimeria), cromatóforos de protistas, plastidios secundarios de organismos microscópicos, incluso algunos ensamblajes “plastid-like” de microbios)
#Descomprimir
gunzip plastid*.gz
#Concatenación en un FASTA único
cat plastid*.fna > plastid_refseq.fna
#Llave genus_species para la lista WCVP filtrada
tail -n +2 lista_neotropical_filtrada.csv \
  | awk -F'|' '{
      genus  = tolower($7)
      species= tolower($9)

      if (genus == "") next

      # Caso 1: especie definida
      if (species != "" && species != "sp." && species != "sp") {
          key = genus"_"species
          taxlevel="species"
      }

      # Caso 2: especie no definida → usar solo género
      else {
          key = genus
          taxlevel="genus"
      }

      print key "\t" $1 "\t" $5 "\t" $7 "\t" $9 "\t" taxlevel
  }' \
  > wcvp_neotrop_keys.tsv
sort -t $'\t' -k1,1 wcvp_neotrop_keys.tsv > wcvp_neotrop_keys.sorted.tsv
#Llaves para los plastomas (FASTA RefSeq)
grep "^>" plastid_refseq.fna > plastid_headers.txt
awk '
{
    sub(/^>/,"",$0)              # Quitar >
    acc = $1                     # Accession
    genus = tolower($2)          # Segundo campo
    species = tolower($3)        # Tercer campo

    # CASO 0: no hay genus → descartar
    if (genus == "" || genus ~ /^[0-9]+$/ ) next

    # CASO 1: especie definida
    if (species != "" && species != "sp" && species != "sp.") {
        key = genus "_" species
        taxlevel = "species"
    }
    # CASO 2: especie ausente o "sp"
    else {
        key = genus
        taxlevel = "genus"
    }

    print key "\t" acc "\t" genus "\t" species "\t" taxlevel
}' plastid_headers.txt \
> plastid_keys.tsv
sort -t $'\t' -k1,1 plastid_keys.tsv > plastid_keys.sorted.tsv
cut -f2 plastid_keys.tsv | cut -c1-2 | sort | uniq -c #N° de accesiones por tipo (NC, NW, etc.) #NC:15225, NW:9
#Hacer el match WCVP ↔ NCBI
join -t $'\t' -1 1 -2 1 \
  wcvp_neotrop_keys.sorted.tsv \
  ncbi/plastid_keys.sorted.tsv \
  > wcvp_plastid_matches.tsv #contiene: key, plant_name_id (WCVP), family (WCVP), genus_wcvp, species_wcvp, taxlevel_wcvp, accession (NCBI, NC_/NZ_), genus_ncbi, species_ncbi, taxlevel_ncbi
cut -f7 wcvp_plastid_matches.tsv \ #N° de accesiones por tipo #NC:3151, NW:0
    | cut -c1-2 \
    | sort \
    | uniq -c
#Generar tabla de cuántos plastomas tiene cada taxón
cut -f1,7 wcvp_plastid_matches.tsv \
    | sort -u \
    | cut -f1 \
    | sort \
    | uniq -c \
    > plastomas_por_taxon.tsv
#Taxones con más de un plastoma
awk '$1 > 1' plastomas_por_taxon.tsv > taxones_con_plastomas_multiples.tsv
#Número de plastomas por familia
awk -F'\t' '{print $3 "\t" $7}' wcvp_plastid_matches.tsv \
  | sort -u \
  | awk -F'\t' '{ counts[$1]++ } END { for (fam in counts) print counts[fam] "\t" fam }' \
  | sort -nr \
  > plastomas_por_familia.tsv
#Número de familias que tienen al menos un plastoma
cut -f2 plastomas_por_familia.tsv | wc -l  #219
#Número de taxones (plant_name_id) por familia
awk -F'\t' '{print $3 "\t" $2}' wcvp_plastid_matches.tsv \
  | sort -u \
  | awk -F'\t' '{tax[$1]++} END {for (fam in tax) print tax[fam] "\t" fam}' \
  | sort -nr \
  > taxones_por_familia.tsv
#Tabla integrada de cobertura por familia (plastomas disponibles vs. taxones neotropicales)
join -t $'\t' -1 2 -2 2 \
  <(sort -t $'\t' -k2,2 plastomas_por_familia.tsv) \
  <(sort -t $'\t' -k2,2 taxones_por_familia.tsv) \
  > plastomas_y_taxones_por_familia.tsv
#Familias que tienen ≥5 plastomas provenientes de géneros distintos
##Extraer familia + género + accesión
awk -F'\t' '{print $3 "\t" $8 "\t" $7}' wcvp_plastid_matches.tsv \
    > fam_genus_acc.tsv
##Contar géneros distintos por familia
cut -f1,2 fam_genus_acc.tsv \
    | sort -u \
    | awk -F'\t' '{count[$1]++} END {for (f in count) print count[f] "\t" f}' \
    | sort -nr \
    > familias_con_generos.tsv
##Filtrar solo familias con ≥5 géneros
awk '$1 >= 5' familias_con_generos.tsv \
    > familias_5_generos.tsv
wc -l familias_5_generos.tsv #N° de familias que tienen ≥5 plastomas provenientes de géneros distintos  #56
##Visualizar y guardar la lista de las primeras 30 familias
column -t familias_5_generos.tsv | head -30  #visualizar
awk '{print $2}' familias_5_generos.tsv | head -30 > familias_seleccionadas.txt #guardar lista
##Asegurar una sola accesión por (familia, género)
sort -u -k1,1 -k2,2 fam_genus_acc.tsv > fam_genus_unique.tsv #Esto deja, por cada combinación family + genus, la primera accesión que aparezca.
sort -k1,1 fam_genus_unique.tsv > fam_genus_unique.sorted.tsv
##Seleccionar hasta 5 plastomas por familia (géneros distintos)
awk 'NR==FNR {sel[$1]=1; next}
     sel[$1] {
         c[$1]++
         if (c[$1] <= 5) print $0
     }' familias_seleccionadas.txt fam_genus_unique.sorted.tsv \
  > plastomas_5porfamilia.tsv
wc -l plastomas_5porfamilia.tsv #Debería salir 150
cut -f1 plastomas_5porfamilia.tsv | sort | uniq -c | sort -nr | head #Para verificar la distribución
#Base FASTA de 150 plastomas
##Sacar la lista de accesiones a un archivo
cut -f3 plastomas_5porfamilia.tsv > acc_plastomas_150.txt  # Columna 3 = accession
head acc_plastomas_150.txt  #revisando
wc -l acc_plastomas_150.txt #debería dar 150
##Extraer solo esas secuencias del FASTA grande (plastid_refseq.fna)
awk 'NR==FNR {
        acc[$1]=1
        next
     }
     /^>/ {
        split($1, a, ">")
        id = a[2]
        keep = (id in acc)
     }
     keep' acc_plastomas_150.txt ncbi/plastid_refseq.fna \
  > plastomas_150_neotrop.fna
grep -c "^>" plastomas_150_neotrop.fna #verificando
head plastomas_150_neotrop.fna #verificando
