#!/bin/bash
# ============================================================
# 02_curacion_plastomas.sh
# Objetivo:
#  - generar llaves taxonómicas
#  - intersectar WCVP ↔ NCBI RefSeq plastid
# ============================================================
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
