#!/bin/bash
# 03_alineamiento.sh
# Alineamiento de los 150 plastomas seleccionados usando MAFFT.
# Produce: plastomas_150_neotrop_aln.fna
## 1. Crear carpeta de trabajo para alineamientos
mkdir -p alineamientos
## 2. Copiar el FASTA de plastomas seleccionados
cd alineamientos
cp ../plastomas_150_neotrop.fna .
## 3. Ejecutar MAFFT sobre los 150 plastomas
#    --auto            : MAFFT elige automáticamente la mejor estrategia
#    --adjustdirection : corrige orientación inversa cuando sea necesario
#    --reorder         : reordena secuencias si es necesario para el alineamiento
#    --thread 8        : usa 8 hilos (ajustar según recursos disponibles)
mafft --auto \
      --adjustdirection \
      --reorder \
      --thread 8 \
      plastomas_150_neotrop.fna \
  > plastomas_150_neotrop_aln.fna
