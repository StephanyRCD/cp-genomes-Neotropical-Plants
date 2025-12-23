# ============================================================
# SELECCIÓN BALANCEADA 10 × 5 (N=50)
# - filtro por n_generos_con_plastoma >= 5 (post-match)
# - Top 10 familias
# - 5 accesiones por familia, 1 por género
# ============================================================

library(readr)
library(dplyr)
library(stringr)

# -------------------------
# (0) Parámetros
# -------------------------
N_FAMILIAS    <- 10
N_POR_FAMILIA <- 5
SEED          <- 123   #hace que el muestreo aleatorio sea reproducible (si se corre el script otra vez, sale lo mismo).
set.seed(SEED)

IN_RANK_WCVP <- "ranking_familias_n_generos_totales_WCVP.csv"
IN_MATCHES   <- "wcvp_plastid_matches.tsv"

OUT_TABLA_MAESTRA <- "tabla_maestra_familias_WCVP_vs_plastoma.csv"
OUT_ELEGIBLES     <- "familias_elegibles_ge10_generos_con_plastoma.csv"
OUT_TOP10         <- "familias_top10_final.csv"

OUT_SEL_CSV <- "seleccion_10x5_accessions.csv"
OUT_SEL_TXT <- "seleccion_10x5_accessions.txt"
OUT_SEL_TSV <- "seleccion_10x5_matches.tsv"

# -------------------------
# (1) Leer ranking pre-match (WCVP)
# -------------------------
rank_wcvp <- read_csv(IN_RANK_WCVP, show_col_types = FALSE) %>%
  mutate(family = str_trim(family))

if (!all(c("family", "n_generos_totales") %in% names(rank_wcvp))) {
  stop("El ranking WCVP debe tener columnas: family, n_generos_totales")
}

# -------------------------
# (2) Leer matches (WCVP × RefSeq) — SIN order
# -------------------------
matches <- read_tsv(IN_MATCHES, col_names = FALSE, show_col_types = FALSE)

# Verificación: deben ser 10 columnas
if (ncol(matches) != 10) {
  stop("Tu wcvp_plastid_matches.tsv no tiene 10 columnas. Tiene: ", ncol(matches),
       ". Ajusta el colnames() según el formato real.")
}

colnames(matches) <- c(
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

matches2 <- matches %>%
  mutate(
    family     = str_trim(family),
    genus_wcvp  = str_trim(genus_wcvp),
    accession  = str_trim(accession),
    is_NC      = str_starts(accession, "NC_")
  ) %>%
  filter(
    !is.na(family), family != "",
    !is.na(genus_wcvp), genus_wcvp != "",
    !is.na(accession), accession != ""
  )

# -------------------------
# (3) Disponibilidad por familia (post-match)
# -------------------------
disp_fam <- matches2 %>%
  distinct(family, genus_wcvp) %>%
  count(family, name = "n_generos_con_plastoma") 


disp_extra <- matches2 %>%
  distinct(family, accession, plant_name_id) %>%
  group_by(family) %>%
  summarise(
    n_plastomas = n_distinct(accession),
    n_especies_con_plastoma = n_distinct(plant_name_id),
    .groups = "drop"
  )

# -------------------------
# (4) Tabla maestra (WCVP + RefSeq)
# -------------------------
tabla_maestra <- rank_wcvp %>%
  left_join(disp_fam, by = "family") %>%
  left_join(disp_extra, by = "family") %>%
  mutate(
    n_generos_con_plastoma = coalesce(n_generos_con_plastoma, 0L),
    n_plastomas = coalesce(n_plastomas, 0L),
    n_especies_con_plastoma = coalesce(n_especies_con_plastoma, 0L)
  ) %>%
  arrange(desc(n_generos_totales), desc(n_generos_con_plastoma))

write_csv(tabla_maestra, OUT_TABLA_MAESTRA)

# -------------------------
# (5) Elegibles: n_generos_con_plastoma >= 5
# -------------------------
elegibles <- tabla_maestra %>%
  filter(n_generos_con_plastoma >= N_POR_FAMILIA) %>%
  arrange(desc(n_generos_totales), desc(n_generos_con_plastoma))

write_csv(elegibles, OUT_ELEGIBLES)

if (nrow(elegibles) < N_FAMILIAS) {
  stop("No hay suficientes familias elegibles con n_generos_con_plastoma ≥ ",
       N_POR_FAMILIA, ". Elegibles: ", nrow(elegibles))
}

# Primero se prioriza familias más diversas en el universo neotropical; entre ellas, se prefiere las que además tienen más genomas disponibles.”
# -------------------------
# (6) Seleccionar 10 familias (robusto): tomar las primeras del ranking elegible
#     que realmente permitan 5 géneros distintos tras priorizar 1/genus
# -------------------------
# Precomputar cuántos géneros únicos tiene cada familia en matches2
gen_por_fam <- matches2 %>%
  distinct(family, genus_wcvp) %>%
  count(family, name = "n_generos_disponibles_en_matches")

elegibles2 <- elegibles %>%
  left_join(gen_por_fam, by = "family") %>%
  mutate(n_generos_disponibles_en_matches = coalesce(n_generos_disponibles_en_matches, 0L)) %>%
  filter(n_generos_disponibles_en_matches >= N_POR_FAMILIA) %>%
  arrange(desc(n_generos_totales), desc(n_generos_con_plastoma))

if (nrow(elegibles2) < N_FAMILIAS) {
  stop("Tras verificar disponibilidad real en matches, no hay suficientes familias con ≥",
       N_POR_FAMILIA, " géneros distintos. Candidatas: ", nrow(elegibles2))
}

top10 <- elegibles2 %>% slice_head(n = N_FAMILIAS)
write_csv(top10, OUT_TOP10)

# -------------------------
# (7) Selección 5 accesiones por familia, garantizando 1 por género
#     (prioriza NC_ dentro de cada género)
# -------------------------
pool <- matches2 %>%
  semi_join(top10 %>% select(family), by = "family")

# 1 accession por (familia, género): elegir la “mejor” priorizando NC_
uno_por_genero <- pool %>%
  arrange(family, genus_wcvp, desc(is_NC), accession) %>%
  group_by(family, genus_wcvp) %>%
  slice_head(n = 1) %>%
  ungroup()

# Selección final: 5 géneros por familia (aleatoria reproducible por SEED)
seleccion_10x5 <- uno_por_genero %>%
  group_by(family) %>%
  slice_sample(n = N_POR_FAMILIA) %>%
  ungroup() %>%
  select(family, genus_wcvp, accession) %>%
  arrange(family, genus_wcvp)

# -------------------------
# (8) Guardar outputs
# -------------------------
write_csv(seleccion_10x5, OUT_SEL_CSV)
write_lines(seleccion_10x5$accession, OUT_SEL_TXT)

matches_sel <- pool %>%
  semi_join(seleccion_10x5, by = c("family", "genus_wcvp", "accession"))

write_tsv(matches_sel, OUT_SEL_TSV)

# Resumen
resumen <- seleccion_10x5 %>%
  summarise(
    n_plastomas = n(),
    n_familias  = n_distinct(family),
    n_generos   = n_distinct(genus_wcvp),
    n_NC        = sum(str_starts(accession, "NC_")),
    .groups = "drop"
  )

print(resumen)
