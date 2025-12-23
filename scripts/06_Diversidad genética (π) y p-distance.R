#!/usr/bin/env Rscript

# ============================================================
# Diversidad genética (π) y p-distance intra/inter por familia
# Control de gaps: pairwise deletion vs complete deletion
# ------------------------------------------------------------
# Inputs:
aln_fasta <- "plastomas_10x5.aln.trimal.fasta"
meta_tsv  <- "seleccion_10x5_matches.tsv"
# ============================================================

suppressPackageStartupMessages({
  library(ape)      # dist.dna
  library(dplyr)
  library(readr)
  library(tidyr)
  library(ggplot2)
})

# ----------------------------
# (1) Leer metadatos
# ----------------------------
meta <- read_tsv(meta_tsv, show_col_types = FALSE) %>%
  mutate(accession = as.character(accession),
         family    = as.character(family)) %>%
  select(accession, family) %>%
  distinct()

# ----------------------------
# (2) Leer alineamiento y mapear familias
# ----------------------------
dna <- read.dna(aln_fasta, format = "fasta")   # DNAbin
acc <- rownames(dna)

meta2 <- meta %>% filter(accession %in% acc)

# Chequeo mínimo: ¿todas las secuencias tienen familia?
missing_meta <- setdiff(acc, meta2$accession)
if (length(missing_meta) > 0) {
  warning("Estas accesiones NO tienen familia en el TSV (se excluirán):\n",
          paste(missing_meta, collapse = ", "))
}

# Mantener solo secuencias con familia
keep_acc <- meta2$accession
dna <- dna[keep_acc, , drop = FALSE]
meta2 <- meta2 %>% filter(accession %in% keep_acc)

# Vector familia alineado al orden de las secuencias en 'dna'
fam_vec <- meta2$family[match(rownames(dna), meta2$accession)]
stopifnot(!any(is.na(fam_vec)))

# ----------------------------
# (3) Funciones auxiliares
# ----------------------------
mean_lower_tri <- function(Dmat) {
  if (nrow(Dmat) < 2) return(NA_real_)
  mean(Dmat[lower.tri(Dmat)], na.rm = TRUE)
}

# Calcula distancias "raw" (= p-distance) con control de gaps:
# - pairwise deletion: pairwise.deletion=TRUE
# - complete deletion: pairwise.deletion=FALSE (sitios con missing en cualquiera se eliminan globalmente)
compute_dist_raw <- function(dna_subset, mode = c("pairwise", "complete")) {
  mode <- match.arg(mode)
  pairwise <- (mode == "pairwise")
  as.matrix(dist.dna(dna_subset, model = "raw", pairwise.deletion = pairwise, as.matrix = TRUE))
}

# ----------------------------
# (4) Intra-familia: π_intra y p-distance_intra
# ----------------------------
families <- sort(unique(fam_vec))

calc_intra <- function(mode = c("pairwise", "complete")) {
  mode <- match.arg(mode)
  
  res <- lapply(families, function(f) {
    idx <- which(fam_vec == f)
    n   <- length(idx)
    
    if (n < 2) {
      return(tibble(family = f, n_seq = n,
                    pi_intra = NA_real_, p_distance_intra = NA_real_))
    }
    
    D <- compute_dist_raw(dna[idx, , drop = FALSE], mode = mode)
    
    # En este contexto: π_intra = promedio de diferencias por sitio entre pares
    # p-distance_intra = distancia p promedio (raw) entre pares
    # (numéricamente iguales si se usa el mismo esquema de deletions)
    m <- mean_lower_tri(D)
    
    tibble(family = f, n_seq = n,
           pi_intra = m,
           p_distance_intra = m)
  })
  
  bind_rows(res) %>%
    arrange(desc(n_seq), family) %>%
    mutate(gap_mode = mode)
}

intra_pairwise <- calc_intra("pairwise")
intra_complete <- calc_intra("complete")

# ----------------------------
# (5) Inter-familia: p-distance_inter por pares de familias
# ----------------------------
calc_inter <- function(mode = c("pairwise", "complete")) {
  mode <- match.arg(mode)
  
  # Distancia global (50×50; aquí es pequeño y eficiente)
  Dall <- compute_dist_raw(dna, mode = mode)
  
  # Tabla de pares familia-familia (i<j)
  fam_levels <- sort(unique(fam_vec))
  pairs <- t(combn(fam_levels, 2))
  pairs_df <- as_tibble(pairs, .name_repair = "minimal") %>%
    setNames(c("family_1", "family_2"))
  
  inter <- pairs_df %>%
    rowwise() %>%
    mutate(
      n1 = sum(fam_vec == family_1),
      n2 = sum(fam_vec == family_2),
      p_distance_inter = {
        i <- which(fam_vec == family_1)
        j <- which(fam_vec == family_2)
        mean(Dall[i, j], na.rm = TRUE)
      }
    ) %>%
    ungroup() %>%
    mutate(gap_mode = mode) %>%
    arrange(desc(p_distance_inter))
  
  inter
}

inter_pairwise <- calc_inter("pairwise")
inter_complete <- calc_inter("complete")

# ----------------------------
# (6) (Opcional) π_total del conjunto completo
# ----------------------------
pi_total <- function(mode = c("pairwise", "complete")) {
  mode <- match.arg(mode)
  Dall <- compute_dist_raw(dna, mode = mode)
  mean_lower_tri(Dall)
}

pi_total_pairwise <- pi_total("pairwise")
pi_total_complete <- pi_total("complete")

# ----------------------------
# (7) Exportar resultados
# ----------------------------
write_csv(intra_pairwise,  "pi_pdist_intra_pairwise.csv")
write_csv(intra_complete,  "pi_pdist_intra_complete.csv")
write_csv(inter_pairwise,  "pdist_inter_pairwise.csv")
write_csv(inter_complete,  "pdist_inter_complete.csv")

summary_tbl <- tibble(
  metric = c("pi_total"),
  pairwise_deletion = c(pi_total_pairwise),
  complete_deletion = c(pi_total_complete)
)
write_csv(summary_tbl, "pi_total_summary.csv")

# ----------------------------
# (8) Gráficos
# ----------------------------

p_intra <- ggplot(intra_pairwise,
                  aes(x = reorder(family, pi_intra),
                      y = pi_intra)) +
  geom_col(width = 0.75, fill = "#4C72B0") +
  coord_flip() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
  geom_text(
    aes(label = sprintf("%.3f", pi_intra)),
    hjust = -0.1,
    size = 3.2
  ) +
  labs(
    x = "Familia",
    y = expression(pi~"(diversidad nucleotídica intra-familia)"),
    title = "Diversidad genética intra-familia (π)"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.title.x = element_text(size = 10, margin = margin(t = 8)),
    axis.title.y = element_text(size = 10, margin = margin(r = 8)),
    axis.line = element_line(color = "black", linewidth = 0.20),
    panel.background = element_rect(fill = "transparent", colour = NA),
    plot.background  = element_rect(fill = "transparent", colour = NA)
  )

ggsave(
  plot = p_intra,
  filename = "diversidad_intra_familia_pi.png",
  width = 18,
  height = 12,
  units = "cm",
  dpi = 300,
  bg = "transparent"
)

p_inter <- ggplot(inter_pairwise,
       aes(x = family_1, y = family_2, fill = p_distance_inter)) +
  geom_tile(color = "white") +
  scale_fill_viridis_c(name = "p-distance") +
  labs(
    x = "Familia",
    y = "Familia",
    title = "Distancia genética promedio entre familias"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1), 
    axis.title.x = element_text(size = 10, margin = margin(t = 8)),
    axis.title.y = element_text(size = 10, margin = margin(r = 8)),
    )

ggsave(
  plot = p_inter,
  filename = "distancia_genetica_inter_familias.png",
  width = 18,
  height = 14,
  units = "cm",
  dpi = 300,
  bg = "white"
)
