library(readr)
library(dplyr)
library(stringr)
library(tibble)
library(tidyr)
library(ggplot2)
library(scales)


ruta <- "D:/MAESTRÍA/Ciclo II/CURSOS/Genómica Evolutiva/PRO/wcvp/"

## 1. Leer tablas
wcvp_names <- read_delim(file.path(ruta, "wcvp_names.csv"),
                         delim = "|",
                         escape_double = FALSE,
                         trim_ws = TRUE)

wcvp_dist <- read_delim(file.path(ruta, "wcvp_distribution.csv"),
                        delim = "|",
                        escape_double = FALSE,
                        trim_ws = TRUE)

## 2. Especies aceptadas
names_accepted <- wcvp_names %>%
  filter(
    taxon_status == "Accepted",
    taxon_rank == "Species"
  )

## 3. Palabras clave geográficas
regiones <- c(
  "Neotropics",
  "Amazon", "Amazon Basin",
  "South America", "Western South America",
  "Northern South America",
  "Central America",
  "Caribbean"
)

paises <- c(
  "Peru", "Ecuador", "Colombia", "Brazil", "Bolivia",
  "Venezuela", "Guyana", "Suriname", "French Guiana",
  "Panama", "Costa Rica", "Nicaragua", "Honduras",
  "Guatemala", "El Salvador", "Belize", "Mexico",
  "Trinidad", "Tobago"
)

patron_geo <- paste(c(regiones, paises), collapse = "|")

## 4. Filtrar distribución
dist_filtrada <- wcvp_dist %>%
  filter(
    str_detect(continent, regex("South America|Central America|Caribbean", ignore_case = TRUE)) |
      str_detect(region,   regex(patron_geo, ignore_case = TRUE)) |
      str_detect(area,     regex(patron_geo, ignore_case = TRUE))
  )

## 5. Unir taxonomía + distribución
dist_simpl <- dist_filtrada %>%
  select(plant_name_id,
         continent, region, area,
         introduced, extinct, location_doubtful)

lista_neotropico <- names_accepted %>%
  inner_join(dist_simpl, by = "plant_name_id") %>%
  distinct(plant_name_id, .keep_all = TRUE)

## 6. Guardar resultado
write_delim(
  lista_neotropico,
  "lista_neotropical_filtrada.csv",
  delim = "|"
)

## 7 Ranking preliminar de familias por n_generos_totales (WCVP)
tabla_generos_totales <- lista_neotropico %>%
  filter(!is.na(family), family != "", !is.na(genus), genus != "") %>%
  distinct(family, genus) %>%
  count(family, name = "n_generos_totales") %>%
  arrange(desc(n_generos_totales))

write_csv(
  tabla_generos_totales,
  "ranking_familias_n_generos_totales_WCVP.csv"
)

# (8) Conteos por etapa (0–4): especies y familias

# ---- Etapa 0: universo WCVP (Species, aceptadas + no aceptadas) ----
n_spec_et0 <- wcvp_names %>%
  filter(taxon_rank == "Species") %>%
  distinct(plant_name_id) %>%
  nrow()

n_fam_et0 <- wcvp_names %>%
  filter(taxon_rank == "Species", !is.na(family), family != "") %>%
  distinct(family) %>%
  nrow()

# ---- Etapa 1: especies con información de distribución (en wcvp_distribution) ----
n_spec_et1 <- wcvp_dist %>%
  distinct(plant_name_id) %>%
  nrow()

n_fam_et1 <- wcvp_dist %>%
  distinct(plant_name_id) %>%
  inner_join(
    wcvp_names %>%
      filter(taxon_rank == "Species") %>%
      select(plant_name_id, family),
    by = "plant_name_id"
  ) %>%
  filter(!is.na(family), family != "") %>%
  distinct(family) %>%
  nrow()

# ---- Etapa 2: especies aceptadas (WCVP) ----
n_spec_et2 <- names_accepted %>%
  distinct(plant_name_id) %>%
  nrow()

n_fam_et2 <- names_accepted %>%
  filter(!is.na(family), family != "") %>%
  distinct(family) %>%
  nrow()

# ---- Etapa 3: aceptadas con información de distribución ----
n_spec_et3 <- names_accepted %>%
  inner_join(wcvp_dist %>% distinct(plant_name_id), by = "plant_name_id") %>%
  distinct(plant_name_id) %>%
  nrow()

n_fam_et3 <- names_accepted %>%
  inner_join(wcvp_dist %>% distinct(plant_name_id), by = "plant_name_id") %>%
  filter(!is.na(family), family != "") %>%
  distinct(family) %>%
  nrow()

# ---- Etapa 4: aceptadas con distribución neotropical (tu lista_neotropico) ----
n_spec_et4 <- lista_neotropico %>%
  distinct(plant_name_id) %>%
  nrow()

n_fam_et4 <- lista_neotropico %>%
  filter(!is.na(family), family != "") %>%
  distinct(family) %>%
  nrow()

# ---- Tabla resumen para guardar ----
tabla_resumen_0_4 <- tibble(
  etapa = 0:4,
  descripcion = c(
    "Universo completo WCVP (Species)",
    "Especies con información de distribución (WCVP)",
    "Especies aceptadas (WCVP)",
    "Aceptadas con información de distribución",
    "Aceptadas con distribución neotropical"
  ),
  n_especies = c(n_spec_et0, n_spec_et1, n_spec_et2, n_spec_et3, n_spec_et4),
  n_familias = c(n_fam_et0, n_fam_et1, n_fam_et2, n_fam_et3, n_fam_et4)
)

write_csv(tabla_resumen_0_4, "tabla_resumen_filtrado_0_4.csv")

# (9) Figuras: reducción por etapas (especies y familias)

# ---- Gráfico 1: especies ----
p_especies_0_4 <- ggplot(tabla_resumen_0_4, aes(x = etapa, y = n_especies, group = 1)) +
  geom_line(linewidth = 1, color = "#E69F00") +
  geom_point(size = 3, color = "#E69F00") +
  geom_text(aes(label = comma(n_especies)), vjust = -0.7, size = 3.5) +
  scale_x_continuous(breaks = 0:4, expand = expansion(mult = c(0.06, 0.06))) +
  scale_y_continuous(labels = comma, expand = expansion(mult = c(0.03, 0.12))) +
  labs(
    x = "Etapas del proceso de filtrado (0–4)",
    y = "Número de especies",
    title = "Reducción progresiva del número de especies"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.title.x = element_text(size = 10, margin = margin(t = 9)),
    axis.title.y = element_text(size = 10, margin = margin(r = 9)),
    axis.line = element_line(color = "black", linewidth = 0.2),
    panel.background = element_rect(fill = "transparent", colour = NA),
    plot.background  = element_rect(fill = "transparent", colour = NA)
  )

ggsave("reduccion_especies_0_4.png", p_especies_0_4, width = 7, height = 4.5, dpi = 300, bg = "transparent")

# ---- Gráfico 2: familias ----
p_familias_0_4 <- ggplot(tabla_resumen_0_4, aes(x = etapa, y = n_familias, group = 1)) +
  geom_line(linewidth = 1, color = "gray50") +
  geom_point(size = 3, color= "gray50") +
  geom_text(aes(label = comma(n_familias)), vjust = -0.7, size = 3.5) +
  scale_x_continuous(breaks = 0:4, expand = expansion(mult = c(0.06, 0.06))) +
  scale_y_continuous(labels = comma, expand = expansion(mult = c(0.03, 0.12))) +
  labs(
    x = "Etapas del proceso de filtrado (0–4)",
    y = "Número de familias",
    title = "Reducción progresiva del número de familias"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.title.x = element_text(size = 10, margin = margin(t = 9)),
    axis.title.y = element_text(size = 10, margin = margin(r = 9)),
    axis.line = element_line(color = "black", linewidth = 0.2),
    panel.background = element_rect(fill = "transparent", colour = NA),
    plot.background  = element_rect(fill = "transparent", colour = NA)
  )

ggsave("reduccion_familias_0_4.png", p_familias_0_4, width = 7, height = 4.5, dpi = 300, bg = "transparent")

