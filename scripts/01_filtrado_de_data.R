library(readr)
library(dplyr)
library(stringr)

ruta <- "D:/MAESTRÍA/Ciclo II/CURSOS/Genómica Evolutiva/DATA/DATA/wcvp"

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
