rm(list = ls())
.rs.restartR(clean = T)
.libPaths()

# renv::install(packages = c("ggprism", "openxlsx", "qs2"), prompt = F)
# renv::install("R.utils", prompt = F)

# Outs directory
outs_dir <- "01_core_mapping_outs"
dir.create(path = outs_dir)

# Need the core to patient map
library(magrittr)
coremap_defacto <- openxlsx::read.xlsx(xlsxFile = "../../sgroi-tnbc-data/sgroi_tma_maps.xlsx", sheet = 6, cols = 9:13)
coremap_defacto$patient <- gsub(pattern = "([a-z])|([A-Z])|( )", replacement = "", x = coremap_defacto$PATH.ID)
coremap_defacto[coremap_defacto$patient == "",]$patient <- "control"
coremap_defacto[grepl(pattern = "^[0-9]", x = coremap_defacto$patient),]$patient %<>% sprintf(fmt = "%02s") %<>% paste("p", ., sep = "")


## tma1 -----------------------------------------------------------------------

tma1 <- data.table::fread("../../sgroi-tnbc-data/tma1/stow/cells.csv.gz") |>
  as.data.frame()
tma1$segmentation_method <- gsub(pattern = "µ", replacement = "u", x = tma1$segmentation_method)
rownames(tma1) <- tma1$cell_id
tma1$run <- 1
tma1$slide <- "tma1"

tma1_cores <- sf::read_sf("../../sgroi-tnbc-data/tma1_cores_BF.geojson") |> dplyr::select(name, geometry)
sf::st_crs(tma1_cores) <- NA

library(ggplot2)
ggplot() +
  geom_sf(data = tma1_cores, mapping = aes(geometry = geometry), fill = "pink", color = "black") +
  scattermore::geom_scattermore(data = tma1, mapping = aes(x = x_centroid/0.2125, y = y_centroid/0.2125)) +
  theme_void()
ggsave(filename = file.path(outs_dir, "tma1_cores.pdf"), width = 6, height = 8)

cell_pts <- sf::st_as_sf(x = tma1 |> dplyr::mutate(X = x_centroid/0.2125, Y = y_centroid/0.2125), coords = c("X", "Y"))
cellstatus <- sf::st_within(x = cell_pts, y = tma1_cores, sparse = T)
cellstatus <- as.vector(cellstatus)
cellstatus[(lengths(cellstatus) == 0)] <- NA
coreids <- plyr::mapvalues(x = unlist(cellstatus), from = rownames(tma1_cores), to = tma1_cores$name)
tma1$core <- coreids

ggplot() +
  geom_sf(data = tma1_cores, mapping = aes(geometry = geometry), fill = NA, color = "black") +
  scattermore::geom_scattermore(data = tma1, mapping = aes(x = x_centroid/0.2125, y = y_centroid/0.2125, color = core)) +
  geom_text(data = tma1 |> dplyr::group_by(core) |> dplyr::summarise(X = median(x_centroid/0.2125)+1500, Y = max(y_centroid/0.2125)+1500),
            mapping = aes(x = X, y = Y, label = core)) +
  theme_void() +
  ggprism::scale_color_prism() +
  Seurat::NoLegend()
ggsave(filename = file.path(outs_dir, "tma1_cores_labeled.pdf"), width = 6, height = 8)

tma1 <- tma1[!is.na(tma1$core),]
tma1$patient <- plyr::mapvalues(x = paste(tma1$slide, tma1$core, sep = ""), from = paste("tma", coremap_defacto$`TMA.#`, coremap_defacto$POS, sep = ""), to = coremap_defacto$patient)

ggplot() +
  geom_sf(data = tma1_cores, mapping = aes(geometry = geometry), fill = NA, color = "black") +
  scattermore::geom_scattermore(data = tma1, mapping = aes(x = x_centroid/0.2125, y = y_centroid/0.2125, color = patient)) +
  geom_text(data = tma1 |> dplyr::group_by(core) |> dplyr::summarise(X = median(x_centroid/0.2125)+1500, Y = max(y_centroid/0.2125)+1500, patient = unique(patient)),
            mapping = aes(x = X, y = Y, label = patient)) +
  theme_void() +
  ggprism::scale_color_prism() +
  Seurat::NoLegend()
ggsave(filename = file.path(outs_dir, "tma1_cores_by_patient.pdf"), width = 6, height = 8)


## tma2 -----------------------------------------------------------------------

tma2 <- data.table::fread("../../sgroi-tnbc-data/tma2/stow/cells.csv.gz") |>
  as.data.frame()
tma2$segmentation_method <- gsub(pattern = "µ", replacement = "u", x = tma2$segmentation_method)
rownames(tma2) <- tma2$cell_id
tma2$run <- 1
tma2$slide <- "tma2"

tma2_cores <- sf::read_sf("../../sgroi-tnbc-data/tma2_cores_BF.geojson") |> dplyr::select(name, geometry)
sf::st_crs(tma2_cores) <- NA

library(ggplot2)
ggplot() +
  geom_sf(data = tma2_cores, mapping = aes(geometry = geometry), fill = "pink", color = "black") +
  scattermore::geom_scattermore(data = tma2, mapping = aes(x = x_centroid/0.2125, y = y_centroid/0.2125)) +
  theme_void()
ggsave(filename = file.path(outs_dir, "tma2_cores.pdf"), width = 6, height = 8)

cell_pts <- sf::st_as_sf(x = tma2 |> dplyr::mutate(X = x_centroid/0.2125, Y = y_centroid/0.2125), coords = c("X", "Y"))
cellstatus <- sf::st_within(x = cell_pts, y = tma2_cores, sparse = T)
cellstatus <- as.vector(cellstatus)
cellstatus[(lengths(cellstatus) == 0)] <- NA
coreids <- plyr::mapvalues(x = unlist(cellstatus), from = rownames(tma2_cores), to = tma2_cores$name)
tma2$core <- coreids

ggplot() +
  geom_sf(data = tma2_cores, mapping = aes(geometry = geometry), fill = NA, color = "black") +
  scattermore::geom_scattermore(data = tma2, mapping = aes(x = x_centroid/0.2125, y = y_centroid/0.2125, color = core)) +
  geom_text(data = tma2 |> dplyr::group_by(core) |> dplyr::summarise(X = median(x_centroid/0.2125)+1500, Y = max(y_centroid/0.2125)+1500),
            mapping = aes(x = X, y = Y, label = core)) +
  theme_void() +
  ggprism::scale_color_prism() +
  Seurat::NoLegend()
ggsave(filename = file.path(outs_dir, "tma2_cores_labeled.pdf"), width = 6, height = 8)

tma2 <- tma2[!is.na(tma2$core),]
tma2$patient <- plyr::mapvalues(x = paste(tma2$slide, tma2$core, sep = ""), from = paste("tma", coremap_defacto$`TMA.#`, coremap_defacto$POS, sep = ""), to = coremap_defacto$patient)

ggplot() +
  geom_sf(data = tma2_cores, mapping = aes(geometry = geometry), fill = NA, color = "black") +
  scattermore::geom_scattermore(data = tma2, mapping = aes(x = x_centroid/0.2125, y = y_centroid/0.2125, color = patient)) +
  geom_text(data = tma2 |> dplyr::group_by(core) |> dplyr::summarise(X = median(x_centroid/0.2125)+1500, Y = max(y_centroid/0.2125)+1500, patient = unique(patient)),
            mapping = aes(x = X, y = Y, label = patient)) +
  theme_void() +
  ggprism::scale_color_prism() +
  Seurat::NoLegend()
ggsave(filename = file.path(outs_dir, "tma2_cores_by_patient.pdf"), width = 6, height = 8)


## tma3 -----------------------------------------------------------------------

tma3 <- data.table::fread("../../sgroi-tnbc-data/tma3/stow/cells.csv.gz") |>
  as.data.frame()
tma3$segmentation_method <- gsub(pattern = "µ", replacement = "u", x = tma3$segmentation_method)
rownames(tma3) <- tma3$cell_id
tma3$run <- 2
tma3$slide <- "tma3"

tma3_cores <- sf::read_sf("../../sgroi-tnbc-data/tma3_cores_BF.geojson") |> dplyr::select(name, geometry)
sf::st_crs(tma3_cores) <- NA

library(ggplot2)
ggplot() +
  geom_sf(data = tma3_cores, mapping = aes(geometry = geometry), fill = "pink", color = "black") +
  scattermore::geom_scattermore(data = tma3, mapping = aes(x = x_centroid/0.2125, y = y_centroid/0.2125)) +
  theme_void()
ggsave(filename = file.path(outs_dir, "tma3_cores.pdf"), width = 6, height = 8)

cell_pts <- sf::st_as_sf(x = tma3 |> dplyr::mutate(X = x_centroid/0.2125, Y = y_centroid/0.2125), coords = c("X", "Y"))
cellstatus <- sf::st_within(x = cell_pts, y = tma3_cores, sparse = T)
cellstatus <- as.vector(cellstatus)
cellstatus[(lengths(cellstatus) == 0)] <- NA
coreids <- plyr::mapvalues(x = unlist(cellstatus), from = rownames(tma3_cores), to = tma3_cores$name)
tma3$core <- coreids

ggplot() +
  geom_sf(data = tma3_cores, mapping = aes(geometry = geometry), fill = NA, color = "black") +
  scattermore::geom_scattermore(data = tma3, mapping = aes(x = x_centroid/0.2125, y = y_centroid/0.2125, color = core)) +
  geom_text(data = tma3 |> dplyr::group_by(core) |> dplyr::summarise(X = median(x_centroid/0.2125)+1500, Y = max(y_centroid/0.2125)+1500),
            mapping = aes(x = X, y = Y, label = core)) +
  theme_void() +
  ggprism::scale_color_prism() +
  Seurat::NoLegend()
ggsave(filename = file.path(outs_dir, "tma3_cores_labeled.pdf"), width = 6, height = 8)

tma3 <- tma3[!is.na(tma3$core),]
tma3$patient <- plyr::mapvalues(x = paste(tma3$slide, tma3$core, sep = ""), from = paste("tma", coremap_defacto$`TMA.#`, coremap_defacto$POS, sep = ""), to = coremap_defacto$patient)

ggplot() +
  geom_sf(data = tma3_cores, mapping = aes(geometry = geometry), fill = NA, color = "black") +
  scattermore::geom_scattermore(data = tma3, mapping = aes(x = x_centroid/0.2125, y = y_centroid/0.2125, color = patient)) +
  geom_text(data = tma3 |> dplyr::group_by(core) |> dplyr::summarise(X = median(x_centroid/0.2125)+1500, Y = max(y_centroid/0.2125)+1500, patient = unique(patient)),
            mapping = aes(x = X, y = Y, label = patient)) +
  theme_void() +
  ggprism::scale_color_prism() +
  Seurat::NoLegend()
ggsave(filename = file.path(outs_dir, "tma3_cores_by_patient.pdf"), width = 6, height = 8)


## tma4 -----------------------------------------------------------------------

tma4 <- data.table::fread("../../sgroi-tnbc-data/tma4/stow/cells.csv.gz") |>
  as.data.frame()
tma4$segmentation_method <- gsub(pattern = "µ", replacement = "u", x = tma4$segmentation_method)
rownames(tma4) <- tma4$cell_id
tma4$run <- 2
tma4$slide <- "tma4"

tma4_cores <- sf::read_sf("../../sgroi-tnbc-data/tma4_cores_BF.geojson") |> dplyr::select(name, geometry)
sf::st_crs(tma4_cores) <- NA

library(ggplot2)
ggplot() +
  geom_sf(data = tma4_cores, mapping = aes(geometry = geometry), fill = "pink", color = "black") +
  scattermore::geom_scattermore(data = tma4, mapping = aes(x = x_centroid/0.2125, y = y_centroid/0.2125)) +
  theme_void()
ggsave(filename = file.path(outs_dir, "tma4_cores.pdf"), width = 6, height = 8)

cell_pts <- sf::st_as_sf(x = tma4 |> dplyr::mutate(X = x_centroid/0.2125, Y = y_centroid/0.2125), coords = c("X", "Y"))
cellstatus <- sf::st_within(x = cell_pts, y = tma4_cores, sparse = T)
cellstatus <- as.vector(cellstatus)
cellstatus[(lengths(cellstatus) == 0)] <- NA
coreids <- plyr::mapvalues(x = unlist(cellstatus), from = rownames(tma4_cores), to = tma4_cores$name)
tma4$core <- coreids

ggplot() +
  geom_sf(data = tma4_cores, mapping = aes(geometry = geometry), fill = NA, color = "black") +
  scattermore::geom_scattermore(data = tma4, mapping = aes(x = x_centroid/0.2125, y = y_centroid/0.2125, color = core)) +
  geom_text(data = tma4 |> dplyr::group_by(core) |> dplyr::summarise(X = median(x_centroid/0.2125)+1500, Y = max(y_centroid/0.2125)+1500),
            mapping = aes(x = X, y = Y, label = core)) +
  theme_void() +
  ggprism::scale_color_prism() +
  Seurat::NoLegend()
ggsave(filename = file.path(outs_dir, "tma4_cores_labeled.pdf"), width = 6, height = 8)

tma4 <- tma4[!is.na(tma4$core),]
tma4$patient <- plyr::mapvalues(x = paste(tma4$slide, tma4$core, sep = ""), from = paste("tma", coremap_defacto$`TMA.#`, coremap_defacto$POS, sep = ""), to = coremap_defacto$patient)

# ** Need a little manual correction here **
# According to my notes in the orientations xlsx file:
# G2 = TNBC 28 and G1 = TNBC 34
tma4[tma4$core == "G2",]$patient <- "p28"
tma4[tma4$core == "G1",]$patient <- "p34"

ggplot() +
  geom_sf(data = tma4_cores, mapping = aes(geometry = geometry), fill = NA, color = "black") +
  scattermore::geom_scattermore(data = tma4, mapping = aes(x = x_centroid/0.2125, y = y_centroid/0.2125, color = patient)) +
  geom_text(data = tma4 |> dplyr::group_by(core) |> dplyr::summarise(X = median(x_centroid/0.2125)+1500, Y = max(y_centroid/0.2125)+1500, patient = unique(patient)),
            mapping = aes(x = X, y = Y, label = patient)) +
  theme_void() +
  ggprism::scale_color_prism() +
  Seurat::NoLegend()
ggsave(filename = file.path(outs_dir, "tma4_cores_by_patient.pdf"), width = 6, height = 8)


## Saving ----------------------------------------------------------------------
qs2::qs_save(
  object = list(
    "tma1" = tma1, 
    "tma2" = tma2, 
    "tma3" = tma3, 
    "tma4" = tma4
  ), 
  file = file.path(outs_dir, "core-mapped_metadata.qs2")
)

