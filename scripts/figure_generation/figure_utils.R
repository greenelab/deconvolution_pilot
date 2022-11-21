library(scales)

alpha_value <- 0.5


theme_set(theme_bw() +
            theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
                  text = element_text(size = 14),
                  strip.background = element_rect(colour = "black",
                                                fill = "white"),
                  plot.title = element_text(hjust = 0.5)
            )
          )

demultiplexing_colors <- c("red" = "#CA0020",
                           "blue" = "#0571B0",
                           "green" = "#4DAC26",
                           "orange" = "#E66101",
                           "gray" = "#999999") #ggsci D3?


colors_samples <- hue_pal()(8) #ggplot defaults

colors_celltypes <- hue_pal()(13) #ggplot defaults

colors_methods <- c("bayesprism" = "#E6AB02", "bisque" = "#66A61E", "cibersortx" = "#E7298A",
                    "epic" = "#1B9E77", "music" = "#7570B3", "nnls" = "#D95F02") #Dark2

colors_bulktypes <- c("rRNA- Chunk" = "#00468B", "rRNA- Dissociated" = "#ED0000", "polyA+ Dissociated" = "#42B540",
                     "even" = "#0099B4", "realistic" = "#925E9F", "sparse" = "#FDAF91", "weighted" = "#AD002A") #lancet

colors_genesets <- c("Adipocytes" = "#FCFC3E", "RBCs" = "#E41A1C", "Endothelial cells" = "#377EB8",
                     "Histones" = "#974EA3", "Other polyA(-)" = "#FF8F00", "MT Genes" = "#4DAF4A", "Other" = "#999999")

heatmap_scale_1d <- scale_fill_gradient(low = "#dcf2fe", high = "#056094")
heatmap_scale_2d <- c("#252974", "#191c4d",
                      colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(20),
                      "#660018", "#80001E")
