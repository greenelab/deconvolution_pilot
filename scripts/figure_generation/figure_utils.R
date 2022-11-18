library(scales)

alpha_value = 0.5
text_size = 18
line_width = 1.5


theme_set(theme_bw() + 
            theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1),
                  text = element_text(size = 14),
                  strip.background=element_rect(colour="black",
                                                fill="white"),
                  plot.title = element_text(hjust = 0.5)
            )
          )

demultiplexing_colors <- c("red"="#CA0020",
                           "blue"="#0571B0",
                           "green"="#4DAC26",
                           "orange"="#E66101",
                           "gray"="#999999") #ggsci D3?


colors_samples <- c("2251"="#E41A1C","2267"="#377EB8","2283"="#4DAF4A","2293"="#974EA3",
                    "2380"="#FF7F00","2428"="#FCFC3E","2467"="#A65628","2497"="#F781BF") #Set1

colors_celltypes <- hue_pal()(13) #ggplot defaults

colors_methods <- c("bayesprism"="#E6AB02", "bisque"="#66A61E", "cibersortx"="#E7298A",
                    "epic"="#1B9E77", "music"="#7570B3", "nnls"="#D95F02") #Dark2

colors_bulktypes <- c("chunk_ribo"="#00468B","dissociated_ribo"="#ED0000","dissociated_polyA"="#42B540",
                     "even"="#0099B4","realistic"="#925E9F","sparse"="#FDAF91","weighted"="#AD002A") #currently lancet, try locuszoom

colors_genesets <- c("Adipocytes" = "#B79F00", "RBCs" = "#F8766D", "Endothelial cells" = "#00BA38",
                     "Histones" = "#00BFC4", "Other polyA(-)" = "#619CFF", "MT Genes" = "#F564E3", "Other" = "#999999")


heatmap_scale_1d <- scale_fill_gradient(low = "#dcf2fe", high = "#056094")
heatmap_scale_2d <- c("#252974", "#191c4d",
                      colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(20),
                      "#660018", "#80001E")
                      



#mypal = pal_npg()(7); show_col(mypal)
#mypal = pal_aaas()(7); show_col(mypal)
#mypal = pal_nejm()(7); show_col(mypal)
#mypal = pal_lancet()(7); show_col(mypal)
#mypal = pal_jama()(7); show_col(mypal) # for samples
#mypal = pal_jco()(7); show_col(mypal)
#mypal = pal_ucscgb()(7); show_col(mypal)
#mypal = pal_d3()(7); show_col(mypal) # for samples?
#mypal = pal_locuszoom()(7); show_col(mypal)
#mypal = pal_igv()(7); show_col(mypal)
#mypal = pal_uchicago()(7); show_col(mypal)
