# One of the paper's reviewers raised a good question: if deconvolution results
# can change based on the size of the reference profile, how big does the
# reference need to be to return good results? So we created several simulated
# reference profiles by proportionally downsampling our single-cell data, down
# to as small as 200 cells. This work eventually became Supplemental Figure 6.

suppressPackageStartupMessages({
  library(data.table)
  library(SingleCellExperiment)
  library(ggplot2)
  library(dplyr)
  library(yaml)
})

params <- read_yaml("../../config.yml")
data_path <- params$data_path
local_data_path <- params$local_data_path
plot_path <- params$plot_path

# Filter down to the methods that use single-cell data, so we can compare the
# effect of difference reference profiles.
single_cell_methods <- c("bayesprism","cibersortx","bisque","music","nnls")

source("evaluation_functions.R")

# Load all deconvolution results into a single dataframe
melted_results <- load_melted_results(reference_comp = TRUE)
melted_results <- subset(melted_results, melted_results$method %in% single_cell_methods)

# Calculate the variance of each cell type x method combo, adding in results
# from a smaller reference profile each time
original_results <- subset(melted_results, melted_results$reference %in% c("genetic", "hashing"))
variance2 <- original_results %>% group_by(method, cell_type, variable, bulk_type) %>%
  summarize(variance = var(proportion))
ggplot(variance2, mapping = aes(x = method, y = log(variance), fill = cell_type)) +
  geom_boxplot() + ggtitle("Just hashing and genetic") + ylim(-50, -3)
variance2 %>% group_by(method) %>% summarize(means = mean(variance)) %>% arrange(means)

reasonable_results <- subset(melted_results, melted_results$reference %in% c("genetic", "hashing","sim2000"))
variance3 <- reasonable_results %>% group_by(method, cell_type, variable, bulk_type) %>%
  summarize(variance = var(proportion))
ggplot(variance3, mapping = aes(x = method, y = log(variance), fill = cell_type)) +
  geom_boxplot() + ggtitle("Just >2000 refs") + ylim(-50, -3)
variance3 %>% group_by(method) %>% summarize(means = mean(variance)) %>% arrange(means)

reasonable_results <- subset(melted_results, melted_results$reference %in% c("genetic", "hashing","sim2000","sim1000"))
variance4 <- reasonable_results %>% group_by(method, cell_type, variable, bulk_type) %>%
  summarize(variance = var(proportion))
ggplot(variance4, mapping = aes(x = method, y = log(variance), fill = cell_type)) +
  geom_boxplot() + ggtitle("Just >1000 refs") + ylim(-50, -3)
variance4 %>% group_by(method) %>% summarize(means = mean(variance)) %>% arrange(means)

reasonable_results <- subset(melted_results, melted_results$reference %in% c("genetic", "hashing","sim2000","sim1000","sim500"))
variance5 <- reasonable_results %>% group_by(method, cell_type, variable, bulk_type) %>%
  summarize(variance = var(proportion))
ggplot(variance5, mapping = aes(x = method, y = log(variance), fill = cell_type)) +
  geom_boxplot() + ggtitle("Just >500 refs") + ylim(-50, -3)
variance5 %>% group_by(method) %>% summarize(means = mean(variance))%>% arrange(means)

variance <- melted_results %>% group_by(method, cell_type, variable, bulk_type) %>%
  summarize(variance = var(proportion))
ggplot(variance, mapping = aes(x = method, y = log(variance), fill = cell_type)) +
  geom_boxplot() + ggtitle("All simulations") + ylim(-50, -3)
variance %>% group_by(method) %>% summarize(means = mean(variance)) %>% arrange(means)

# If we treat the results from the biggest reference profile (single cells that
# were assigned by genetic demultiplexing) as default, then subtract the smaller
# profile results from it, do we see any cell type trends? Doesn't look like it.
genetic_results <- subset(melted_results, reference=="genetic")
genetic_results$reference <- NULL; setnames(genetic_results, "proportion", "genetic_proportion")
other_results <- subset(melted_results, reference!="genetic")
deconvolution <- full_join(genetic_results, other_results)

deconvolution$prop_diff <- deconvolution$proportion - deconvolution$genetic_proportion 

ggplot(subset(deconvolution, reference=="sim200")) + geom_boxplot(mapping = aes(x = cell_type, y= prop_diff, fill = method, color = method)) +
  geom_boxplot(mapping = aes(x = cell_type, y = proportion, fill = method), outlier.colour = NA) +
  scale_color_manual(name = "Method", values = colors_methods) +
  scale_fill_manual(name = "Method", values = colors_methods) +
  ylab("Estimated proportion - single cell proportion") +
  theme(axis.title.x = element_blank()) + ggtitle("sim200")


# A quick wrapper function to plot side-by-side pie charts of the different cell
# type estimates for any two reference profiles for a single sample, bulk type,
# and deconvolution method.
plot_any_two <- function(sample, bulk_type_name, reference1, reference2, methodname) {
  subset1 <- subset(melted_results, variable==sample & bulk_type==bulk_type_name &
                      reference==reference1 & method==methodname)
  subset2 <- subset(melted_results, variable==sample & bulk_type==bulk_type_name &
                      reference==reference2 & method==methodname)
  
  gA <- ggplot(subset1, aes(x="", y=proportion, fill=cell_type)) + 
    geom_bar(stat="identity", width = 1) + coord_polar("y", start = 0) +
    scale_fill_manual(name = "Cell type", values = colors_celltypes) +
    theme_void() + ggtitle(reference1) + theme(plot.title = element_text(hjust = 0.5))
  
  gB <- ggplot(subset2, aes(x="", y=proportion, fill=cell_type)) + 
    geom_bar(stat="identity", width = 1) + coord_polar("y", start = 0) +
    scale_fill_manual(name = "Cell type", values = colors_celltypes) +
    theme_void() + ggtitle(reference2) + theme(plot.title = element_text(hjust = 0.5))
  
  
  gA + gB + plot_layout(guides = "collect") + plot_annotation(
    paste(sample, methodname, bulk_type_name, sep = " ")
  )
}

plot_any_two("2267", "dissociated_ribo", "genetic", "sim200", "bayesprism")
plot_any_two("2497", "dissociated_ribo", "genetic", "sim200", "bayesprism")
plot_any_two("2251", "dissociated_polyA", "genetic", "sim200", "bayesprism")

plot_any_two("2251_sim_1", "sparse", "genetic", "sim200", "bayesprism")
plot_any_two("2380_sim_30", "sparse", "genetic", "sim200", "bayesprism")

plot_any_two("2267", "dissociated_ribo", "genetic", "sim200", "cibersortx")
plot_any_two("2497", "dissociated_ribo", "genetic", "sim200", "cibersortx")
plot_any_two("2251", "dissociated_polyA", "genetic", "sim200", "cibersortx")

plot_any_two("2251_sim_1", "sparse", "genetic", "sim200", "cibersortx")
plot_any_two("2380_sim_30", "sparse", "genetic", "sim200", "cibersortx")

plot_any_two("2267", "dissociated_ribo", "genetic", "sim200", "music")
plot_any_two("2497", "dissociated_ribo", "genetic", "sim200", "music")
plot_any_two("2251", "dissociated_polyA", "genetic", "sim200", "music")

plot_any_two("2251_sim_1", "sparse", "genetic", "sim200", "music")
plot_any_two("2380_sim_30", "sparse", "genetic", "sim200", "music")