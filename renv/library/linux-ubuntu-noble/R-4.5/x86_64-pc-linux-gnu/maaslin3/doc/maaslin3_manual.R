## ----eval=FALSE, cache = FALSE------------------------------------------------
# if (!require("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# 
# BiocManager::install("biobakery/maaslin3")

## ----eval=TRUE, cache = FALSE, echo=FALSE-------------------------------------
for (lib in c('maaslin3', 'dplyr', 'ggplot2', 'knitr')) {
    suppressPackageStartupMessages(require(lib, character.only = TRUE))
}

## ----echo = TRUE, results = 'hide', warning = FALSE, cache = FALSE------------
# Read abundance table
taxa_table_name <- system.file("extdata", "HMP2_taxonomy.tsv",
                                package = "maaslin3")
taxa_table <- read.csv(taxa_table_name, sep = '\t', row.names = 1)

# Read metadata table
metadata_name <- system.file("extdata", "HMP2_metadata.tsv",
                                package = "maaslin3")
metadata <- read.csv(metadata_name, sep = '\t', row.names = 1)

metadata$diagnosis <-
    factor(metadata$diagnosis, levels = c('nonIBD', 'UC', 'CD'))
metadata$dysbiosis_state <-
    factor(metadata$dysbiosis_state, levels = c('none', 'dysbiosis_UC',
                                                'dysbiosis_CD'))
metadata$antibiotics <-
    factor(metadata$antibiotics, levels = c('No', 'Yes'))

# Fit models
set.seed(1)
fit_out <- maaslin3(input_data = taxa_table,
                    input_metadata = metadata,
                    output = 'hmp2_output',
                    formula = '~ diagnosis + dysbiosis_state +
                    antibiotics + age + reads',
                    normalization = 'TSS',
                    transform = 'LOG',
                    augment = TRUE,
                    standardize = TRUE,
                    max_significance = 0.1,
                    median_comparison_abundance = TRUE,
                    median_comparison_prevalence = FALSE,
                    max_pngs = 20,
                    cores = 1)

## ----eval=FALSE---------------------------------------------------------------
# maaslin_log_arguments(input_data = taxa_table,
#                     input_metadata = metadata,
#                     output = 'hmp2_output',
#                     formula = '~ diagnosis + dysbiosis_state +
#                     antibiotics + age + reads',
#                     normalization = 'TSS',
#                     transform = 'LOG',
#                     augment = TRUE,
#                     standardize = TRUE,
#                     max_significance = 0.1,
#                     median_comparison_abundance = TRUE,
#                     median_comparison_prevalence = FALSE,
#                     max_pngs = 20,
#                     cores = 1)
# 
# read_data_list <- maaslin_read_data(taxa_table,
#                                     metadata)
# read_data_list <- maaslin_reorder_data(read_data_list$data,
#                                     read_data_list$metadata)
# 
# data <- read_data_list$data
# metadata <- read_data_list$metadata
# 
# formulas <- maaslin_check_formula(data,
#                                 metadata,
#                                 '~ diagnosis + dysbiosis_state +
#                                     antibiotics + age + reads')
# 
# formula <- formulas$formula
# random_effects_formula <- formulas$random_effects_formula
# 
# normalized_data = maaslin_normalize(data,
#                                     output = 'hmp2_output',
#                                     normalization = 'TSS')
# 
# filtered_data = maaslin_filter(normalized_data,
#                             output = 'hmp2_output')
# 
# transformed_data = maaslin_transform(filtered_data,
#                                     output = 'hmp2_output',
#                                     transform = 'LOG')
# 
# standardized_metadata = maaslin_process_metadata(metadata,
#                                                 formula,
#                                                 standardize = TRUE)
# 
# maaslin_results = maaslin_fit(filtered_data,
#                             transformed_data,
#                             standardized_metadata,
#                             formula,
#                             random_effects_formula,
#                             data = data,
#                             median_comparison_abundance = TRUE,
#                             median_comparison_prevalence = FALSE,
#                             cores = 1)
# 
# maaslin_write_results(output = 'hmp2_output',
#                     maaslin_results$fit_data_abundance,
#                     maaslin_results$fit_data_prevalence,
#                     random_effects_formula,
#                     max_significance = 0.1)
# 
# # Invisible to avoid dumping plots
# invisible(maaslin_plot_results(output = 'hmp2_output',
#                             transformed_data,
#                             metadata,
#                             maaslin_results$fit_data_abundance,
#                             maaslin_results$fit_data_prevalence,
#                             normalization = "TSS",
#                             transform = "LOG",
#                             median_comparison_abundance = TRUE,
#                             median_comparison_prevalence = FALSE,
#                             max_significance = 0.1,
#                             max_pngs = 20))

## ----echo = TRUE, results = 'hide', warning = FALSE, cache = FALSE------------
# This section is necessary for updating
# the summary plot and the association plots

# Rename results file with clean titles
all_results <- read.csv('hmp2_output/all_results.tsv', sep='\t')
all_results <- all_results %>%
    mutate(metadata = case_when(metadata == 'age' ~ 'Age',
                                metadata == 'antibiotics' ~ 'Abx',
                                metadata == 'diagnosis' ~ 'Diagnosis',
                                metadata == 'dysbiosis_state' ~ 'Dysbiosis',
                                metadata == 'reads' ~ 'Read depth'),
        value = case_when(value == 'dysbiosis_CD' ~ 'CD',
                            value == 'dysbiosis_UC' ~ 'UC',
                            value == 'Yes' ~ 'Used', # Antibiotics
                            value == 'age' ~ 'Age',
                            value == 'reads' ~ 'Read depth',
                            TRUE ~ value),
        feature = gsub('_', ' ', feature) %>%
            gsub(pattern = 'sp ', replacement = 'sp. '))

# Write results
write.table(all_results, 'hmp2_output/all_results.tsv', sep='\t')

# Set the new heatmap and coefficient plot variables and order them
heatmap_vars = c('Dysbiosis UC', 'Diagnosis UC',
                'Abx Used', 'Age', 'Read depth')
coef_plot_vars = c('Dysbiosis CD', 'Diagnosis CD')

# This section is necessary for updating the association plots
colnames(taxa_table) <- gsub('_', ' ', colnames(taxa_table)) %>%
    gsub(pattern = 'sp ', replacement = 'sp. ')

# Rename the features in the norm transformed data file
data_transformed <-
    read.csv('hmp2_output/features/data_transformed.tsv', sep='\t')
colnames(data_transformed) <-
    gsub('_', ' ', colnames(data_transformed)) %>%
    gsub(pattern = 'sp ', replacement = 'sp. ')
write.table(data_transformed,
            'hmp2_output/features/data_transformed.tsv',
            sep='\t', row.names = FALSE)

# Rename the metadata like in the outputs table
colnames(metadata) <-
    case_when(colnames(metadata) == 'age' ~ 'Age',
            colnames(metadata) == 'antibiotics' ~ 'Abx',
            colnames(metadata) == 'diagnosis' ~ 'Diagnosis',
            colnames(metadata) == 'dysbiosis_state' ~ 'Dysbiosis',
            colnames(metadata) == 'reads' ~ 'Read depth',
            TRUE ~ colnames(metadata))
metadata <- metadata %>%
    mutate(Dysbiosis = case_when(Dysbiosis == 'dysbiosis_UC' ~ 'UC',
                                Dysbiosis == 'dysbiosis_CD' ~ 'CD',
                                Dysbiosis == 'none' ~ 'None') %>%
            factor(levels = c('None', 'UC', 'CD')),
        Abx = case_when(Abx == 'Yes' ~ 'Used',
                        Abx == 'No' ~ 'Not used') %>%
            factor(levels = c('Not used', 'Used')),
        Diagnosis = case_when(Diagnosis == 'nonIBD' ~ 'non-IBD',
                                TRUE ~ Diagnosis) %>%
            factor(levels = c('non-IBD', 'UC', 'CD')))

# Recreate the plots
scatter_plots <- maaslin_plot_results_from_output(
    output = 'hmp2_output',
    metadata,
    normalization = "TSS",
    transform = "LOG",
    median_comparison_abundance = TRUE,
    median_comparison_prevalence = FALSE,
    max_significance = 0.1,
    max_pngs = 20)

## ----out.width='100%', echo=FALSE, cache = FALSE, include=FALSE---------------
knitr::include_graphics("hmp2_output/figures/summary_plot.png")

## ----echo=FALSE, fig.show='hold', cache = FALSE, include=FALSE, eval=FALSE----
# prefix <- "hmp2_output/figures/association_plots"
# plot_vec <- c("/Age/linear/Age_Enterocloster clostridioformis_linear.png",
#             "/Dysbiosis/linear/Dysbiosis_Escherichia coli_linear.png",
#             "/Age/logistic/Age_Bifidobacterium longum_logistic.png",
#             paste0("/Dysbiosis/logistic/",
#                 "Dysbiosis_Faecalibacterium prausnitzii_logistic.png"
#                 ))
# knitr::include_graphics(paste0(prefix, plot_vec))

## ----echo = TRUE, results = 'hide', warning = FALSE, cache = FALSE------------
contrast_mat <- matrix(c(1, 1, 0, 0, 0, 0, 1, 1), 
    ncol = 4, nrow = 2, byrow = TRUE)
    
colnames(contrast_mat) <- c("diagnosisUC",
                            "diagnosisCD",
                            "dysbiosis_statedysbiosis_UC",
                            "dysbiosis_statedysbiosis_CD")
                            
rownames(contrast_mat) <- c("diagnosis_test", "dysbiosis_test")

contrast_mat
                            
maaslin_contrast_test(maaslin3_fit = fit_out,
                        contrast_mat = contrast_mat)

## -----------------------------------------------------------------------------
sessionInfo()

## -----------------------------------------------------------------------------
# Clean-up generated outputs
unlink('hmp2_output', recursive = TRUE)

