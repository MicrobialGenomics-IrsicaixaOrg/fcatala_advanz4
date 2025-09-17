## ----eval=FALSE, cache = FALSE------------------------------------------------
# sessionInfo()

## ----eval=FALSE, cache = FALSE------------------------------------------------
# if (!require("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# 
# BiocManager::install("biobakery/maaslin3")

## ----eval=TRUE, cache = FALSE, echo=FALSE-------------------------------------
for (lib in c('maaslin3', 'dplyr', 'ggplot2', 'knitr', 'kableExtra')) {
    suppressPackageStartupMessages(require(lib, character.only = TRUE))
}

## ----cache = FALSE------------------------------------------------------------
# Read abundance table
taxa_table_name <- system.file("extdata", "HMP2_taxonomy.tsv",
                                package = "maaslin3")
taxa_table <- read.csv(taxa_table_name, sep = '\t', row.names = 1)

# Read metadata table
metadata_name <- system.file("extdata", "HMP2_metadata.tsv",
                            package = "maaslin3")
metadata <- read.csv(metadata_name, sep = '\t', row.names = 1)

# Factor the categorical variables to test IBD against healthy controls
metadata$diagnosis <-
    factor(metadata$diagnosis, levels = c('nonIBD', 'UC', 'CD'))
metadata$dysbiosis_state <-
    factor(metadata$dysbiosis_state, levels =
                c('none', 'dysbiosis_UC', 'dysbiosis_CD'))
metadata$antibiotics <-
    factor(metadata$antibiotics, levels = c('No', 'Yes'))

taxa_table[1:5, 1:5]
metadata[1:5, 1:5]

## ----echo = TRUE, results = 'hide', warning = FALSE, cache = FALSE------------
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
                    max_pngs = 10,
                    cores = 1)

## ----echo = TRUE, results = 'hide', warning = FALSE, eval = FALSE-------------
# se <- SummarizedExperiment(
#     assays = list(taxa_table = t(taxa_table)),
#     colData = metadata
# )
# 
# fit_out <- maaslin3(input_data = se,
#                     output = 'hmp2_output',
#                     formula = '~ diagnosis + dysbiosis_state +
#                         antibiotics + age + reads',
#                     normalization = 'TSS',
#                     transform = 'LOG',
#                     augment = TRUE,
#                     standardize = TRUE,
#                     max_significance = 0.1,
#                     median_comparison_abundance = TRUE,
#                     median_comparison_prevalence = FALSE,
#                     max_pngs = 10,
#                     cores = 1)

## ----echo = FALSE, cache = FALSE----------------------------------------------
signif_results <- read.csv('hmp2_output/significant_results.tsv',
sep='\t')
head(signif_results, 20) %>%
    dplyr::mutate_if(is.numeric, .funs = function(x){(
        as.character(signif(x, 3)))}) %>%
    knitr::kable() %>%
    kableExtra::kable_styling("striped", position = 'center') %>%
    kableExtra::scroll_box(width = "800px", height = "400px")

## ----out.width='100%', echo=FALSE, cache = FALSE, include=FALSE, eval=FALSE----
# # Rename summary plot to avoid knitting issues later
# quiet_out <- file.rename('hmp2_output/figures/summary_plot.png',
#                         'hmp2_output/figures/summary_plot_first.png')
# 
# knitr::include_graphics("hmp2_output/figures/summary_plot_first.png")

## ----echo=FALSE, fig.show='hold',include=FALSE, eval=FALSE--------------------
# prefix <- "hmp2_output/figures/association_plots"
# plot_vec <-
#     c("/age/linear/age_Enterocloster_clostridioformis_linear.png",
#     "/dysbiosis_state/linear/dysbiosis_state_Escherichia_coli_linear.png",
#     "/age/logistic/age_Bifidobacterium_longum_logistic.png",
#     paste0("/dysbiosis_state/logistic/",
#         "dysbiosis_state_Faecalibacterium_prausnitzii_logistic.png"
#         ))
# knitr::include_graphics(paste0(prefix, plot_vec))

## ----echo = TRUE, results = 'hide', warning = FALSE, cache = FALSE------------
# This section is necessary for updating the
# summary plot and the association plots

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
taxa_table_copy <- taxa_table
colnames(taxa_table_copy) <- gsub('_', ' ', colnames(taxa_table_copy)) %>%
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
metadata_copy <- metadata
colnames(metadata_copy) <-
    case_when(colnames(metadata_copy) == 'age' ~ 'Age',
            colnames(metadata_copy) == 'antibiotics' ~ 'Abx',
            colnames(metadata_copy) == 'diagnosis' ~ 'Diagnosis',
            colnames(metadata_copy) == 'dysbiosis_state' ~ 'Dysbiosis',
            colnames(metadata_copy) == 'reads' ~ 'Read depth',
            TRUE ~ colnames(metadata_copy))
metadata_copy <- metadata_copy %>%
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
    metadata = metadata_copy,
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

## ----cache = FALSE------------------------------------------------------------
# Abundance table
taxa_table_name <- system.file("extdata", "abundance_spike_in_ex.tsv",
                            package = "maaslin3")
spike_in_taxa_table <- read.csv(taxa_table_name, sep = '\t', row.names = 1)

# Metadata table
metadata_name <- system.file("extdata", "metadata_spike_in_ex.tsv",
                            package = "maaslin3")
spike_in_metadata <- read.csv(metadata_name, sep = '\t', row.names = 1)
for (col in c('Metadata_1', 'Metadata_2', 'Metadata_5')) {
    spike_in_metadata[,col] <- factor(spike_in_metadata[,col])
}

# Spike-in table
unscaled_name <- system.file("extdata",
"scaling_factors_spike_in_ex.tsv",
                            package = "maaslin3")
spike_in_unscaled <- read.csv(unscaled_name, sep = '\t', row.names = 1)

spike_in_taxa_table[c(1:5, 101),1:5]
spike_in_metadata[1:5,]
spike_in_unscaled[1:5, , drop=FALSE]

## ----echo = TRUE, results = 'hide', warning = FALSE, cache = FALSE------------
fit_out <- maaslin3(
    input_data = spike_in_taxa_table,
    input_metadata = spike_in_metadata,
    output = 'spike_in_demo',
    formula = '~ Metadata_1 + Metadata_2 + Metadata_3 +
        Metadata_4 + Metadata_5',
    normalization = 'TSS',
    transform = 'LOG',
    median_comparison_abundance = FALSE,
    unscaled_abundance = spike_in_unscaled)

## ----cache = FALSE------------------------------------------------------------
rownames(fit_out$fit_data_abundance$results) <- NULL
head(fit_out$fit_data_abundance$results, 20) %>%
    dplyr::mutate_if(is.numeric, .funs =
                            function(x){(as.character(signif(x, 3)))}) %>%
    knitr::kable() %>%
    kableExtra::kable_styling("striped", position = 'center') %>%
    kableExtra::scroll_box(width = "800px", height = "400px")

## ----cache = FALSE------------------------------------------------------------
# Abundance table
taxa_table_name <- system.file("extdata", "abundance_total_ex.tsv",
    package = "maaslin3")
total_scaling_taxa_table <- read.csv(taxa_table_name, sep = '\t', row.names = 1)

# Metadata table
metadata_name <- system.file("extdata", "metadata_total_ex.tsv",
                                package = "maaslin3")
total_scaling_metadata <- read.csv(metadata_name, sep = '\t', row.names = 1)
for (col in c('Metadata_1', 'Metadata_3', 'Metadata_5')) {
    spike_in_metadata[,col] <- factor(spike_in_metadata[,col])
}

# Total abundance table
unscaled_name <- system.file("extdata", "scaling_factors_total_ex.tsv",
                                package = "maaslin3")
total_scaling_unscaled <- read.csv(unscaled_name, sep = '\t', row.names = 1)

total_scaling_taxa_table[1:5, 1:5]
total_scaling_metadata[1:5,]
total_scaling_unscaled[1:5, , drop=FALSE]

## ----echo = TRUE, results = 'hide', warning = FALSE, cache = FALSE------------
fit_out <- maaslin3(
    input_data = total_scaling_taxa_table,
    input_metadata = total_scaling_metadata,
    output = 'total_scaling_demo',
    formula = '~ Metadata_1 + Metadata_2 + Metadata_3 +
        Metadata_4 + Metadata_5',
    normalization = 'TSS',
    transform = 'LOG',
    median_comparison_abundance = FALSE,
    unscaled_abundance = total_scaling_unscaled)

## ----cache = FALSE------------------------------------------------------------
rownames(fit_out$fit_data_abundance$results) <- NULL
head(fit_out$fit_data_abundance$results, n = 20) %>%
    dplyr::mutate_if(is.numeric,
                    .funs = function(x){(as.character(signif(x, 3)))}) %>%
    knitr::kable() %>%
    kableExtra::kable_styling("striped", position = 'center') %>%
    kableExtra::scroll_box(width = "800px", height = "400px")

## ----echo = TRUE, results = 'hide', warning = FALSE, messages = FALSE---------
# Subset to only CD cases for time; taxa are subset automatically
fit_out <- maaslin3(
    input_data = taxa_table,
    input_metadata = metadata[metadata$diagnosis == 'CD',],
    output = 'random_effects_output',
    formula = '~ dysbiosis_state + antibiotics +
        age + reads + (1|participant_id)',
    plot_summary_plot = FALSE,
    plot_associations = FALSE)

## ----echo = FALSE, cache = FALSE----------------------------------------------
signif_results <-
read.csv('random_effects_output/significant_results.tsv',
                            sep='\t')
head(signif_results, n = 20) %>%
    dplyr::mutate_if(is.numeric, .funs =
                        function(x){(as.character(signif(x, 3)))}) %>%
    knitr::kable() %>%
    kableExtra::kable_styling("striped", position = 'center') %>%
    kableExtra::scroll_box(width = "800px", height = "400px")

## ----echo = TRUE, results = 'hide', warning = FALSE, cache = FALSE------------
metadata <- metadata %>%
    mutate(dysbiosis_general = ifelse(dysbiosis_state != 'none',
                                    'dysbiosis', 'none')) %>%
    mutate(dysbiosis_general = factor(dysbiosis_general, levels =
                                        c('none', 'dysbiosis')))

fit_out <- maaslin3(
    input_data = taxa_table,
    input_metadata = metadata,
    output = 'interaction_output',
    formula = '~ diagnosis + diagnosis:dysbiosis_general +
        antibiotics + age + reads')

## ----cache = FALSE------------------------------------------------------------
full_results <- rbind(fit_out$fit_data_abundance$results,
                    fit_out$fit_data_prevalence$results)
full_results <- full_results %>%
    dplyr::arrange(qval_joint) %>%
    dplyr::filter(metadata == "diagnosis")
rownames(full_results) <- NULL
head(full_results, n = 20) %>%
    dplyr::mutate_if(is.numeric,
                    .funs = function(x){(as.character(signif(x, 3)))}) %>%
    knitr::kable() %>%
    kableExtra::kable_styling("striped", position = 'center') %>%
    kableExtra::scroll_box(width = "800px", height = "400px")

## ----echo = TRUE, results = 'hide', warning = FALSE, cache = FALSE------------
# Put the red meat consumption responses in order
metadata <- metadata %>%
    mutate(red_meat = ifelse(
        red_meat == 'No, I did not consume these products in the last 7 days',
                            'Not in the last 7 days',
                            red_meat) %>%
            factor(levels = c('Not in the last 7 days',
                                'Within the past 4 to 7 days',
                                'Within the past 2 to 3 days',
                                'Yesterday, 1 to 2 times',
                                'Yesterday, 3 or more times'))
    )

# Create the model with only non-IBD subjects
fit_out <- maaslin3(
    input_data = taxa_table,
    input_metadata = metadata[metadata$diagnosis == 'nonIBD',],
    output = 'ordered_outputs',
    formula = '~ ordered(red_meat) + antibiotics + age + reads',
    plot_summary_plot = TRUE,
    plot_associations = TRUE,
    heatmap_vars = c('red_meat Within the past 4 to 7 days',
                    'red_meat Within the past 2 to 3 days',
                    'red_meat Yesterday, 1 to 2 times',
                    'red_meat Yesterday, 3 or more times'),
    max_pngs = 30)

## ----out.width='100%', echo=FALSE, cache = FALSE, include=FALSE---------------
knitr::include_graphics("ordered_outputs/figures/summary_plot.png")

## ----echo = FALSE, cache = FALSE----------------------------------------------
full_results <- rbind(fit_out$fit_data_abundance$results,
                    fit_out$fit_data_prevalence$results)
full_results <- full_results %>%
    dplyr::filter(metadata %in% c("red_meat") &
                    feature == 'Alistipes_shahii' &
                    model == 'logistic') %>%
    dplyr::mutate(
        value =
            factor(
                value,
                levels =
                    c('No, I did not consume these products in the last 7 days',
                                                'Within the past 4 to 7 days',
                                                'Within the past 2 to 3 days',
                                                'Yesterday, 1 to 2 times',
                                            'Yesterday, 3 or more times'))) %>%
    dplyr::arrange(value)
rownames(full_results) <- NULL
head(full_results, n = 20) %>%
    dplyr::mutate_if(is.numeric,
                    .funs = function(x){(as.character(signif(x, 3)))}) %>%
    knitr::kable() %>%
    kableExtra::kable_styling("striped", position = 'center') %>%
    kableExtra::scroll_box(width = "800px", height = "200px")

## ----echo = T, results = 'hide', warning = FALSE, cache = FALSE---------------
fit_out <- maaslin3(
    input_data = taxa_table,
    input_metadata = metadata[metadata$diagnosis == 'nonIBD',],
    output = 'group_outputs',
    formula = '~ group(red_meat) + antibiotics + age + reads',
    plot_summary_plot = TRUE,
    plot_associations = TRUE,
    heatmap_vars = c('red_meat Within the past 4 to 7 days',
                    'red_meat Within the past 2 to 3 days',
                    'red_meat Yesterday, 1 to 2 times',
                    'red_meat Yesterday, 3 or more times'),
    max_pngs = 200)

## ----echo = FALSE, cache = FALSE----------------------------------------------
full_results <- rbind(fit_out$fit_data_abundance$results,
                    fit_out$fit_data_prevalence$results)
full_results <- full_results %>%
    dplyr::filter(metadata %in% c("red_meat") &
                    feature == 'Alistipes_shahii' &
                    model == 'logistic') %>%
    dplyr::mutate(
        value =
            factor(
                value,
                levels =
                    c('No, I did not consume these products in the last 7 days',
                                                'Within the past 4 to 7 days',
                                                'Within the past 2 to 3 days',
                                                'Yesterday, 1 to 2 times',
                                            'Yesterday, 3 or more times'))) %>%
    dplyr::arrange(value)
rownames(full_results) <- NULL
head(full_results, n = 20) %>%
    dplyr::mutate_if(is.numeric,
                    .funs = function(x){(as.character(signif(x, 3)))}) %>%
    knitr::kable() %>%
    kableExtra::kable_styling("striped", position = 'center') %>%
    kableExtra::scroll_box(width = "800px", height = "100px")

## ----echo = TRUE, results = 'hide', warning = FALSE, cache = FALSE------------
fit_out <- maaslin3(input_data = taxa_table,
                    input_metadata = metadata,
                    output = 'contrast_test_output',
                    formula = '~ diagnosis + dysbiosis_state +
                        antibiotics + age + reads',
                    plot_summary_plot = FALSE,
                    plot_associations = FALSE,
                    cores = 1)

## ----echo = TRUE, warning = FALSE, cache = FALSE------------------------------
contrast_mat <- matrix(c(1, -1, 0, 0, 0, 0, 1, -1), 
    ncol = 4, nrow = 2, byrow = TRUE)
    
colnames(contrast_mat) <- c("diagnosisUC",
                            "diagnosisCD",
                            "dysbiosis_statedysbiosis_UC",
                            "dysbiosis_statedysbiosis_CD")
                            
rownames(contrast_mat) <- c("diagnosis_test", "dysbiosis_test")

contrast_mat

contrast_out <- maaslin_contrast_test(maaslin3_fit = fit_out,
                        contrast_mat = contrast_mat)

head(contrast_out$fit_data_abundance$results, n = 20) %>%
    knitr::kable() %>%
    kableExtra::kable_styling("striped", position = 'center') %>%
    kableExtra::scroll_box(width = "800px", height = "400px")

## -----------------------------------------------------------------------------
sessionInfo()

## -----------------------------------------------------------------------------
# Clean-up
unlink('hmp2_output', recursive = TRUE)
unlink('spike_in_demo', recursive = TRUE)
unlink('total_scaling_demo', recursive = TRUE)
unlink('random_effects_output', recursive = TRUE)
unlink('interaction_output', recursive = TRUE)
unlink('ordered_outputs', recursive = TRUE)
unlink('group_outputs', recursive = TRUE)
unlink('contrast_test_output', recursive = TRUE)

