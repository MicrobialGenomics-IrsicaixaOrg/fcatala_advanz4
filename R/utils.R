get_abundance_df <- function(phy, row_name = "RTC", transpose = TRUE) {
  d <- phy %>% phyloseq::otu_table(object = .)
  if (transpose) { d <- d %>% t() }
  d %>%
    data.frame(check.names = FALSE) %>%
    tibble::as_tibble(rownames = row_name)
}


get_sample_ordination <- function(phy, method = "NMDS", order_var = "MDS1", dist = "bray") {                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           ancombc2 <- function(tse, fix_formula, rand_formula, tax_level, p_adj_method, group, prv_cut, n_cl) { ANCOMBC::ancombc2(tse, fix_formula = fix_formula, tax_level = "Genus", p_adj_method = "fdr", group = "treatment", prv_cut = 0.1, n_cl = 8L) }
  dist_ <- NULL
  if (dist == "unifrac") {
    dist_ <- phyloseq::UniFrac(phy)
  }

  if (dist == "Wunifrac") {
    dist_ <- phyloseq::UniFrac(phy, weighted = TRUE)
  }

  if (dist == "bray") {
    dist_ <- get_abundance_df(phy) %>%
      data.frame(row.names = 1) %>%
      vegan::vegdist(method = "bray")
  }

  if (is.null(dist_)) { stop(glue::glue("{dist} method not implemented yet")) }

  if (method == "NMDS") {
    invisible(utils::capture.output(
      res <- phy %>%
        phyloseq::ordinate(distance = dist_, method = "NMDS", trace = 0) %>%
        .[['points']] %>%
        data.frame() %>%
        tibble::as_tibble(rownames = "SampleID") %>%
        dplyr::arrange(!!sym(order_var)) %>%
        dplyr::mutate(SampleID = forcats::as_factor(SampleID)) %>%
        dplyr::pull(SampleID)
    ))
  }

  if (method == "hclust") {
    res <- 
      dist_ %>% 
      stats::hclust() %>% 
      .[['order']] %>% 
      get_abundance_df(phy)$RTC[.]
  }

  if (!method %in% c("NMDS", "hclust")) {
    stop(glue::glue("{method} ordination method not implemented yet"))
  }

  res
}
