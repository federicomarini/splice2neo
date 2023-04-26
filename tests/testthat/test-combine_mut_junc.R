test_that("combine_mut_junc works", {

  # get spliceAI junctions
  spliceai_file <- system.file("extdata", "spliceai_output.vcf", package = "splice2neo")

  spliceai_annot_df <- parse_spliceai(spliceai_file) %>%
    format_spliceai() %>%
    annotate_mut_effect(toy_transcripts, toy_transcripts_gr)

  # get pangolin junctions
  pangolin_file <- system.file("extdata", "spliceai_output.pangolin.vcf", package = "splice2neo")

  pangolin_annot_df <- parse_pangolin(pangolin_file) %>%
    format_pangolin() %>%
    annotate_mut_effect(toy_transcripts, toy_transcripts_gr)

  # get mmsplice junctions
  mmsplice_file <- system.file("extdata", "mmsplice_pred.csv", package = "splice2neo")
  mmsplice_annot_df <- parse_mmsplice(mmsplice_file) %>%
    annotate_mmsplice(toy_transcripts)

  mmsplice_annot_df$pathogenicity <- "xx"
  mmsplice_annot_df$effect <- "xx"
  mmsplice_annot_df$efficiency <- "xx"
  # junc_annot <- annotate_mmsplice(mmsplice_df, toy_transcripts)

  junc_data_list = list(
    "spliceai" = spliceai_annot_df,
    "pangolin" = pangolin_annot_df,
    "mmsplice" = mmsplice_annot_df
  )

  df_comb <- combine_mut_junc(junc_data_list)

  expect_true(nrow(df_comb) >= 1)

})

test_that("combine_mut_junc gives unique output", {

  # an (ugly) example
  pangolin_file <- system.file("extdata", "spliceai_output.pangolin.vcf", package = "splice2neo")
  df_pango <- parse_pangolin(pangolin_file) %>% rowwise() %>% filter(!grepl(",", ALT)) %>% ungroup() %>% dplyr::slice(1:4)

  # modify for a match with spliceai
  df_pango$CHROM[1] <- "chr4"
  df_pango$POS[1] <- "53751995"
  df_pango$REF[1] <- "G"
  df_pango$ALT[1] <- "A"
  df_pango$increase_pos[1] <- 39
  df_pango$increase_score[1] <- 0.04
  df_pango$decrease_score[1] <- 0.05
  df_pango$decrease_pos[1] <- 38
  df_pango <- format_pangolin(df_pango)
  df_pango_annot <- annotate_mut_effect(df_pango, toy_transcripts, toy_transcripts_gr)

  spliceai_file <- system.file("extdata", "spliceai_output.vcf", package = "splice2neo")
  df_spliceai <- parse_spliceai(spliceai_file) %>% rowwise() %>% filter(!grepl(",", ALT)) %>% ungroup() %>% dplyr::slice(1:4)
  # modify for match with pangolin
  df_spliceai$CHROM[1] <- "chr4"
  df_spliceai$POS[1] <- "53751995"
  df_spliceai$REF[1] <- "G"
  df_spliceai$ALT[1] <- "A"
  df_spliceai$DS_AG[1] <- 0.04
  df_spliceai$DS_AL[1] <- 0.01
  df_spliceai$DS_DG[1] <- 0.05
  df_spliceai$DS_DL[1] <- 0.05
  df_spliceai$DP_AG[1] <- -14
  df_spliceai$DP_AL[1] <- 38
  df_spliceai$DP_DG[1] <- 39
  df_spliceai$DP_DL[1] <- -9
  df_spliceai <- format_spliceai(df_spliceai)
  df_spliceai_annot <- annotate_mut_effect(df_spliceai, toy_transcripts, toy_transcripts_gr)

  in_file <- system.file("extdata", "mmsplice_pred.csv", package = "splice2neo")
  df_mmsplcie <- parse_mmsplice(in_file) %>% filter(transcript_id %in% names(toy_transcripts))
  df_mmsplcie_annot <- annotate_mmsplice(df_mmsplcie[1:10,], toy_transcripts)

  df_combined <- combine_mut_junc(list(
    pangolin = df_pango_annot,
    spliceai = df_spliceai_annot,
    mmsplice =df_mmsplcie_annot))

  df_combined %>%
    distinct(mut_id, tx_id, junc_id)

  df_combined %>%
    filter(junc_id == "chr4:53752033-53752034:-", tx_id == "ENST00000388940.8_2", mut_id == "chr4_53751995_G_A" ) %>%
    dplyr::select(mut_id, tx_id, junc_id, pangolin_effect, spliceai_effect)


  expect_equal(nrow(df_combined), nrow(df_combined %>% distinct(mut_id, tx_id, junc_id)))

})

test_that("combine_mut_junc gives unique output with minimal input", {

  df_pango_annot <- tribble(
    ~mut_id,               ~tx_id,             ~junc_id,                    ~effect,
    "chr2_179415988_C_CA", "ENST00000419746",  "chr2:179415962-179447299:+", "DG",
    "chr2_179415988_C_CA", "ENST00000419746",  "chr2:179415962-179447299:+", "AG",
  )
  df_spliceai_annot <- tribble(
    ~mut_id,               ~tx_id,             ~junc_id,                    ~effect,
    "chr2_179415988_C_CA", "ENST00000419746",  "chr2:179415962-179447299:+", "DG"
  )


  df_combined <- combine_mut_junc(list(
    pangolin = df_pango_annot,
    spliceai = df_spliceai_annot))

})
