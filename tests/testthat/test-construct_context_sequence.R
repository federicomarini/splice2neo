test_that("construct_mutated_range works for A3SS event on positive strand", {
  junc_id <- "chr2_158156177_158157190_+"
  junc_df <- tibble(junc_id = junc_id) %>%
    separate(
      junc_id,
      sep = "_",
      into = c("chr", "pos1", "pos2", "strand"),
      remove = FALSE
    ) %>%
    mutate(pos1 = as.integer(pos1), pos2 = as.integer(pos2))
  junc_pos1 <-
    GenomicRanges::GRanges(junc_df$chr, junc_df$pos1, strand = junc_df$strand)
  junc_pos2 <-
    GenomicRanges::GRanges(junc_df$chr, junc_df$pos2, strand = junc_df$strand)
  transcript_id <- "ENST00000259056"
  # genomic ranges of wt transcript
  wt_transcript_range <- toy_transcripts[[transcript_id]]

  mutated_range <- construct_mutated_range(
    wt_transcript_range = wt_transcript_range,
    strand_direction = "+",
    junc_pos1 = junc_pos1,
    junc_pos2 = junc_pos2
  )
  expect_true(mutated_range@ranges@start[7] == junc_pos2@ranges@start)
  expect_true(mutated_range@ranges@start[6] +  mutated_range@ranges@width[6] - 1 == junc_pos1@ranges@start)
  expect_true(length(mutated_range) == length(wt_transcript_range))

})

test_that("construct_mutated_range works for exon-skipping event on negative strand", {
  junc_id <- "chr2_152727125_152728931_-"
  junc_df <- tibble(junc_id = junc_id) %>%
    separate(
      junc_id,
      sep = "_",
      into = c("chr", "pos1", "pos2", "strand"),
      remove = FALSE
    ) %>%
    mutate(pos1 = as.integer(pos1), pos2 = as.integer(pos2))
  junc_pos1 <-
    GenomicRanges::GRanges(junc_df$chr, junc_df$pos1, strand = junc_df$strand)
  junc_pos2 <-
    GenomicRanges::GRanges(junc_df$chr, junc_df$pos2, strand = junc_df$strand)
  transcript_id <- "ENST00000637224"
  # genomic ranges of wt transcript
  wt_transcript_range <- toy_transcripts[[transcript_id]]

  mutated_range <- construct_mutated_range(
    wt_transcript_range = wt_transcript_range,
    strand_direction = "-",
    junc_pos1 = junc_pos1,
    junc_pos2 = junc_pos2
  )
  expect_true(mutated_range@ranges@start[6] == junc_pos2@ranges@start)
  expect_true(mutated_range@ranges@start[7] +  mutated_range@ranges@width[7] - 1 == junc_pos1@ranges@start)
  expect_true(length(mutated_range) == length(wt_transcript_range) - 1)

})


test_that("juncid2context works", {
  junc_id <- "chr2_158156177_158157190_+"
  requireNamespace("BSgenome.Hsapiens.UCSC.hg19", quietly = TRUE)
  genome <-  BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19
  res <- juncid2context(junc_id =junc_id, tx_id = "ENST00000259056",
                 transcript_db = toy_transcripts,
                 genome_db = genome,
                 window_size = 200)
  expect_true(nchar(res) == 405)
})
