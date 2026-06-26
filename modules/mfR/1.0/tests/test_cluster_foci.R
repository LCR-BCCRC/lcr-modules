# Run from modules/mutation_foci/1.0/ :  Rscript tests/test_cluster_foci.R
source("src/R/cluster_foci.R")   # defines cluster_one_chromosome, runs nothing

run <- function(maf, hclust_method = "centroid")
  cluster_one_chromosome(maf, pos_col = "Start_Position",
                         dist_method = "euclidean",
                         hclust_method = hclust_method,
                         h_min = 1L, h_max = 100L)

## Case 1 — two tight clumps + two stragglers
maf1 <- data.frame(
  Chromosome     = "1",
  Start_Position = c(1000, 1004, 1009,     # clump A (spread 9 bp)
                     50000, 50006, 50011,  # clump B (49 kb away)
                     200000,               # straggler
                     999999),              # straggler
  Tumor_Sample_Barcode = paste0("S", 1:8)
)
r1 <- run(maf1)
print(r1$maf[, c("Start_Position", "group")])
g <- r1$maf$group
stopifnot(
  nrow(r1$maf) == 8,                  # left_join kept every row
  !anyNA(g),                          # everyone got a focus
  length(unique(g[1:3])) == 1,        # clump A is one focus
  length(unique(g[4:6])) == 1,        # clump B is one focus
  g[1] != g[4],                       # A and B differ
  g[7] != g[1], g[7] != g[4],         # straggler stands alone
  length(unique(g)) == 4              # 4 foci total
)
cat("Case 1 OK (best_h =", r1$best_h, ")\n")

## Case 2 — same position, different samples -> same focus
maf2 <- data.frame(
  Chromosome = "1",
  Start_Position = c(1000, 1000, 1000, 5000),
  Tumor_Sample_Barcode = paste0("S", 1:4)
)
r2 <- run(maf2)
stopifnot(length(unique(r2$maf$group[1:3])) == 1,
          r2$maf$group[1] != r2$maf$group[4])
cat("Case 2 OK\n")

## Case 3 — degenerate chromosomes
r0 <- run(maf1[0, ]); stopifnot(nrow(r0$maf) == 0, is.na(r0$best_h))
r_one <- run(maf1[1, ]); stopifnot(r_one$maf$group == 1L, is.na(r_one$best_h))
cat("Case 3 (empty + single position) OK\n")

## Case 4 — nothing within h_max -> all singletons, no foci
maf4 <- data.frame(
  Chromosome = "1",
  Start_Position = c(1000, 5000, 9000, 13000),   # all >100 bp apart
  Tumor_Sample_Barcode = paste0("S", 1:4)
)
stopifnot(length(unique(run(maf4)$maf$group)) == 4)
cat("Case 4 OK\n\nAll cases passed.\n")