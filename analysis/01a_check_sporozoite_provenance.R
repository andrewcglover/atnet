# Provenance check for the sporozoite summary used by the EIP fit.
#
# The summary SPZ_pooled_summary_APR26.csv was assembled by hand in Excel. This
# script reconstructs the assay groups from the workbooks with the package
# functions and reports, for each row of the summary, which group or groups it
# reproduces. It writes nothing and changes nothing.
#
# Run with the working directory at the root of the package.

devtools::load_all(quiet = TRUE)

raw_dir <- "data_private/archive/raw/2026-05_exp_data"
summary_path <- file.path(raw_dir, "SPZ_pooled_summary_APR26.csv")
workbooks <- c(
  individual = "200mgPCLnetdata_individualreps.xlsx",
  pooled_reference = "200mgPCLnetdata_pooledreps_reference.xlsx",
  pooled = "200mgPCLnetdata_pooledreps.xlsx"
)

# Pair each control sub-group with the treated sub-group in the same position,
# for the same exposure time and number of blood meals.
paired_groups <- function(path, book) {
  tab <- tabulate_assay_workbook(path)
  if (!nrow(tab)) {
    return(NULL)
  }
  key <- interaction(tab$sheet, tab$exposure_h, tab$n_bloodmeals, tab$subgroup,
                     drop = TRUE)
  out <- lapply(split(tab, key), function(g) {
    ctl <- g[g$arm == "control", ]
    trt <- g[g$arm == "treated", ]
    if (nrow(ctl) != 1L || nrow(trt) != 1L) {
      return(NULL)
    }
    data.frame(
      book = book, sheet = ctl$sheet, subgroup = ctl$subgroup,
      exposure_h = ctl$exposure_h, n_bloodmeals = ctl$n_bloodmeals,
      n_CTL = ctl$n, pos_CTL = ctl$n_positive,
      n_200mg = trt$n, pos_200mg = trt$n_positive,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, out)
}

candidates <- do.call(rbind, lapply(names(workbooks), function(book) {
  paired_groups(file.path(raw_dir, workbooks[[book]]), book)
}))

target <- read.csv(summary_path)
counts <- c("n_CTL", "pos_CTL", "n_200mg", "pos_200mg")

message(sprintf("%d candidate groups; %d summary rows\n", nrow(candidates),
                nrow(target)))

unmatched <- integer(0)
for (i in seq_len(nrow(target))) {
  row <- target[i, ]
  hit <- which(apply(candidates[counts], 1L, function(z) {
    all(z == unlist(row[counts]))
  }))
  label <- sprintf("row %2d | sheet %s, %gh, %dBF, day %d | ", i,
                   row$pooled_sheet, row$post_time, row$BF, row$day)
  if (!length(hit)) {
    unmatched <- c(unmatched, i)
    message(label, "NO SOURCE FOUND")
  } else {
    message(label, paste(sprintf("%s[%s #%d]", candidates$book[hit],
                                 candidates$sheet[hit], candidates$subgroup[hit]),
                         collapse = "; "))
  }
}

message(sprintf("\n%d of %d rows reproduced.", nrow(target) - length(unmatched),
                nrow(target)))
if (length(unmatched)) {
  message("Unreproduced rows: ", paste(unmatched, collapse = ", "),
          ".\nSee data_private/DATA_SOURCES.md for the corrections that ",
          "explain them.")
}
