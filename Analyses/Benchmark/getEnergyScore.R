#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(stringr)
})

msg_stop <- function(...) { message(...); quit(save="no", status=1) }

# =========================
#  E-statistics helpers
# =========================
rho_XY <- function(x, y) {
  K <- length(x); if (length(y) != K) stop("x,y length mismatch")
  # mean pairwise L2 distances between X and Y (标量等价于绝对值)
  as.numeric( sum(abs(outer(x, y, "-"))) / (K^2) )
}

rho_within <- function(v) {
  K <- length(v); if (K < 2) return(0)
  idx <- utils::combn(K, 2)
  as.numeric( 2 * sum(abs(v[idx[1,]] - v[idx[2,]])) / (K*(K-1)) )
}

energy_stat <- function(x, y) {
  rXY <- rho_XY(x, y)
  rX  <- rho_within(x)
  rY  <- rho_within(y)
  list(D = 2*rXY - rX - rY, rho_XY = rXY, rho_X = rX, rho_Y = rY)
}

# =========================
#  CLI options
# =========================
option_list <- list(
  make_option(c("-i","--input"),  type="character", help="scRBP TRS 结果表 (CSV/TSV，支持 .gz)"),
  make_option(c("-t","--target"), type="character", help="目标细胞类型（与数据中 cell_type 列匹配）"),
  make_option(c("-s","--score_col"), type="character", default=NA,
              help="用于计算的分数列（默认自动在 TRS_mm/TRS/z_TRS 中选择）"),
  make_option(c("--celltype_col"), type="character", default="cell_type",
              help="细胞类型列名 [默认: %default]"),
  make_option(c("--regulon_col"),  type="character", default="RBP",
              help="regulon/RBP 列名 [默认: %default]"),
  make_option(c("-b","--background"), type="character", default="others",
              help="背景细胞类型集合：'others'=除目标外全部；或逗号分隔列表，如 'T_cells,NK_cells'"),
  make_option(c("-o","--out"), type="character", default=NULL,
              help="输出前缀（写出 <prefix>.summary.csv 与 <prefix>.per_regulon.csv）"),
  make_option(c("--zscore_within_regulon"), action="store_true", default=FALSE,
              help="对每个 regulon 跨细胞类型做 z-score 再计算 E-statistics")
)
opt <- parse_args(OptionParser(option_list=option_list))

if (is.null(opt$input) || is.null(opt$target)) {
  msg_stop("Usage:\n  Rscript getEnergyScore.R --input file.csv --target 'Monocytes' [--score_col TRS_mm]")
}
if (!file.exists(opt$input)) msg_stop("找不到输入文件: ", opt$input)

# =========================
#  Robust reader (CSV 强制, 失败则回退)
# =========================
fread_smart <- function(fp) {
  # 先取首行用于回退时拆列名
  first <- readLines(fp, n=1, warn=FALSE, encoding="UTF-8")
  first <- sub("^\ufeff", "", first, perl=TRUE)  # 去 BOM

  # 尝试：强制逗号，清理空白
  dt <- tryCatch(
    data.table::fread(
      fp, sep=",", header=TRUE, data.table=TRUE,
      strip.white=TRUE, quote='"', fill=Inf, showProgress=FALSE
    ),
    error = function(e) NULL
  )
  if (!is.null(dt) && ncol(dt) > 1) return(dt)

  # 回退1：用首行手动拆分列名再读
  cn <- tryCatch(strsplit(first, ",", fixed=TRUE)[[1]], error=function(e) NULL)
  if (!is.null(cn)) {
    cn <- trimws(sub("^\ufeff","", cn, perl=TRUE))
    dt2 <- tryCatch(
      data.table::fread(
        fp, sep=",", header=TRUE, data.table=TRUE,
        strip.white=TRUE, quote="", fill=TRUE, col.names=cn, showProgress=FALSE
      ),
      error=function(e) NULL
    )
    if (!is.null(dt2) && ncol(dt2) > 1) return(dt2)
  }

  # 回退2：若是 .tsv/.gz.tsv 之类
  if (grepl("\\.tsv(\\.gz)?$", fp, ignore.case=TRUE)) {
    dt3 <- tryCatch(
      data.table::fread(fp, sep="\t", header=TRUE, data.table=TRUE,
                        strip.white=TRUE, quote="", fill=TRUE, showProgress=FALSE),
      error=function(e) NULL
    )
    if (!is.null(dt3)) return(dt3)
  }

  stop("无法可靠读取该文件，请检查分隔符/首行。")
}

# 读取数据（这一步此前你漏掉了）
dt <- tryCatch(fread_smart(opt$input), error=function(e) msg_stop("读取文件失败: ", e$message))

# =========================
#  列名清洗/匹配
# =========================
normalize_names <- function(nms) {
  nms <- iconv(nms, from="", to="UTF-8")
  nms <- sub("^\ufeff", "", nms, perl=TRUE)  # 去 BOM
  nms <- trimws(nms)
  nms <- gsub("\\s+", "_", nms)
  nms
}
names(dt) <- normalize_names(names(dt))

match_col <- function(want, pool) {
  if (is.na(want)) return(NA_character_)
  want_norm <- tolower(normalize_names(want))
  pool_lc   <- tolower(pool)
  hit <- which(pool_lc == want_norm)
  if (length(hit) == 1) pool[hit] else NA_character_
}

ct_col <- match_col(opt$celltype_col, names(dt))
rg_col <- match_col(opt$regulon_col,  names(dt))
if (is.na(ct_col)) {
  message("可用列名：", paste(names(dt), collapse=", "))
  msg_stop("未找到细胞类型列: ", opt$celltype_col)
}
if (is.na(rg_col)) {
  message("可用列名：", paste(names(dt), collapse=", "))
  msg_stop("未找到 regulon 列: ", opt$regulon_col)
}

score_col <- opt$score_col
if (is.na(score_col)) {
  cand <- c("TRS_mm","TRS","z_TRS","trs","TRS_mean","TRS_z")
  score_col <- cand[cand %in% names(dt)][1]
  if (is.na(score_col)) {
    message("可用列名：", paste(names(dt), collapse=", "))
    msg_stop("未找到用于计算的分数列，请用 --score_col 指定。")
  }
} else {
  score_col <- match_col(score_col, names(dt))
  if (is.na(score_col)) {
    message("可用列名：", paste(names(dt), collapse=", "))
    msg_stop("分数列不存在: ", opt$score_col)
  }
}

# =========================
#  背景细胞类型集合
# =========================
all_ct <- sort(unique(dt[[ct_col]]))
if (!(opt$target %in% all_ct)) {
  msg_stop("目标细胞类型不存在: ", opt$target, "；数据中的 cell_type 有：", paste(all_ct, collapse=", "))
}
bg_ct <- if (identical(opt$background, "others")) {
  setdiff(all_ct, opt$target)
} else {
  z <- trimws(unlist(strsplit(opt$background, ",")))
  miss <- setdiff(z, all_ct)
  if (length(miss)) msg_stop("背景细胞类型不存在: ", paste(miss, collapse=", "))
  z
}
if (!length(bg_ct)) msg_stop("背景集合为空。")

# =========================
#  可选：每个 regulon 跨细胞类型 z-score
# =========================
if (isTRUE(opt$zscore_within_regulon)) {
  dt[, (score_col) := {
    v <- get(score_col)
    (v - mean(v, na.rm=TRUE)) / sd(v, na.rm=TRUE)
  }, by = rg_col]
}

# =========================
#  组装 X/Y 并计算 E-statistics
# =========================
# 目标细胞类型：每个 regulon 取均值
x_dt <- dt[get(ct_col) == opt$target,
           .(x = mean(get(score_col), na.rm=TRUE)), by = c(rg_col)]

# 背景：指定背景细胞类型集合的均值
y_dt <- dt[get(ct_col) %in% bg_ct,
           .(y = mean(get(score_col), na.rm=TRUE)), by = c(rg_col)]

xy <- merge(x_dt, y_dt, by = rg_col)
xy <- xy[is.finite(x) & is.finite(y)]
K  <- nrow(xy)
if (K < 2) msg_stop("可用于计算的 regulon 数不足 (K<2)。")

es <- energy_stat(xy$x, xy$y)

summary_dt <- data.table(
  target_cell_type = opt$target,
  background       = if (identical(opt$background, "others")) "others" else paste(bg_ct, collapse="|"),
  score_col        = score_col,
  K_regulons       = K,
  rho_XY           = es$rho_XY,
  rho_X            = es$rho_X,
  rho_Y            = es$rho_Y,
  E_statistic_D    = es$D
)

per_reg_dt <- xy[, .(regulon = get(rg_col), x, y,
                     diff = x - y, abs_diff = abs(x - y))]

cat("==== E-statistics ====\n")
print(summary_dt)

if (!is.null(opt$out)) {
  fwrite(summary_dt, paste0(opt$out, ".summary.csv"))
  fwrite(per_reg_dt[order(-abs_diff)], paste0(opt$out, ".per_regulon.csv"))
  cat("\nWritten:\n  ", paste0(opt$out, ".summary.csv"), "\n  ",
      paste0(opt$out, ".per_regulon.csv"), "\n", sep="")
}

