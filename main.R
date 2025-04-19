# メインスクリプト
# このファイルでは、ネットワークメタアナリシスとメタ回帰分析のワークフローを実行します

# 環境セットアップの確認
if (!file.exists("setup.R")) {
  stop("setup.Rファイルが見つかりません。正しいディレクトリにいるか確認してください。")
}
source("setup.R")

# ライブラリとプロジェクトソースの読み込み
message("必要なライブラリとプロジェクトソースを読み込んでいます...")
library(here)  # パス管理用
library(yaml)  # 設定読み込み用
library(logger) # 構造化ログ用
library(fs)    # ファイル操作用

# 関数の読み込み
source(here("R", "01_utils.R"))
source(here("R", "02_data_prep.R"))
source(here("R", "03_models.R"))
source(here("R", "04_visualization.R"))

# 設定の読み込み
config <- yaml::read_yaml(here("config", "config.yml"))
message("設定読み込み完了")

# タイムスタンプディレクトリ作成
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
results_dir <- path(here("results"), timestamp)
dir_create(results_dir, recurse = TRUE)
message(paste("結果ディレクトリを作成しました:", results_dir))

# 最新結果へのシンボリックリンク更新
latest_link <- path(here("results"), "latest")
if (dir_exists(latest_link)) dir_delete(latest_link)
if (Sys.info()["sysname"] == "Windows") {
  # Windowsではシンボリックリンクを作成せず、バッチファイルを作成
  batch_content <- paste0("@echo off\ncd ", normalizePath(results_dir))
  writeLines(batch_content, path(here("results"), "latest.bat"))
  message("Windowsでは'latest.bat'ファイルを作成しました。")
} else {
  # Unix系OSではシンボリックリンクを作成
  link_create(results_dir, latest_link)
  message("最新結果へのシンボリックリンクを更新しました。")
}

# ロガー設定
log_appender(appender_file(path(results_dir, "run.log")))
log_threshold(INFO)
log_info("分析開始")

# 設定ファイルをコピー（記録用）
file_copy(here("config", "config.yml"), path(results_dir, "config_used.yml"))

# 乱数シードの設定（再現性のため）
set.seed(config$seed)
log_info(paste("乱数シード設定:", config$seed))

# データ準備
log_info("データ準備開始")
data <- prepare_data(
  file_path = here("data", "raw", config$data$filename),
  markers = config$data$markers
)
saveRDS(data, path(results_dir, "prepared_data.rds"))
log_info("データ準備完了")

# NMAモデル実行
log_info("ネットワークメタアナリシス実行開始")
nma_results <- run_nma_model(
  data = data,
  iter = config$nma$iter,
  chains = config$nma$chains,
  cores = config$nma$cores,
  adapt_delta = config$nma$adapt_delta,
  stepsize = config$nma$stepsize,
  max_treedepth = config$nma$max_treedepth,
  results_dir = results_dir
)
saveRDS(nma_results, path(results_dir, "nma_results.rds"))

# NMA結果のCSV出力
write.csv(nma_results$mu, path(results_dir, "nma_sensitivity_specificity.csv"))
write.csv(nma_results$dor, path(results_dir, "nma_dor.csv"))
write.csv(nma_results$sindex, path(results_dir, "nma_sindex.csv"))
log_info("ネットワークメタアナリシス実行完了")

# メタ回帰モデル実行（設定で有効になっている場合）
regression_results <- NULL
if (config$meta_regression$run) {
  log_info("メタ回帰モデル実行開始")
  regression_results <- run_regression_model(
    data = data,
    covariate = config$meta_regression$covariate,
    iter = config$meta_regression$iter,
    chains = config$meta_regression$chains,
    cores = config$meta_regression$cores,
    adapt_delta = config$meta_regression$adapt_delta,
    stepsize = config$meta_regression$stepsize,
    max_treedepth = config$meta_regression$max_treedepth,
    results_dir = results_dir
  )
  saveRDS(regression_results, path(results_dir, "regression_results.rds"))
  
  # メタ回帰結果のCSV出力
  write.csv(regression_results$mu, path(results_dir, "regression_sensitivity_specificity.csv"))
  write.csv(regression_results$dor, path(results_dir, "regression_dor.csv"))
  write.csv(regression_results$sindex, path(results_dir, "regression_sindex.csv"))
  write.csv(regression_results$beta1, path(results_dir, "regression_coefficients.csv"))
  log_info("メタ回帰モデル実行完了")
}

# 結果の可視化
log_info("可視化作成開始")
plots <- create_visualizations(
  nma_results = nma_results, 
  regression_results = regression_results, 
  reference_marker = config$data$reference_marker
)

# プロットを保存
for (i in seq_along(plots)) {
  plot_name <- names(plots)[i]
  plot_path <- path(results_dir, paste0(plot_name, ".png"))
  ggsave(plot_path, plots[[i]], width = 10, height = 8, dpi = 300)
  log_info(paste("プロット保存:", plot_name))
}
log_info("可視化作成完了")

# 実行情報を記録
run_info <- list(
  timestamp = Sys.time(),
  r_version = R.version.string,
  platform = Sys.info()["sysname"],
  packages = sapply(required_packages, function(pkg) packageVersion(pkg)),
  config = config,
  execution_time = list(
    nma = nma_results$execution_info$elapsed_time
  )
)

if (!is.null(regression_results)) {
  run_info$execution_time$regression <- regression_results$execution_info$elapsed_time
}

# JSONとRDSで保存
jsonlite::write_json(run_info, path(results_dir, "run_info.json"), pretty = TRUE, auto_unbox = TRUE)
saveRDS(run_info, path(results_dir, "run_info.rds"))
log_info("実行情報を記録しました")

# 古い結果の削除（設定で有効になっている場合）
if (config$results$prune_old_results) {
  log_info("古い結果を削除しています...")
  pruned_count <- prune_old_results(
    results_path = here("results"), 
    keep_last = config$results$keep_last_n_results
  )
  log_info(paste("削除された古い結果ディレクトリ数:", pruned_count))
}

log_info("分析完了")
message("分析が完了しました。結果は以下のディレクトリに保存されています:")
message(results_dir)
