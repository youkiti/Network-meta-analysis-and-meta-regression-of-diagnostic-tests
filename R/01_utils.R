# ユーティリティ関数
# このファイルには、固定長文字列処理や結果抽出などの共通ユーティリティ関数が含まれます

#' 文字列配列を固定長に変換する関数
#'
#' @param string.arr 文字列配列
#' @return 固定長に変換された文字列配列
fixlen.str <- function(string.arr = NULL) {
  max.len <- max(nchar(as.character(string.arr)))
  rule <- paste('% -', as.character(max.len), 's', sep = '')
  string.arr.fix <- sprintf(rule, string.arr)
  return(string.arr.fix)
}

#' 95%信頼区間の上限を取得する関数
#'
#' @param x 数値ベクトル
#' @return 95%信頼区間の上限値
get.ub.95 <- function(x) {
  y <- sort(x)
  len <- length(x)
  return(y[round(len * 0.975)])
}

#' 95%信頼区間の下限を取得する関数
#'
#' @param x 数値ベクトル
#' @return 95%信頼区間の下限値
get.lb.95 <- function(x) {
  y <- sort(x)
  len <- length(x)
  return(y[round(len * 0.025)])
}

#' Stan modelの結果を抽出して保存する関数
#'
#' @param fit.data Stanモデルのfit結果
#' @param filename.save 保存するファイル名
#' @return 抽出されたシミュレーションデータ
fit.extract <- function(fit.data = NULL, filename.save = NULL) {
  results.sum <- summary(fit.data)[[1]]
  results.simdata <- fit.data@sim$samples[[1]]
  save(results.sum, results.simdata, file = filename.save)
  return(results.simdata)
}

#' 古い結果ディレクトリを削除する関数
#'
#' @param results_path 結果ディレクトリのパス
#' @param keep_last 保持する最新の結果ディレクトリ数
#' @return 削除したディレクトリ数
prune_old_results <- function(results_path, keep_last = 5) {
  # タイムスタンプディレクトリのみをリスト
  dirs <- list.dirs(results_path, recursive = FALSE)
  # "latest"などの特殊ディレクトリを除外
  pattern <- "^\\d{8}_\\d{6}$"
  timestamp_dirs <- dirs[grepl(pattern, basename(dirs))]
  
  # 日付順にソート
  timestamp_dirs <- sort(timestamp_dirs, decreasing = TRUE)
  
  # 保持する数を超えた分を削除
  if (length(timestamp_dirs) > keep_last) {
    dirs_to_remove <- timestamp_dirs[(keep_last + 1):length(timestamp_dirs)]
    for (dir in dirs_to_remove) {
      unlink(dir, recursive = TRUE)
    }
    return(length(dirs_to_remove))
  }
  
  return(0)
}
