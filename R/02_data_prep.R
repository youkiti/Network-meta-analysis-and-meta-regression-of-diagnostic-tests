# データ準備関数
# このファイルには、生データを読み込み、分析用に前処理するための関数が含まれます

#' ネットワークメタアナリシスのためのデータを準備する関数
#'
#' @param file_path データファイルのパス
#' @param markers 対象のマーカー名のベクトル
#' @return データと前処理されたパラメータのリスト
prepare_data <- function(file_path, markers) {
  # 元のNMA_code.rのgetpara関数と同様にシンプルに読み込み
  message(paste("Reading data from:", file_path))
  data.mk <- read.csv(file_path)
  message("Data file read successfully")
  
  names.study <- unique(fixlen.str(data.mk$Study))
  hash.study <- 1:length(names.study)
  names(hash.study) <- names.study
  
  hash.mk <- 1:length(markers)
  names(hash.mk) <- markers
  
  para <- data.frame(
    study.names = fixlen.str(data.mk$Study),
    Study = array(hash.study[as.character(fixlen.str(data.mk$Study))]),
    TP = round(data.mk$Case * data.mk$Sensitivity),
    TN = round(data.mk$Control * data.mk$Specificity),
    Dis = data.mk$Case,
    NDis = data.mk$Control,
    Test = array(hash.mk[as.character(data.mk$Marker)])
  )
  
  # NMAモデル用のパラメータリストを作成
  para.as.ls.nma <- list(
    N = nrow(para),
    Nt = length(unique(para$Test)),
    Ns = length(unique(para$Study)),
    TP = para$TP,
    Dis = para$Dis,
    TN = para$TN,
    NDis = para$NDis,
    Study = as.numeric(para$Study),
    Test = para$Test
  )
  
  # メタ回帰モデル用のパラメータリストを作成
  para.as.ls.1mk <- list(
    N = nrow(para),
    Nt = length(unique(para$Test)),
    Ns = length(unique(para$Study)),
    TP = para$TP,
    Dis = para$Dis,
    TN = para$TN,
    NDis = para$NDis,
    Study = as.numeric(para$Study),
    Test = para$Test,
    Covar1 = data.mk$is.prevalence05
  )
  
  return(list(
    data = data.mk,
    para = para,
    para.as.ls.nma = para.as.ls.nma,
    para.as.ls.1mk = para.as.ls.1mk,
    markers = markers[unique(para$Test)]
  ))
}
