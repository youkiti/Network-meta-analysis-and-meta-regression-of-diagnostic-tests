# 可視化関数
# このファイルには、結果の可視化関連の関数が含まれます

#' 感度・特異度のフォレストプロットを作成する関数
#'
#' @param results get_nma_results関数で得られた結果のリスト
#' @param title プロットのタイトル
#' @param reference_marker リファレンスマーカー名
#' @return ggplotオブジェクト
create_forest_plot <- function(results, title = "Forest Plot of Sensitivity and Specificity", reference_marker = "qSOFA") {
  # マーカー情報を取得
  markers <- results$marker
  
  # 結果のデータフレームを作成
  df <- data.frame(
    Marker = rep(markers, 2),
    Type = c(rep("Sensitivity", length(markers)), rep("Specificity", length(markers))),
    Mean = c(results$mu$mean_se, results$mu$mean_sp),
    Lower = c(results$mu$lb_se, results$mu$lb_sp),
    Upper = c(results$mu$ub_se, results$mu$ub_sp),
    stringsAsFactors = FALSE
  )
  
  # リファレンスマーカーに印をつける
  df$Reference <- ifelse(df$Marker == reference_marker, "Yes", "No")
  
  # マーカーを結果の順にソート
  df$Marker <- factor(df$Marker, levels = rev(markers))
  
  # プロット作成
  plot <- ggplot2::ggplot(df, ggplot2::aes(x = Marker, y = Mean, color = Type)) +
    ggplot2::geom_point(position = ggplot2::position_dodge(width = 0.5), 
               size = 3, 
               ggplot2::aes(shape = Reference)) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = Lower, ymax = Upper), 
                  width = 0.25,
                  position = ggplot2::position_dodge(width = 0.5)) +
    ggplot2::coord_flip() +
    ggplot2::scale_color_manual(values = c("Sensitivity" = "blue", "Specificity" = "red")) +
    ggplot2::scale_shape_manual(values = c("Yes" = 18, "No" = 16)) +
    ggplot2::labs(title = title,
         y = "Proportion",
         x = "Biomarker") +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "bottom")
  
  return(plot)
}

#' DORのプロットを作成する関数
#'
#' @param results get_nma_results関数で得られた結果のリスト
#' @param title プロットのタイトル
#' @param reference_marker リファレンスマーカー名
#' @return ggplotオブジェクト
create_dor_plot <- function(results, title = "Diagnostic Odds Ratio (DOR)", reference_marker = "qSOFA") {
  # マーカー情報を取得
  markers <- results$marker
  
  # 結果のデータフレームを作成
  df <- data.frame(
    Marker = markers,
    DOR = results$dor$mean,
    Lower = results$dor$lb,
    Upper = results$dor$ub,
    stringsAsFactors = FALSE
  )
  
  # リファレンスマーカーに印をつける
  df$Reference <- ifelse(df$Marker == reference_marker, "Yes", "No")
  
  # DORの値でソート
  df <- df[order(df$DOR, decreasing = TRUE), ]
  df$Marker <- factor(df$Marker, levels = df$Marker)
  
  # プロット作成
  plot <- ggplot2::ggplot(df, ggplot2::aes(x = Marker, y = DOR)) +
    ggplot2::geom_point(ggplot2::aes(shape = Reference), size = 3, color = "darkblue") +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = Lower, ymax = Upper), width = 0.25, color = "darkblue") +
    ggplot2::coord_flip() +
    ggplot2::scale_shape_manual(values = c("Yes" = 18, "No" = 16)) +
    ggplot2::labs(title = title,
         y = "Diagnostic Odds Ratio",
         x = "Biomarker") +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "bottom")
  
  return(plot)
}

#' サマリーインデックスのプロットを作成する関数
#'
#' @param results get_nma_results関数で得られた結果のリスト
#' @param title プロットのタイトル
#' @param reference_marker リファレンスマーカー名
#' @return ggplotオブジェクト
create_sindex_plot <- function(results, title = "Summary Index", reference_marker = "qSOFA") {
  # マーカー情報を取得
  markers <- results$marker
  
  # 結果のデータフレームを作成
  df <- data.frame(
    Marker = rep(markers, 3),
    Type = c(rep("Overall", length(markers)), rep("Sensitivity", length(markers)), rep("Specificity", length(markers))),
    Index = c(results$sindex$mean, results$sindex_se$mean, results$sindex_sp$mean),
    Lower = c(results$sindex$lb, results$sindex_se$lb, results$sindex_sp$lb),
    Upper = c(results$sindex$ub, results$sindex_se$ub, results$sindex_sp$ub),
    stringsAsFactors = FALSE
  )
  
  # リファレンスマーカーに印をつける
  df$Reference <- ifelse(df$Marker == reference_marker, "Yes", "No")
  
  # マーカーを結果の順にソート（Overallの値で）
  overall_index <- df[df$Type == "Overall", c("Marker", "Index")]
  overall_index <- overall_index[order(overall_index$Index, decreasing = TRUE), ]
  df$Marker <- factor(df$Marker, levels = overall_index$Marker)
  
  # プロット作成
  plot <- ggplot2::ggplot(df, ggplot2::aes(x = Marker, y = Index, color = Type)) +
    ggplot2::geom_point(position = ggplot2::position_dodge(width = 0.5), 
               size = 3, 
               ggplot2::aes(shape = Reference)) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = Lower, ymax = Upper), 
                  width = 0.25,
                  position = ggplot2::position_dodge(width = 0.5)) +
    ggplot2::coord_flip() +
    ggplot2::scale_color_manual(values = c("Overall" = "purple", "Sensitivity" = "blue", "Specificity" = "red")) +
    ggplot2::scale_shape_manual(values = c("Yes" = 18, "No" = 16)) +
    ggplot2::labs(title = title,
         y = "Superiority Index",
         x = "Biomarker") +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "bottom")
  
  return(plot)
}

#' 可視化関数のラッパー：すべてのプロットを作成
#'
#' @param nma_results get_nma_results関数で得られたNMA結果のリスト
#' @param regression_results get_regression_results関数で得られたメタ回帰結果のリスト（オプション）
#' @param reference_marker リファレンスマーカー名
#' @return プロットのリスト
create_visualizations <- function(nma_results, regression_results = NULL, reference_marker = "qSOFA") {
  # プロットのリストを初期化
  plots <- list()
  
  # NMAの結果をプロット
  plots$nma_forest <- create_forest_plot(
    nma_results, 
    title = "Forest Plot of Sensitivity and Specificity (NMA)",
    reference_marker = reference_marker
  )
  
  plots$nma_dor <- create_dor_plot(
    nma_results, 
    title = "Diagnostic Odds Ratio (NMA)",
    reference_marker = reference_marker
  )
  
  plots$nma_sindex <- create_sindex_plot(
    nma_results, 
    title = "Superiority Index (NMA)",
    reference_marker = reference_marker
  )
  
  # メタ回帰の結果をプロット（もしあれば）
  if (!is.null(regression_results)) {
    plots$reg_forest <- create_forest_plot(
      regression_results, 
      title = paste0("Forest Plot of Sensitivity and Specificity (Regression: ", 
                    regression_results$execution_info$covariate, ")"),
      reference_marker = reference_marker
    )
    
    plots$reg_dor <- create_dor_plot(
      regression_results, 
      title = paste0("Diagnostic Odds Ratio (Regression: ", 
                   regression_results$execution_info$covariate, ")"),
      reference_marker = reference_marker
    )
    
    plots$reg_sindex <- create_sindex_plot(
      regression_results, 
      title = paste0("Superiority Index (Regression: ", 
                    regression_results$execution_info$covariate, ")"),
      reference_marker = reference_marker
    )
    
    # 回帰係数のプロット（beta1）を追加
    plots$beta_coef <- create_beta_coefficient_plot(
      regression_results,
      title = paste0("Regression Coefficients (", 
                    regression_results$execution_info$covariate, ")"),
      reference_marker = reference_marker
    )
  }
  
  return(plots)
}

#' メタ回帰係数のプロットを作成する関数
#'
#' @param results get_regression_results関数で得られた結果のリスト
#' @param title プロットのタイトル
#' @param reference_marker リファレンスマーカー名
#' @return ggplotオブジェクト
create_beta_coefficient_plot <- function(results, title = "Regression Coefficients", reference_marker = "qSOFA") {
  # マーカー情報を取得
  markers <- results$marker
  
  # 結果のデータフレームを作成
  df <- data.frame(
    Marker = rep(markers, 2),
    Type = c(rep("Sensitivity", length(markers)), rep("Specificity", length(markers))),
    Coefficient = c(results$beta1$mean_se, results$beta1$mean_sp),
    Lower = c(results$beta1$lb_se, results$beta1$lb_sp),
    Upper = c(results$beta1$ub_se, results$beta1$ub_sp),
    stringsAsFactors = FALSE
  )
  
  # リファレンスマーカーに印をつける
  df$Reference <- ifelse(df$Marker == reference_marker, "Yes", "No")
  
  # マーカーを結果の順にソート
  df$Marker <- factor(df$Marker, levels = rev(markers))
  
  # 0の基準線を示す
  plot <- ggplot2::ggplot(df, ggplot2::aes(x = Marker, y = Coefficient, color = Type)) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
    ggplot2::geom_point(position = ggplot2::position_dodge(width = 0.5), 
               size = 3, 
               ggplot2::aes(shape = Reference)) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = Lower, ymax = Upper), 
                  width = 0.25,
                  position = ggplot2::position_dodge(width = 0.5)) +
    ggplot2::coord_flip() +
    ggplot2::scale_color_manual(values = c("Sensitivity" = "blue", "Specificity" = "red")) +
    ggplot2::scale_shape_manual(values = c("Yes" = 18, "No" = 16)) +
    ggplot2::labs(title = title,
         y = "Regression Coefficient",
         x = "Biomarker") +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "bottom")
  
  return(plot)
}
