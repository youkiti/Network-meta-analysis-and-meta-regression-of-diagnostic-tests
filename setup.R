# 環境セットアップスクリプト
# このファイルでは、必要なパッケージのインストールと環境の準備を行います

# 必要パッケージのリスト
required_packages <- c(
  # 基本パッケージ
  "rstan", "bayesplot", "ggplot2", "doParallel", "devtools",
  # パス管理・設定管理
  "here", "yaml", 
  # ログ・ファイル操作
  "logger", "fs", 
  # その他のユーティリティ
  "jsonlite", "digest"
)

# パッケージのインストールチェック
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message(paste("Installing package:", pkg))
    install.packages(pkg)
  }
  library(pkg, character.only = TRUE)
}

# 環境変数の設定
Sys.setenv(PATH = paste("C:/rtools40/mingw64/bin/", Sys.getenv("PATH"), sep = ";"))
Sys.setenv(BINPREF = "C:/rtools40/mingw64/bin/")

# Stan設定
message("Setting up rstan options...")
rstan::rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# プロジェクトのディレクトリ構造を確認
required_dirs <- c(
  "R", "config", "data/raw", "data/processed", "results"
)

for (dir in required_dirs) {
  if (!dir.exists(dir)) {
    message(paste("Creating directory:", dir))
    dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  }
}

# 日本語ロケールの設定
Sys.setlocale("LC_ALL", "Japanese")

message("環境のセットアップが完了しました。")
message("続いてmain.Rを実行してください。")
