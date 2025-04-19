# 環境セットアップスクリプト
# このファイルでは、必要なパッケージのインストールと環境の準備を行います

# CRANミラーの設定
options(repos = c(CRAN = "https://cloud.r-project.org"))

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

# Rtools44の設定 - 新しいパス構造に対応
rtools_home <- "C:/rtools44"
Sys.setenv(RTOOLS44_HOME = rtools_home)

# Windows上のパスは\\を使う必要があるため、こちらのようにパスを設定
rtools_path <- paste0(rtools_home, "\\x86_64-w64-mingw32.static.posix\\bin;", 
                     rtools_home, "\\usr\\bin")
current_path <- Sys.getenv("PATH")
Sys.setenv(PATH = paste(rtools_path, current_path, sep = ";"))
Sys.setenv(BINPREF = "C:/rtools44/x86_64-w64-mingw32.static.posix/bin/")

# Rtoolsが正しく設定されているか確認
message("Checking Rtools configuration...")
if (!requireNamespace("pkgbuild", quietly = TRUE)) {
  message("Installing pkgbuild package...")
  install.packages("pkgbuild")
}
library(pkgbuild)

# パスの確認
message("PATH environment variable:")
message(Sys.getenv("PATH"))

# g++が見つかるか確認
message("g++ location:")
message(Sys.which("g++"))

# g++のバージョン確認
message("g++ version:")
system("g++ --version", intern = TRUE)

# Rtoolsの検出確認
message("Checking if Rtools is found:")
has_rtools <- pkgbuild::has_rtools()
message(paste("has_rtools() returns:", has_rtools))

# Rcppなどの重要パッケージを再インストール
message("Reinstalling critical C++ packages...")
install.packages(c("Rcpp", "RcppEigen", "BH", "StanHeaders", "rstan"), 
                 repos = "https://cloud.r-project.org",
                 dependencies = TRUE)

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

# ロケール設定
# Windows環境ではUTF-8を使用
Sys.setlocale("LC_ALL", "English_United States.1252")

message("Setup completed successfully.")
message("Please run main.R next.")
