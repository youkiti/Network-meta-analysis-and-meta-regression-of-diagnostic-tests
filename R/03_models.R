# モデル定義と実行関数
# このファイルにはStanモデルコードとモデル実行関数が含まれます

# ネットワークメタアナリシスのStanコード
STAN_CODE_NMA <- '
data{
     int N;  //number of comparison??? - 121
     int Nt; //number of test - 10
     int Ns; //number of study - 72
     int TP[N];
     int Dis[N];  //diseased
     int TN[N];
     int NDis[N]; //non-diseased
     int Study[N];
     int Test[N];
}
parameters{
           matrix[2, Nt] logitmu;
           vector[Ns] nu[2];
           matrix[Ns, Nt] delta[2];
           vector<lower=0>[Nt] tau[2]; //*
           vector<lower=0>[2] sigmab; 
           real<lower=-1, upper=1> rho;
}
transformed parameters{
                       matrix[Ns, 2] p_i[Nt];
                       matrix[2, Nt] MU;
                       matrix[2, Nt] RR;
                       matrix[2, Nt] OR;
                       vector[Nt] DOR;
                       vector[Nt] S;
                       vector[Nt] S_se;
                       vector[Nt] S_sp;
                       matrix[Nt, Nt] A;
                       matrix[Nt, Nt] B;
                       matrix[Nt, Nt] C;
                       matrix[Nt, Nt] A_se;
                       matrix[Nt, Nt] B_se;
                       matrix[Nt, Nt] C_se;
                       matrix[Nt, Nt] A_sp;
                       matrix[Nt, Nt] B_sp;
                       matrix[Nt, Nt] C_sp;
                       vector<lower=0>[Nt] tausq[2];
                       vector<lower=0>[2] sigmabsq;
                       matrix[Nt, Nt] sigmasq[2];
                       matrix[Nt, Nt] rhow[2];

    for (i in 1:Ns){
        for (j in 1:2){
            for (k in 1:Nt)
                p_i[k][i,j] = inv_logit(logitmu[j,k] +  nu[j][i] + delta[j][i,k]);
        }
    }
 
    for (j in 1:2){
        for (k in 1:Nt){
                         MU[j,k] = mean(col(p_i[k], j));
        }
        tausq[j] = (tau[j]).*(tau[j]);
    }

    for (j in 1:2){
        for (k in 1:Nt){
                         RR[j, k] = MU[j, k]/MU[j, 1]; 
                         OR[j, k] = (MU[j, k]*(1 - MU[j, 1]))/(MU[j, 1]*(1 - MU[j, k]));
         }
    }

    for (l in 1:Nt){
                     DOR[l] = (MU[1, l]*MU[2, l])/((1 - MU[1, l])*(1 - MU[2, l]));

        for(m in 1:Nt){
                        A[l, m] = if_else((MU[1, l] > MU[1, m]) && (MU[2, l] > MU[2, m]), 1, 0);
                        B[l, m] = if_else((MU[1, l] < MU[1, m]) && (MU[2, l] < MU[2, m]), 1, 0);
                        C[l, m] = if_else((MU[1, l] == MU[1, m]) && (MU[2, l] == MU[2, m]), 1, 0);

                        A_se[l, m] = if_else((MU[1, l] > MU[1, m]), 1, 0);
                        B_se[l, m] = if_else((MU[1, l] < MU[1, m]), 1, 0);
                        C_se[l, m] = if_else((MU[1, l] == MU[1, m]), 1, 0);

                        A_sp[l, m] = if_else((MU[2, l] > MU[2, m]), 1, 0);
                        B_sp[l, m] = if_else((MU[2, l] < MU[2, m]), 1, 0);
                        C_sp[l, m] = if_else((MU[2, l] == MU[2, m]), 1, 0);
        }

        S[l] = (2*sum(row(A, l)) + sum(row(C, l)))/(2*sum(row(B, l)) + sum(row(C, l)));
        S_se[l] = (2*sum(row(A_se, l)) + sum(row(C_se, l)))/(2*sum(row(B_se, l)) + sum(row(C_se, l)));
        S_sp[l] = (2*sum(row(A_sp, l)) + sum(row(C_sp, l)))/(2*sum(row(B_sp, l)) + sum(row(C_sp, l)));
    }
    
    sigmabsq = (sigmab).*(sigmab);

    for (j in 1:2){
        for (k in 1:Nt){
            for (l in 1:Nt){
                             sigmasq[j][k,l] = (sigmabsq[j] + tausq[j][k])*((sigmabsq[j] + tausq[j][l]));
                             rhow[j][k,l] = sigmabsq[j]/sqrt(sigmasq[j][k,l]);
            }
        }
    }

}
model{
	   //Priors
       for (j in 1:2){
                       logitmu[j] ~ normal(0, 5);
		               tau[j] ~ cauchy(0, 2.5);
       }

         sigmab ~ cauchy(0, 2.5);
	     rho ~ uniform(-1, 1);
         nu[2] ~ normal(0, sigmab[2]);
      
       for (i in 1:Ns){
                       nu[1][i] ~ normal((sigmab[1]/sigmab[2])*rho*nu[2][i], sqrt(sigmabsq[1]*(1 - (rho*rho))));
          for (j in 1:2){
              for (k in 1:Nt)
                              delta[j][i,k] ~ normal(0, tau[j][k]);
        }
    }

    for (n in 1:N){
                    TP[n] ~ binomial(Dis[n], p_i[Test[n]][Study[n], 1]);
                    TN[n] ~ binomial(NDis[n], p_i[Test[n]][Study[n], 2]);
    }

}
generated quantities{
    
    vector[2*N] loglik;

    for (n in 1:N)
                   loglik[n] = binomial_lpmf(TN[n]| NDis[n], p_i[Test[n]][Study[n], 1]);

    for (n in (N+1):(2*N))
                   loglik[n] = binomial_lpmf(TN[n-N]| NDis[n-N], p_i[Test[n-N]][Study[n-N], 2]);

}
'

# メタ回帰モデルのStanコード
STAN_CODE_1MK <- '
data{
      int N;         //number of comparison (n=121)
      int Nt;        //number of test (n=7)
      int Ns;        //number of study (n=107)
      int TP[N];
      int Dis[N];    //diseased
      int TN[N];
      int NDis[N];   //non-diseased
      int Study[N];
      int Test[N];
      int Covar1[N]; //one additional variable
    }

parameters{
            matrix[2,Nt] logitmu;
            vector[Ns] nu[2];
            matrix[Ns,Nt] delta[2];
            vector<lower=0>[Nt] tau[2];
            vector<lower=0>[2] sigmab; 
            real<lower=-1,upper=1> rho;
            matrix[2,Nt] beta1;
        }

transformed parameters{
                        matrix[Ns,2] p_i[Nt];
                        matrix[2,Nt] MU;
                        matrix[2,Nt] RR;
                        matrix[2,Nt] OR;
                        vector[Nt] DOR;
                        vector[Nt] S;
                        vector[Nt] S_se;
                        vector[Nt] S_sp;
                        matrix[Nt,Nt] A;
                        matrix[Nt,Nt] B;
                        matrix[Nt,Nt] C;
                        matrix[Nt, Nt] A_se;
                        matrix[Nt, Nt] B_se;
                        matrix[Nt, Nt] C_se;
                        matrix[Nt, Nt] A_sp;
                        matrix[Nt, Nt] B_sp;
                        matrix[Nt, Nt] C_sp;

                        vector<lower=0>[Nt] tausq[2];
                        vector<lower=0>[2] sigmabsq;

                        matrix[Nt,Nt] sigmasq[2];
                        matrix[Nt,Nt] rhow[2];

                for (i in 1:Ns){
                     for (j in 1:2){
                          for (k in 1:Nt){
                                           p_i[k][i,j]=inv_logit(logitmu[j,k]+nu[j][i]+delta[j][i,k]+beta1[j,k]*Covar1[i]);
                           }
                      }
                 }

                for (j in 1:2){
                     for (k in 1:Nt){
                                      MU[j,k]=mean(col(p_i[k],j));
                     }
                      tausq[j]=(tau[j]).*(tau[j]);
                }
                  
                for (j in 1:2){
                     for (k in 1:Nt){
                                      RR[j,k]=MU[j,k]/MU[j,1]; 
                                      OR[j,k]=(MU[j,k]*(1-MU[j,1]))/(MU[j,1]*(1-MU[j,k]));
                       }
                 }

                for (l in 1:Nt){
                                 DOR[l]=(MU[1,l]*MU[2,l])/((1-MU[1,l])*(1-MU[2,l]));
                     for(m in 1:Nt){
                                     A[l,m]=if_else((MU[1,l]>MU[1,m]) && (MU[2,l]>MU[2,m]),1,0);
                                     B[l,m]=if_else((MU[1,l]<MU[1,m]) && (MU[2,l]<MU[2,m]),1,0);
                                     C[l,m]=if_else((MU[1,l]==MU[1,m]) && (MU[2,l]==MU[2,m]),1,0);
                                     
                                     A_se[l, m] = if_else((MU[1, l] > MU[1, m]), 1, 0);
                                     B_se[l, m] = if_else((MU[1, l] < MU[1, m]), 1, 0);
                                     C_se[l, m] = if_else((MU[1, l] == MU[1, m]), 1, 0);
 
                                     A_sp[l, m] = if_else((MU[2, l] > MU[2, m]), 1, 0);
                                     B_sp[l, m] = if_else((MU[2, l] < MU[2, m]), 1, 0);
                                     C_sp[l, m] = if_else((MU[2, l] == MU[2, m]), 1, 0);
                      }
                        S[l]=(2*sum(row(A,l))+sum(row(C,l)))/(2*sum(row(B,l))+sum(row(C,l)));
                        S_se[l] = (2*sum(row(A_se, l)) + sum(row(C_se, l)))/(2*sum(row(B_se, l)) + sum(row(C_se, l)));
                        S_sp[l] = (2*sum(row(A_sp, l)) + sum(row(C_sp, l)))/(2*sum(row(B_sp, l)) + sum(row(C_sp, l)));
                }
        
                 sigmabsq=(sigmab).*(sigmab);

                for (j in 1:2){
                     for (k in 1:Nt){
                          for (l in 1:Nt){
                                           sigmasq[j][k,l]=(sigmabsq[j]+tausq[j][k])*((sigmabsq[j]+tausq[j][l]));
                                           rhow[j][k,l]=sigmabsq[j]/sqrt(sigmasq[j][k,l]);
                            }
                       }
                 }
            }

                model{
                      //Priors
                      for (j in 1:2){
                                      logitmu[j]~normal(0,5);
                                      tau[j]~cauchy(0,2.5);
                                      
                           for (k in 1:Nt){
                                            beta1[j,k]~normal(0,5);
                           }            
                      }
                       sigmab~cauchy(0, 2.5);
                       rho~uniform(-1, 1);
                       nu[2]~normal(0,sigmab[2]);

                      for (i in 1:Ns){
                                       nu[1][i]~normal((sigmab[1]/sigmab[2])*rho*nu[2][i],sqrt(sigmabsq[1]*(1-(rho*rho))));
                           for (j in 1:2){
                               for (k in 1:Nt){
                                                delta[j][i,k]~normal(0,tau[j][k]);
                              }
                         }
                      }

                     for (n in 1:N){
                                     TP[n]~binomial(Dis[n],p_i[Test[n]][Study[n],1]);
                                     TN[n]~binomial(NDis[n],p_i[Test[n]][Study[n],2]);
                      }

                    }

                generated quantities{
                                      vector[2*N] loglik;

                                      for (n in 1:N){
                                                      loglik[n]=binomial_lpmf(TN[n]|NDis[n],p_i[Test[n]][Study[n],1]);
                                       }
                                      for (n in (N+1):(2*N)){
                                                              loglik[n]=binomial_lpmf(TN[n-N]|NDis[n-N],p_i[Test[n-N]][Study[n-N],2]);
                                       }
                }
'

#' ネットワークメタアナリシスモデルを実行する関数
#'
#' @param data prepare_data関数で準備されたデータリスト
#' @param iter MCMCの反復回数
#' @param chains MCMCのチェーン数
#' @param cores 使用するCPUコア数
#' @param adapt_delta 適応パラメータ
#' @param stepsize ステップサイズ
#' @param max_treedepth 最大ツリー深度
#' @param results_dir 結果を保存するディレクトリ
#' @return 分析結果のリスト
run_nma_model <- function(data, iter = 20000, chains = 4, cores = 4, adapt_delta = 0.99, 
                          stepsize = 0.01, max_treedepth = 17, results_dir = NULL) {
  
  # モデルをコンパイル
  stan.model.nma <- rstan::stan(
    model_code = STAN_CODE_NMA,
    data = data$para.as.ls.nma, 
    iter = 1, 
    warmup = 0, 
    chains = 2
  )
  
  # モデルを実行
  start_time <- Sys.time()
  fit.nma <- rstan::stan(
    fit = stan.model.nma,
    data = data$para.as.ls.nma,
    iter = iter,
    chains = chains,
    cores = cores,
    control = list(
      adapt_delta = adapt_delta,
      stepsize = stepsize,
      max_treedepth = max_treedepth
    )
  )
  end_time <- Sys.time()
  elapsed_time <- end_time - start_time
  
  # 結果を保存
  if (!is.null(results_dir)) {
    filename.simdata.save <- file.path(results_dir, 'fit_nma_simdata.RData')
  } else {
    filename.simdata.save <- 'fit_nma_simdata.RData'
  }
  
  results.nma <- fit.extract(fit.data = fit.nma, filename.save = filename.simdata.save)
  
  # 結果を抽出
  results.nma.95 <- get_nma_results(
    data.ls = list(
      mk = data$markers,
      data.nma = results.nma
    ),
    mk.names = data$markers,
    conf = 95
  )
  
  # 実行情報を追加
  results.nma.95$execution_info <- list(
    elapsed_time = elapsed_time,
    iter = iter,
    chains = chains,
    cores = cores,
    adapt_delta = adapt_delta,
    stepsize = stepsize,
    max_treedepth = max_treedepth
  )
  
  return(results.nma.95)
}

#' メタ回帰モデルを実行する関数
#'
#' @param data prepare_data関数で準備されたデータリスト
#' @param covariate 使用する共変量変数名
#' @param iter MCMCの反復回数
#' @param chains MCMCのチェーン数
#' @param cores 使用するCPUコア数
#' @param adapt_delta 適応パラメータ
#' @param stepsize ステップサイズ
#' @param max_treedepth 最大ツリー深度
#' @param results_dir 結果を保存するディレクトリ
#' @return 分析結果のリスト
run_regression_model <- function(data, covariate = "is.prevalence05", iter = 20000, chains = 4, 
                                cores = 4, adapt_delta = 0.99, stepsize = 0.01, 
                                max_treedepth = 17, results_dir = NULL) {
  
  # モデルをコンパイル
  stan.model.1mk <- rstan::stan(
    model_code = STAN_CODE_1MK,
    data = data$para.as.ls.1mk, 
    iter = 1, 
    warmup = 0, 
    chains = 1
  )
  
  # モデルを実行
  start_time <- Sys.time()
  fit.1mk <- rstan::stan(
    fit = stan.model.1mk,
    data = data$para.as.ls.1mk,
    iter = iter,
    chains = chains,
    cores = cores,
    control = list(
      adapt_delta = adapt_delta,
      stepsize = stepsize,
      max_treedepth = max_treedepth
    )
  )
  end_time <- Sys.time()
  elapsed_time <- end_time - start_time
  
  # 結果を保存
  if (!is.null(results_dir)) {
    filename.simdata.save <- file.path(results_dir, 'fit_1mk_simdata.RData')
  } else {
    filename.simdata.save <- 'fit_1mk_simdata.RData'
  }
  
  results.1mk <- fit.extract(fit.data = fit.1mk, filename.save = filename.simdata.save)
  
  # 結果を抽出
  results.1mk.95 <- get_regression_results(
    data.ls = list(
      mk = data$markers,
      data.nma = results.1mk
    ),
    mk.names = data$markers,
    conf = 95
  )
  
  # 実行情報を追加
  results.1mk.95$execution_info <- list(
    covariate = covariate,
    elapsed_time = elapsed_time,
    iter = iter,
    chains = chains,
    cores = cores,
    adapt_delta = adapt_delta,
    stepsize = stepsize,
    max_treedepth = max_treedepth
  )
  
  return(results.1mk.95)
}

#' ネットワークメタアナリシスの結果を抽出する関数
#'
#' @param data.ls データリスト
#' @param mk.names マーカー名のベクトル
#' @param n.study 研究数
#' @param conf 信頼区間（%）
#' @return 結果のリスト
get_nma_results <- function(data.ls = NULL, mk.names = NULL, n.study = 1, conf = 95) {
  data.nma <- data.ls$data.nma
  
  func.ub <- NULL
  func.lb <- NULL
  
  if (conf == 95) {
    func.ub <- get.ub.95
    func.lb <- get.lb.95
  } else if (conf == 90) {
    func.ub <- get.ub.90
    func.lb <- get.lb.90
  }
  
  data.nma.muse <- data.nma[paste('MU[', 1, ',', 1:length(data.ls$mk), ']', sep = '')]
  data.nma.musp <- data.nma[paste('MU[', 2, ',', 1:length(data.ls$mk), ']', sep = '')]
  data.nma.rrse <- data.nma[paste('RR[', 1, ',', 1:length(data.ls$mk), ']', sep = '')]
  data.nma.rrsp <- data.nma[paste('RR[', 2, ',', 1:length(data.ls$mk), ']', sep = '')]
  data.nma.orse <- data.nma[paste('OR[', 1, ',', 1:length(data.ls$mk), ']', sep = '')]
  data.nma.orsp <- data.nma[paste('OR[', 2, ',', 1:length(data.ls$mk), ']', sep = '')]
  data.nma.dor <- data.nma[paste('DOR[', 1:length(data.ls$mk), ']', sep = '')]
  data.nma.sindex <- data.nma[paste('S[', 1:length(data.ls$mk), ']', sep = '')]
  data.nma.sindex_se <- data.nma[paste('S_se[', 1:length(data.ls$mk), ']', sep = '')]
  data.nma.sindex_sp <- data.nma[paste('S_sp[', 1:length(data.ls$mk), ']', sep = '')]
  
  mu.info <- data.frame(
    mean_sp = array(unlist(lapply(data.nma.musp, mean))),
    ub_sp = array(unlist(lapply(data.nma.musp, func.ub))),
    lb_sp = array(unlist(lapply(data.nma.musp, func.lb))),
    mean_se = array(unlist(lapply(data.nma.muse, mean))),
    ub_se = array(unlist(lapply(data.nma.muse, func.ub))),
    lb_se = array(unlist(lapply(data.nma.muse, func.lb)))
  )
  
  rr.info <- data.frame(
    mean_sp = array(unlist(lapply(data.nma.rrsp, mean))),
    ub_sp = array(unlist(lapply(data.nma.rrsp, func.ub))),
    lb_sp = array(unlist(lapply(data.nma.rrsp, func.lb))),
    mean_se = array(unlist(lapply(data.nma.rrse, mean))),
    ub_se = array(unlist(lapply(data.nma.rrse, func.ub))),
    lb_se = array(unlist(lapply(data.nma.rrse, func.lb)))
  )
  
  or.info <- data.frame(
    mean_sp = array(unlist(lapply(data.nma.orsp, mean))),
    ub_sp = array(unlist(lapply(data.nma.orsp, func.ub))),
    lb_sp = array(unlist(lapply(data.nma.orsp, func.lb))),
    mean_se = array(unlist(lapply(data.nma.orse, mean))),
    ub_se = array(unlist(lapply(data.nma.orse, func.ub))),
    lb_se = array(unlist(lapply(data.nma.orse, func.lb)))
  )
  
  dor.info <- data.frame(
    mean = array(unlist(lapply(data.nma.dor, mean))),
    ub = array(unlist(lapply(data.nma.dor, func.ub))),
    lb = array(unlist(lapply(data.nma.dor, func.lb)))
  )
  
  sindex.info <- data.frame(
    mean = array(unlist(lapply(data.nma.sindex, mean))),
    ub = array(unlist(lapply(data.nma.sindex, func.ub))),
    lb = array(unlist(lapply(data.nma.sindex, func.lb)))
  )
  
  sindex_se.info <- data.frame(
    mean = array(unlist(lapply(data.nma.sindex_se, mean))),
    ub = array(unlist(lapply(data.nma.sindex_se, func.ub))),
    lb = array(unlist(lapply(data.nma.sindex_se, func.lb)))
  )
  
  sindex_sp.info <- data.frame(
    mean = array(unlist(lapply(data.nma.sindex_sp, mean))),
    ub = array(unlist(lapply(data.nma.sindex_sp, func.ub))),
    lb = array(unlist(lapply(data.nma.sindex_sp, func.lb)))
  )
  
  # 行名を設定（マーカー名）
  row.names(mu.info) <- mk.names[mk.names %in% data.ls$mk]
  row.names(rr.info) <- mk.names[mk.names %in% data.ls$mk]
  row.names(or.info) <- mk.names[mk.names %in% data.ls$mk]
  row.names(dor.info) <- mk.names[mk.names %in% data.ls$mk]
  row.names(sindex.info) <- mk.names[mk.names %in% data.ls$mk]
  row.names(sindex_se.info) <- mk.names[mk.names %in% data.ls$mk]
  row.names(sindex_sp.info) <- mk.names[mk.names %in% data.ls$mk]
  
  results.nma <- list(
    marker = mk.names[mk.names %in% data.ls$mk],
    mu = mu.info,
    rr = rr.info,
    or = or.info,
    dor = dor.info,
    sindex = sindex.info,
    sindex_se = sindex_se.info,
    sindex_sp = sindex_sp.info
  )
  
  return(results.nma)
}

#' メタ回帰分析の結果を抽出する関数
#'
#' @param data.ls データリスト
#' @param mk.names マーカー名のベクトル
#' @param n.study 研究数
#' @param conf 信頼区間（%）
#' @return 結果のリスト
get_regression_results <- function(data.ls = NULL, mk.names = NULL, n.study = 1, conf = 95) {
  data.nma <- data.ls$data.nma
  
  func.ub <- NULL
  func.lb <- NULL
  
  if (conf == 95) {
    func.ub <- get.ub.95
    func.lb <- get.lb.95
  } else if (conf == 90) {
    func.ub <- get.ub.90
    func.lb <- get.lb.90
  }
  
  data.nma.muse <- data.nma[paste('MU[', 1, ',', 1:length(data.ls$mk), ']', sep = '')]
  data.nma.musp <- data.nma[paste('MU[', 2, ',', 1:length(data.ls$mk), ']', sep = '')]
  data.nma.rrse <- data.nma[paste('RR[', 1, ',', 1:length(data.ls$mk), ']', sep = '')]
  data.nma.rrsp <- data.nma[paste('RR[', 2, ',', 1:length(data.ls$mk), ']', sep = '')]
  data.nma.orse <- data.nma[paste('OR[', 1, ',', 1:length(data.ls$mk), ']', sep = '')]
  data.nma.orsp <- data.nma[paste('OR[', 2, ',', 1:length(data.ls$mk), ']', sep = '')]
  data.nma.dor <- data.nma[paste('DOR[', 1:length(data.ls$mk), ']', sep = '')]
  data.nma.sindex <- data.nma[paste('S[', 1:length(data.ls$mk), ']', sep = '')]
  data.nma.sindex_se <- data.nma[paste('S_se[', 1:length(data.ls$mk), ']', sep = '')]
  data.nma.sindex_sp <- data.nma[paste('S_sp[', 1:length(data.ls$mk), ']', sep = '')]
  data.nma.beta1se <- data.nma[paste('beta1[1,', 1:length(data.ls$mk), ']', sep = '')]
  data.nma.beta1sp <- data.nma[paste('beta1[2,', 1:length(data.ls$mk), ']', sep = '')]
  
  mu.info <- data.frame(
    mean_sp = array(unlist(lapply(data.nma.musp, mean))),
    ub_sp = array(unlist(lapply(data.nma.musp, func.ub))),
    lb_sp = array(unlist(lapply(data.nma.musp, func.lb))),
    mean_se = array(unlist(lapply(data.nma.muse, mean))),
    ub_se = array(unlist(lapply(data.nma.muse, func.ub))),
    lb_se = array(unlist(lapply(data.nma.muse, func.lb)))
  )
  
  rr.info <- data.frame(
    mean_sp = array(unlist(lapply(data.nma.rrsp, mean))),
    ub_sp = array(unlist(lapply(data.nma.rrsp, func.ub))),
    lb_sp = array(unlist(lapply(data.nma.rrsp, func.lb))),
    mean_se = array(unlist(lapply(data.nma.rrse, mean))),
    ub_se = array(unlist(lapply(data.nma.rrse, func.ub))),
    lb_se = array(unlist(lapply(data.nma.rrse, func.lb)))
  )
  
  or.info <- data.frame(
    mean_sp = array(unlist(lapply(data.nma.orsp, mean))),
    ub_sp = array(unlist(lapply(data.nma.orsp, func.ub))),
    lb_sp = array(unlist(lapply(data.nma.orsp, func.lb))),
    mean_se = array(unlist(lapply(data.nma.orse, mean))),
    ub_se = array(unlist(lapply(data.nma.orse, func.ub))),
    lb_se = array(unlist(lapply(data.nma.orse, func.lb)))
  )
  
  dor.info <- data.frame(
    mean = array(unlist(lapply(data.nma.dor, mean))),
    ub = array(unlist(lapply(data.nma.dor, func.ub))),
    lb = array(unlist(lapply(data.nma.dor, func.lb)))
  )
  
  sindex.info <- data.frame(
    mean = array(unlist(lapply(data.nma.sindex, mean))),
    ub = array(unlist(lapply(data.nma.sindex, func.ub))),
    lb = array(unlist(lapply(data.nma.sindex, func.lb)))
  )
  
  sindex_se.info <- data.frame(
    mean = array(unlist(lapply(data.nma.sindex_se, mean))),
    ub = array(unlist(lapply(data.nma.sindex_se, func.ub))),
    lb = array(unlist(lapply(data.nma.sindex_se, func.lb)))
  )
  
  sindex_sp.info <- data.frame(
    mean = array(unlist(lapply(data.nma.sindex_sp, mean))),
    ub = array(unlist(lapply(data.nma.sindex_sp, func.ub))),
    lb = array(unlist(lapply(data.nma.sindex_sp, func.lb)))
  )
  
  beta1.info <- data.frame(
    mean_sp = array(unlist(lapply(data.nma.beta1sp, mean))),
    ub_sp = array(unlist(lapply(data.nma.beta1sp, func.ub))),
    lb_sp = array(unlist(lapply(data.nma.beta1sp, func.lb))),
    mean_se = array(unlist(lapply(data.nma.beta1se, mean))),
    ub_se = array(unlist(lapply(data.nma.beta1se, func.ub))),
    lb_se = array(unlist(lapply(data.nma.beta1se, func.lb)))
  )
  
  # 行名を設定（マーカー名）
  row.names(mu.info) <- mk.names[mk.names %in% data.ls$mk]
  row.names(rr.info) <- mk.names[mk.names %in% data.ls$mk]
  row.names(or.info) <- mk.names[mk.names %in% data.ls$mk]
  row.names(dor.info) <- mk.names[mk.names %in% data.ls$mk]
  row.names(sindex.info) <- mk.names[mk.names %in% data.ls$mk]
  row.names(sindex_se.info) <- mk.names[mk.names %in% data.ls$mk]
  row.names(sindex_sp.info) <- mk.names[mk.names %in% data.ls$mk]
  row.names(beta1.info) <- mk.names[mk.names %in% data.ls$mk]
  
  results.nma <- list(
    marker = mk.names[mk.names %in% data.ls$mk],
    mu = mu.info,
    rr = rr.info,
    or = or.info,
    dor = dor.info,
    sindex = sindex.info,
    sindex_se = sindex_se.info,
    sindex_sp = sindex_sp.info,
    beta1 = beta1.info
  )
  
  return(results.nma)
}
