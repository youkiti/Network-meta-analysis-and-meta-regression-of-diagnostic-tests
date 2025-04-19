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
