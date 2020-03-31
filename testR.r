
#======================== library ==============================

library(rgl)
library(dplyr)
library(outliers)
library(JGR)
library(rJava)
library(caret)
library(doParallel)
library(gtools)
library(glmnet)
library(BBmisc) #normalizeを使うため
library(stringr)
library(Metrics)

#======================= read data ============================

data_ready = function()
{
  
  data = read.csv("C:\\ICF_AutoCapsule_disabled\\R\\testdata_ball_GaussNoise_HighLevel_ver2.csv")
  #data = filter(data, Y_L_noNoise<5)
  #data = filter(data, Y_P_NoNoise<70)
  #data = filter(data, Y_P_NoNoise>5)
  
  Xdata = data[,1:7]
  Ydata = data[,18]
  XYdata = cbind(Xdata,Ydata)
  write.csv(XYdata, "XYdata.csv")
  
  #filterしたい時はここをいじる（異常値が有ったりする行を消す）
  XYdata = filter(XYdata, Ydata<0.5)
  XYdata = filter(XYdata, Ydata>0)
  
  #----------------- 各descriptorの名前(Symbol)作成 X1, X2 ,,,--------
  
  number_X = ncol(XYdata)
  Xsymbol = c()
  
  for (ii in 1:number_X)
  {
    Xsymbol_each = paste("X", ii, sep="")
    Xsymbol = c(Xsymbol, Xsymbol_each)
  }
  
  Xsymbol_and_DescriptorName = cbind(colnames(XYdata), Xsymbol)
  
  Xsymbol_Xdata_original = XYdata
  colnames(Xsymbol_Xdata_original) = Xsymbol
  return(Xsymbol_Xdata_original)
  
}


#========================= ここからしばらくはfuntion作成 =================

#===================== 積の項を作る ===============================

make_product = function(Xsymbol_Xdata)
{
  
  X_symbol_name = colnames(Xsymbol_Xdata)
  
  product_term_name = combinations(n=ncol(Xsymbol_Xdata), r=2, X_symbol_name)

  
  
  product_term = matrix(0, ncol=1, nrow=nrow(Xsymbol_Xdata))
  
  for (nn in 1:nrow(product_term_name))
  {
    product = Xsymbol_Xdata[ ,product_term_name[nn,1]]*Xsymbol_Xdata[ ,product_term_name[nn,2]]
    
    product_term = cbind(product_term, product)
  
  }
  
  product_term = product_term[ ,-1]
  colnames(product_term)=paste("(", product_term_name[ ,1], product_term_name[, 2],")", sep ="")
  
  return(product_term)
  
}

#======================== 商の項を作る ================================

make_quotient = function(Xsymbol_Xdata)
{
  
  X_symbol_name = colnames(Xsymbol_Xdata)
  
  quotient_term_name = permutations(n=ncol(Xsymbol_Xdata), r=2, X_symbol_name, repeats.allowed = F)
  quotient_term = matrix(0, ncol=1, nrow=nrow(Xsymbol_Xdata))
  
  for (nn in 1:nrow(quotient_term_name))
  {
    
    quotient = Xsymbol_Xdata[ ,quotient_term_name[nn,1]] / Xsymbol_Xdata[ ,quotient_term_name[nn,2]]
    quotient_term = cbind(quotient_term, quotient)
    
  }
  
  quotient_term = quotient_term[ ,-1]
  colnames(quotient_term)=paste("(", quotient_term_name[ ,1],"/", quotient_term_name[, 2], ")" , sep ="")
  return(quotient_term)
  
}
#======================== 二乗の項を作る =================================

make_square = function(Xsymbol_Xdata)
{
  
  X_symbol_name = colnames(Xsymbol_Xdata)
  
  square_term_name = cbind(X_symbol_name, X_symbol_name)
  
  square_term = matrix(0, ncol=1, nrow=nrow(Xsymbol_Xdata))
  
  for (nn in 1:nrow(square_term_name))
  {
    
    square = Xsymbol_Xdata[ ,square_term_name[nn,1]]*Xsymbol_Xdata[ ,square_term_name[nn,2]]
    square_term = cbind(square_term, square)
    
  }
  
  square_term = square_term[ ,-1]
  colnames(square_term)=paste("(", square_term_name[ ,1], square_term_name[, 2],")", sep ="")
  
  return(square_term)
  
}

#=========================== 指数の項を作る ================================
#あとで
#===========================　和の項を作る =================================

make_sum = function(Xsymbol_Xdata)
{
  
  X_symbol_name = colnames(Xsymbol_Xdata)
  
  sum_term_name = combinations(n=ncol(Xsymbol_Xdata), r=2, X_symbol_name)
  
  sum_term = matrix(0, ncol=1, nrow=nrow(Xsymbol_Xdata))
  
  for (nn in 1:nrow(sum_term_name))
  {
    
    sum1 = Xsymbol_Xdata[ ,sum_term_name[nn,1]]+Xsymbol_Xdata[ ,sum_term_name[nn,2]]
    sum_term = cbind(sum_term, sum1)
    
  }
  
  sum_term = sum_term[ ,-1]
  colnames(sum_term)=paste("(", sum_term_name[ ,1], "+" ,sum_term_name[, 2],")", sep ="")
  
  return(sum_term)
  
}




#=========================== 差の項を作る =================================

make_difference = function(Xsymbol_Xdata)
{
  
  X_symbol_name = colnames(Xsymbol_Xdata)
  
  difference_term_name = permutations(n=ncol(Xsymbol_Xdata), r=2, X_symbol_name, repeats.allowed = F)
  
  difference_term = matrix(0, ncol=1, nrow=nrow(Xsymbol_Xdata))
  
  for (nn in 1:nrow(difference_term_name))
  {
    
    difference = Xsymbol_Xdata[ ,difference_term_name[nn,1]] - Xsymbol_Xdata[ ,difference_term_name[nn,2]]
    difference_term = cbind(difference_term, difference)
    
  }
  
  difference_term = difference_term[ ,-1]
  colnames(difference_term)=paste("(", difference_term_name[ ,1],"-", difference_term_name[, 2], ")" , sep ="")
  
  return(difference_term)
  
}

#============================ penalty factorを作る =======================

make_penalty_factor = function(new_Xdata, Xsymbol_Xdata_original)
{
  X_symbol_name = colnames(new_Xdata)
  X_symbol_name_removed = str_extract_all(X_symbol_name, pattern="X", simplify=FALSE)
  penalty_factor = c()
  print(length(X_symbol_name_removed))
  for (ii in 1:length(X_symbol_name_removed))
  {
    penalty_factor = c(penalty_factor, length(X_symbol_name_removed[[ii]]))
  }

  #ある一定以上の次数のpenaltyを極端に大きくする。とりあえず、Xsymbol_Xdata_originalの項数+2を閾値としてデフォルトにしている。
  #the penalty factors areinternally rescaled to sum to nvars
  #つまり各penalty.factorは内部で規格化される、規格化されたpenalty.factor = （各penalty.factorの値 * 項の数 / 全部のpenalty.factorの値の和）
  #penalty_factor = replace(penalty_factor, (penalty_factor > ncol(Xsymbol_Xdata_original)+1),100)
  #penalty_factor = replace(penalty_factor, (penalty_factor <= ncol(Xsymbol_Xdata_original)+1),1)
  #もしくは、自分で指定する（今は4）
  penalty_factor = replace(penalty_factor, (penalty_factor > 4),100)
  penalty_factor = replace(penalty_factor, (penalty_factor <= 4),1)

  return(penalty_factor)
  
}

#========================== 標準化した回帰係数を計算するfunction =================

make_std_coefficient = function(model, Xdata)
{
  cf = coef(model, s="lambda.min")[ ,1]
  sds = sapply(as.data.frame(Xdata), sd)
  mus = sapply(as.data.frame(Xdata), mean)
  cf_std <- cf * c(1, c(sds[names(cf)][-1]))
  #cf_std[1] <- cf[1] + sum(cf[-1] * mus[names(cf)][-1]) #interceptアリの時はこれも入れる
  
  return(cf_std)
}

#========================== 標準化した回帰係数を計算するfunction =================

make_std_coefficient_intercept = function(model, Xdata)
{
  cf = coef(model, s="lambda.min")[ ,1]
  sds = sapply(as.data.frame(Xdata), sd)
  mus = sapply(as.data.frame(Xdata), mean)
  cf_std <- cf * c(1, c(sds[names(cf)][-1]))
  cf_std[1] <- cf[1] + sum(cf[-1] * mus[names(cf)][-1]) #interceptアリの時はこれも入れる
  
  return(cf_std)
}

X_Y = data_ready()

#データ数をちょっと減らす　2000くらいにする
X_Y = X_Y[1:5000, ]

#Xsymbol_Xdata_originalがある状態から始める
#Xsymbol_Xdata = Xsymbol_Xdata_original

nnn = 1 #ランダムにに実行する回数（nfoldを少しずつ変えながら）
#pmax_vector = c(100, 80, 60, 40, 20, 10, 10) #pmaxの値は少しずつ減らしていく。ここで指定
#pmax_vector = c(300, 300, 300, 300) #pmaxの値は少しずつ減らしていく。ここで指定。reccurentする回数分の要素がある

rrr = ncol(X_Y)  #reccurentする回数。ここでは、最初の説明変数の数に設定してみる
#rrr = 3

FirstLassoCutoff = 30 #最初の変数選択での変数の数の上限
FinalCutoff = 5　#最後に出てくる式の変数の上限

max_pmax = 1000

object_record = c()
error_record = c()
#term_record = matrix(0, nrow=1, ncol=max_pmax)
#coefficient_record = matrix(0, nrow=1, ncol=max_pmax)
summary_data = matrix(0, nrow=1, ncol=max_pmax)

i_error = 10000  #誤差の初期値　大きめに取っておけば問題ない

#X_YのデータをXとYに分ける。それを網羅的にやる
#for (pp in 1:ncol(X_Y))　#全部試すときはこっち
for (pp in ncol(X_Y):ncol(X_Y)) #X5をYにする時だけ
{
  Ydata = X_Y[ ,pp,drop=FALSE]
  Ydata = normalize(Ydata, method="range", range=c(0,1))
  Xsymbol_Xdata = X_Y[, -pp , drop=FALSE]
  Xsymbol_Xdata_original = Xsymbol_Xdata
  
  for (nn in 1:nnn)
  {
    
    for (ii in 1:rrr) #lambdaの下限を指定してだんだん選ばれる項の数を減らす 
    {
      
      #積と商と和と差と二乗と指数の項を作成する
      product_term = make_product(Xsymbol_Xdata)
      quotient_term = make_quotient(Xsymbol_Xdata)
      square_term = make_square(Xsymbol_Xdata)
      sum_term = make_sum(Xsymbol_Xdata)
      difference_term = make_difference(Xsymbol_Xdata)
      new_Xdata = cbind(Xsymbol_Xdata, product_term, quotient_term, square_term, sum_term, difference_term)
      
      #NAやNaNやInfがたまに出てくるので、それが含まれる列を除去
      new_Xdata[!grepl("Inf", new_Xdata)]
      new_Xdata[!grepl("NA", new_Xdata)]
      new_Xdata[!grepl("NaN", new_Xdata)]
      
      #全く同じ項が選ばれることがあるので、その場合は片方を除去
      new_Xdata = t(new_Xdata)
      new_Xdata = unique(new_Xdata)
      new_Xdata = t(new_Xdata)
      
      #履歴をプリント
      moji = paste("pp", pp, "nn", nn, "ii", ii, sep=" ")
      print(moji) 
      moji = paste("number of Xsymbol_Xdata", ncol(Xsymbol_Xdata), sep=" ")
      print(moji)
      moji = paste("number of new_Xdata", ncol(new_Xdata), sep=" ")
      print(moji)

      #penalty.factorを作る。次数（項の名前に入っているXの数）が大きい項は選ばれにくくなるようにする。
      penalty_factor = make_penalty_factor(new_Xdata, Xsymbol_Xdata_original)
      
      cl = makeCluster(4)  #10コアで並列化開始
      registerDoParallel(cl)  
      #cross validation
      fitLassoCV1 = cv.glmnet(x = as.matrix(new_Xdata), y = as.matrix(Ydata), family ="gaussian", alpha = 1
                              ,nfolds = 5, parallel=TRUE, standardize = TRUE, penalty.factor=penalty_factor
                              ,thresh = 1E-7, maxit = 10^5, nlambda = 1000, intercept=FALSE, lambda.min.ratio = 0.00000001
                              ,lambda = 2^(-40:5) )
      stopCluster(cl)  #並列化おしまい
      print(fitLassoCV1)
      coefficient = coef(fitLassoCV1, s="lambda.min") #たくさん選ばれるように1seは使わない
      coefficient = as.matrix(coefficient)
      coefficient = coefficient[coefficient !=0,]
      #print(coefficient)
      used_term_name = names(coefficient)
      #print(used_term_name)
    }
  }
}

