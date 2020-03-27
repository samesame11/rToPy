
#======================== library ==============================

library(rgl)
library(dplyr)
library(outliers)
library(JGR)
library(rJava)
library(caret)
library(doParallel)
library(gtools)
library(dplyr)
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

#========================================================================
#====== non-linear reccurent adaptive interpretable regression (NRAI regression) ======

#------------------ input parameter -------------------------------------

#読み込むファイル名とかは上の方をいじってね.目的変数は決めておらず、フレキシブルに網羅的に試す
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
			coefficient = coef(fitLassoCV1, s="lambda.min") #たくさん選ばれるように1seは使わない
			coefficient = as.matrix(coefficient)
			coefficient = coefficient[coefficient !=0,]
			used_term_name = names(coefficient)
			#used_term_name = used_term_name[-which(used_term_name %in% "(Intercept)")] #intercept入れてるときはこれも入れる

			used_term_data = as.matrix(new_Xdata[ ,used_term_name])
			colnames(used_term_data) = used_term_name ##

			#標準化係数の計算
			std_coef = make_std_coefficient(fitLassoCV1, new_Xdata)
			#標準化係数の上位いくつかのみを選ぶ
			abs_std_coef = abs(std_coef)
			abs_std_coef = abs_std_coef[abs_std_coef !=0]
			sort_abs_std_coef = sort(abs_std_coef, decreasing=T)
			sort_coef_name = names(sort_abs_std_coef)
			top_coef_name = head(sort_coef_name, FirstLassoCutoff)
			top_coef = std_coef[top_coef_name]
			std_selected_term_data =  as.matrix(new_Xdata[ ,top_coef_name])

			#標準化係数でピックアップしたデータのみを使用して再度LASSO
			#cl = makeCluster(4)  #10コアで並列化開始
			#registerDoParallel(cl)  
			#cross validation
			#fitLassoCV1 = cv.glmnet(x = as.matrix(std_selected_term_data), y = as.matrix(Ydata), family ="gaussian", alpha = 1
    	   	      #    	    ,nfolds = 10, parallel=TRUE, standardize = TRUE, penalty.factor=penalty_factor
    	           	#          ,thresh = 1E-7, maxit = 10^5, nlambda = 1000, intercept=FALSE, lambda.min.ratio = 0.00000001
			#          ,lambda = 2^(-40:5) )
			#stopCluster(cl)  #並列化おしまい

			#coefficient = coef(fitLassoCV1, s="lambda.min") #たくさん選ばれるように1seは使わない
			#coefficient = as.matrix(coefficient)
			#coefficient = coefficient[coefficient !=0,]
			#used_term_name = names(coefficient)
			#used_term_name = used_term_name[-which(used_term_name %in% "(Intercept)")] #intercept入れてるときはこれも入れる
			#used_term_data = as.matrix(std_selected_term_data[ ,used_term_name])
			#used_term_data = as.matrix(new_Xdata[ ,used_term_name])
			#colnames(used_term_data) = used_term_name ##

			#履歴をプリント
			moji = paste("number of selected data with 1st LASSO", length(coefficient), sep=" ")
			print(moji)
			moji = paste("number of selected data with 2nd std shrink", ncol(std_selected_term_data), sep=" ")
			print(moji)
			moji = paste("error", min(fitLassoCV1$cvm), sep=" ")
			print(moji)

			Xsymbol_Xdata = cbind(Xsymbol_Xdata_original, std_selected_term_data)
	
			#各iterationで選ばれた項の記録を残す
			used_term_name_record = used_term_name
			if (is.null(used_term_name_record)==TRUE){
				used_term_name_record = c(NA)}
			length(used_term_name_record) = max_pmax
			#term_record = rbind(term_record, used_term_name_record)

			#さらにstd_coefのsortで選ばれた項の記録を残す
			selected_coef_name_record = top_coef_name
			if (is.null(selected_coef_name_record)==TRUE){
				selected_coef_name_record = c(NA)}
			length(selected_coef_name_record) = max_pmax
			#term_record = rbind(term_record, used_term_name_record)

			#各coefficientの履歴も残す
			coefficient_one_record = coefficient
			coefficient_one_record = as.vector(coefficient_one_record)
			if (is.null(coefficient_one_record)==TRUE){
				coefficient_one_record = c(NA)}
			#coefficient_one_record = c("-", coefficient_one_record) #intercept入れてるときはこっち
			coefficient_one_record = c("-","-", coefficient_one_record)　#intercept入ってないときはこっち
			length(coefficient_one_record) = max_pmax			

			#各std_coefの履歴も残す
			std_coefficient_one_record = top_coef
			std_coefficient_one_record = as.vector(std_coefficient_one_record)
			if (is.null(std_coefficient_one_record)==TRUE){
				std_coefficient_one_record = c(NA)}
			#coefficient_one_record = c("-", coefficient_one_record) #intercept入れてるときはこっち
			std_coefficient_one_record = c("-","-", std_coefficient_one_record)　#intercept入ってないときはこっち
			length(std_coefficient_one_record) = max_pmax			


			#errorの履歴を残す
			error_one = min(fitLassoCV1$cvm)

			#目的変数の履歴も残す
			object_one = colnames(Ydata)

			#used_term_name, coettifient, error(fitLassoCV1$cvm), object_record(colnames(Ydata))を合併した行列を作成
			summary_data_EOU = c(error_one, object_one, used_term_name_record)
			length(summary_data_EOU) = max_pmax
			summary_data_EOstdU = c("-", "-", top_coef_name)
			length(summary_data_EOstdU) = max_pmax


			summary_data = rbind(summary_data, summary_data_EOU)
			summary_data = rbind(summary_data, coefficient_one_record)
			summary_data = rbind(summary_data, summary_data_EOstdU)
			summary_data = rbind(summary_data, std_coefficient_one_record)

			#間をあけてわかりやすくする
			aida = rep("-", max_pmax)
			summary_data = rbind(summary_data, aida)


			#各iterationでのエラー（cvm）の記録を残す
			#error_record = c(error_record, min(fitLassoCV1$cvm))

			#各iterationでの目的変数の名前の履歴を残す
			#object_record = c(object_record, colnames(Ydata))



			#新しいXsymbol_XdataにXsymbol_Xdata_originalのモノが含まれているとあとでエラーになるので、2重のモノは消す
			X_symbol_name_original = colnames(Xsymbol_Xdata_original)
			Xsymbol_Xdata_add = Xsymbol_Xdata[, -which(colnames(Xsymbol_Xdata) %in% X_symbol_name_original), drop=FALSE]
			#項が1個しか残らなかった場合は、なぜかXsymbol_Xdata_addの列名が消えてしまうのでdrop=FALSEを付けた
			Xsymbol_Xdata = cbind(Xsymbol_Xdata_original, Xsymbol_Xdata_add)

			#一番cv誤差が小さかったモデルを記録
			if (min(fitLassoCV1$cvm) < i_error){
				i_error = min(fitLassoCV1$cvm)
				best_model = fitLassoCV1
				best_model_x_used = new_Xdata
				best_model_y_used = Ydata
				best_model_x_selected = std_selected_term_data
				}
		}
	}
}


#一番良かったモデルの標準化係数を計算
std_coef_final = make_std_coefficient(best_model, best_model_x_selected)

#標準化係数の上位いくつかのみを選ぶ
abs_std_coef = abs(std_coef_final)
abs_std_coef = abs_std_coef[abs_std_coef !=0]
sort_abs_std_coef = sort(abs_std_coef, decreasing=T)
sort_coef_name = names(sort_abs_std_coef)
top_coef_name = head(sort_coef_name, FinalCutoff)　
top_coef = std_coef_final[top_coef_name]
std_selected_term_data =  as.matrix(best_model_x_used[ ,top_coef_name])

#一番良かったモデルの標準化係数上位数個のもので線形回帰
#同じパッケージでやりたかったので、λを極力小さいものにした
#最後はinterceptをTRUEにして、切片有にした
cl = makeCluster(4)
registerDoParallel(cl)  
#cross validation
FinalLinearModel = cv.glmnet(x = as.matrix(std_selected_term_data), y = as.matrix(best_model_y_used), family ="gaussian", alpha = 1
   	          	    ,nfolds = 5, parallel=TRUE, standardize = TRUE, ,thresh = 1E-7
			    ,maxit = 10^5, nlambda = 100, intercept=TRUE, lambda.min.ratio = 0.00000000001
		          ,lambda = c(0.000000000000000001,5000))
stopCluster(cl)  #並列化おしまい

#最終モデルの係数
coefficient_final = coef(FinalLinearModel, s="lambda.min")
coefficient_final = as.matrix(coefficient_final)
coefficient_final = coefficient_final[coefficient_final !=0,]
used_term_name_final = names(coefficient_final)
final_error = min(FinalLinearModel$cvm)

#最終モデルの標準化した係数
std_coef_finalfinal = make_std_coefficient_intercept(FinalLinearModel, std_selected_term_data)
std_coef_finalfinal_name = names(std_coef_finalfinal)



#summary.csvファイルの作成
summary_data = summary_data[-1,]
final_coefficient_name = c(final_error, used_term_name_final)
length(final_coefficient_name) = max_pmax
final_coefficient = c("-", coefficient_final)
length(final_coefficient) = max_pmax
std_coef_finalfinal = c("-", std_coef_finalfinal)
length(std_coef_finalfinal)=max_pmax
std_coef_finalfinal_name = c("-", std_coef_finalfinal_name)
length(std_coef_finalfinal_name) = max_pmax

summary_data = rbind(summary_data, final_coefficient_name)
summary_data = rbind(summary_data, final_coefficient)
summary_data = rbind(summary_data, std_coef_finalfinal_name)
summary_data = rbind(summary_data, std_coef_finalfinal)

colnames(summary_data) = c("error", "object/intercept", rep("X", max_pmax-2))
write.csv(summary_data, "summary.csv")


#プロット
predictY = predict(FinalLinearModel, newx=as.matrix(std_selected_term_data))
plot(predictY, as.matrix(best_model_y_used), xlim=c(0,1), ylim=c(0,1))































#------------------- 比較のためのただのLASSO ----------------------------

cl = makeCluster(10)  #10コアで並列化開始
registerDoParallel(cl)  
#cross validation
fitLassoCV1 = cv.glmnet(x = as.matrix(Xsymbol_Xdata_original), y = as.matrix(Ydata), family ="gaussian", alpha = 1
          	          ,nfolds = 50, pmax = 5, standardize = TRUE)
stopCluster(cl)  #並列化おしまい
predictY = predict(fitLassoCV1, as.matrix(Xsymbol_Xdata_original))
plot(as.matrix(predictY), as.matrix(Ydata), xlim=c(0,1), ylim=c(0,1))
coef(fitLassoCV1, s="lambda.min")
min(fitLassoCV1$cvm)

#------------------- 比較のためのただのNN ----------------------------

XXdata = Xsymbol_Xdata_original[,1:7]
YYdata = Ydata[,1]

cl = makeCluster(10)  #10コアで並列化開始
registerDoParallel(cl) 
 
ffitcontrol = trainControl(method = "repeatedcv", number = 4, repeats = 1)
ttgrid = expand.grid(layer1 = (5:10), layer2 = (5:10), layer3 = (0:10))
model.tune.NN = train(XXdata, YYdata, method = "mlpML", trControl = ffitcontrol,
                   tuneGrid = ttgrid, metric = "Rsquare", preProc =c("center","scale"),maxit = 1000)
model.tune.NN

stopCluster(cl)  #並列化おしまい

plot(YYdata, predict(model.tune.NN, XXdata), xlim=c(0,1),ylim=c(0,1))
#Mean Square Errorに統一する
mse(predict(model.tune.NN, XXdata),YYdata )




















#同じ列を探してそれを消す
testdata = new_Xdata

testdata = t(testdata)
testdata = unique(testdata)
testdata = t(testdata)
colnames(testdata)





Xsymbol_Xdata = X_Y[, -3 , drop=FALSE]












cl = makeCluster(4)  #10コアで並列化開始
registerDoParallel(cl)  

fitLassoCV1 = cv.glmnet(x = as.matrix(new_Xdata), y = as.matrix(Ydata), family ="gaussian", alpha = 1
    	   	          ,nfolds = 10, parallel=TRUE, standardize = TRUE,penalty.factor=penalty_factor
    	           	    ,thresh = 1E-10, maxit = 10^5, nlambda = 100, intercept=TRUE, lambda.min.ratio = 0.00000000001
			    ,lambda = 2^(-80:-5))

stopCluster(cl) 

plot(fitLassoCV1)












data = read.csv("testdata_VanDerWaals.csv")
#data = filter(data, Y_L_noNoise<5)
#data = filter(data, Y_P_NoNoise<70)
#data = filter(data, Y_P_NoNoise>5)

Xdata = data[,1:5]
Ydata = data[,6]
XYdata = cbind(Xdata,Ydata)

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










test = used_term_name_record

if (is.null(test)==TRUE){
	test = c(NA)}











X_symbol_name = colnames(new_Xdata)
X_symbol_name_removed = str_extract_all(X_symbol_name, pattern="X", simplify=FALSE)

penalty_factor = c()
for (ii in 1:length(X_symbol_name_removed))
{
penalty_factor = c(penalty_factor, length(X_symbol_name_removed[[ii]]))
}




charmatch(X_symbol_name, "X")























predictYY = (Xdata[ ,3]*0.08314*Xdata[ ,2])/(data[ ,1] - data[ ,5]*data[ ,2]) - (data[,4]*data[,2]*data[,2])/(data[,1]*data[,1])

plot(predictYY, Ydata)













#積と商と和と差と二乗と指数の項を作成する
product_term = make_product(Xsymbol_Xdata)
quotient_term = make_quotient(Xsymbol_Xdata)
square_term = make_square(Xsymbol_Xdata)
sum_term = make_sum(Xsymbol_Xdata)
difference_term = make_difference(Xsymbol_Xdata)

new_Xdata = cbind(Xsymbol_Xdata_original, product_term, quotient_term, square_term, sum_term, difference_term)

cl = makeCluster(6)  #4コアで並列化開始
registerDoParallel(cl)  
#cross validation
fitLassoCV1 = cv.glmnet(x = as.matrix(new_Xdata), y = as.matrix(Ydata), family ="gaussian", alpha = 1
                         ,nfolds = 100, pmax = i_pmax)
stopCluster(cl)  #並列化おしまい

coefficient = coef(fitLassoCV1, s="lambda.min")
coefficient = as.matrix(coefficient)
coefficient = coefficient[coefficient !=0,]
used_term_name = names(coefficient)
used_term_name = used_term_name[-which(used_term_name %in% "(Intercept)")]
#used_term_name = used_term_name[-1]

used_term_data = new_Xdata[ ,used_term_name]
Xsymbol_Xdata = cbind(Xsymbol_Xdata_original, used_term_data)

#各iterationで選ばれた項の記録を残す
used_term_name_record = used_term_name
length(used_term_name_record) = i_pmax
term_record = rbind(term_record, used_term_name_record) 



















#Xsymbol_Xdataがある状態から始める
X_symbol_name = colnames(Xsymbol_Xdata)

#積と商と和と差と二乗と指数の項を作成する
product_term = make_product(Xsymbol_Xdata)
quotient_term = make_quotient(Xsymbol_Xdata)
square_term = make_square(Xsymbol_Xdata)
sum_term = make_sum(Xsymbol_Xdata)
difference_term = make_difference(Xsymbol_Xdata)

new_Xdata = cbind(Xsymbol_Xdata, product_term, quotient_term, square_term, sum_term, difference_term)

cl = makeCluster(4)  #4コアで並列化開始
registerDoParallel(cl)  
#cross validation
fitLassoCV1 = cv.glmnet(x = as.matrix(new_Xdata), y = as.matrix(Ydata), family ="gaussian", alpha = 1
                         ,nfolds = 100, pmax = 15)
stopCluster(cl)  #並列化おしまい

coefficient = coef(fitLassoCV1, s="lambda.min")
coefficient = as.matrix(coefficient)
coefficient = coefficient[coefficient !=0,]
used_term_name = names(coefficient)
used_term_name = used_term_name[-1]

used_term_data = new_Xdata[ ,used_term_name]















pred_fitLassoCV1 <- predict(fitLassoCV1, newx=as.matrix(new_Xdata))

plot(pred_fitLassoCV1, Ydata, xlim=c(0,5), ylim=c(0,5))






































dday = "20190528" #virtual screeningを実行する日

#========================= data read ===========================

filename = paste("data_for_study/",dday, "summary.csv", sep="")
data = read.csv(filename)

#========================= データの前処理 ===========================

#------------------------- 組成情報が無いものを除去 --------------------

data = filter(data, composition1>-100 | composition2>-100 | composition3>-100)

#---------- 20181205Fe100_Si320 295 1.0Paは組成情報が変なので除去 ----

data = filter(data, filename != "20181205 Fe110_Si320 295sec 1.0Pa")

#-------- pad_idが207のものを取り除く。昔は最後のサンプルを２重で測定してたから -----------------

data = filter(data, pad_id != 207)

#---------- thermopower_uV（熱電電圧）が閾値以上であるものを除去 ----------------------
#------(今のところMAXがFeAlTeの1.2*10^-5uVくらい。ただしannealサンプルは除く) --------------

data = filter(data, thermopower_uV　< 1.5*10^-5)

#---------- thermopower_uV（熱電電圧）が閾値以下であるものを除去 ----------------------

data = filter(data, thermopower_uV　> -10^-6)

#---------- resist_ohm(電気抵抗)が閾値以上であるものを除去 ----------------------

data = filter(data, resist_ohm　< 5000)

#---------- rhickness_nm（膜厚）がマイナスになっているものは除去 --------------------------

data = filter(data, thickness_nm　> 0)

#=============== magpie用のインプットデータ作成 ====================================

magpie_input = paste(data$element1, data$composition1, data$element2, data$composition2, data$element3, data$composition3, sep=",")
magpie_input = gsub("NA", "", magpie_input)
magpie_input = gsub(",,,", ",", magpie_input)
magpie_input = paste(magpie_input, ", 0 0", sep="")
#magpie_input = c("Name Thermopower Stability //head", magpie_input)

magpie_input_filename = paste("magpie_input_", dday, ".txt", sep="")
write.table(magpie_input, magpie_input_filename, quote=FALSE, row.names=FALSE, col.names=TRUE)

#各行の後ろに"0 0"が付いているのはmagpieの使用によるものなので仕方ない
#できた.txtファイルをmagpieに入れるとmagpie descriptrの完成。主導でやってもいいが、以下を実行しても作れる

#----------------------------- magpie実行 ----------------------------------

order1 = paste("copy /Y ", magpie_input_filename, " C:\\", "ICF_Autocapsule_Disabled\\magpie\\datasets\\", "magpie_input.txt", sep="")
shell(order1)

system("java -jar C:\\ICF_Autocapsule_Disabled\\magpie\\dist\\Magpie.jar C:\\ICF_Autocapsule_Disabled\\magpie\\magpie_input.in")

order2 = paste("copy /Y ", "C:\\ICF_Autocapsule_Disabled\\magpie\\magpie_input.csv", " C:\\ICF_Autocapsule_Disabled\\R\\virtual_screening\\magpie_descriptor\\magpie_descriptor_", dday, ".csv", sep="")
shell(order2)

#これをやると、magpie_descriptorフォルダの下に、magpieのデータファイル（csv）ができる

#==================== データの読み込みとdescriptorの前処理 ==========================

#--------------------------- 学習データの読み込みと作成 ----------------------------

filename_study_x = paste("magpie_descriptor\\magpie_descriptor_", dday, ".csv", sep="")
study_x = read.csv(filename_study_x)
study_y = data$thermopower_uV

#--------------------------- 予測データの読み込みと作成 ----------------------------
#ここはとても重いので、一度やったらやめたほうが良い

Fe_magpie = read.csv("input_for_prediction\\Fe_totaldata.csv")
Co_magpie = read.csv("input_for_prediction\\Co_totaldata.csv")
Ni_magpie = read.csv("input_for_prediction\\Ni_totaldata.csv")

Fe_comp =  read.csv("input_for_prediction\\head_Fe_totaldata.txt")
Co_comp =  read.csv("input_for_prediction\\head_Co_totaldata.txt")
Ni_comp =  read.csv("input_for_prediction\\head_Ni_totaldata.txt")

predict_x_magpie = rbind(Fe_magpie, Co_magpie, Ni_magpie)
predict_x_comp = rbind(Fe_comp, Co_comp, Ni_comp)

#------------------ descriptorの前処理その０　膜厚データを追加する -------------------------

study_x = cbind(study_x, data$thickness_nm)
names(study_x)[ which( names(study_x)=="data$thickness_nm" ) ] <- "thickness_predict"

thickness_predict = rep(100:100, nrow(predict_x_magpie))
predict_x_magpie = cbind(predict_x_magpie, thickness_predict)

#------------------ descriptorの前処理その1　descriptorから文字列（Class）の列を除去 ----------------

study_x = study_x[ , colnames(study_x) != "Class"]
predict_x_magpie = predict_x_magpie[ , colnames(predict_x_magpie) != "Class"]

#---------- descriptorの前処理その2　全部同じ数のdescriptor(相関行列でNAと表示)は除去 -----------------

pearson_matrix = cor(study_x, method = "pearson")
exclude = pearson_matrix[!complete.cases(pearson_matrix[ ,2]),]
exclude_name = rownames(exclude)

study_x = study_x[ , !colnames(study_x) %in% exclude_name]
predict_x_magpie = predict_x_magpie[ , !colnames(predict_x_magpie) %in% exclude_name]

#----------------- descriptorの前処理その3　相関係数が1もしくは-1のものは(最初の方を)除去 ------------------------

study_x2 = study_x
exclude_name2 = c()
exclude_name_summary = c()
pearson_matrix2 = cor(study_x2, method = "pearson")
colname = colnames(pearson_matrix2)
exclude_name2 = c()

#まずは相関係数が1のやつについてやる


repeat{

	for (ii in colname)
	{
		if (table(abs(pearson_matrix2[,ii]) ==1)[[2]] > 1) {
			exclude_name2 = c(exclude_name2, ii) }
	}

	exclude_name_summary = c(exclude_name_summary, exclude_name2[1])
	study_x2 = study_x2[ , colnames(study_x2) != exclude_name2[1]]
	pearson_matrix2 = cor(study_x2, method = "pearson")
	colname = colnames(pearson_matrix2)

	if(length(exclude_name2)==0) { break }
	exclude_name2 = c()

}

study_x = study_x[ , !colnames(study_x) %in% exclude_name_summary]
predict_x_magpie = predict_x_magpie[ , !colnames(predict_x_magpie) %in% exclude_name_summary]

#---------------------- descriptorの前処理その4 正規化 ----------------------------------------------- 

#descriptorの標準化
XX = rbind(study_x, predict_x_magpie)
std_XX = scale(XX)

std_study_x = std_XX[1:nrow(study_x),]
std_predict_x_magpie = std_XX[(nrow(study_x)+1):nrow(XX),]

#出力（thermopwoer_uV）の0-1正規化（正規化の方がデータを戻しやすそうだから）
n_study_y = (study_y - min(study_y))/(max(study_y)-min(study_y))

#------------------------ データの保存 ---------------------------------------------------------------

#write.csv(std_study_x, "study_x_20190528.csv", row.names = FALSE)
#write.csv(n_study_y, "study_y_20190528.csv", row.names = FALSE)

#=========================== モデル作成と予測 ===========================================================

#----------------------------- LASSOモデル作成と予測 ------------------------------------------------------

cl = makeCluster(10)  #10コアで並列化開始
registerDoParallel(cl)  

#caretでcross validation
fitcontrol = trainControl(method = "repeatedcv", number = 20, repeats = 1)
tgrid = expand.grid(lambda = 2^(-25:-15))
model_tune_LASSO = train(std_study_x, n_study_y, method = "rqlasso", trControl = fitcontrol,
                   tuneGrid = tgrid)
#lambda=4.76*10^-7くらいがよいらしい


stopCluster(cl)  #並列化おしまい

model_tune_LASSO
plot(n_study_y, predict(model_tune_LASSO, std_study_x), xlim=c(-0.1,1),ylim=c(-0.1,1))


#　LASSOモデルで予測

LASSO_predict = predict(model_tune_LASSO, std_predict_x_magpie)

#　予測値でソート 

sort_matrix = cbind(predict_x_comp, LASSO_predict)
sort_matrix = sort_matrix[order(sort_matrix$LASSO_predict, decreasing=T), ]
filemameLASSO = paste("predict_LASSO_", dday, ".csv", sep="")
write.csv(sort_matrix, filemameLASSO)

#これで予測値のアウトプットファイルができる

#----------------------------- NNモデル作成と予測 ------------------------------------------------------

cl = makeCluster(10)  #10コアで並列化開始
registerDoParallel(cl)  

#caretでcross validation

ffitcontrol = trainControl(method = "repeatedcv", number = 20, repeats = 1)
ttgrid = expand.grid(decay = 2^(-18:-1), size = (2:10)*2)
model.tune.NN = train(std_study_x, n_study_y, method = "nnet", trControl = ffitcontrol,
                   tuneGrid = ttgrid, metric = "Rsquared")
model.tune.NN

#decay=7.629*10^-6, size=6辺りが良いらしい

stopCluster(cl)  #並列化おしまい

plot(n_study_y, predict(model.tune.NN, std_study_x), xlim=c(-0.1,1),ylim=c(-0.1,1))

#　NNモデルで予測

NN_predict = predict(model.tune.NN, std_predict_x_magpie)

sort_matrix = cbind(predict_x_comp, NN_predict)
sort_matrix = sort_matrix[order(sort_matrix$NN_predict, decreasing=T), ]
filemameLASSO = paste("NN_predict_", dday, ".csv", sep="")
write.csv(sort_matrix, filemameLASSO)














ave(model_tune_SVM, file = "modeltuneSVM3.dat")
#load("modeltuneSVM.dat")

#opt_sigma = model_tune_SVM$bestTune$sigma
#opt_C = model_tune_SVM$bestTune$C

#------------------------------ LASSO prediction ----------------------------------

input_matrix = as.matrix(input)
output_matrix = as.matrix(output)
input_for_predict_matrix = as.matrix(input_for_predict)

#SVM_reg = ksvm(input_matrix, output_matrix, kernel ="rbfdot", kpar = list(sigma = opt_sigma),
#C = opt_C, epsilon = 0.1, prob.model = FALSE,
#class.weights = NULL, cross = 0, fit = TRUE, cache = 40,
#tol = 0.001, shrinking = TRUE)

SVM_reg = ksvm(input_matrix, output_matrix, kernel ="rbfdot", kpar = list(sigma = 0.02),
C = 2, epsilon = 0.1, prob.model = FALSE,
class.weights = NULL, cross = 0, fit = TRUE, cache = 40,
tol = 0.001, shrinking = TRUE)

output_predict = predict(SVM_reg, input_for_predict_matrix)

write.table(output_predict,"prediction_SVM.txt", quote=F, col.names=F)

































































plot(study_x[ ,2], std_study_x[ ,2])











study_x2 = study_x2[ , colnames(study_x2) != exclude_name2[1]]






for (ii in colname)
{
	if ((table(pearson_matrix2[,ii] ==1)[[2]] > 1) || (table(pearson_matrix2[,ii] ==-1)[[2]] > 0)){
		exclude_name2 = c(exclude_name2, ii) 
													}
}









if ((table(pearson_matrix2[,ii] ==1)[[2]] > 1) || (table(pearson_matrix2[,ii] ==-1)[[2]] > 0)){}




table(study_x$NComp ==1)[2] > 1


table(pearson_matrix2[,"NComp"] == 1)[[2]]
table(pearson_matrix2[ ,2] == 1)


for (ii in colname)
{
if (table(pearson_matrix2[,ii] == 1)[[2]] > 1){}
}







for (ii in colname)
{
	if (table(study_xii ==1)[2] > 1) exclude_name2 = c(exclude_name, ii)
		
}


if (count(pearson_matrix2[ ,2] ==1) > 1){exclude_name2 = c(exclude_name, ii)}


table(pearson_matrix2[ ,2] ==1)[2]

if (table(study_x[ ,2] ==1)[2] > 1)
 {
exclude_name2 = c(exclude_name2, "test")
}



Fe_ver1 = read.csv("input_for_prediction\\magpie_data_for_prediction_Fe_first_half.csv")
Fe_ver2 = read.csv("input_for_prediction\\magpie_data_for_prediction_Fe_second_half.csv")
Co_ver1 = read.csv("input_for_prediction\\magpie_data_for_prediction_Co_first_half.csv")
Co_ver2 = read.csv("input_for_prediction\\magpie_data_for_prediction_Co_second_half.csv")
Ni_ver1 = read.csv("input_for_prediction\\magpie_data_for_prediction_Ni_first_half.csv")
Ni_ver2 = read.csv("input_for_prediction\\magpie_data_for_prediction_Ni_second_half.csv")

#magpie_descriptorの最後の項目"Class"は取り除く
Fe_ver1 = Fe_ver1[ , colnames(Fe_ver1) != "Class"]
Fe_ver2 = Fe_ver2[ , colnames(Fe_ver2) != "Class"]
Co_ver1 = Co_ver1[ , colnames(Co_ver1) != "Class"]
Co_ver2 = Co_ver2[ , colnames(Co_ver2) != "Class"]
Ni_ver1 = Ni_ver1[ , colnames(Ni_ver1) != "Class"]
Ni_ver2 = Ni_ver2[ , colnames(Ni_ver2) != "Class"]

Fe_comp_ver1 = read.csv("input_for_prediction\\head_data_for_prediction_Fe_first_half.txt", header = F)
Fe_comp_ver2 = read.csv("input_for_prediction\\head_data_for_prediction_Fe_second_half.txt", header = F)
Co_comp_ver1 = read.csv("input_for_prediction\\head_data_for_prediction_Co_first_half.txt", header = F)
Co_comp_ver2 = read.csv("input_for_prediction\\head_data_for_prediction_Co_second_half.txt", header = F)
Ni_comp_ver1 = read.csv("input_for_prediction\\head_data_for_prediction_Ni_first_half.txt", header = F)
Ni_comp_ver2 = read.csv("input_for_prediction\\head_data_for_prediction_Ni_second_half.txt", header = F)

predict_x_magpie = rbind(Fe_ver1, Fe_ver2, Co_ver1, Co_ver2, Ni_ver1, Ni_ver2)
predict_x_comp = rbind(Fe_comp_ver1, Fe_comp_ver2, Co_comp_ver1, Co_comp_ver2, Ni_comp_ver1, Ni_comp_ver2)









shell("dir")







barplot(data$thickness_nm)

barplot(data$thermopower_uV)

barplot(data$resist_ohm)

write.csv(data, "test.csv")















#========================= FeCoIr ternary =======================

data = read.csv("FeCoIr_ternary_SL9.csv")

comp_Fe = data[ ,2]*100
comp_Co = data[ ,4]*100
comp_Ir = data[ ,6]*100
total_spin = data[ ,19]*1000

triangle_X=100-comp_Fe / 2-comp_Co
triangle_Y=comp_Fe * 0.866

triangle_X = (triangle_X - min(triangle_X))/(max(triangle_X)-min(triangle_X))*100
triangle_Y = (triangle_Y - min(triangle_Y))/(max(triangle_Y)-min(triangle_Y))*86.5

scale_max = 2.41 *1000
scale_min = 1.97 * 1000

total_spin = (total_spin - scale_min)/(scale_max-scale_min)
total_spin = round(total_spin*1000)

total_spin_col =  abs((total_spin - 1000))+1

#spinV[spinV < 1] = 1

#col_start = 0.4
#col_end = 1

col_start = 0
col_end = 0.8

col_number = 1001

spheres3d(triangle_X, triangle_Y, 0, col = rainbow(col_number, start=col_start, end=col_end)[total_spin_col],
          radius = 2.4    )

z=5
nlwd = 0.5

rgl.viewpoint( theta = 0, phi = 0, fov = 0, zoom = 1)
par3d("windowRect" = c(137, 10, 744, 544))
segments3d(matrix(c(0, 0, z, 50, 86.6, z), nrow=2, ncol=3, byrow=T),lwd=nlwd)  #
segments3d(matrix(c(0, 0, z, 100, 0, z), nrow=2, ncol=3, byrow=T),lwd=nlwd)  #
segments3d(matrix(c(100, 0, z, 50, 86.6, z), nrow=2, ncol=3, byrow=T),lwd=nlwd)  #
segments3d(matrix(c(33.3, 0, z, 66.6, 57.67, z), nrow=2, ncol=3, byrow=T), lwd=nlwd)　#
segments3d(matrix(c(66.6, 0, z, 83.27, 28.87, z), nrow=2, ncol=3, byrow=T), lwd=nlwd) #
segments3d(matrix(c(33.3, 0, z, 16.67, 28.87, z), nrow=2, ncol=3, byrow=T), lwd=nlwd) #
segments3d(matrix(c(66.6, 0, z, 33.3, 57.67, z), nrow=2, ncol=3, byrow=T), lwd=nlwd) #
segments3d(matrix(c(16.67, 28.87, z, 83.27, 28.87, z), nrow=2, ncol=3, byrow=T), lwd=nlwd) #
segments3d(matrix(c(33.3, 57.67, z, 66.6, 57.67, z), nrow=2, ncol=3, byrow=T), lwd=nlwd) #

rgl.snapshot(".png")








#========================= Fe75C025IrPt ternary =======================

data = read.csv("Fe75Co25IrPt_removed.csv")

comp_Fe = data[ ,2]*100
comp_Co = data[ ,3]*100
comp_Ir = data[ ,5]*100
comp_Pt = data[ ,7]*100
total_spin = data[ ,8]*1000

triangle_X=100-comp_Ir / 2-comp_Pt
triangle_Y=comp_Ir * 0.866

triangle_X = (triangle_X - min(triangle_X))/(max(triangle_X)-min(triangle_X))*100
triangle_Y = (triangle_Y - min(triangle_Y))/(max(triangle_Y)-min(triangle_Y))*86.5

total_spin = (total_spin - min(total_spin))/(max(total_spin)-min(total_spin))
total_spin = round(total_spin*1000)+1

total_spin_col =  abs((total_spin - max(total_spin)))+1

#spinV[spinV < 1] = 1

#col_start = 0.4
#col_end = 1

col_start = 0
col_end = 0.8

col_number = max(total_spin)

spheres3d(triangle_X, triangle_Y, 0, col = rainbow(col_number, start=col_start, end=col_end)[total_spin_col],
          radius = 2.4    )
#spheres3d(triangle_X, triangle_Y, 0, col = heat.colors(col_number)[total_spin],
#          radius = 2    )

#spheres3d(triangle_X, triangle_Y, 0, col = rainbow(col_number, start=0, end=1)[total_spin],
#          radius = 2    )




rgl.viewpoint( theta = 0, phi = 0, fov = 0, zoom = 1)
par3d("windowRect" = c(137, 10, 744, 544))
segments3d(matrix(c(0, 0, 0, 50, 86.6, 0), nrow=2, ncol=3, byrow=T))
segments3d(matrix(c(0, 0, 0, 100, 0, 0), nrow=2, ncol=3, byrow=T))
segments3d(matrix(c(100, 0, 0, 50, 86.6, 0), nrow=2, ncol=3, byrow=T))
segments3d(matrix(c(20, 0, 0, 60, 69.28, 0), nrow=2, ncol=3, byrow=T), lwd=0.1)
segments3d(matrix(c(40, 0, 0, 70, 51.96, 0), nrow=2, ncol=3, byrow=T), lwd=0.1)
segments3d(matrix(c(60, 0, 0, 80, 34.64, 0), nrow=2, ncol=3, byrow=T), lwd=0.1)
segments3d(matrix(c(80, 0, 0, 90, 17.32, 0), nrow=2, ncol=3, byrow=T), lwd=0.1)
segments3d(matrix(c(20, 0, 0, 10, 17.32, 0), nrow=2, ncol=3, byrow=T), lwd=0.1)
segments3d(matrix(c(40, 0, 0, 20, 34.64, 0), nrow=2, ncol=3, byrow=T), lwd=0.1)
segments3d(matrix(c(60, 0, 0, 30, 51.96, 0), nrow=2, ncol=3, byrow=T), lwd=0.1)
segments3d(matrix(c(80, 0, 0, 40, 69.28, 0), nrow=2, ncol=3, byrow=T), lwd=0.1)
segments3d(matrix(c(10, 17.32, 0, 90, 17.32, 0), nrow=2, ncol=3, byrow=T), lwd=0.1)
segments3d(matrix(c(20, 34.64, 0, 80, 34.64, 0), nrow=2, ncol=3, byrow=T), lwd=0.1)
segments3d(matrix(c(30, 51.96, 0, 70, 51.96, 0), nrow=2, ncol=3, byrow=T), lwd=0.1)
segments3d(matrix(c(40, 69.28, 0, 60, 69.28, 0), nrow=2, ncol=3, byrow=T), lwd=0.1)
text3d(c(-5, -5, 0), texts=composition_name2, cex=1.5)
text3d(c(105, -5, 0), texts=composition_name3, cex=1.5)
text3d(c(50, 91.6, 0), texts=composition_name1, cex=1.5)

rgl.snapshot("volt.png")

#-------------------color bar ------------------------------

X = numeric(100)
Y = seq(1,100, by=1)

plot(X, Y, pch=15, col=rainbow(100, start=col_start, end=col_end)[Y])



#========================= FeCoIr ternary =======================

data = read.csv("FeCoIr_ternary_SL9.csv")

comp_Fe = data[ ,2]*100
comp_Co = data[ ,4]*100
comp_Ir = data[ ,6]*100
total_spin = data[ ,19]*1000

triangle_X=100-comp_Fe / 2-comp_Co
triangle_Y=comp_Fe * 0.866

triangle_X = (triangle_X - min(triangle_X))/(max(triangle_X)-min(triangle_X))*100
triangle_Y = (triangle_Y - min(triangle_Y))/(max(triangle_Y)-min(triangle_Y))*86.5

scale_max = 2.4 *1000
scale_min = 2.0 * 1000

total_spin = (total_spin - scale_min)/(scale_max-scale_min)
total_spin = round(total_spin*1000)+1

total_spin_col =  abs((total_spin - max(total_spin)))+1

#spinV[spinV < 1] = 1

#col_start = 0.4
#col_end = 1

col_start = 0
col_end = 0.8

col_number = 1001

spheres3d(triangle_X, triangle_Y, 0, col = rainbow(col_number, start=col_start, end=col_end)[total_spin_col],
          radius = 2.4    )

z=5
nlwd = 0.5

rgl.viewpoint( theta = 0, phi = 0, fov = 0, zoom = 1)
par3d("windowRect" = c(137, 10, 744, 544))
segments3d(matrix(c(0, 0, z, 50, 86.6, z), nrow=2, ncol=3, byrow=T),lwd=nlwd)  #
segments3d(matrix(c(0, 0, z, 100, 0, z), nrow=2, ncol=3, byrow=T),lwd=nlwd)  #
segments3d(matrix(c(100, 0, z, 50, 86.6, z), nrow=2, ncol=3, byrow=T),lwd=nlwd)  #
segments3d(matrix(c(33.3, 0, z, 66.6, 57.67, z), nrow=2, ncol=3, byrow=T), lwd=nlwd)　#
segments3d(matrix(c(66.6, 0, z, 83.27, 28.87, z), nrow=2, ncol=3, byrow=T), lwd=nlwd) #
segments3d(matrix(c(33.3, 0, z, 16.67, 28.87, z), nrow=2, ncol=3, byrow=T), lwd=nlwd) #
segments3d(matrix(c(66.6, 0, z, 33.3, 57.67, z), nrow=2, ncol=3, byrow=T), lwd=nlwd) #
segments3d(matrix(c(16.67, 28.87, z, 83.27, 28.87, z), nrow=2, ncol=3, byrow=T), lwd=nlwd) #
segments3d(matrix(c(33.3, 57.67, z, 66.6, 57.67, z), nrow=2, ncol=3, byrow=T), lwd=nlwd) #

rgl.snapshot(".png")





















#------------------FeCoIr_asdepo-----------------

filenamelist = list.files("C:/ICF_Autocapsule_Disabled/R/FeCoX_combiVSM/NEC/NEC/20181218_FeCoIr_asdepo", full.names=T)
aveNo = 10
FeCoIr_asdepo = numeric(20)

for (ii in 1:20)
{
	data = read.table(filenamelist[ii], skip=41, sep=",")
	
	average1 = mean(data[1:aveNo,3])
	average2 = mean(data[(752-aveNo):752,3])*(-1)
	average3 = mean(data[752:(752+aveNo),3])*(-1)
	average4 = mean(data[(1504-aveNo):1504,3])
	average = mean(average1, average2, average3, average4)

	FeCoIr_asdepo[ii] = average
}

#------------------FeCoIr_asdepo-----------------

filenamelist = list.files("C:/ICF_Autocapsule_Disabled/R/FeCoX_combiVSM/NEC/NEC/20181219_FeCoNi_asdepo", full.names=T)
aveNo = 10
FeCoNi_asdepo = numeric(20)

for (ii in 1:20)
{
	data = read.table(filenamelist[ii], skip=41, sep=",")
	
	average1 = mean(data[1:aveNo,3])
	average2 = mean(data[(752-aveNo):752,3])*(-1)
	average3 = mean(data[752:(752+aveNo),3])*(-1)
	average4 = mean(data[(1504-aveNo):1504,3])
	average = mean(average1, average2, average3, average4)

	FeCoNi_asdepo[ii] = average
}

#------------------FeCoPt_asdepo-----------------

filenamelist = list.files("C:/ICF_Autocapsule_Disabled/R/FeCoX_combiVSM/NEC/NEC/20181219_FeCoPt_asdepo", full.names=T)
aveNo = 10
FeCoPt_asdepo = numeric(20)

for (ii in 1:20)
{
	data = read.table(filenamelist[ii], skip=41, sep=",")
	
	average1 = mean(data[1:aveNo,3])
	average2 = mean(data[(752-aveNo):752,3])*(-1)
	average3 = mean(data[752:(752+aveNo),3])*(-1)
	average4 = mean(data[(1504-aveNo):1504,3])
	average = mean(average1, average2, average3, average4)

	FeCoPt_asdepo[ii] = average
}


plot(FeCoPt_asdepo)


#------------------FeCoIr_600deg-----------------

filenamelist = list.files("C:/ICF_Autocapsule_Disabled/R/FeCoX_combiVSM/NEC/NEC/20181220_FeCoIr_600deg", full.names=T)
aveNo = 10
FeCoIr_600deg = numeric(20)

for (ii in 1:20)
{
	data = read.table(filenamelist[ii], skip=41, sep=",")
	
	average1 = mean(data[1:aveNo,3])
	average2 = mean(data[(752-aveNo):752,3])*(-1)
	average3 = mean(data[752:(752+aveNo),3])*(-1)
	average4 = mean(data[(1504-aveNo):1504,3])
	average = mean(average1, average2, average3, average4)

	FeCoIr_600deg[ii] = average
}

plot(FeCoIr_600deg)

#------------------FeCoPt_600deg-----------------

filenamelist = list.files("C:/ICF_Autocapsule_Disabled/R/FeCoX_combiVSM/NEC/NEC/20181220_FeCoPt_600deg", full.names=T)
aveNo = 10
FeCoPt_600deg = numeric(20)

for (ii in 1:20)
{
	data = read.table(filenamelist[ii], skip=41, sep=",")
	
	average1 = mean(data[1:aveNo,3])
	average2 = mean(data[(752-aveNo):752,3])*(-1)
	average3 = mean(data[752:(752+aveNo),3])*(-1)
	average4 = mean(data[(1504-aveNo):1504,3])
	average = mean(average1, average2, average3, average4)

	FeCoPt_600deg[ii] = average
}


#------------------FeCoNi_600deg-----------------

filenamelist = list.files("C:/ICF_Autocapsule_Disabled/R/FeCoX_combiVSM/NEC/NEC/20181221_FeCoNi_600deg", full.names=T)
aveNo = 10
FeCoNi_600deg = numeric(20)

for (ii in 1:20)
{
	data = read.table(filenamelist[ii], skip=41, sep=",")
	
	average1 = mean(data[1:aveNo,3])
	average2 = mean(data[(752-aveNo):752,3])*(-1)
	average3 = mean(data[752:(752+aveNo),3])*(-1)
	average4 = mean(data[(1504-aveNo):1504,3])
	average = mean(average1, average2, average3, average4)

	FeCoNi_600deg[ii] = average
}

plot(FeCoNi_600deg)

#--------------------   graph raw     ---------------------

composition = read.csv("composition_Pt_Ir_Ni.csv")
Ptcomposition = composition[ ,2]
Ircomposition = composition[ ,3]
Nicomposition = composition[ ,4]

#1は端っこなので外す
#as-depo sample

plot(Ptcomposition, FeCoPt_asdepo[2:20], col ="red", pch = 16, cex = 2, xlim=c(0,8), ylim=c(0.00040, 0.00073))
par(new=T)
plot(Ircomposition, FeCoIr_asdepo[2:20], col ="blue", pch = 16, cex = 2, xlim=c(0,8), ylim=c(0.00040, 0.00073))
par(new=T)
plot(Nicomposition, FeCoNi_asdepo[2:20], col ="black", pch = 16, cex = 2, xlim=c(0,8), ylim=c(0.00040, 0.00073))

#600deg samples

plot(Ptcomposition, FeCoPt_600deg[2:20], col ="red", pch = 16, cex = 2, xlim=c(0,8), ylim=c(0.00040, 0.00073))
par(new=T)
plot(Ircomposition, FeCoIr_600deg[2:20], col ="blue", pch = 16, cex = 2, xlim=c(0,8), ylim=c(0.00040, 0.00073))
par(new=T)
plot(Nicomposition, FeCoNi_600deg[2:20], col ="black", pch = 16, cex = 2, xlim=c(0,8), ylim=c(0.00040, 0.00073))


#--------------------   graph wmu/cc     ---------------------

sampleX = 0.80 #mm
sampleY = 3.80 #mm
thickness = 100 #nm
area = sampleX * sampleY * thickness * 10^-9


composition = read.csv("composition_Pt_Ir_Ni.csv")
Ptcomposition = composition[ ,2]
Ircomposition = composition[ ,3]
Nicomposition = composition[ ,4]

#1は端っこなので外す
#as-depo sample

plot(Ptcomposition, FeCoPt_asdepo[2:20]/area, col ="red", pch = 16, cex = 2, xlim=c(0,8), ylim=c(1400, 2300))
par(new=T)
plot(Ircomposition, FeCoIr_asdepo[2:20]/area, col ="blue", pch = 16, cex = 2, xlim=c(0,8), ylim=c(1400, 2300))
par(new=T)
plot(Nicomposition, FeCoNi_asdepo[2:20]/area, col ="black", pch = 16, cex = 2, xlim=c(0,8), ylim=c(1400, 2300))

#600deg samples

plot(Ptcomposition, FeCoPt_600deg[2:20]/area, col ="red", pch = 16, cex = 2, xlim=c(0,8), ylim=c(1400, 2300))
par(new=T)
plot(Ircomposition, FeCoIr_600deg[2:20]/area, col ="blue", pch = 16, cex = 2, xlim=c(0,8), ylim=c(1400, 2300))
par(new=T)
plot(Nicomposition, FeCoNi_600deg[2:20]/area, col ="black", pch = 16, cex = 2, xlim=c(0,8), ylim=c(1400, 2300))


FeCoPt_asdepo[1]/area












#-------------------- composition interpolation spline ----------

rawdata = read.csv("composition_ver2.csv")
position = read.csv("position.csv", head = F)

spPt = splinefun(rawdata[,1], rawdata[,2])

Ptcomposition = numeric(20)

for (ii in  0:20)
{
	Ptcomposition[ii] = spPt(0.5+ii) 
}

Ptcomposition

plot(Ptcomposition)

plot(rawdata[ ,1], rawdata[, 2], xlim=c(0,25),ylim=c(0,9))
par(new=T)
plot(1:20, Ptcomposition, xlim=c(0,25), ylim=c(0,9))
write.csv(Ptcomposition, "Ptcomposition.csv")

-----

rawdata = read.csv("composition_ver2.csv")
position = read.csv("position.csv", head = F)

spIr = splinefun(rawdata[,1], rawdata[,3])

Ircomposition = numeric(20)

for (ii in  0:20)
{
	Ircomposition[ii] = spIr(0.5+ii) 
}

Ircomposition

plot(Ircomposition)

plot(rawdata[ ,1], rawdata[, 3], xlim=c(0,25),ylim=c(0,9))
par(new=T)
plot(1:20, Ircomposition, xlim=c(0,25), ylim=c(0,9))
write.csv(Ircomposition, "Ircomposition.csv")

-----

rawdata = read.csv("composition_ver2.csv")
position = read.csv("position.csv", head = F)

spNi = splinefun(rawdata[,1], rawdata[,4])

Nicomposition = numeric(20)

for (ii in  0:20)
{
	Nicomposition[ii] = spNi(0.5+ii) 
}

Nicomposition

plot(Nicomposition)

plot(rawdata[ ,1], rawdata[, 4], xlim=c(0,25),ylim=c(0,9))
par(new=T)
plot(1:20, Nicomposition, xlim=c(0,25), ylim=c(0,9))
write.csv(Nicomposition, "Nicomposition.csv")









plot(FeCoPt_asdepo, ylim = c(0.0005, 0.0007), col=1)
par(new=T)
plot(FeCoIr_asdepo, ylim = c(0.0005, 0.0007), col=2)
par(new=T)
plot(FeCoNi_asdepo, ylim = c(0.0005, 0.0007), col=3)












#-------------------- FePtEr data -----------------------------

conbi_file = read.csv("numeric_FeCoPt_FeCoIr.csv")

#=============== plot for spinV ============================

triangle_X=100-conbi_file[,1] / 2-conbi_file[,2]
triangle_Y=conbi_file[,1] * 0.866

col_start = 0.6
col_end = 1


spheres3d(triangle_X, triangle_Y, 0, col = 1,
          radius = 0.1    )
rgl.viewpoint( theta = 0, phi = 0, fov = 0, zoom = 1)
par3d("windowRect" = c(137, 10, 744, 544))
segments3d(matrix(c(0, 0, 0, 50, 86.6, 0), nrow=2, ncol=3, byrow=T))
segments3d(matrix(c(0, 0, 0, 100, 0, 0), nrow=2, ncol=3, byrow=T))
segments3d(matrix(c(100, 0, 0, 50, 86.6, 0), nrow=2, ncol=3, byrow=T))
segments3d(matrix(c(20, 0, 0, 60, 69.28, 0), nrow=2, ncol=3, byrow=T), lwd=0.1)
segments3d(matrix(c(40, 0, 0, 70, 51.96, 0), nrow=2, ncol=3, byrow=T), lwd=0.1)
segments3d(matrix(c(60, 0, 0, 80, 34.64, 0), nrow=2, ncol=3, byrow=T), lwd=0.1)
segments3d(matrix(c(80, 0, 0, 90, 17.32, 0), nrow=2, ncol=3, byrow=T), lwd=0.1)
segments3d(matrix(c(20, 0, 0, 10, 17.32, 0), nrow=2, ncol=3, byrow=T), lwd=0.1)
segments3d(matrix(c(40, 0, 0, 20, 34.64, 0), nrow=2, ncol=3, byrow=T), lwd=0.1)
segments3d(matrix(c(60, 0, 0, 30, 51.96, 0), nrow=2, ncol=3, byrow=T), lwd=0.1)
segments3d(matrix(c(80, 0, 0, 40, 69.28, 0), nrow=2, ncol=3, byrow=T), lwd=0.1)
segments3d(matrix(c(10, 17.32, 0, 90, 17.32, 0), nrow=2, ncol=3, byrow=T), lwd=0.1)
segments3d(matrix(c(20, 34.64, 0, 80, 34.64, 0), nrow=2, ncol=3, byrow=T), lwd=0.1)
segments3d(matrix(c(30, 51.96, 0, 70, 51.96, 0), nrow=2, ncol=3, byrow=T), lwd=0.1)
segments3d(matrix(c(40, 69.28, 0, 60, 69.28, 0), nrow=2, ncol=3, byrow=T), lwd=0.1)
text3d(c(-5, -5, 0), texts="Pt", cex=1.5)
text3d(c(105, -5, 0), texts="Co", cex=1.5)
text3d(c(50, 91.6, 0), texts="Fe", cex=1.5)

rgl.snapshot("volt.png")

xx = which.max(spinV)
conbi_file[xx, ]

#=============== plot for spinP ============================


triangle_X=100-conbi_file[,3] / 2-conbi_file[,4]
triangle_Y=conbi_file[,3] * 0.866

spinP = conbi_file$spinP
spinP = (spinP*100)
spinP = round(spinP)
resist = round(conbi_file$resist)
thickness = conbi_file$thickness

spinP[spinP < 1] = 1

col_start = 0.6
col_end = 1

#color_cut_off = mean(spinP)
#spinP[spinP < color_cut_off] = color_cut_off
#spinP = spinP - color_cut_off + 1
col_number = max(spinP)

spheres3d(triangle_X, triangle_Y, 0, col = rainbow(col_number, start=col_start, end=col_end)[spinP],
          radius = 2    )
rgl.viewpoint( theta = 0, phi = 0, fov = 0, zoom = 1)
par3d("windowRect" = c(137, 10, 744, 544))
segments3d(matrix(c(0, 0, 0, 50, 86.6, 0), nrow=2, ncol=3, byrow=T))
segments3d(matrix(c(0, 0, 0, 100, 0, 0), nrow=2, ncol=3, byrow=T))
segments3d(matrix(c(100, 0, 0, 50, 86.6, 0), nrow=2, ncol=3, byrow=T))
segments3d(matrix(c(20, 0, 0, 60, 69.28, 0), nrow=2, ncol=3, byrow=T), lwd=0.1)
segments3d(matrix(c(40, 0, 0, 70, 51.96, 0), nrow=2, ncol=3, byrow=T), lwd=0.1)
segments3d(matrix(c(60, 0, 0, 80, 34.64, 0), nrow=2, ncol=3, byrow=T), lwd=0.1)
segments3d(matrix(c(80, 0, 0, 90, 17.32, 0), nrow=2, ncol=3, byrow=T), lwd=0.1)
segments3d(matrix(c(20, 0, 0, 10, 17.32, 0), nrow=2, ncol=3, byrow=T), lwd=0.1)
segments3d(matrix(c(40, 0, 0, 20, 34.64, 0), nrow=2, ncol=3, byrow=T), lwd=0.1)
segments3d(matrix(c(60, 0, 0, 30, 51.96, 0), nrow=2, ncol=3, byrow=T), lwd=0.1)
segments3d(matrix(c(80, 0, 0, 40, 69.28, 0), nrow=2, ncol=3, byrow=T), lwd=0.1)
segments3d(matrix(c(10, 17.32, 0, 90, 17.32, 0), nrow=2, ncol=3, byrow=T), lwd=0.1)
segments3d(matrix(c(20, 34.64, 0, 80, 34.64, 0), nrow=2, ncol=3, byrow=T), lwd=0.1)
segments3d(matrix(c(30, 51.96, 0, 70, 51.96, 0), nrow=2, ncol=3, byrow=T), lwd=0.1)
segments3d(matrix(c(40, 69.28, 0, 60, 69.28, 0), nrow=2, ncol=3, byrow=T), lwd=0.1)
text3d(c(-5, -5, 0), texts=composition_name2, cex=1.5)
text3d(c(105, -5, 0), texts=composition_name3, cex=1.5)
text3d(c(50, 91.6, 0), texts=composition_name1, cex=1.5)
rgl.snapshot("power.png")

yy = which.max(spinP)
conbi_file[yy, ]

#=============== plot for resist ============================

triangle_X=100-conbi_file[,3] / 2-conbi_file[,4]
triangle_Y=conbi_file[,3] * 0.866

resist = round(conbi_file$resist * 10)

col_start = 0.6
col_end = 1
col_number = max(resist)

spheres3d(triangle_X, triangle_Y, 0, col = rainbow(col_number, start=col_start, end=col_end)[resist],
          radius = 2    )
rgl.viewpoint( theta = 0, phi = 0, fov = 0, zoom = 1)
par3d("windowRect" = c(137, 10, 744, 544))
segments3d(matrix(c(0, 0, 0, 50, 86.6, 0), nrow=2, ncol=3, byrow=T))
segments3d(matrix(c(0, 0, 0, 100, 0, 0), nrow=2, ncol=3, byrow=T))
segments3d(matrix(c(100, 0, 0, 50, 86.6, 0), nrow=2, ncol=3, byrow=T))
segments3d(matrix(c(20, 0, 0, 60, 69.28, 0), nrow=2, ncol=3, byrow=T), lwd=0.1)
segments3d(matrix(c(40, 0, 0, 70, 51.96, 0), nrow=2, ncol=3, byrow=T), lwd=0.1)
segments3d(matrix(c(60, 0, 0, 80, 34.64, 0), nrow=2, ncol=3, byrow=T), lwd=0.1)
segments3d(matrix(c(80, 0, 0, 90, 17.32, 0), nrow=2, ncol=3, byrow=T), lwd=0.1)
segments3d(matrix(c(20, 0, 0, 10, 17.32, 0), nrow=2, ncol=3, byrow=T), lwd=0.1)
segments3d(matrix(c(40, 0, 0, 20, 34.64, 0), nrow=2, ncol=3, byrow=T), lwd=0.1)
segments3d(matrix(c(60, 0, 0, 30, 51.96, 0), nrow=2, ncol=3, byrow=T), lwd=0.1)
segments3d(matrix(c(80, 0, 0, 40, 69.28, 0), nrow=2, ncol=3, byrow=T), lwd=0.1)
segments3d(matrix(c(10, 17.32, 0, 90, 17.32, 0), nrow=2, ncol=3, byrow=T), lwd=0.1)
segments3d(matrix(c(20, 34.64, 0, 80, 34.64, 0), nrow=2, ncol=3, byrow=T), lwd=0.1)
segments3d(matrix(c(30, 51.96, 0, 70, 51.96, 0), nrow=2, ncol=3, byrow=T), lwd=0.1)
segments3d(matrix(c(40, 69.28, 0, 60, 69.28, 0), nrow=2, ncol=3, byrow=T), lwd=0.1)
text3d(c(-5, -5, 0), texts=composition_name2, cex=1.5)
text3d(c(105, -5, 0), texts=composition_name3, cex=1.5)
text3d(c(50, 91.6, 0), texts=composition_name1, cex=1.5)
rgl.snapshot("resist.png")

zz = which.max(resist)
conbi_file[zz, ]

#=============== plot for thickness ============================

triangle_X=100-conbi_file$Fe/2-conbi_file$Pt
triangle_Y=conbi_file$Fe * 0.866

thickness = conbi_file$thickness

col_start = 0.6
col_end = 1
dev.new()

col_number = max(thickness)

spheres3d(triangle_X, triangle_Y, 0, col = rainbow(col_number, start=col_start, end=col_end)[thickness],
          radius = 3    )
rgl.viewpoint( theta = 0, phi = 0, fov = 0, zoom = 1)
par3d("windowRect" = c(137, 10, 744, 544))
segments3d(matrix(c(0, 0, 0, 50, 86.6, 0), nrow=2, ncol=3, byrow=T))
segments3d(matrix(c(0, 0, 0, 100, 0, 0), nrow=2, ncol=3, byrow=T))
segments3d(matrix(c(100, 0, 0, 50, 86.6, 0), nrow=2, ncol=3, byrow=T))
segments3d(matrix(c(20, 0, 0, 60, 69.28, 0), nrow=2, ncol=3, byrow=T), lwd=0.1)
segments3d(matrix(c(40, 0, 0, 70, 51.96, 0), nrow=2, ncol=3, byrow=T), lwd=0.1)
segments3d(matrix(c(60, 0, 0, 80, 34.64, 0), nrow=2, ncol=3, byrow=T), lwd=0.1)
segments3d(matrix(c(80, 0, 0, 90, 17.32, 0), nrow=2, ncol=3, byrow=T), lwd=0.1)
segments3d(matrix(c(20, 0, 0, 10, 17.32, 0), nrow=2, ncol=3, byrow=T), lwd=0.1)
segments3d(matrix(c(40, 0, 0, 20, 34.64, 0), nrow=2, ncol=3, byrow=T), lwd=0.1)
segments3d(matrix(c(60, 0, 0, 30, 51.96, 0), nrow=2, ncol=3, byrow=T), lwd=0.1)
segments3d(matrix(c(80, 0, 0, 40, 69.28, 0), nrow=2, ncol=3, byrow=T), lwd=0.1)
segments3d(matrix(c(10, 17.32, 0, 90, 17.32, 0), nrow=2, ncol=3, byrow=T), lwd=0.1)
segments3d(matrix(c(20, 34.64, 0, 80, 34.64, 0), nrow=2, ncol=3, byrow=T), lwd=0.1)
segments3d(matrix(c(30, 51.96, 0, 70, 51.96, 0), nrow=2, ncol=3, byrow=T), lwd=0.1)
segments3d(matrix(c(40, 69.28, 0, 60, 69.28, 0), nrow=2, ncol=3, byrow=T), lwd=0.1)
text3d(c(-5, -5, 0), texts=composition_name2, cex=1.5)
text3d(c(105, -5, 0), texts=composition_name3, cex=1.5)
text3d(c(50, 91.6, 0), texts=composition_name1, cex=1.5)
























#=============== plot for spinV ============================

triangle_X=100-conbi_file[,3] / 2-conbi_file[,4]
triangle_Y=conbi_file[,3] * 0.866

spinV = round(conbi_file$spinV * 100000000)
resist = round(conbi_file$resist)
thickness = conbi_file$thickness

spinV[spinV < 1] = 1

col_start = 0.6
col_end = 1
dev.new()

col_number = max(spinV)


plot(triangle_X, triangle_Y, col = rainbow(col_number, start=col_start, end=col_end)[spinV]
     ,pch=19, cex=1.5, xlim=c(-10, 110), ylim=c(-10, 100),ann=F)
segments(0, 0, 50, 86.6)
segments(0, 0, 100, 0)
segments(100, 0, 50, 86.6)
text(0, 0, composition_name2, cex=2.5)
text(100, 0, composition_name3, cex=2.5)
text(50, 86.6, composition_name1, cex=2.5)
arrows(0, 20, 25, 60)
segments(20, 0, 60, 69.28, lty="dotted", lwd=0.5)
segments(40, 0, 70, 51.96, lty="dotted", lwd=0.5)
segments(60, 0, 80, 34.64, lty="dotted", lwd=0.5)
segments(80, 0, 90, 17.32, lty="dotted", lwd=0.5)
segments(20, 0, 10, 17.32, lty="dotted", lwd=0.5)
segments(40, 0, 20, 34.64, lty="dotted", lwd=0.5)
segments(60, 0, 30, 51.96, lty="dotted", lwd=0.5)
segments(80, 0, 40, 69.28, lty="dotted", lwd=0.5)
segments(10, 17.32, 90, 17.32, lty="dotted", lwd=0.5)
segments(20, 34.64, 80, 34.64, lty="dotted", lwd=0.5)
segments(30, 51.96, 70, 51.96, lty="dotted", lwd=0.5)
segments(40, 69.28, 60, 69.28, lty="dotted", lwd=0.5)

xx = which.max(spinV)
conbi_file[xx, ]

#=============== plot for spinP ============================


triangle_X=100-conbi_file[,3] / 2-conbi_file[,4]
triangle_Y=conbi_file[,3] * 0.866

spinP = round(conbi_file$spinP * 100)
resist = round(conbi_file$resist)
thickness = conbi_file$thickness

spinP[spinP < 1] = 1

col_start = 0.6
col_end = 1
dev.new()

#color_cut_off = mean(spinP)
#spinP[spinP < color_cut_off] = color_cut_off
#spinP = spinP - color_cut_off + 1
col_number = max(spinP)

plot(triangle_X, triangle_Y, col = rainbow(col_number, start=col_start, end=col_end)[spinP]
     ,pch=19, cex=1.5, xlim=c(-10, 110), ylim=c(-10, 100),ann=F)


segments(0, 0, 50, 86.6)
segments(0, 0, 100, 0)
segments(100, 0, 50, 86.6)
text(0, 0, composition_name2, cex=2.5)
text(100, 0, composition_name3, cex=2.5)
text(50, 86.6, composition_name1, cex=2.5)
arrows(0, 20, 25, 60)
segments(20, 0, 60, 69.28, lty="dotted", lwd=0.5)
segments(40, 0, 70, 51.96, lty="dotted", lwd=0.5)
segments(60, 0, 80, 34.64, lty="dotted", lwd=0.5)
segments(80, 0, 90, 17.32, lty="dotted", lwd=0.5)
segments(20, 0, 10, 17.32, lty="dotted", lwd=0.5)
segments(40, 0, 20, 34.64, lty="dotted", lwd=0.5)
segments(60, 0, 30, 51.96, lty="dotted", lwd=0.5)
segments(80, 0, 40, 69.28, lty="dotted", lwd=0.5)
segments(10, 17.32, 90, 17.32, lty="dotted", lwd=0.5)
segments(20, 34.64, 80, 34.64, lty="dotted", lwd=0.5)
segments(30, 51.96, 70, 51.96, lty="dotted", lwd=0.5)
segments(40, 69.28, 60, 69.28, lty="dotted", lwd=0.5)

yy = which.max(spinP)
conbi_file[yy, ]

#=============== plot for resist ============================

triangle_X=100-conbi_file[,3] / 2-conbi_file[,4]
triangle_Y=conbi_file[,3] * 0.866


resist = round(conbi_file$resist * 10)


col_start = 0.6
col_end = 1
dev.new()

col_number = max(resist)

plot(triangle_X, triangle_Y, col = rainbow(col_number, start=col_start, end=col_end)[resist]
     ,pch=19, cex=1.5, xlim=c(-10, 110), ylim=c(-10, 100),ann=F)
segments(0, 0, 50, 86.6)
segments(0, 0, 100, 0)
segments(100, 0, 50, 86.6)
text(0, 0, composition_name2, cex=2.5)
text(100, 0, composition_name3, cex=2.5)
text(50, 86.6, composition_name1, cex=2.5)
arrows(0, 20, 25, 60)
segments(20, 0, 60, 69.28, lty="dotted", lwd=0.5)
segments(40, 0, 70, 51.96, lty="dotted", lwd=0.5)
segments(60, 0, 80, 34.64, lty="dotted", lwd=0.5)
segments(80, 0, 90, 17.32, lty="dotted", lwd=0.5)
segments(20, 0, 10, 17.32, lty="dotted", lwd=0.5)
segments(40, 0, 20, 34.64, lty="dotted", lwd=0.5)
segments(60, 0, 30, 51.96, lty="dotted", lwd=0.5)
segments(80, 0, 40, 69.28, lty="dotted", lwd=0.5)
segments(10, 17.32, 90, 17.32, lty="dotted", lwd=0.5)
segments(20, 34.64, 80, 34.64, lty="dotted", lwd=0.5)
segments(30, 51.96, 70, 51.96, lty="dotted", lwd=0.5)
segments(40, 69.28, 60, 69.28, lty="dotted", lwd=0.5)

zz = which.max(resist)
conbi_file[zz, ]

#=============== plot for thickness ============================

triangle_X=100-conbi_file$Fe/2-conbi_file$Pt
triangle_Y=conbi_file$Fe * 0.866

thickness = conbi_file$thickness

col_start = 0.6
col_end = 1
dev.new()

col_number = max(thickness)

plot(triangle_X, triangle_Y, col = rainbow(col_number, start=col_start, end=col_end)[thickness]
     ,pch=19, cex=1.5, xlim=c(-10, 110), ylim=c(-10, 100),ann=F)
segments(0, 0, 50, 86.6)
segments(0, 0, 100, 0)
segments(100, 0, 50, 86.6)
text(0, 0, "Pt", cex=2.5)
text(100, 0, "Co", cex=2.5)
text(50, 86.6, "Fe", cex=2.5)
arrows(0, 20, 25, 60)
segments(20, 0, 60, 69.28, lty="dotted", lwd=0.5)
segments(40, 0, 70, 51.96, lty="dotted", lwd=0.5)
segments(60, 0, 80, 34.64, lty="dotted", lwd=0.5)
segments(80, 0, 90, 17.32, lty="dotted", lwd=0.5)
segments(20, 0, 10, 17.32, lty="dotted", lwd=0.5)
segments(40, 0, 20, 34.64, lty="dotted", lwd=0.5)
segments(60, 0, 30, 51.96, lty="dotted", lwd=0.5)
segments(80, 0, 40, 69.28, lty="dotted", lwd=0.5)
segments(10, 17.32, 90, 17.32, lty="dotted", lwd=0.5)
segments(20, 34.64, 80, 34.64, lty="dotted", lwd=0.5)
segments(30, 51.96, 70, 51.96, lty="dotted", lwd=0.5)
segments(40, 69.28, 60, 69.28, lty="dotted", lwd=0.5)












