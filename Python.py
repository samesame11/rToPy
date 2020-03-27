#======================== library ==============================
import pandas as pd
from itertools import combinations 
from itertools import permutations 
import numpy as np
#======================= read data ============================
def data_ready(): 
    data = pd.read_csv("C:\\ICF_AutoCapsule_disabled\\R\\testdata_ball_GaussNoise_HighLevel_ver2.csv")
    # select data
    Xdata = data.iloc[:,0:7] 
    Ydata = data.iloc[:,17]

    #find name of column
    nameofcolumnYdata = data.columns.values.tolist()
    
    #combind column with name column
    Xdata = pd.DataFrame(Xdata)
    Xdata[nameofcolumnYdata[17]] = Ydata
    XYdata = Xdata

    #filter
    XYdata = XYdata[XYdata.iloc[:, 7] <0.5]
    XYdata = XYdata[XYdata.iloc[:, 7] > 0]
    
#----------------- 各descriptorの名前(Symbol)作成 X1, X2 ,,,--------
    #number of columns
    number_X = len(XYdata.columns.values.tolist())
    # print(number_X)
    Xsymbol = []

    for ii in range(1,number_X+1):
        Xsymbol_each = "X"+str(ii)
        Xsymbol.append(Xsymbol_each)
    # print(Xsymbol)
    
    #combination
    nameofcolumnXYdata = XYdata.columns.values.tolist()
    dictOfNameXYdata = pd.DataFrame(nameofcolumnXYdata) 
    dictOfNameXYdata[''] = Xsymbol
    Xsymbol_and_DescriptorName = dictOfNameXYdata
    # print(Xsymbol_and_DescriptorName)

    #create column name
    Xsymbol_Xdata_original = XYdata
    for name in nameofcolumnXYdata:
        index = nameofcolumnXYdata.index(name)
        Xsymbol_Xdata_original = Xsymbol_Xdata_original.rename(columns ={name:Xsymbol[index]})
    # print(Xsymbol_Xdata_original)   
    return Xsymbol_Xdata_original
#========================= ここからしばらくはfuntion作成 =================

#===================== 積の項を作る ===============================

def make_product(Xsymbol_Xdata):
    #find name of variable
    X_symbol_name = Xsymbol_Xdata.columns.values.tolist()
    product_term_name = list(combinations(X_symbol_name, 2)) 
    product_term_name = np.array(product_term_name)
    # create matrix
    row = len(Xsymbol_Xdata)
    cloumn = 1
    product_term = [0] * row
    for x in range (row):
        product_term[x] = [0] * cloumn
    #for 0 - production_term_name

    count_row = len(product_term_name) 

    for nn in range(0,count_row+1):
        product = Xsymbol_Xdata.iloc[:, Xsymbol_Xdata.columns.get_loc(str(product_term_name[nn,0]))]*Xsymbol_Xdata.iloc[:, Xsymbol_Xdata.columns.get_loc(str(product_term_name[nn,1]))]
        #combination product_term and product
        product_term = np.insert(product_term, [nn+1], product, axis=1)
        print(product_term)

    print(product_term)
    #create name column
    nameofproduct_term = product_term.columns.values.tolist()
    product_term = product_term.rename(columns ={nameofproduct_term:"("+product_term_name[ :,0]+product_term_name[:, 1]+")"})
    return(product_term)
#======================== 商の項を作る ================================
def make_quotient(Xsymbol_Xdata):
   #find name of variable
    X_symbol_name = Xsymbol_Xdata.columns.values.tolist()
    quotient_term_name = permutations(X_symbol_name,2)
    quotient_term = []
        # For input 
    for i in range(Xsymbol_Xdata.shape[0]):          # A for loop for row entries 
        a =[] 
        # A for loop for column entries 
        a.append(0) 
        quotient_term.append(a) 
    
    #FIXME: For printing the matrix was created for check if evrythinhs is okay gonna delete it
    for i in range(Xsymbol_Xdata.shape[0]):  
        print(quotient_term[i][1], end = " ") 

    count_row = quotient_term_name.shape[0] 
    for nn in range(1,count_row+1):
        quotient = Xsymbol_Xdata[: ,quotient_term_name[nn,1]] / Xsymbol_Xdata[: ,quotient_term_name[nn,2]]
        #combination quotient_term and quotient
        dictOfquotient_term = pd.DataFrame(quotient_term) 
        dictOfquotient_term[''] = quotient
        quotient_term = dictOfquotient_term

    quotient_term = quotient_term[ :,-1]
    #create name column
    nameofQuotient_term = quotient_term.columns.values.tolist()
    quotient_term = quotient_term.rename(columns ={nameofQuotient_term:"("+ quotient_term_name[: ,1]+"/"+ quotient_term_name[:, 2]+")"})
    return(quotient_term)
#======================== 二乗の項を作る =================================
def make_square(Xsymbol_Xdata):
    X_symbol_name = Xsymbol_Xdata.columns.values.tolist()
    dictOfNameXsymbol = pd.DataFrame(X_symbol_name) 
    dictOfNameXsymbol[''] = Xsymbol_Xdata
    square_term_name = dictOfNameXsymbol
    square_term = [] 
    # For input 
    for i in range(Xsymbol_Xdata.shape[0]):          # A for loop for row entries 
        a =[] 
        # A for loop for column entries 
        a.append(0) 
        square_term.append(a) 
    
   #FIXME: For printing the matrix was created for check if evrythinhs is okay gonna delete it
    for i in range(Xsymbol_Xdata.shape[0]):  
        print(square_term[i][1], end = " ") 

    count_row = square_term_name.shape[0]
    for nn in range(1,count_row+1):
        square = Xsymbol_Xdata[ :,square_term_name[nn,1]]*Xsymbol_Xdata[: ,square_term_name[nn,2]]
        dictOfNamesquare_term = pd.DataFrame(square_term) 
        dictOfNamesquare_term[''] = square
        square_term = dictOfNamesquare_term

    square_term = square_term[ :,-1]
    #create name column
    nameofsquare_term = square_term.columns.values.tolist()
    square_term = product_term.rename(columns ={nameofsquare_term:"("+square_term_name[:,1]+square_term_name[:, 2]+")"})
    return(square_term)
#=========================== 指数の項を作る ================================
#あとで
#===========================　和の項を作る =================================
def make_sum(Xsymbol_Xdata):
    X_symbol_name = Xsymbol_Xdata.columns.values.tolist()
    sum_term_name = combinations(X_symbol_name,2)
    sum_term = [] 
    # For input 
    for i in range(Xsymbol_Xdata.shape[0]):          # A for loop for row entries 
        a =[] 
        # A for loop for column entries 
        a.append(0) 
        sum_term.append(a) 
    
    #FIXME: For printing the matrix was created for check if evrythinhs is okay gonna delete it
    for i in range(Xsymbol_Xdata.shape[0]):  
        print(sum_term[i][1], end = " ") 

    count_row = sum_term_name.shape[0]
    for nn in range(1,count_row+1):
        sum1 = Xsymbol_Xdata[: ,sum_term_name[nn,1]]+Xsymbol_Xdata[ :,sum_term_name[nn,2]]
        dictOfNamesum_term = pd.DataFrame(sum_term) 
        dictOfNamesum_term[''] = sum1
        square_term = dictOfNamesum_term
    sum_term = sum_term[: ,-1]
    #create name column
    nameofsum_term = sum_term.columns.values.tolist()
    sum_term = sum_term.rename(columns ={nameofsum_term:"("+sum_term_name[ :,1]+"+" +sum_term_name[:, 2]+")"})
    return(sum_term)
#=========================== 差の項を作る =================================
def make_difference(Xsymbol_Xdata):
    X_symbol_name = Xsymbol_Xdata.columns.values.tolist()
    difference_term_name = permutations(X_symbol_name,2)
    difference_term = []

    # For input 
    for i in range(Xsymbol_Xdata.shape[0]):          # A for loop for row entries 
        a =[] 
        # A for loop for column entries 
        a.append(0) 
        difference_term.append(a) 
    
    #FIXME: For printing the matrix was created for check if evrythinhs is okay gonna delete it
    for i in range(Xsymbol_Xdata.shape[0]):  
        print(difference_term[i][1], end = " ") 

    count_row = difference_term_name.shape[0]
    for nn in range(1,count_row+1):
        difference = Xsymbol_Xdata[ :,difference_term_name[nn,1]] - Xsymbol_Xdata[ :,difference_term_name[nn,2]]
        dictOfNamedifference_term = pd.DataFrame(difference_term) 
        dictOfNamedifference_term[''] = difference
        square_term = dictOfNamedifference_term

    difference_term = difference_term[ :,-1]
    #create name column
    nameofdifference_term = difference_term.columns.values.tolist()
    product_term = difference_term.rename(columns ={nameofdifference_term:"("+difference_term_name[: ,1]+"-"+ difference_term_name[:, 2]+ ")"})
    return(difference_term)
#============================ penalty factorを作る =======================
def make_penalty_factor(new_Xdata, Xsymbol_Xdata_original):
    X_symbol_name = Xsymbol_Xdata.columns.values.tolist()
    s = pd.Series(X_symbol_name)
    #FIXME: this function not sure 
    X_symbol_name_removed = s.str.extractall("X")
    penalty_factor =[]

    count_row = X_symbol_name_removed.shape[0]
    for ii in range(1,count_row+1):
        penalty_factor.append(len(X_symbol_name_removed[ii])) 
    #ある一定以上の次数のpenaltyを極端に大きくする。とりあえず、Xsymbol_Xdata_originalの項数+2を閾値としてデフォルトにしている。
    #the penalty factors areinternally rescaled to sum to nvars
    #つまり各penalty.factorは内部で規格化される、規格化されたpenalty.factor = （各penalty.factorの値 * 項の数 / 全部のpenalty.factorの値の和）
    #penalty_factor = replace(penalty_factor, (penalty_factor > ncol(Xsymbol_Xdata_original)+1),100)
    #penalty_factor = replace(penalty_factor, (penalty_factor <= ncol(Xsymbol_Xdata_original)+1),1)
    #もしくは、自分で指定する（今は4）
    #FIXME: this function not sure 
    penalty_factor = penalty_factor.replace(penalty_factor> 4,100)
    penalty_factor = penalty_factor.replace(penalty_factor<= 4,1)
    return(penalty_factor)
#========================== 標準化した回帰係数を計算するfunction =================
def make_std_coefficient(model, Xdata):
    #TODO: variable
    #interceptアリの時はこれも入れる
	return(cf_std)

def make_std_coefficient_intercept(model, Xdata):
    #TODO: variable
    #interceptアリの時はこれも入れる
	return(cf_std)
# normolize for python
def nom(X, x_min, x_max):
            nom = (X-Ydata.min(axis=0))*(x_max-x_min)
            denom = X.max(axis=0) - X.min(axis=0)
            denom[denom==0] = 1
            return x_min + nom/denom 
#========================================================================
#====== non-linear reccurent adaptive interpretable regression (NRAI regression) ======

#------------------ input parameter -------------------------------------

if __name__ == "__main__":
    #読み込むファイル名とかは上の方をいじってね.目的変数は決めておらず、フレキシブルに網羅的に試す
    X_Y = data_ready()
    #データ数をちょっと減らす　2000くらいにする
    X_Y = X_Y.iloc[0:5000,]
    #Xsymbol_Xdata_originalがある状態から始める
    #Xsymbol_Xdata = Xsymbol_Xdata_original

    nnn = 1 
    #ランダムにに実行する回数（nfoldを少しずつ変えながら）
    #pmax_vector = c(100, 80, 60, 40, 20, 10, 10) #pmaxの値は少しずつ減らしていく。ここで指定
    #pmax_vector = c(300, 300, 300, 300) #pmaxの値は少しずつ減らしていく。ここで指定。reccurentする回数分の要素がある
    rrr = X_Y.shape[1]  
    #reccurentする回数。ここでは、最初の説明変数の数に設定してみる
    #rrr = 3

    FirstLassoCutoff = 30 
    #最初の変数選択での変数の数の上限
    FinalCutoff = 5 
    #最後に出てくる式の変数の上限
    max_pmax = 1000

    object_record = []
    error_record = []

    n = 1
    m = 1000
    summary_data = [0] * n
    for x in range (n):
        summary_data[x] = [0] * m

    i_error = 10000  
    #誤差の初期値　大きめに取っておけば問題ない
    # #X_YのデータをXとYに分ける。それを網羅的にやる
    # #for (pp in 1:ncol(X_Y))　#全部試すときはこっち
    # #X5をYにする時だけ
    # print(X_Y)
    for pp in range(X_Y.shape[1]-1,X_Y.shape[1]+1): 
        #TODO:
        Ydata = X_Y.iloc[:,pp]
        Ydata = pd.DataFrame(Ydata)
        Ydata = nom(Ydata, 0, 1)    
        # in R use [,-pp] is mean except -pp
        Xsymbol_Xdata = X_Y.iloc[:,0:pp]
        Xsymbol_Xdata_original = Xsymbol_Xdata
        for nn in range(0,nnn+1):
            for ii in range(0,rrr+1):
                #積と商と和と差と二乗と指数の項を作成する
                product_term = make_product(Xsymbol_Xdata)
                print(product_term)
                #quotient_term = make_quotient(Xsymbol_Xdata)
                # print(quotient_term)
                #square_term = make_square(Xsymbol_Xdata)
                # print(square_term)
                #sum_term = make_sum(Xsymbol_Xdata)
                # print(sum_term)
                #difference_term = make_difference(Xsymbol_Xdata)
                # print(difference_term)
                #new_Xdata = []
                
    #             #TODO: newdata


    #             #NAやNaNやInfがたまに出てくるので、それが含まれる列を除去
    #             #TODO: cut NA NaN Inf

    #             #全く同じ項が選ばれることがあるので、その場合は片方を除去
    #             #TODO: delete duplicate

    #             #履歴をプリント
    #             #TODO: record history

    #             #penalty.factorを作る。次数（項の名前に入っているXの数）が大きい項は選ばれにくくなるようにする。
    #             #TODO: from 321 - 349

    #             #標準化係数でピックアップしたデータのみを使用して再度LASSO
    #             #cl = makeCluster(4)  #10コアで並列化開始
    #             #registerDoParallel(cl)  
    #             #cross validation
    #             #fitLassoCV1 = cv.glmnet(x = as.matrix(std_selected_term_data), y = as.matrix(Ydata), family ="gaussian", alpha = 1
    #                 #    	    ,nfolds = 10, parallel=TRUE, standardize = TRUE, penalty.factor=penalty_factor
    #                     #          ,thresh = 1E-7, maxit = 10^5, nlambda = 1000, intercept=FALSE, lambda.min.ratio = 0.00000001
    #             #          ,lambda = 2^(-40:5) )
    #             #stopCluster(cl)  #並列化おしまい

    #             #coefficient = coef(fitLassoCV1, s="lambda.min") #たくさん選ばれるように1seは使わない
    #             #coefficient = as.matrix(coefficient)
    #             #coefficient = coefficient[coefficient !=0,]
    #             #used_term_name = names(coefficient)
    #             #used_term_name = used_term_name[-which(used_term_name %in% "(Intercept)")] #intercept入れてるときはこれも入れる
    #             #used_term_data = as.matrix(std_selected_term_data[ ,used_term_name])
    #             #used_term_data = as.matrix(new_Xdata[ ,used_term_name])
    #             #colnames(used_term_data) = used_term_name ##

    #             #履歴をプリント
    #             #TODO: record history
    #             #TODO: from 371 - 460

    #             #一番cv誤差が小さかったモデルを記録
                
    # #一番良かったモデルの標準化係数を計算
    # std_coef_final = make_std_coefficient(best_model, best_model_x_selected)