#======================== library ==============================
import pandas as pd
from itertools import combinations 
from itertools import permutations 
import numpy as np

import glmnet_python
from glmnet import glmnet

import scipy, importlib, pprint, matplotlib.pyplot as plt, warnings
from glmnet import glmnet; from glmnetPlot import glmnetPlot
from glmnetPrint import glmnetPrint; from glmnetCoef import glmnetCoef; from glmnetPredict import glmnetPredict
from cvglmnet import cvglmnet; from cvglmnetCoef import cvglmnetCoef
from cvglmnetPlot import cvglmnetPlot; from cvglmnetPredict import cvglmnetPredict
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

    #combination [maxtrix table]
    nameofcolumnXYdata = XYdata.columns.values.tolist()
    dictOfNameXYdata = pd.DataFrame(nameofcolumnXYdata) 
    dictOfNameXYdata[''] = Xsymbol
    Xsymbol_and_DescriptorName = dictOfNameXYdata

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

    count_row = len(product_term_name) 

    for nn in range(0, count_row):
        product = Xsymbol_Xdata.iloc[:, Xsymbol_Xdata.columns.get_loc(str(product_term_name[nn,0]))]*Xsymbol_Xdata.iloc[:, Xsymbol_Xdata.columns.get_loc(str(product_term_name[nn,1]))]
        #combination product_term and product [array]
        a = np.array(product_term)
        b = np.array(product)
        product_term = np.concatenate((a,b[:,None]),axis=1)
    
    product_term = product_term[ : ,1:]
    #create name column
    names = [] 
    reshaped = product_term.reshape(len(Xsymbol_Xdata),count_row)
    for n in range(0, count_row):
        name = "".join(str(x) for x in product_term_name[n])
        names.append("("+name+")")
    product_term =  pd.DataFrame(reshaped, columns=names)
    return(product_term)
#======================== 商の項を作る ================================
def make_quotient(Xsymbol_Xdata):
   #find name of variable
    X_symbol_name = Xsymbol_Xdata.columns.values.tolist()
    quotient_term_name = list(permutations(X_symbol_name, 2)) 
    quotient_term_name = np.array(quotient_term_name)
    # create matrix
    row = len(Xsymbol_Xdata)
    cloumn = 1
    quotient_term = [0] * row
    for x in range (row):
        quotient_term[x] = [0] * cloumn
    
    count_row = len(quotient_term_name) 
     
    for nn in range(0, count_row):
        quotient = Xsymbol_Xdata.iloc[:, Xsymbol_Xdata.columns.get_loc(str(quotient_term_name[nn,0]))]/Xsymbol_Xdata.iloc[:, Xsymbol_Xdata.columns.get_loc(str(quotient_term_name[nn,1]))]
        #combination product_term and product
        a = np.array(quotient_term)
        b = np.array(quotient)
        quotient_term = np.concatenate((a,b[:,None]),axis=1)
    

    quotient_term = quotient_term[ : ,1:]
    #create name column
    names = [] 
    reshaped = quotient_term.reshape(len(Xsymbol_Xdata),count_row)
    for n in range(0, count_row):
        name = "/".join(str(x) for x in quotient_term_name[n])
        names.append("("+name+")")
    quotient_term =  pd.DataFrame(reshaped, columns=names)
    return(quotient_term)
#======================== 二乗の項を作る =================================
def make_square(Xsymbol_Xdata):
    X_symbol_name = Xsymbol_Xdata.columns.values.tolist()
    #combination
    dictOfNameX_symbol_name = pd.DataFrame(X_symbol_name) 
    dictOfNameX_symbol_name[''] = X_symbol_name
    square_term_name = np.array(dictOfNameX_symbol_name)

    # create matrix
    row = len(Xsymbol_Xdata)
    cloumn = 1
    square_term = [0] * row
    for x in range (row):
        square_term[x] = [0] * cloumn
    
    count_row = len(square_term_name) 
    
    for nn in range(0, count_row):
        square = Xsymbol_Xdata.iloc[:, Xsymbol_Xdata.columns.get_loc(str(square_term_name[nn,0]))]*Xsymbol_Xdata.iloc[:, Xsymbol_Xdata.columns.get_loc(str(square_term_name[nn,1]))]
        #combination product_term and product
        a = np.array(square_term)
        b = np.array(square)
        square_term = np.concatenate((a,b[:,None]),axis=1)

    square_term = square_term[ : ,1:]
    #create name column
    names = [] 
    reshaped = square_term.reshape(len(Xsymbol_Xdata),count_row)
    for n in range(0, count_row):
        name = "".join(str(x) for x in square_term_name[n])
        names.append("("+name+")")
    square_term =  pd.DataFrame(reshaped, columns=names)
    return(square_term)
#=========================== 指数の項を作る ================================
#あとで
#===========================　和の項を作る =================================
def make_sum(Xsymbol_Xdata):
    #find name of variable
    X_symbol_name = Xsymbol_Xdata.columns.values.tolist()
    sum_term_name = list(combinations(X_symbol_name, 2)) 
    sum_term_name = np.array(sum_term_name)
    # create matrix
    row = len(Xsymbol_Xdata)
    cloumn = 1
    sum_term = [0] * row
    for x in range (row):
        sum_term[x] = [0] * cloumn

    count_row = len(sum_term_name) 

    for nn in range(0, count_row):
        sum1 = Xsymbol_Xdata.iloc[:, Xsymbol_Xdata.columns.get_loc(str(sum_term_name[nn,0]))]+Xsymbol_Xdata.iloc[:, Xsymbol_Xdata.columns.get_loc(str(sum_term_name[nn,1]))]
        #combination product_term and product
        a = np.array(sum_term)
        b = np.array(sum1)
        sum_term = np.concatenate((a,b[:,None]),axis=1)
    
    sum_term = sum_term[ : ,1:]
    #create name column
    names = [] 
    reshaped = sum_term.reshape(len(Xsymbol_Xdata),count_row)
    for n in range(0, count_row):
        name = "+".join(str(x) for x in sum_term_name[n])
        names.append("("+name+")")
    sum_term =  pd.DataFrame(reshaped, columns=names)
    return(sum_term)
#=========================== 差の項を作る =================================
def make_difference(Xsymbol_Xdata):
    #find name of variable
    X_symbol_name = Xsymbol_Xdata.columns.values.tolist()
    difference_term_name = list(permutations(X_symbol_name, 2)) 
    difference_term_name = np.array(difference_term_name)
    # create matrix
    row = len(Xsymbol_Xdata)
    cloumn = 1
    difference_term = [0] * row
    for x in range (row):
        difference_term[x] = [0] * cloumn

    count_row = len(difference_term_name) 

    for nn in range(0, count_row):
        difference = Xsymbol_Xdata.iloc[:, Xsymbol_Xdata.columns.get_loc(str(difference_term_name[nn,0]))]-Xsymbol_Xdata.iloc[:, Xsymbol_Xdata.columns.get_loc(str(difference_term_name[nn,1]))]
        #combination product_term and product
        a = np.array(difference_term)
        b = np.array(difference)
        difference_term = np.concatenate((a,b[:,None]),axis=1)
    
    difference_term = difference_term[ : ,1:]
    #create name column
    names = [] 
    reshaped = difference_term.reshape(len(Xsymbol_Xdata),count_row)
    for n in range(0, count_row):
        name = "-".join(str(x) for x in difference_term_name[n])
        names.append("("+name+")")
    difference_term =  pd.DataFrame(reshaped, columns=names)
    return(difference_term)
#============================ penalty factorを作る =======================
def make_penalty_factor(new_Xdata, Xsymbol_Xdata_original):
    X_symbol_name = new_Xdata.columns.values.tolist()
    penalty_factor =[]

    penalty_factor = np.char.count(X_symbol_name, sub ='X')  

    #ある一定以上の次数のpenaltyを極端に大きくする。とりあえず、Xsymbol_Xdata_originalの項数+2を閾値としてデフォルトにしている。
    #the penalty factors areinternally rescaled to sum to nvars
    #つまり各penalty.factorは内部で規格化される、規格化されたpenalty.factor = （各penalty.factorの値 * 項の数 / 全部のpenalty.factorの値の和）
    #penalty_factor = replace(penalty_factor, (penalty_factor > ncol(Xsymbol_Xdata_original)+1),100)
    #penalty_factor = replace(penalty_factor, (penalty_factor <= ncol(Xsymbol_Xdata_original)+1),1)
    #もしくは、自分で指定する（今は4）
    penalty_factor = np.array(penalty_factor)
    penalty_factor = np.where(penalty_factor > 4, 100, penalty_factor) 
    penalty_factor = np.where(penalty_factor <= 4, 1, penalty_factor) 
   
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
def nb(y=1930):
    debut = 1816
    MatDFemale = matrix(D . Female, nrow=111)
    colnames(MatDFemale) .set(range((debut + 0), 198))
    cly = range((y - debut + 1), 111)
    deces = diag(MatDFemale[:, cly[set(cly) & set(range(1, 199))]])
    return tuple(B . Female[B . Year == y], deces)
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
    for pp in range(X_Y.shape[1]-1,X_Y.shape[1]): 
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
                quotient_term = make_quotient(Xsymbol_Xdata)
                square_term = make_square(Xsymbol_Xdata)
                sum_term = make_sum(Xsymbol_Xdata)
                difference_term = make_difference(Xsymbol_Xdata)
                
                #name of each columns
                nameofXsymbol_Xdata = Xsymbol_Xdata.columns.values.tolist()
                lenghtOfXsymbol_Xdata = len(Xsymbol_Xdata.columns.values.tolist())
                nameofproduct_term = product_term.columns.values.tolist()
                nameofquotient_term = quotient_term.columns.values.tolist()
                nameofsquare_term = square_term.columns.values.tolist()
                nameofsum_term = sum_term.columns.values.tolist()
                nameofdifference_term = difference_term.columns.values.tolist()

                # table to array
                arrayXsymbol_Xdata = np.array(Xsymbol_Xdata)
                product_term = np.array(product_term)
                quotient_term = np.array(quotient_term)
                square_term = np.array(square_term)
                sum_term = np.array(sum_term)
                difference_term = np.array(difference_term)

                # add product_term
                for each in nameofproduct_term:
                    nameofXsymbol_Xdata.append(each)

                nametotal = nameofXsymbol_Xdata
                new_Xdata = np.concatenate((arrayXsymbol_Xdata,product_term),axis=1)
                reshaped = new_Xdata.reshape(len(arrayXsymbol_Xdata),len(nametotal))
                new_Xdata =  pd.DataFrame(reshaped, columns=nametotal)

                # 1st add quotient
                for each in nameofquotient_term:
                    nameofXsymbol_Xdata.append(each)

                nametotal = nameofXsymbol_Xdata
                new_Xdata = np.concatenate((new_Xdata,quotient_term),axis=1)
                reshaped = new_Xdata.reshape(len(product_term),len(nametotal))
                new_Xdata =  pd.DataFrame(reshaped, columns=nametotal)

                #2nd add square
                for each in nameofsquare_term:
                    nameofXsymbol_Xdata.append(each)

                nametotal = nameofXsymbol_Xdata
                new_Xdata = np.concatenate((new_Xdata,square_term),axis=1)
                reshaped =  new_Xdata.reshape(len(new_Xdata),len(nametotal))
                new_Xdata =  pd.DataFrame(reshaped, columns=nametotal)

                #3rd add square
                for each in nameofsum_term:
                    nameofXsymbol_Xdata.append(each)

                nametotal = nameofXsymbol_Xdata
                new_Xdata = np.concatenate((new_Xdata,sum_term),axis=1)
                reshaped =  new_Xdata.reshape(len(new_Xdata),len(nametotal))
                new_Xdata =  pd.DataFrame(reshaped, columns=nametotal) 

                #4th add difference_term
                for each in nameofdifference_term:
                    nameofXsymbol_Xdata.append(each)

                nametotal = nameofXsymbol_Xdata
                new_Xdata = np.concatenate((new_Xdata,difference_term),axis=1)
                reshaped =  new_Xdata.reshape(len(new_Xdata),len(nametotal))
                new_Xdata =  pd.DataFrame(reshaped, columns=nametotal) 
              
                #find string and delete column
                new_Xdata.drop(columns=new_Xdata.columns[(new_Xdata == 'Inf').any()])
                new_Xdata.drop(columns=new_Xdata.columns[(new_Xdata == 'NA').any()])
                new_Xdata.drop(columns=new_Xdata.columns[(new_Xdata == 'NaN').any()])

                #TODO: for unique

                moji = "pp "+str(pp+1)+" nn "+str(nn)+" ii "+str(ii)
                print(moji)
                moji = "number of Xsymbol_Xdata "+ str(lenghtOfXsymbol_Xdata)
                print(moji)
                moji = "number of new_Xdata "+ str(len(nametotal))
                print(moji)

                penalty_factor = make_penalty_factor(new_Xdata, Xsymbol_Xdata_original)

                #標準化係数でピックアップしたデータのみを使用して再度LASSO
               
                #cross validation
                fitLassoCV1 = glmnet(x = np.array(new_Xdata), y = np.array(Ydata), family = 'gaussian', 
                                     alpha = 1, nlambda = 1000)
                glmnetPrint(fit)

