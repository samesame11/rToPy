
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
library(BBmisc) #normalize���g������
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
  
  #filter���������͂�����������i�ُ�l���L�����肷��s�������j
  XYdata = filter(XYdata, Ydata<0.5)
  XYdata = filter(XYdata, Ydata>0)
  
  #----------------- �edescriptor�̖��O(Symbol)�쐬 X1, X2 ,,,--------
  
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


#========================= �������炵�΂炭��funtion�쐬 =================

#===================== �ς̍������ ===============================

make_product = function(Xsymbol_Xdata)
{
  
  X_symbol_name = colnames(Xsymbol_Xdata)
  
  product_term_name = combinations(n=ncol(Xsymbol_Xdata), r=2, X_symbol_name)

  
  
  product_term = matrix(0, ncol=1, nrow=nrow(Xsymbol_Xdata))
  
  for (nn in 1:nrow(product_term_name))
  {
    product = Xsymbol_Xdata[ ,product_term_name[nn,1]]*Xsymbol_Xdata[ ,product_term_name[nn,2]]
    
    product_term = cbind(product_term, product)
    print(product_term)
  
  }
  
  product_term = product_term[ ,-1]
  colnames(product_term)=paste("(", product_term_name[ ,1], product_term_name[, 2],")", sep ="")
  
  return(product_term)
  
}

#======================== ���̍������ ================================

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
#======================== ���̍������ =================================

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

#=========================== �w���̍������ ================================
#���Ƃ�
#===========================�@�a�̍������ =================================

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




#=========================== ���̍������ =================================

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

#============================ penalty factor����� =======================

make_penalty_factor = function(new_Xdata, Xsymbol_Xdata_original)
{
  X_symbol_name = colnames(new_Xdata)
  X_symbol_name_removed = str_extract_all(X_symbol_name, pattern="X", simplify=FALSE)
  
  penalty_factor = c()
  for (ii in 1:length(X_symbol_name_removed))
  {
    penalty_factor = c(penalty_factor, length(X_symbol_name_removed[[ii]]))
  }
  
  
  #������ȏ�̎�����penalty���ɒ[�ɑ傫������B�Ƃ肠�����AXsymbol_Xdata_original�̍���+2��臒l�Ƃ��ăf�t�H���g�ɂ��Ă���B
  #the penalty factors areinternally rescaled to sum to nvars
  #�܂�epenalty.factor�͓����ŋK�i�������A�K�i�����ꂽpenalty.factor = �i�epenalty.factor�̒l * ���̐� / �S����penalty.factor�̒l�̘a�j
  #penalty_factor = replace(penalty_factor, (penalty_factor > ncol(Xsymbol_Xdata_original)+1),100)
  #penalty_factor = replace(penalty_factor, (penalty_factor <= ncol(Xsymbol_Xdata_original)+1),1)
  #�������́A�����Ŏw�肷��i����4�j
  penalty_factor = replace(penalty_factor, (penalty_factor > 4),100)
  penalty_factor = replace(penalty_factor, (penalty_factor <= 4),1)
  
  return(penalty_factor)
  
}

#========================== �W����������A�W�����v�Z����function =================

make_std_coefficient = function(model, Xdata)
{
  cf = coef(model, s="lambda.min")[ ,1]
  sds = sapply(as.data.frame(Xdata), sd)
  mus = sapply(as.data.frame(Xdata), mean)
  cf_std <- cf * c(1, c(sds[names(cf)][-1]))
  #cf_std[1] <- cf[1] + sum(cf[-1] * mus[names(cf)][-1]) #intercept�A���̎��͂���������
  
  return(cf_std)
}

#========================== �W����������A�W�����v�Z����function =================

make_std_coefficient_intercept = function(model, Xdata)
{
  cf = coef(model, s="lambda.min")[ ,1]
  sds = sapply(as.data.frame(Xdata), sd)
  mus = sapply(as.data.frame(Xdata), mean)
  cf_std <- cf * c(1, c(sds[names(cf)][-1]))
  cf_std[1] <- cf[1] + sum(cf[-1] * mus[names(cf)][-1]) #intercept�A���̎��͂���������
  
  return(cf_std)
}

X_Y = data_ready()

#�f�[�^����������ƌ��炷�@2000���炢�ɂ���
X_Y = X_Y[1:5000, ]

#Xsymbol_Xdata_original�������Ԃ���n�߂�
#Xsymbol_Xdata = Xsymbol_Xdata_original

nnn = 1 #�����_���ɂɎ��s����񐔁infold���������ς��Ȃ���j
#pmax_vector = c(100, 80, 60, 40, 20, 10, 10) #pmax�̒l�͏��������炵�Ă����B�����Ŏw��
#pmax_vector = c(300, 300, 300, 300) #pmax�̒l�͏��������炵�Ă����B�����Ŏw��Breccurent����񐔕��̗v�f������

rrr = ncol(X_Y)  #reccurent����񐔁B�����ł́A�ŏ��̐����ϐ��̐��ɐݒ肵�Ă݂�
#rrr = 3

FirstLassoCutoff = 30 #�ŏ��̕ϐ��I���ł̕ϐ��̐��̏��
FinalCutoff = 5�@#�Ō�ɏo�Ă��鎮�̕ϐ��̏��

max_pmax = 1000

object_record = c()
error_record = c()
#term_record = matrix(0, nrow=1, ncol=max_pmax)
#coefficient_record = matrix(0, nrow=1, ncol=max_pmax)
summary_data = matrix(0, nrow=1, ncol=max_pmax)

i_error = 10000  #�덷�̏����l�@�傫�߂Ɏ���Ă����Ζ��Ȃ�
#X_Y�̃f�[�^��X��Y�ɕ�����B�����ԗ��I�ɂ��

Ydata = X_Y[ ,8,drop=FALSE]
Ydata = normalize(Ydata, method="range", range=c(0,1))
Xsymbol_Xdata = X_Y[, -8 , drop=FALSE]
Xsymbol_Xdata_original = Xsymbol_Xdata
  
product_term = make_product(Xsymbol_Xdata)
