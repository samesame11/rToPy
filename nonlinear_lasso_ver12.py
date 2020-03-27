#======================== library ==============================
import csv 

#======================= read data ============================
def data_ready():    
    with open('C:\\ICF_AutoCapsule_disabled\\R\\test.csv','r') as file:
        data = csv.reader(file)
        for each_data in data:
            Xdata = []
            Ydata = []
            XYdata = [] 
            #set Xdata
            Xdata = each_data[0:3]
            #set Ydata
            Ydata = each_data[2]
            #combine Xdata and Ydata into XYdata
            Xdata.append(Ydata)
            XYdata = Xdata
            print(XYdata)
 #Todo : filter --- python read every line in files    
        #count number of column    
        number_X = len(XYdata)
        print(number_X)
        Xsymbol =[]
        #set variable of Xsymbol
        for ii in range(1,number_X+1):
            Xsymbol_each = "X"+str(ii)
            Xsymbol.append(Xsymbol_each)
        print(Xsymbol)  
        





data_ready()