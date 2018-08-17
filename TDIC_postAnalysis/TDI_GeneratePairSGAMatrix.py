import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

class Timer:
    def __init__(self):
        self.start = time.time()

    def restart(self):
        self.start = time.time()

    def get_time_hhmmss(self):
        end = time.time()
        m, s = divmod(end - self.start, 60)
        h, m = divmod(m, 60)
        time_str = "%02d:%02d:%02d" % (h, m, s)
        return time_str
		
def generatePairSGAs(FILE_SGA, filterthd,  FILE_SGAwPairs, FILE_SGAcountHist, FIlE_SGApairInfo):
    #df = pd.read_csv("C:/Users/XIM33/Documents/NetBeansProjects/DataSource/PanCancer13tts.SGAmatrix.4TCI.csv", index_col = 0)
    df = pd.read_csv(FILE_SGA, index_col = 0)
    # df = df.astype(int)
    # add SGA pairs to df_append dataframe
    cols = list(df)
    df_append=pd.DataFrame(index = df.index)
    for i in range (0,len(cols)-1):
        for j in range(i+1,len(cols)):
            df_append[cols[i]+"-"+cols[j]] = df[cols[i]] * df[cols[j]]

    #fileter the pair by sum >= 10
    cols_append = list(df_append)
    for col in cols_append:
         if df_append[col].sum() < filterthd :
             del df_append[col]

    print "There are "  + str(len(list(df_append))) + " pairs of SGAs."


    df_all = pd.concat([df,df_append], axis = 1)
    df_all.to_csv(FILE_SGAwPairs)
    count_ori = df.sum()
    count_all = df_all.sum()

    with open(FIlE_SGApairInfo, 'w') as f:
        f.write("There are " + str(len(list(df_append))) + " pairs of SGAs.\n")
    with open(FIlE_SGApairInfo, 'a') as f:
        count_all.to_csv(f)


    bins = np.linspace(min(count_ori.min(),count_all.min()), count_ori.max(), 1001)
    plt.hist(count_ori, bins, alpha=0.5, label='ori')
    plt.hist(count_all, bins, alpha=0.5, label='all')
    plt.savefig(FILE_SGAcountHist)
    plt.show()

FILE_SGA = "/home/xim33/DataSource/PanCancer13tts.SGAmatrix.4TCI.1217.csv"
FILE_SGAwPairs = "/home/xim33/Python_Projects/TDI_ResultsAnalysis/output/tGtMwPairSGAs.csv"
FILE_SGAcountHist = "/home/xim33/Python_Projects/TDI_ResultsAnalysis/output/SGAcountHist.png"
FIlE_SGApairInfo = "/home/xim33/Python_Projects/TDI_ResultsAnalysis/output/SGApairInfo.csv"
filterthd = 10
generatePairSGAs(FILE_SGA, filterthd,  FILE_SGAwPairs, FILE_SGAcountHist, FIlE_SGApairInfo)
my_timer = Timer()
generatePairSGAs(FILE_SGA, filterthd,  FILE_SGAwPairs, FILE_SGAcountHist, FIlE_SGApairInfo)
time_hhmmss = my_timer.get_time_hhmmss()
with open(FIlE_SGApairInfo, 'a') as f:
    f.write("Total time elapsed: %s\n" % time_hhmmss)