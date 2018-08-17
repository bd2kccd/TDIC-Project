import pandas as pd
import sys

def initGlobalDriver(file_GtM, file_GeM, file_dict, file_globalDriver):

    df_sga = pd.read_csv(file_GtM, index_col=0)
    df_sga.insert(0,'A0',1)
    df_gd = pd.read_csv(file_GeM, index_col=0)
    df_dict = pd.read_csv(file_dict, header = None, names=('a','b'))

    #convert df_dict to a dictionary
    d = dict(zip(df_dict['a'], df_dict['b']))

    # #create an empty dataframe with df_sga.index and d.keys columns
    # df_gd = pd.DataFrame(index = df_sga.index, columns=d.keys())
    # in order to get the same deg sequence of DEGmatrix, the above mathod cannot be used

    #iterate d, and replace df_gd column with df_sga column
    for key, value in d.iteritems():
        df_gd[key] = df_sga[value]

    #write back df_gd to globaldriver file
    df_gd.to_csv(file_globalDriver)

# file_GtM = sys.argv[1]
# file_GeM = sys.argv[2]
# file_dict = sys.argv[3]
# file_globalDriver = sys.argv[4]

# file_GtM ="C:/Users/XIM33/Documents/NetBeansProjects/DataSource/tGtM.csv"
# file_GeM ="C:/Users/XIM33/Documents/NetBeansProjects/DataSource/tGeM.csv"
# file_dict = "C:/Users/XIM33/Documents/NetBeansProjects/DataSource/tdictNew.csv"
# file_globalDriver = "C:/Users/XIM33/Documents/TDIC_iteration/DataSource/tGdM.csv"
# initGlobalDriver(file_GtM, file_GeM, file_dict, file_globalDriver)
