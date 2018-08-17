# Public functions included:
#
# 1. extractTriplets_real(PATH_TDI,PATH_triplets)

# 2. extractTriplets_rand(PATH_TDI_rand,PATH_triplets_rand)

# 3. extractTriplets_singleSGA(PATH_TDI_singleSGA,PATH_triplets_singleSGA)

# 4. determinePosteriorThd_parallel(PATH_randTriplets_in,FILE_SGAMatrix_in,FILE_threshold_out, sig=0.05)

# 5. filterTriplets(FILEorPATH_in, FILE_threshold, FILEorPATH_out)

# 6. determinePosteriorThd (PATH_randTriplets_in,FILE_SGAMatrix_in,FILE_threshold_out,start=1, end=999, sig=0.05)

if 'os' not in dir():
    import os
if 'time' not in dir():
    import time
if 'pd' not in dir():
    import pandas as pd
if 'mp' not in dir():
    import multiprocessing as mp
if 'np' not in dir():
    import numpy as np
if 're' not in dir():
    import re

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


def getTopGtpostprob_parallel(PATH_TDI, FILE):
    df = pd.read_csv(PATH_TDI + FILE, index_col=0)
    tumor = FILE[:-4]  # remove '.csv'
    TSDPs = []
    for col in df:
        topIdx = df[col].idxmax()
        SGA = topIdx
        DEG = col
        postprob = df[col].loc[topIdx]
        TSDP = [tumor, SGA, DEG, postprob]
        TSDPs.append(TSDP)
    return TSDPs


def outputTriplets(triplets, PATHNAME_triplet):
    df_triplets = pd.DataFrame(data=triplets,
                               columns=['patient_name', 'cause_gene_name', 'result_gene_name', 'posterior'])
    df_triplets.to_csv(PATHNAME_triplet, index=False)

def extractTriplets_parallel(PATH_TDI,FILE_triplets):

    num_cores = mp.cpu_count()
    pool = mp.Pool()  # defalut is num_cores ==   pool = mp.Pool(num_cores)

    #    triplets = [pool.apply_async( getTopGtpostprob_parallel, args=(PATH_TDI,FILE) ) for FILE in os.listdir(PATH_TDI))]
    results = []
    for FILE in os.listdir(PATH_TDI):
        results.append(pool.apply_async(getTopGtpostprob_parallel, args=(PATH_TDI, FILE)))

    triplets = []
    for result in results:
        triplet = result.get()
        triplets += triplet

    outputTriplets(triplets, FILE_triplets)

def extractTriplets_real(path_in, path_out):
    my_timer = Timer()
    file_out = ""
    if path_in[-1] != "\\" and path_in[-1] != "/":
        path_in += "/"

    if path_out[-4:] == ".csv":
        file_out = path_out
    else:
        if path_out[-1] != "\\" and path_out[-1] != "/":
            path_out += "/"
        file_out = path_out+"triplet.csv"

    print "Extracting triplets from " + path_in + "..."
    extractTriplets_parallel(path_in,file_out)
    print "Triplet has been extracted successfully! "


    time_hhmmss = my_timer.get_time_hhmmss()
    print("Total time elapsed: %s\n" % time_hhmmss)

def extractTriplets_rand(path_in, path_out, start=1, end=999):
    my_timer = Timer()

    if path_in[-1] != "\\" and path_in[-1] != "/":
        path_in += "/"
    if path_out[-1] != "\\" and path_out[-1] != "/":
        path_out += "/"

    dirList = os.listdir(path_in)
 
    dirList_sorted = sorted(dirList, key=lambda x: (int(re.sub('\D', '', x)), x))
    count = 0
    for item in dirList_sorted:
        if os.path.isfile(path_in + item ): #given the specific rand path in path_in
            #get the rand number
            randnum = int(re.sub('\D', '', path_in[path_in.rfind('/')+1:]))
            print "Extracting triplets from " + path_in + "..."
            extractTriplets_parallel(path_in, path_out+"triplet.csv" + str(randnum))
            print "triplet.csv" + str(randnum) + " has been successfully extracted to " + path_out
            break
        else:
            count += 1
            if count < start:
                continue
            if count > end:
                break
            print "Extracting triplets from " + path_in +  item + "..."
            extractTriplets_parallel(path_in + item + "/", path_out + "triplet.rand." + str(count) +".csv")
            print "triplet.rand." + str(count) +".csv" + " has been successfully extracted to " +  path_out + "\n"

    time_hhmmss = my_timer.get_time_hhmmss()
    print("Total time elapsed: %s\n" % time_hhmmss)

def extractTriplets_singleSGA(path_in, path_out, SGAs=""):
    my_timer = Timer()

    if path_in[-1] != "\\" and path_in[-1] != "/":
        path_in += "/"
    if path_out[-1] != "\\" and path_out[-1] != "/":
        path_out += "/"

    if SGAs == '':
        for item in os.listdir(path_in):
            print "Extracting triplets from " + path_in + item + "..."
            extractTriplets_parallel(path_in + item + "/", path_out + "triplet." + item + ".csv")
            print item + ".csv" + " has been successfully extracted to " + path_out + "\n"
    else: #SGAs is a one or multiple SGAs separated by comma
        l_SGAs = SGAs.split(',')
        for item in l_SGAs:
            print "Extracting triplets from " + path_in + item + "..."
            extractTriplets_parallel(path_in + item + "/", path_out + "triplet." + item + ".csv")
            print item + ".csv" + " has been successfully extracted to " + path_out + "\n"

    time_hhmmss = my_timer.get_time_hhmmss()
    print("Total time elapsed: %s\n" % time_hhmmss)

def getHistCount_parallel(PATH_randTriplets_in,file,bins):
    if PATH_randTriplets_in[-1] != "\\" and PATH_randTriplets_in[-1] != "/":
        PATH_randTriplets_in += "/"
    df = pd.read_csv(PATH_randTriplets_in + file)
    dict_SGAhist = {}
    for SGA, group in df.groupby("cause_gene_name")['posterior']:
        count, bins_pp = np.histogram(group,bins)
        dict_SGAhist[SGA]=count
    df_hist = pd.DataFrame(dict_SGAhist)#rows are bins, columns are SGA
    return df_hist

def determinePosteriorThd_parallel(PATH_randTriplets_in,FILE_SGAMatrix_in,FILE_threshold_out, sig=0.05):
    my_timer = Timer()
    print "Processing ...."    
    bins = np.linspace(0, 1, num=1001)
    pool = mp.Pool()  # defalut is num_cores ==   pool = mp.Pool(num_cores)
    results = []
    for file in os.listdir(PATH_randTriplets_in):
        results.append(pool.apply_async(getHistCount_parallel, args=(PATH_randTriplets_in,file,bins)))
    
    df_hist_total = pd.DataFrame()
    for result in results:
        df_hist = result.get()
        df_hist_total = df_hist_total.add(df_hist, fill_value=0)

    # normalize the count of each column(SGA)
    df_norm = df_hist_total / df_hist_total.sum()
    # calculate the CDF of each column(SGA)
    df_CDF = df_norm.cumsum()
    # set index of df_CDF  to the smaller side of bin, so we choose bins[:-1]
    df_CDF.index = bins[:-1]
    # get the index of threshold for each column(SGA)
    #   this is a seriers with index of SGAs and value of bins
    SGA_thd = abs(df_CDF - (1 - sig)).idxmin()
    SGA_thd = SGA_thd.round(3)
    # Get the total SGA names of SGAMatrix
    #   read first line of SGAMatrix to get SGA name list
    if "csv" not in dir():
        import csv
    with open(FILE_SGAMatrix_in, 'r') as f:
        reader = csv.reader(f)
        SGAs = next(reader)[1:]  # line start from comma
    # convert SGAs list to pandas dataframe
    df_SGAthd = pd.DataFrame(SGAs, columns=['SGA'])
    df_SGAthd['PostProbcutoff'] = df_SGAthd['SGA'].map(SGA_thd).fillna(1)

    df_SGAthd.to_csv(FILE_threshold_out, index=False)
    print "Process finished successfully"
    time_hhmmss = my_timer.get_time_hhmmss()
    print("Total time elapsed: %s\n" % time_hhmmss)
    
def determinePosteriorThd(PATH_randTriplets_in,FILE_SGAMatrix_in,FILE_threshold_out,start=1, end=999, sig=0.05):
    my_timer = Timer()
    print "Processing ...."
    bins = np.linspace(0,1,num=1001)

    df_hist_total = pd.DataFrame()
    files = os.listdir(PATH_randTriplets_in)
    # sort file name by it's digital number
    files_sorted = sorted(files, key=lambda x: (int(re.sub('\D', '', x)), x))  # replace all non-digits to empty, then convert to int
    filecount = 0
    for file in files_sorted:
        filecount += 1
        if filecount < start:
            continue
        if filecount > end:
            break
        df = pd.read_csv(PATH_randTriplets_in + file)

        dict_SGAhist = {}
        for SGA, group in df.groupby("cause_gene_name")['posterior']:
            count, bins_pp = np.histogram(group,bins)
            dict_SGAhist[SGA]=count

        df_hist = pd.DataFrame(dict_SGAhist)#rows are bins, columns are SGA
        #add dataframe of each cell by column(SGA), fill the difference with 0
        df_hist_total = df_hist_total.add(df_hist, fill_value=0)

    #normalize the count of each column(SGA)
    df_norm = df_hist_total/df_hist_total.sum()
    #calculate the CDF of each column(SGA)
    df_CDF = df_norm.cumsum()
    #set index of df_CDF  to the smaller side of bin, so we choose bins[:-1]
    df_CDF.index = bins[:-1]
    #get the index of threshold for each column(SGA)
    #   this is a seriers with index of SGAs and value of bins
    SGA_thd = abs(df_CDF - (1-sig)).idxmin()
    SGA_thd = SGA_thd.round(3)
    #Get the total SGA names of SGAMatrix
    #   read first line of SGAMatrix to get SGA name list
    with open(FILE_SGAMatrix_in, 'r') as f:
        reader = csv.reader(f)
        SGAs = next(reader)[1:]  #line start from comma
    #convert SGAs list to pandas dataframe
    df_SGAthd = pd.DataFrame(SGAs,columns=['SGA'])
    df_SGAthd['PostProbcutoff']=df_SGAthd['SGA'].map(SGA_thd).fillna(1)

    df_SGAthd.to_csv(FILE_threshold_out,index=False)

    print "Process finished successfully"
    time_hhmmss = my_timer.get_time_hhmmss()
    print("Total time elapsed: %s\n" % time_hhmmss)

def filterTriplets(FILEorPATH_in, FILE_threshold, FILEorPATH_out ):
    # if os.path.isfile(FILEorPATH_in) != os.path.isfile(FILEorPATH_out) :
        # print "The file path Input and output are not consistant"
        # exit(1)
    if os.path.isfile(FILEorPATH_in):
        FILE_triplets = FILEorPATH_in
        FILE_fileterdTriplets = FILEorPATH_out
        df = pd.read_csv(FILE_triplets)
        df_posterior = pd.read_csv(FILE_threshold)
        tmp = dict(zip(df_posterior[df_posterior.columns[0]], df_posterior[df_posterior.columns[1]]))
        df['a'] = df.cause_gene_name.map(tmp)
        # cutoff by posterior
        df_updated_cutoff = df.loc[df['posterior'] >= df['a']]
        df_updated_cutoff.to_csv(FILE_fileterdTriplets)
        return df_updated_cutoff
    

# if __name__ == '__main__':
    # sig = 0.05
    # PATH_randTriplets_in = "/mnt/Data_2TB/TDIReults_Analysis/PanCancer13tts_20171217/ExtractInfo/triplets/100DEGfreqRands/"
    # FILE_SGAMatrix_in = "/home/xim33/DataSource/PanCancer13tts.SGAmatrix.4TCI.csv"
    # FILE_threshold_out = "/mnt/Data_2TB/TDIReults_Analysis/PanCancer13tts_20171217/ExtractInfo/postprobthd.p="+str(sig)+".csv"
    # start = 1  # the start random experiement
    # end = 3 # the end random experiement
    # determinePosteriorThd(PATH_randTriplets_in,FILE_SGAMatrix_in,FILE_threshold_out)


# PATH_TDI = "/home/xim33/TDIC_01_CPU_-fGtCM_SGAprior/OV/ori"
# PATH_triplets = "/home/xim33/Python_Projects/TDI_ResultsAnalysis/OV/triplet/real"
# extractTriplets_real(PATH_TDI,PATH_triplets)

# PATH_TDI = "/home/xim33/TDIC_01_CPU_-fGtCM_SGAprior/OV/rand"
# PATH_triplets = "/home/xim33/Python_Projects/TDI_ResultsAnalysis/OV/triplet/rand"
# extractTriplets_rand(PATH_TDI,PATH_triplets)

# PATH_TDI = "/home/xim33/TDIC_01_CPU_-fGtCM_SGAprior/OV/randbyfreq"
# PATH_triplets = "/home/xim33/Python_Projects/TDI_ResultsAnalysis/OV/triplet/randbyfreq"
# extractTriplets_rand(PATH_TDI,PATH_triplets)

# sig = 0.05
# PATH_randTriplets_in = "/home/xim33/Python_Projects/TDI_ResultsAnalysis/OV/triplet/randbyfreq"
# FILE_SGAMatrix_in = "/home/xim33/DataSource/PanCancer13tts.SGAmatrix.4TCI.1217_OV.csv"
# FILE_threshold_out = "/home/xim33/Python_Projects/TDI_ResultsAnalysis/OV/postprobthd.p="+str(sig)+".csv"
# determinePosteriorThd_parallel(PATH_randTriplets_in,FILE_SGAMatrix_in,FILE_threshold_out,sig)
    
FILEorPATH_in = "/home/xim33/Python_Projects/TDI_ResultsAnalysis/OV/Triplet_OV.csv"
FILE_threshold = "/home/xim33/Python_Projects/TDI_ResultsAnalysis/OV/postprobthd.p=0.05.csv"
FILEorPATH_out = "/home/xim33/Python_Projects/TDI_ResultsAnalysis/OV/FilteredTriplet_OV.csv"
filterTriplets(FILEorPATH_in, FILE_threshold, FILEorPATH_out)
