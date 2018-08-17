if 'os' not in dir():
    import os
if 'pd' not in dir():
    import pandas as pd
if 'np' not in dir():
    import numpy as np
if 'mp' not in dir():
    import multiprocessing as mp
if 'time' not in dir():
    import time


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


def getHistCount_parallel(PATH_randTriplets_in, file, bins):
    if PATH_randTriplets_in[-1] != "\\" and PATH_randTriplets_in[-1] != "/":
        PATH_randTriplets_in += "/"
    df = pd.read_csv(PATH_randTriplets_in + file)
    dict_SGAhist = {}
    for SGA, group in df.groupby("cause_gene_name")['posterior']:
        count, bins_pp = np.histogram(group, bins)
        dict_SGAhist[SGA] = count
    df_hist = pd.DataFrame(dict_SGAhist)  # rows are bins, columns are SGA
    return df_hist


def determinePosteriorThd_parallel(PATH_randTriplets_in, FILE_SGAMatrix_in, FILE_threshold_out, sig=0.05):
    my_timer = Timer()
    print "Processing ...."
    bins = np.linspace(0, 1, num=1001)
    pool = mp.Pool()  # defalut is num_cores ==   pool = mp.Pool(num_cores)
    results = []
    for file in os.listdir(PATH_randTriplets_in):
        results.append(pool.apply_async(getHistCount_parallel, args=(PATH_randTriplets_in, file, bins)))

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


# if __name__ == '__main__':
#     sig = 0.05
#     PATH_randTriplets = "/mnt/Data_2TB/TDIReults_Analysis/PanCancer13tts_20171217/ExtractInfo/triplets/100DEGfreqRands/"
#     FILE_SGAMatrix = "/home/xim33/DataSource/PanCancer13tts.SGAmatrix.4TCI.csv"
#     FILE_threshold = "/mnt/Data_2TB/TDIReults_Analysis/PanCancer13tts_20171217/ExtractInfo/postprobthd.p=0.05.csv"
#
#     determinePosteriorThd_parallel(PATH_randTriplets, FILE_SGAMatrix, FILE_threshold)

PATH_TripletsRandbyfreq = "/home/xim33/DataSource/PanCancer13tts_pairExp_extraction/triplets_randbyfreq/"
FILE_SGA = "/home/xim33/DataSource/pairExperiement/PanCancer13tts.SGAwPairMatrix.csv"
sig = 0.05
FILE_Threshold = "/home/xim33/DataSource/PanCancer13tts_pairExp_extraction/postprobthd.p=" + str(sig) + ".csv"
determinePosteriorThd_parallel(PATH_TripletsRandbyfreq,FILE_SGA,FILE_Threshold, sig)
