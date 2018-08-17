#from CommonFuctions import *
import os
import pandas as pd
import numpy as np
import multiprocessing as mp


def getHistCount(SGA,group,bins,dict_SGAhist):

    count, bins_pp = np.histogram(group, bins)
    dict_SGAhist[SGA] = count

def determineCutoffThd_parallel(PATH_randTriplets_in, FILE_SGAMatrix_in, FILE_threshold_out, start, end):
    df_hist_total = pd.DataFrame()
    bins = np.linspace(0, 1, num=11)
    count = 0
    for file in os.listdir(PATH_randTriplets_in):
        count += 1
        if count < start:
            continue
        if count > end:
            break
        df = pd.read_csv(PATH_randTriplets_in + file)

        manager = mp.Manager()
        dict_SGAhist = manager.dict()
        jobs = []
        for SGA, group in df.groupby("cause_gene_name")['posterior']:
            p = mp.Process(target=getHistCount, args=(SGA, group, bins, dict_SGAhist))
            jobs.append(p)
            p.start()
        for proc in jobs:
            proc.join()

        df_hist = pd.DataFrame(dict_SGAhist)  # rows are bins, columns are SGA
        # add dataframe of each cell by column(SGA), fill the difference with 0
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


if __name__ == '__main__':

    PATH_randTriplets_in = "C:/Users/XIM33/Documents/NetBeansProjects/TDIC_01_CPU_-fGtCM_SGAprior/randTriplets/"
    FILE_SGAMatrix_in = "C:/Users/XIM33/Documents/NetBeansProjects/DataSource/tGtM.csv"
    FILE_threshold_out = "C:/Users/XIM33/Documents/NetBeansProjects/TDIC_01_CPU_-fGtCM_SGAprior/output/threshold.csv"
    start = 1  #start random experiement
    end = 3 #end random experiement
    sig = 0.05
    determineCutoffThd_parallel(PATH_randTriplets_in, FILE_SGAMatrix_in, FILE_threshold_out, start, end)

