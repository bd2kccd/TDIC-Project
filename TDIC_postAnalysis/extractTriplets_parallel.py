if 'os' not in dir():
    import os
if 'time' not in dir():
    import time
if 'pd' not in dir():
    import pandas as pd
if 'mp' not in dir():
    import multiprocessing as mp


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


def extractTriplets_parallel(PATH_TDI, FILE_triplets):
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

    if path_out[-4:] == ".csv":
        file_out = path_out
        extractTriplets_parallel(path_in, file_out)
        print "Triplets " + file_out + " has been successfully extracted ! \n"
    else:

        if path_in[-1] != "\\" and path_in[-1] != "/":
            path_in += "/"
        if path_out[-1] != "\\" and path_out[-1] != "/":
            path_out += "/"

        print "Extracting triplets from " + path_in + "..."
        extractTriplets_parallel(path_in, path_out + "triplets.csv")
        print "Triplets has been successfully extracted ! \n"

    time_hhmmss = my_timer.get_time_hhmmss()
    print("Total time elapsed: %s\n" % time_hhmmss)


def extractTriplets_rand(path_in, path_out, start=1, end=999):
    my_timer = Timer()
    if os.path.isfile(path_out):
        file_out = path_out
        extractTriplets_parallel(path_in, file_out)
    else:
        if path_in[-1] != "\\" and path_in[-1] != "/":
            path_in += "/"
        if path_out[-1] != "\\" and path_out[-1] != "/":
            path_out += "/"

        dirList = os.listdir(path_in)
        if 're' not in dir():
            import re
        # dirList_sorted = sorted(dirList, key=lambda x: (int(re.sub('\D', '', x)), x))

        for item in dirList:
            if os.path.isfile(path_in + item):  # given the specific rand path in path_in
                if os.path.isfile(path_out):
                    extractTriplets_parallel(path_in, path_out)
                    break
                else:
                    print "Please specify the output triplets file name instead of path"
                    # get the rand number
                    randnum = int(re.sub('\D', '', path_in[path_in.rfind('/') + 1:]))
                    print "Extracting triplets from " + path_in + "..."
                    extractTriplets_parallel(path_in, path_out + "triplet.csv" + str(randnum))
                    print "triplet.csv" + str(randnum) + " has been successfully extracted to " + path_out
                    break
            else:
                randnuminItem = int(re.sub('\D', '',item))
                if randnuminItem in range(start, end+1):
                    print "Extracting triplets from " + path_in + item + "..."
                    extractTriplets_parallel(path_in + item + "/", path_out + "triplet.rand." + str(randnuminItem) + ".csv")
                    print "triplet.rand." + str(randnuminItem) + ".csv" + " has been successfully extracted to " + path_out + "\n"
                else:
                    continue
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
    else:  # SGAs is a one or multiple SGAs separated by comma
        l_SGAs = SGAs.split(',')
        for item in l_SGAs:
            print "Extracting triplets from " + path_in + item + "..."
            extractTriplets_parallel(path_in + item + "/", path_out + "triplet." + item + ".csv")
            print item + ".csv" + " has been successfully extracted to " + path_out + "\n"

    time_hhmmss = my_timer.get_time_hhmmss()
    print("Total time elapsed: %s\n" % time_hhmmss)


######################################################################################################
# Main program
#

# if __name__ == "__main__":
#     PATH_TDIrandbyfreq = "C:/Users/XIM33/Documents/NetBeansProjects/TDIC_01_CPU_Tumors_SGAprior/output_rand"
#     PATH_TripletsRandbyfreq = "C:/Users/XIM33/Documents/NetBeansProjects/TDIC_01_CPU_Tumors_SGAprior/triplets"
#     extractTriplets_rand(PATH_TDIrandbyfreq, PATH_TripletsRandbyfreq, 4, 5)
# from postAnalysisDataPreparision import *
#PATH_TDIrandbyfreq = "/home/xim33/TDIC_01_CPU_-fGtCM_SGAprior/pairexp_randbyfreq/"
#PATH_TripletsRandbyfreq = "/home/xim33/DataSource/PanCancer13tts_pairExp_extraction/triplets_randbyfreq"
# PATH_TDIrand = "/home/xim33/TDIC_01_CPU_-fGtCM_SGAprior/pairexp_rand/"
# PATH_TripletsRand = "/home/xim33/DataSource/PanCancer13tts_pairExp_extraction/"

#extractTriplets_rand(PATH_TDIrandbyfreq, PATH_TripletsRandbyfreq, 4, 5)
PATH_TDI = "/home/xim33/TDIC_01_CPU_-fGtCM_SGAprior/output/"
PATH_triplets = "/home/xim33/TDIC_01_CPU_-fGtCM_SGAprior/triplet/triplet.csv"
extractTriplets_real(PATH_TDI,PATH_triplets)