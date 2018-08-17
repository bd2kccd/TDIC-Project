import os
import time
import pandas as pd



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

def getTopGtpostprob(tumor,df):
    TSDPs = []
    for col in df:
        topIdx = df[col].idxmax()
        SGA = topIdx
        DEG = col
        postprob = df[col].loc[topIdx]
        TSDP = [tumor, SGA, DEG, postprob]
        TSDPs.append(TSDP)
    # df_triplet = pd.DataFrame(SDPs,columns=['SGA','DEG','postprob'])
    #df_triplet.insert(0,"Tumor",tumor)
    return TSDPs

def extractTriplets(PATH_TDI,PATHNAME_triplet,fstart,fend):
    triplets = []
    fcount = 0
    for FILE in os.listdir(PATH_TDI):
        fcount += 1
        if fcount < fstart:
            continue
        if fcount > fend:
            break

        df_postprob = pd.read_csv(PATH_TDI+FILE, index_col=0)
        tumor = FILE[:-4]
        triplet = getTopGtpostprob(tumor,df_postprob)
        triplets += triplet
    return triplets


def outputTriplets(triplets,FILE_triplet):
    df_triplets = pd.DataFrame(data=triplets, columns=['Tumor', 'SGA', 'DEG', 'postprob'])
    df_triplets.to_csv(FILE_triplet,index=False)




######################################################################################################
# Main program
#
PATH_TDI = "C:/Users/XIM33/Documents/NetBeansProjects/TDIC_01_CPU_Tumors_SGAprior/output/"
PATHNAME_triplet = "C:/Users/XIM33/Documents/NetBeansProjects/TDIC_01_CPU_Tumors_SGAprior/triplets.csv"

try:
    fstart = int(raw_input('Enter the number of the start file: '))
except:
    fstart = 0
try:
    fend = int(raw_input('Enter the number of the end file: '))
except:
    fend = 99999

my_timer = Timer()

triplets = extractTriplets(PATH_TDI,PATHNAME_triplet,fstart,fend)
outputTriplets(triplets, PATHNAME_triplet)
time_hhmmss = my_timer.get_time_hhmmss()
print("Time elapsed: %s\n" % time_hhmmss)
#triplets = extractTriplets(PATH_TDI,PATHNAME_triplet,fstart,fend)
# if __name__ == "__main__":
#     PATH_TDI = "C:/Users/XIM33/Documents/NetBeansProjects/TDIC_01_CPU_Tumors_SGAprior/output/"
#     PATHNAME_triplet = "C:/Users/XIM33/Documents/NetBeansProjects/TDIC_01_CPU_Tumors_SGAprior/triplets.csv"
#
#     my_timer = Timer()
#
#     num_cores = multiprocessing.cpu_count()
#     triplets = Parallel(n_jobs=num_cores)(delayed(getTopGtpostprob_parallel)(PATH_TDI,FILE) for FILE in os.listdir(PATH_TDI))
#
#     outputTriplets(triplets,PATHNAME_triplet)
#
#     time_hhmmss = my_timer.get_time_hhmmss()
#     print("Time elapsed: %s\n" % time_hhmmss)