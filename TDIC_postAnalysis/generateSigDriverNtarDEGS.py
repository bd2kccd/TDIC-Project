import pandas as pd
import csv
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

def generateSigDriverNtarDEGS(FILE_filterdTriplets, FILE_SGA, PATH_output, FILE_driverTriplet, FILE_sigDrivers="sigDrivers.csv", FILE_sigTarDEGs = 'sigTarDEGs.csv'):
    my_timer = Timer()
    df = pd.read_csv(FILE_filterdTriplets)
    # DEG/tumor/SGA >= 5
    driverTriplet = df.groupby(['patient_name', 'cause_gene_name']).filter(lambda x: len(x) >= 5)
    driverTriplet.to_csv(FILE_driverTriplet, index=False)
    # tumor perDriver>= 30
    sigDriverTripletsTemp = driverTriplet.groupby('cause_gene_name').filter(lambda x: len(x) >= 30)

    df_SGA = pd.read_csv(FILE_SGA, index_col= 0)
    SGAs = list(df_SGA)

    #Drivercallrate >= 50%
    sigDrivers = {}
    count_driver = sigDriverTripletsTemp.groupby('cause_gene_name').size()
    count_SGA = df_SGA.sum()
    for sga, count in count_driver.iteritems():
        if sga in SGAs:
            rate = count_driver[sga] / float(count_SGA[sga])
            if rate >= 0.5:
                sigDrivers[sga] = rate

    sigDriverList = list(sigDrivers.keys())
    sigDriverTriplets = sigDriverTripletsTemp[sigDriverTripletsTemp['cause_gene_name'].isin(sigDriverList)]

    #get tarDEGs
    # sigDriverTriplets50 = sigDriverTriplets.groupby('cause_gene_name').filter(lambda x: len(x) >= 50)
    # sigDriver50list = sigDriverTriplets50['cause_gene_name'].unique()
    #??????

    # sigTarDEGs = {}
    # for SGA in sigDriverList:
        # SGAtriplets = sigDriverTriplets[sigDriverTriplets['cause_gene_name']== SGA]
        # tumor_count = len(SGAtriplets['patient_name'].unique())
        # tarDEGs = []
        # count_DEG = SGAtriplets.groupby('result_gene_name').size()
        # for DEG, count in count_DEG.items():
            # DEGrate = count/float(tumor_count)
            # if DEGrate >= 0.2:
                # tarDEGs.append(DEG)
        # if tarDEGs != []:
            # sigTarDEGs[SGA] = tarDEGs

    SGA_grouped =sigDriverTriplets.groupby('cause_gene_name')
    sigTarDEGs = {}
    for SGA, group in SGA_grouped:
        tumor_count = len(group['patient_name'].unique())
        tarDEGs = []
        for DEG, count in group.groupby('result_gene_name').size().items():
            DEGrate = count/float(tumor_count)
            if DEGrate >= 0.2:
                tarDEGs.append(DEG)
        if tarDEGs != []:
            sigTarDEGs[SGA] = tarDEGs

    #output sigDriver, sigTarDEGs
    if PATH_output[-1] != "\\" and PATH_output[-1] != "/":
        PATH_output += "/"
    with open(PATH_output + FILE_sigDrivers, 'w') as fdriver:
        for key in sigDrivers:
            fdriver.write(key + "," + str(sigDrivers[key]) + "\n")

    with open(PATH_output + FILE_sigTarDEGs, 'w') as ftarDEGs:
        for key in sigTarDEGs.keys():
            ftarDEGs.write(key + ":")
            for item in sigTarDEGs[key]:
                ftarDEGs.write(","+ item)
            ftarDEGs.write("\n")

    print "Process finished successfully"
    time_hhmmss = my_timer.get_time_hhmmss()
    print("Total time elapsed: %s\n" % time_hhmmss)

PATH_output = "/home/xim33/Python_Projects/TDI_ResultsAnalysis/OV"
FILE_filterdTriplets = "/home/xim33/Python_Projects/TDI_ResultsAnalysis/OV/FilteredTriplet_OV.csv"
FILE_SGA = "/home/xim33/DataSource/PanCancer13tts.SGAmatrix.4TCI.1217_OV.csv"
FILE_driverTriplet = "/home/xim33/Python_Projects/TDI_ResultsAnalysis/OV/driverTriplet_OV.csv"
generateSigDriverNtarDEGS(FILE_filterdTriplets, FILE_SGA, PATH_output, FILE_driverTriplet,FILE_sigDrivers="sigDrivers_OV.csv", FILE_sigTarDEGs = 'sigTarDEGs_OV.csv')

