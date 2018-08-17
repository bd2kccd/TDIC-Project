import updateGlobalDriverMatrix
import initGlobalDriverMatrix
import initDEGposterirsAndDEGtoSGAs
import extractTriplets
# import subprocess
import os
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


my_timer = Timer()

# file_GtM ="./DataSource/tGtM.csv"
# file_GeM ="./DataSource/tGeM.csv"
# file_dict = "./DataSource/tdict.csv"
# file_gprior = "./DataSource/tGlobalPrior.csv"
# file_prior = "./DataSource/tprior.csv"
# file_triplet_prefix = "./DataSource/tTriplets"
# path_TDI = "./tTDICoutput"

# file_DEGposteriors = "./DataSource/tDEGposteriors.csv"
# file_DEGtoSGAs = "./DataSource/tDEGtoSGAs.csv"
# file_GdM = "./DataSource/tGdM.csv"


# file_GtM ="/home/xim33/DataSource/PanCancer13tts.SGAmatrix.4TCI.1217.csv"
# file_GeM ="/home/xim33/DataSource/PanCancer13tts.DEGmatrix.4TCI.1217.csv"
# file_dict = "/home/xim33/DataSource/PanCancer13tts.dict.1217.csv"
# file_gprior = "/home/xim33/DataSource/strcprior.poplv.20171217.csv"
# file_prior = "/home/xim33/DataSource/strcprior.tumorspecific.20171217.csv"
# file_triplet_prefix = "./DataSource/Triplets"
# path_TDI = "./TDICoutput"

# file_DEGposteriors = "./DataSource/DEGposteriors.csv"
# file_DEGtoSGAs = "./DataSource/DEGtoSGAs.csv"
# file_GdM = "./DataSource/GdM.csv"


file_GtM ="/home/xim33/DataSource/PanCancer13tts/PanCancer13tts.SGAmatrix.4TCI.csv"
file_GeM ="/home/xim33/DataSource/PanCancer13tts/PanCancer13tts.DEGmatrix.4TCI.csv"
file_dict = "/home/xim33/DataSource/PanCancer13tts/dict.pvmin=0.005.csv"
file_gprior = "/home/xim33/DataSource/PanCancer13tts/strcpriors.poplv.pvmin=0.005.csv"
file_prior = "/home/xim33/DataSource/PanCancer13tts/strcpriors.tumorspecific.pvmin=0.005.csv"

path_TDI = "./PanCancer13tts/TDIoutput_pvmin=0.005"
path_Results = "./PanCancer13tts/Results_pvmin=0.005"

file_triplet_prefix = path_Results + "/Triplets"
file_DEGposteriors = path_Results +  "/DEGposteriors.csv"
file_DEGtoSGAs = path_Results + "/DEGtoSGAs.csv"
file_GdM = path_Results + "/GdM.csv"

if not os.path.exists(path_TDI):
    os.makedirs(path_TDI)

if not os.path.exists(path_Results):
    os.makedirs(path_Results)
	
# generate global driver
os.system("./TDIC_GD_exeOMP -p " + file_gprior + " -f " + file_GtM + " -d " + file_GeM + " -o " + file_dict);

#convert global driver list ot globle driver matrix
# os.system("initGlobalDriverMatrix.py " + file_GtM + " " + file_GeM + " " + file_dict + " " + file_GdM) #does not work
initGlobalDriverMatrix.initGlobalDriver(file_GtM, file_GeM, file_dict, file_GdM)

# init DEGposteriors and DEGtoSGAs
initDEGposterirsAndDEGtoSGAs.initDEGposteriorAndSGAs(file_GeM, file_DEGposteriors, file_DEGtoSGAs)

count = 0
while True:
    count += 1
    
    ## Run TDI
    os.system("./PanCanTDICexeOOMP -p "+ file_prior +  " -f " + file_GtM + " -d " + file_GeM + " -g " + file_GdM + " -o " + path_TDI)

    ##extract triplets
    file_triplet = file_triplet_prefix+"_"+str(count)+".csv"
    extractTriplets.extractTriplets_real(path_TDI, file_triplet)

    #update GlobalDriverMatrix
    print "This is the " + str(count) + "th run"
    converged = updateGlobalDriverMatrix.updateGlobalDriver(count, file_triplet,file_GtM,file_DEGposteriors,file_DEGtoSGAs,file_GdM)
    if converged or count == 5:
        print "\nTotal times of iterations is: " + str(count)
        break

time_hhmmss = my_timer.get_time_hhmmss()
print("Final time elapsed: %s\n" % time_hhmmss)
