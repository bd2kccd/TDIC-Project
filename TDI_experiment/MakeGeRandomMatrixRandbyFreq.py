import csv, random,sys

def readinMatrix(file_in):
    f_in = open(file_in, "r")

    count = -1
    cols = []
    tumors = []
    for row in csv.reader(f_in):
        count += 1
        if( count == 0) :#title
            cols = [0] * len(row)
        else:
            tumors.append(row[0])
            for i in range(1,len(row)): #i=0 it is column of tumor name
                cols[i] += int(row[i])
    f_in.close()


    #remove dulicated number by convert list to set
    colsSet = set(cols[1:])
    return tumors, colsSet


def make_random_cols(setCols, totalNum):
    map_freqToListRandNums = {}
    for freq in setCols:
        randSet = set()
        while True:
            randNum = random.randrange(1, totalNum + 1)
            randSet.add(randNum)
            if len(randSet) >= freq:
                break
        randList = list(randSet)
        map_freqToListRandNums[freq]=randList
    return map_freqToListRandNums

def outputMatrix(file_out, tumors, m_freqToRands):
    f_out = open(file_out, "w")
    #output the first title line
    for key in m_freqToRands.keys():
        f_out.write(",freq"+ str(key))
    f_out.write("\n")
    #output the random matrix
    count = 0
    for tumor in tumors:
        f_out.write(tumor)
        count += 1
        for key in m_freqToRands.keys():
            rands = m_freqToRands[key]
            if count in rands:
                f_out.write("," + '1')
            else:
                f_out.write("," + '0')
        f_out.write("\n")
    f_out.close()


file_in = sys.argv[1]
file_outPrefix = sys.argv[2]
file_out_range1 = int(sys.argv[3])	
file_out_range2 = int(sys.argv[4])+1
l_tumors, set_freqs = readinMatrix(file_in)	
for i in range(file_out_range1,file_out_range2):
    file_out = file_outPrefix+"_randbyfreq"+str(i)+".csv"
    map_freqToRands = make_random_cols(set_freqs, len(l_tumors))
    outputMatrix(file_out, l_tumors, map_freqToRands)
	
	
	
# file_in = "../DataSource/PanCancer13tts.DEGmatrix.twostates.csv"

# l_tumors, set_freqs = readinMatrix(file_in)
# for i in range(1,6):
    # file_out = "./DataSource/PanCancer13tts.DEGmatrix.twostates_randbyfreq"+str(i)+".csv"
    # map_freqToRands = make_random_cols(set_freqs, len(l_tumors))
    # outputMatrix(file_out, l_tumors, map_freqToRands)
