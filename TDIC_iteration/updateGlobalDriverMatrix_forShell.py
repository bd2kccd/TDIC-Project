import csv
import pandas as pd
import sys


def updateGlobalDriver(count, file_triplet,file_SGA,file_DEGposteriors,file_DEGtoSGAs,file_GdMatrix):
    #read in SGAMatrix
    df_SGA = pd.read_csv(file_SGA, index_col=0)
    df_SGA.insert(0,'A0',1)
    #read in GdMatrix
    df_GdMatrix = pd.read_csv(file_GdMatrix, index_col=0)
    #read in triplet
    df_TSDPtriplet = pd.read_csv(file_triplet)
    # readinDEGposteriorFile into DEGtoPosteriors dict
    mapDEGposteriors = {}
    with open(file_DEGposteriors, 'r') as fin1:
        reader = csv.reader(fin1)
        for row in reader:
            DEG = row[0]
            posteriors = [float(i) for i in row[1:] ]
            mapDEGposteriors[DEG]= posteriors

    # readinDEGtoSGAs file into mapDEGtoSGAs dict
    mapDEGtoSGAs = {}
    with open(file_DEGtoSGAs, 'r') as fin2:
        reader = csv.reader(fin2)
        for row in reader:
            DEG = row[0]
            SGAs = row[1:]
            mapDEGtoSGAs[DEG]= SGAs


    allConverged = True

    grouped = df_TSDPtriplet.groupby('result_gene_name')
    for DEG, group in grouped:
        PosteriorsList = mapDEGposteriors[DEG]
        if not PosteriorsList:
            lastPosterior = 0
        else:
            lastPosterior = PosteriorsList[-1]
        #if posterior is converged, do nothing
        if (lastPosterior == 1 or lastPosterior == 0.99999) :
            continue #goto the next DEG
        #calculate the ave posterior of this DEG
        currentPosterior = group['posterior'].mean()
        diff = currentPosterior - lastPosterior
        if diff > 0.05: #not converged
            allConverged = False

            #update mapDEGPosteriors
            PosteriorsList.append(currentPosterior)
            mapDEGposteriors[DEG] = PosteriorsList

            #update mapDEGtoSGAs

            ##get the most n frequent SGA for each DEG
            SGAlist = []
            tops = group['cause_gene_name'].value_counts().head(count)
            ## added to SGAlist
            for idx in range(0,tops.size):  #tops.size maybe less than count, use count may cause index overflow
                SGAlist.append(tops.index[idx])
            ## update mapDEGtoSGAs
            mapDEGtoSGAs[DEG] = SGAlist

            #update DEGglobalDriverMatrix

            ## create an 0 series with the same index of df_SGA
            sumGt = pd.Series(data=0, index=df_SGA.index)
            ## add all df_SGA item (or operation
            for item in SGAlist:
                sumGt += df_SGA[item]
            ## update df_GdMatrix
            df_GdMatrix[DEG]=sumGt
        elif diff > 0:#converged
            #append 1 to posteriorList
            PosteriorsList.append(1)
            mapDEGposteriors[DEG] = PosteriorsList
        else:
            PosteriorsList.append(0.99999)
            mapDEGposteriors[DEG] = PosteriorsList
    # write back file_DEGposterior
    with open(file_DEGposteriors, 'w') as fout1:
        for key, value in mapDEGposteriors.iteritems():
            fout1.write(key)
            for item in value:
                fout1.write(","+str(item))
            fout1.write("\n")

    # write back file_DEGtoSGAs
    with open(file_DEGtoSGAs, 'w') as fout2:
        for key, value in mapDEGtoSGAs.iteritems():
            fout2.write(key)
            for item in value:
                fout2.write("," + item)
            fout2.write("\n")

    # write back GdMatrix
    df_GdMatrix.to_csv(file_GdMatrix)

    if allConverged:
        exit(1)
    else:
        exit(0)


count = sys.argv[1]
file_triplet = sys.argv[2]
file_SGA = sys.argv[3]
file_DEGposteriors = sys.argv[4]
file_DEGtoSGAs = sys.argv[5]
file_GdMatrix = sys.argv[6]

# file_triplet = "C:/Users/XIM33/Documents/TDIC_iteration/DataSource/tTriplets.csv"
# file_SGA = "C:/Users/XIM33/Documents/NetBeansProjects/DataSource/tGtM.csv"
# file_DEGposteriors = "C:/Users/XIM33/Documents/TDIC_iteration/DataSource/tDEGposteriors.csv"
# file_DEGtoSGAs = "C:/Users/XIM33/Documents/TDIC_iteration/DataSource/tDEGtoSGAs.csv"
# file_GdMatrix = "C:/Users/XIM33/Documents/TDIC_iteration/DataSource/tGdM.csv"
#
updateGlobalDriver(count,file_triplet,file_SGA,file_DEGposteriors,file_DEGtoSGAs,file_GdMatrix)