import csv

# FILE_DEGMatrix = "C:/Users/XIM33/Documents/NetBeansProjects/DataSource/tGeM.csv"
# FILE_DEGposteriors = "C:/Users/XIM33/Documents/TDIC_iteration/DataSource/tDEGposteriors.csv"
# FILE_DEGtoSGAs = "C:/Users/XIM33/Documents/TDIC_iteration/DataSource/tDEGtoSGAs.csv"
def initDEGposteriorAndSGAs(FILE_DEGMatrix,FILE_DEGposteriors,FILE_DEGtoSGAs):
	#read the first line of DEGmatrix to get the DEG names
	with open(FILE_DEGMatrix, 'r') as fin:
		reader = csv.reader(fin)
		DEGs = next(reader)[1:]

	#write to DEGposteriors files
	with open(FILE_DEGposteriors, 'w') as fout1:
		for item in DEGs:
			fout1.write(item+"\n")

	#write to DEGposteriors files
	with open(FILE_DEGtoSGAs, 'w') as fout2:
		for item in DEGs:
			fout2.write(item+"\n")
