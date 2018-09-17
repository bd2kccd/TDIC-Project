fileNameIn = "./DataSource/fges.txt"
fileNameOut = "./DataSource/fgesEdges.csv"
with open(fileNameIn) as fin:
    content = fin.readlines()

with open(fileNameOut, 'w') as fout:

    for line in content:
        if "-->" not in line:
            continue
        pos=line.find('.')
        line = line[pos+2:]
        line = line.replace(" --> ",",")
        fout.write(line)