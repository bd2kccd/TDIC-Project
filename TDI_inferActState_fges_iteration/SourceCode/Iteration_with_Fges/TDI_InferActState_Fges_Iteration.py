import subprocess,csv

prevProb = 0
count = 0
while True:
    count += 1
    subprocess.call("java -jar causal-cmd-0.3.2-jar-with-dependencies.jar --dataset ./DataSource/combinedMatrix.csv --out ./DataSource --prefix fges --data-type discrete --delimiter comma --knowledge ./DataSource/knowledge.txt --algorithm fges --score bdeu --samplePrior 1.0 --structurePrior 1.0 --symmetricFirstStep --maxDegree 100 --bootstrapSampleSize 0 --bootstrapEnsemble 1 --skip-latest",shell=True)

    subprocess.call("python convertFgesEdgesToInferActStateEdges.py", shell=True)

    subprocess.call("inferinnetwork_tdi.exe -m ./DataSource/combinedMatrix.csv -s ./DataSource/SGAstate.csv -e ./DataSource/fgesEdges.csv -o ./DataSource -x 2",shell=True)

    with open('./DataSource/jointProbsAll.csv') as f:
        reader = csv.reader(f)
        line = next(reader)
    currProb = float(line[-1])
    if abs(currProb - prevProb) < 0.01 or count >= 20:
        break
    else:
        prevProb = currProb