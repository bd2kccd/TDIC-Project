#!/bin/bash

SGAori=./DataSource/SGAmatrix.csv
DEGori=./DataSource/DEGmatrix.csv
globalPrior=./DataSource/strcprior.poplv.20171217.csv
prior=./DataSource/strcprior.tumorspecific.20171217.csv

DEGrandOutPrefix=./DataSource/DEGrand/DEG
dictPrefix=./DataSource/dict

DEGrandRange1=1
DEGrandRange2=3

DEGrandbyfreqRange1=1
DEGrandbyfreqRange2=5

mkdir -p ./DataSource/DEGrand
python Make_random_matrix.py $DEGori $DEGrandOutPrefix $DEGrandRange1 $DEGrandRange2

python MakeGeRandomMatrixRandbyFreq.py $DEGori $DEGrandOutPrefix $DEGrandbyfreqRange1 $DEGrandbyfreqRange2;


mkdir -p ./DataSource
./TDIC_GD_exeOMP -p $globalPrior -f $SGAori -d $DEGori -o "$dictPrefix.csv";

for value in $(seq $DEGrandRange1 $DEGrandRange2)
do
echo $DEGrandOutPrefix"_rand"$value".csv" 
./TDIC_GD_exeOMP -p $globalPrior -f $SGAori -d $DEGrandOutPrefix"_rand"$value".csv" -o $dictPrefix"_rand"$value".csv";
done

for value in $( seq $DEGrandbyfreqRange1 $DEGrandbyfreqRange2)
do
./TDIC_GD_exeOMP -p $globalPrior -f $SGAori -d $DEGrandOutPrefix"_randbyfreq"$value".csv" -o $dictPrefix"_randbyfreq"$value".csv";
done

mkdir -p ./Results/TDI
./PanCanTDICexeOOMP -p $prior -f $SGAori -d $DEGori -g $dictPrefix".csv" -o ./Results/TDI

mkdir -p ./Results/TDI_rand
for value in $(seq $DEGrandRange1 $DEGrandRange2)
do
mkdir -p ./Results/TDI_rand/"rand"$value
nohup ./PanCanTDICexeOOMP -p $prior -f $SGAori -d $DEGrandOutPrefix"_rand"$value".csv"  -g $dictPrefix"_rand"$value".csv" -o "./Results/TDI_rand/rand"$value
done

mkdir -p ./Results/TDI_randbyfreq
for value in $(seq $DEGrandbyfreqRange1 $DEGrandbyfreqRange2)
do
mkdir -p ./Results/TDI_randbyfreq/"randbyfreq"$value
nohup ./PanCanTDICexeOOMP -p $prior -f $SGAori -d $DEGrandOutPrefix"_randbyfreq"$value".csv"  -g $dictPrefix"_randbyfreq"$value".csv" -o "./Results/TDI_randbyfreq/randbyfreq"$value
done

mkdir -p ./Results/Triplets_randbyfreq

python postAnalysisDataPreparation.py $SGAori;

python generateSigDriverNtarDEGS.py $SGAori;






















# nohup python MakeGeRandomMatrixRandbyFreq.py;

# nohup ./TDIC_GD_exeOMP -p ../DataSource/strcprior.poplv.20171217.csv -f ../DataSource/PanCancer13tts.SGAmatrix.4TCI.1217.csv -d ../DataSource/PanCancer13tts.DEGmatrix.twostates.csv -o ./DataSource/PanCancer13tts.dict.4TCI.1217.csv;

# nohup ./TDIC_GD_exeOMP -p ../DataSource/strcprior.poplv.20171217.csv -f ../DataSource/PanCancer13tts.SGAmatrix.4TCI.1217.csv -d ./DataSource/PanCancer13tts.DEGmatrix.twostates_rand1.csv -o ./DataSource/PanCancer13tts.dict.4TCI.1217_rand1.csv;
# nohup ./TDIC_GD_exeOMP -p ../DataSource/strcprior.poplv.20171217.csv -f ../DataSource/PanCancer13tts.SGAmatrix.4TCI.1217.csv -d ./DataSource/PanCancer13tts.DEGmatrix.twostates_rand2.csv -o ./DataSource/PanCancer13tts.dict.4TCI.1217_rand2.csv;

# nohup ./TDIC_GD_exeOMP -p ../DataSource/strcprior.poplv.20171217.csv -f ../DataSource/PanCancer13tts.SGAmatrix.4TCI.1217.csv -d ./DataSource/PanCancer13tts.DEGmatrix.twostates_randbyfreq1.csv -o ./DataSource/PanCancer13tts.dict.4TCI.1217_randbyfreq1.csv;
# nohup ./TDIC_GD_exeOMP -p ../DataSource/strcprior.poplv.20171217.csv -f ../DataSource/PanCancer13tts.SGAmatrix.4TCI.1217.csv -d ./DataSource/PanCancer13tts.DEGmatrix.twostates_randbyfreq2.csv -o ./DataSource/PanCancer13tts.dict.4TCI.1217_randbyfreq2.csv;
# nohup ./TDIC_GD_exeOMP -p ../DataSource/strcprior.poplv.20171217.csv -f ../DataSource/PanCancer13tts.SGAmatrix.4TCI.1217.csv -d ./DataSource/PanCancer13tts.DEGmatrix.twostates_randbyfreq3.csv -o ./DataSource/PanCancer13tts.dict.4TCI.1217_randbyfreq3.csv;
# nohup ./TDIC_GD_exeOMP -p ../DataSource/strcprior.poplv.20171217.csv -f ../DataSource/PanCancer13tts.SGAmatrix.4TCI.1217.csv -d ./DataSource/PanCancer13tts.DEGmatrix.twostates_randbyfreq4.csv -o ./DataSource/PanCancer13tts.dict.4TCI.1217_randbyfreq4.csv;
# nohup ./TDIC_GD_exeOMP -p ../DataSource/strcprior.poplv.20171217.csv -f ../DataSource/PanCancer13tts.SGAmatrix.4TCI.1217.csv -d ./DataSource/PanCancer13tts.DEGmatrix.twostates_randbyfreq5.csv -o ./DataSource/PanCancer13tts.dict.4TCI.1217_randbyfreq5.csv;

# nohup ./PanCanTDICexeOOMP -p ../DataSource/strcprior.tumorspecific.20171217.csv -f ../DataSource/PanCancer13tts.SGAmatrix.4TCI.1217.csv -d ../DataSource/PanCancer13tts.DEGmatrix.twostates.csv  -g ./DataSource/PanCancer13tts.dict.4TCI.1217.csv -o ./Results/TDI

# nohup ./PanCanTDICexeOOMP -p ../DataSource/strcprior.tumorspecific.20171217.csv -f ../DataSource/PanCancer13tts.SGAmatrix.4TCI.1217.csv -d ./DataSource/PanCancer13tts.DEGmatrix.twostates_rand1.csv  -g ./DataSource/PanCancer13tts.dict.4TCI.1217_rand1.csv -o ./Results/TDI_rand/rand1;
# nohup ./PanCanTDICexeOOMP -p ../DataSource/strcprior.tumorspecific.20171217.csv -f ../DataSource/PanCancer13tts.SGAmatrix.4TCI.1217.csv -d ./DataSource/PanCancer13tts.DEGmatrix.twostates_rand2.csv  -g ./DataSource/PanCancer13tts.dict.4TCI.1217_rand2.csv -o ./Results/TDI_rand/rand2;

# nohup ./PanCanTDICexeOOMP -p ../DataSource/strcprior.tumorspecific.20171217.csv -f ../DataSource/PanCancer13tts.SGAmatrix.4TCI.1217.csv -d ./DataSource/PanCancer13tts.DEGmatrix.twostates_randbyfreq1.csv  -g ./DataSource/PanCancer13tts.dict.4TCI.1217_randbyfreq1.csv -o ./Results/TDI_randbyfreq/randbyfreq1;
# nohup ./PanCanTDICexeOOMP -p ../DataSource/strcprior.tumorspecific.20171217.csv -f ../DataSource/PanCancer13tts.SGAmatrix.4TCI.1217.csv -d ./DataSource/PanCancer13tts.DEGmatrix.twostates_randbyfreq2.csv  -g ./DataSource/PanCancer13tts.dict.4TCI.1217_randbyfreq2.csv -o ./Results/TDI_randbyfreq/randbyfreq2;
# nohup ./PanCanTDICexeOOMP -p ../DataSource/strcprior.tumorspecific.20171217.csv -f ../DataSource/PanCancer13tts.SGAmatrix.4TCI.1217.csv -d ./DataSource/PanCancer13tts.DEGmatrix.twostates_randbyfreq3.csv  -g ./DataSource/PanCancer13tts.dict.4TCI.1217_randbyfreq3.csv -o ./Results/TDI_randbyfreq/randbyfreq3;
# nohup ./PanCanTDICexeOOMP -p ../DataSource/strcprior.tumorspecific.20171217.csv -f ../DataSource/PanCancer13tts.SGAmatrix.4TCI.1217.csv -d ./DataSource/PanCancer13tts.DEGmatrix.twostates_randbyfreq4.csv  -g ./DataSource/PanCancer13tts.dict.4TCI.1217_randbyfreq4.csv -o ./Results/TDI_randbyfreq/randbyfreq4;
# nohup ./PanCanTDICexeOOMP -p ../DataSource/strcprior.tumorspecific.20171217.csv -f ../DataSource/PanCancer13tts.SGAmatrix.4TCI.1217.csv -d ./DataSource/PanCancer13tts.DEGmatrix.twostates_randbyfreq5.csv  -g ./DataSource/PanCancer13tts.dict.4TCI.1217_randbyfreq5.csv -o ./Results/TDI_randbyfreq/randbyfreq5;

# nohup python postAnalysisDataPreparation.py;

# nohup python generateSigDriverNtarDEGS.py;
