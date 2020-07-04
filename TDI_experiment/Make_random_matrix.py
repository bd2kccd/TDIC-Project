#!/usr/bin/python
import os, csv, sys, time, random, math

def get_random_number(p_candidates,p_keepNo):
   soluTp = {}
   while len(soluTp)< p_keepNo:
      rn = random.randrange(len(p_candidates))
      soluTp[p_candidates[rn]] = 1
   return soluTp


def readMatrix(p_input,p_output):
   #read a 0/1 matrix, generate a new random matrix
   # such that the number of 1 in the new matrix is
   # the same with the old matrix.
   
   # print p_input
   # print p_output
   L_rets = {}; L_genes = []
   ff = open(p_input, 'r')
   count = 0
   L_title = []
   for ln in ff:
      count += 1
      ln1 = ln.split('\n')
      ln2 = ln1[0].split(',')
      if count == 1:
         L_title = ln2
         for idx in range(1, len(ln2)):
            tu = L_title[idx]
            L_rets[tu] = 0
         # print '--- column number:', len(ln2)
      else:
         gene = ln2[0]
         L_genes.append(gene)
         for idx in range(1, len(ln2)):
            va = int(float(ln2[idx]))
            tu = L_title[idx]
            if not va==0:
               L_rets[tu] += 1
   ff.close()
   # print '--- row number',(len(L_genes)+1)
   L_random = {}
   for idx in range(1, len(L_title)):
      tu = L_title[idx]
      va = L_rets[tu]
      tRand = get_random_number(L_genes, va)
      L_random[tu] = tRand
   f1 = open(p_output,'w')
   output = ','.join(L_title)
   output += '\n'
   f1.write(output)
   for gene in L_genes:
      outTp = [gene]
      for idx in range(1, len(L_title)):
         tu = L_title[idx]
         if gene in L_random[tu]:
            outTp.append('1')
         else:
            outTp.append('0')

      output = ','.join(outTp)+'\n'
      f1.write(output)
   f1.close()
   
   print ('------- End --------')

#-------------------
# f_input = '../DataSource/PanCancer13tts.DEGmatrix.4TCI.1217_OV.csv'
# f_output = '../DataSource/PanCancer13tts.DEGmatrix.4TCI.1217_OV_rand1_ori.csv'
f_input = sys.argv[1]
f_outputPrefix = sys.argv[2]
f_output_range1 = int(sys.argv[3])
f_output_range2 = int(sys.argv[4])+1
# print f_input
# print f_outputPrefix
# print f_output_range1
# print f_output_range2
for i in range(f_output_range1, f_output_range2):
    f_output = f_outputPrefix+"_rand"+str(i)+".csv"
    readMatrix(f_input, f_output)
