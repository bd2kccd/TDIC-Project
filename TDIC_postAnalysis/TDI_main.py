import ExtractTriplets_parallel
if __name__ == "__main__":
    PATH_TDI = "C:/Users/XIM33/Documents/NetBeansProjects/TDIC_01_CPU_Tumors_SGAprior/output/"
    PATHNAME_triplet = "C:/Users/XIM33/Documents/NetBeansProjects/TDIC_01_CPU_Tumors_SGAprior/triplets.csv"
    ExtractTriplets_parallel.ExtractTriplet_parallel(PATH_TDI, PATHNAME_triplet)

