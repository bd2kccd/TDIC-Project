import  pandas as pd

file_sga = "./DataSource/InitialDriverstate.csv"
file_deg = "./DataSource/DEGmatrix.csv"
file_combined = "./DataSource/combinedMatrix.csv"
df_sga = pd.read_csv(file_sga, index_col=0)
df_sga.rename(columns=lambda x: x+"_sga", inplace=True)
df_deg = pd.read_csv(file_deg, index_col=0)
df_deg.rename(columns=lambda x: x+"_deg", inplace=True)
df_combined = pd.concat([df_sga, df_deg], axis=1)
df_combined.to_csv(file_combined, index=False)
