import pandas as pd

def csv2excel(filename):

    df = pd.read_csv(filename)
    df.to_excel(filename.replace(".csv",".xlsx"))