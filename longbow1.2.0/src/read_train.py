import pandas as pd

def read_train_file(filename : str) -> tuple:

    df = pd.read_csv(filename)
    df_cleaned = df.dropna()
    drop_colums = ["Sample", "Config"] + ["R" + str(i) for i in range(0, 91)] + ["id"] + ["Q0"] + ["K" + str(i) for i in range(0, 91)]

    df_final = df_cleaned.drop(columns = drop_colums, axis = 1)
    
    y = df_final["Label"].astype(int)
    X = df_final.drop("Label", axis = 1)
    X = X.values

    return (X, y)

"""
if __name__ == "__main__":
    test = read_train_file("/home/mencius/test/longbow1.2.0/model/train_guppy.csv")
    for i in test[0]:
        print(i)
"""
