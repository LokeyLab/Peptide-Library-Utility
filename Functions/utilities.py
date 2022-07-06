# External
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import multivariate_normal


def ReadFileToDataframe(f):

    df = pd.read_csv(f)

    return df


def CyclicPeptidesToDataframe(cyclicPeptides):
    columns = ["Name", "Exact Mass", "TPSA", "ALogP", "Predicted PappE-6", "SMILES String"]

    df = pd.DataFrame(columns=columns)

    for i in range(0,len(cyclicPeptides)):
        df.loc[len(df.index)] = { "Name" : cyclicPeptides[i].name,
                                  "Exact Mass" : cyclicPeptides[i].exactMass,
                                  "TPSA" : cyclicPeptides[i].TPSA,
                                  "ALogP" : cyclicPeptides[i].ALogP,
                                  "Predicted PappE-6": cyclicPeptides[i].predictedPapp,
                                  "SMILES String" : cyclicPeptides[i].smilesString }

    return df


def PrintDataframe(df):
    print("\n" + df.to_string())


def DataframeToCSV(df):
    df.to_csv("Output/output.csv")

def PlotDataframe(df):
    fig = plt.figure(figsize=(12,8))
    ax = fig.add_subplot(projection='3d')
    x = df["Exact Mass"]
    y = df["TPSA"]
    z = df["ALogP"]
    p = (x + y + z) # df["Predicted PappE-6"]

    print(x.max())
    print(y.min())
    print(z.min())

    ax.scatter(x, y, z, c = p, cmap="viridis", edgecolors='black', s=20)
    ax.plot(x, z, "r+", zdir='y', zs=y.max())
    ax.plot(y, z, "g+", zdir='x', zs=x.min())
    ax.plot(x, y, "b+", zdir='z', zs=z.min())
    ax.set_xlabel('Exact Mass')
    ax.set_ylabel('TPSA')
    ax.set_zlabel('ALogP')

    plt.show()

