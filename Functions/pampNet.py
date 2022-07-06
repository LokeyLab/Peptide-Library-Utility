from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import statistics

import tensorflow as tf
from tensorflow import keras
from tensorflow.keras import layers

from sklearn.preprocessing import LabelEncoder
from sklearn.preprocessing import StandardScaler
from sklearn.neural_network import MLPRegressor
from sklearn import tree
from sklearn.model_selection import train_test_split



np.set_printoptions(precision=3, suppress=True)

def ReadTrainingData():
    df = pd.read_csv("TrainingData/hex.csv")
    X = df.iloc[0:, 10]
    y = df.iloc[0:, 14]

    # print(df.to_string())

    return X, y

def PreprocessData(X, y):
    df = pd.DataFrame(columns=["Exact Mass", "TPSA", "ALogP", "PappE-6"])

    indexCount = X.size
    index = 0
    for x in range(0, indexCount):
        mol = Chem.MolFromSmiles(X[x])

        molH = Chem.AddHs(mol)

        AllChem.EmbedMolecule(molH)

        AllChem.MMFFOptimizeMolecule(molH)


        """
        atomList = Chem.MolToMolBlock(molH).split("\n")
    
        newAtomList = []
    
        for i in range(4, len(atomList)):
            atom = atomList[i].split(" ")
    
            tempAtom = []
    
            for j in atom:
                if j.strip():
                    tempAtom.append(j)
    
            lb = LabelEncoder()
    
            atom = tempAtom[0:4]
            isAtom = atom[-1].isalpha()
    
            if isAtom:
                newAtomList.append(atom)
            else:
                break
    
        newAtomList = np.ravel(newAtomList)
        """

        exactMass = Descriptors.ExactMolWt(molH)
        TPSA = Descriptors.TPSA(molH)
        aLogP = Descriptors.MolLogP(molH)

        data = {"Exact Mass" : exactMass, "TPSA" : TPSA, "ALogP" : aLogP, "PappE-6" : y[x]}
        df.loc[len(df.index)] = data

        index+=1

        print(str(index) + "/" + str(indexCount))

    print(df)
    df.to_csv("TrainingData/hexData.csv")

    return df

def Train(data):
    # df = pd.read_csv(data)
    df = pd.read_csv("TrainingData/hex.csv")

    X = df.loc[:, ["Mass", "molAlogP"]]
    y = df.loc[:, "PappE-6"]

    X = np.nan_to_num(np.array(X))
    y = np.nan_to_num(np.array(y))

    X_train = X[:901]
    y_train = y[:901]
    X_test = X[901:]
    y_test = y[901:]

    print(len(X_train), len(X_test))
    """
    pampNetModel = tf.keras.Sequential([
        layers.Dense(len(X)),
        layers.Dense(len(X)),
        layers.Dense(len(X)),
        layers.Dense(len(X)),
        layers.Dense(len(X)),
                                        ])

    pampNetModel.compile(loss = tf.keras.losses.MeanSquaredError(),
                         optimizer = tf.keras.optimizers.Adam())

    pampNetModel.fit(X, y, epochs=10)"""
    """
    scaler = StandardScaler()
    
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.1, random_state=42)

    print("X_train shape = " + str(X_train.shape))
    print("y_train shape = " + str(y_train.shape))
    print("X_test shape = " + str(X_test.shape))
    print("y_test shape = " + str(y_test.shape))

    X_train.fillna(0, inplace=True)
    X_test.fillna(0, inplace=True)
    y_train.fillna(0, inplace=True)
    y_test.fillna(0, inplace=True)

    scaler.fit(X_train)
    X_train = scaler.transform(X_train)
    X_test = scaler.transform(X_test)

    y_train = y_train.values.ravel()
    y_test = y_test.values.ravel()

    regr = MLPRegressor(solver="adam", activation="relu", hidden_layer_sizes=(200,200),
                        learning_rate="invscaling", learning_rate_init=1e-5, max_iter=10000, random_state = 42)

    print("\nFitting...\n")
    regr = regr.fit(X_train, y_train)

    print("Predicting...\n")
    y_pred = regr.predict(X_test)
    difSum = 0

    print("Evaluating Difference...\n")
    for i in range(0, len(y_pred)):
        # print("Predicted: " + str(y_pred[i]) +
        #      "    True: " + str(y_test[i]) +
        #      "    Difference: " + str(y_pred[i] - y_test[i]))

        difSum += y_pred[i] - y_test[i]

    print("Mean Difference = " + str(difSum / len(y_pred)))
    print("SD Difference = " + str(statistics.stdev(y_pred - y_test)))

    print("\nEvaluating Model...")
    print('\nTest R^2 Score : %.3f' % regr.score(X_test, y_test))
    print('Training R^2 Score : %.3f' % regr.score(X_train, y_train))
    print("Loss : ", regr.loss_)
    # print("Number of Coefs : ", len(regr.coefs_))
    # print("Number of Intercepts : ", len(regr.intercepts_))
    print("Number of Iterations for Which Estimator Ran : ", regr.n_iter_)
    # print("Name of Output Layer Activation Function : ", regr.out_activation_)

    """

def PlotMol():
    """
    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(projection='3d')

    xs = df["X"].tolist()
    ys = df["Y"].tolist()
    zs = df["Z"].tolist()
    atoms = df["Atom"].tolist()

    for i in range(0, len(xs)):
        if atoms[i] == "H":
            ax.scatter(float(xs[i]), float(ys[i]), float(zs[i]),
                       color = "white", marker = "o", edgecolors='black', label = "Hydrogen")
        if atoms[i] == "C":
            ax.scatter(float(xs[i]), float(ys[i]), float(zs[i]),
                       color = "black", marker = "o", label = "Carbon")
        if atoms[i] == "N":
            ax.scatter(float(xs[i]), float(ys[i]), float(zs[i]),
                       color = "blue", marker = "o", label = "Nitrogen")
        if atoms[i] == "O":
            ax.scatter(float(xs[i]), float(ys[i]), float(zs[i]),
                       color = "red", marker = "o", label = "Oxygen")
        if atoms[i] == "S":
            ax.scatter(float(xs[i]), float(ys[i]), float(zs[i]),
                       color = "yellow", marker = "o", label = "Sulfur")

    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')

    plt.show()
    """
