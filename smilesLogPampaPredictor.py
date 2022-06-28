from pandas import read_csv
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from sklearn.preprocessing import LabelEncoder

from sklearn.model_selection import train_test_split
from sklearn import tree
from sklearn.neural_network import MLPRegressor
from sklearn.metrics import accuracy_score

from rdkit.Chem import MolFromSmiles
from rdkit.Chem import Descriptors


def TrainModel():
    df = read_csv("hex.csv")

    X = df.iloc[0:, 10]
    y = df.iloc[0:, 14]

    longestSmiles = len(max(df.iloc[1:, 10], key = len)) + 100

    print("Longest SMILES string = " + str(longestSmiles))

    lb = LabelEncoder()

    for i in range(0, X.size):
        temp = list(X.iloc[i])
        temp = lb.fit_transform(temp)
        if len(temp) < longestSmiles:
            padSize = longestSmiles - temp.shape[0]
            temp = np.pad(temp, (0, padSize), 'constant')

        X[i] = temp

    X = np.array(X.tolist())
    y = np.array(y.tolist())

    X = pd.DataFrame(X)
    y = pd.DataFrame(y)

    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size = 0.1, random_state = 42)
    # X_test, X_valid, y_test, y_valid = train_test_split(X_test, y_test, test_size=0.5, random_state=42)

    print("X_train shape = " + str(X_train.shape))
    print("y_train shape = " + str(y_train.shape))
    print("X_test shape = " + str(X_test.shape))
    print("y_test shape = " + str(y_test.shape))

    y_train.fillna(0, inplace=True)
    y_test.fillna(0, inplace=True)
    # y_valid.fillna(0, inplace=True)

    mlp = MLPRegressor(solver='adam', alpha=1, activation = 'relu', max_iter = 2000,
                       hidden_layer_sizes = 200, random_state = 42)

    fig, ax = plt.subplots(figsize=(6, 4))

    # lcValid = mlp.fit(X_valid, y_valid.values.ravel())
    # plt.plot(lcValid.loss_curve_, label = "validation")

    mlp = mlp.fit(X_train, y_train.values.ravel())
    # plt.plot(mlp.loss_curve_, label = "train")

    """
    plt.legend()

    ax.set_xlabel('Number of iterations')
    ax.set_ylabel('Loss')
    plt.show()
    plt.close()
    """
    y_pred = mlp.predict(X_test)

    print('Test R^2 Score : %.3f' % mlp.score(X_test, y_test))
    print('Training R^2 Score : %.3f' % mlp.score(X_train, y_train))
    print("Loss : ", mlp.loss_)
    print("Number of Coefs : ", len(mlp.coefs_))
    print("Number of Intercepts : ", len(mlp.intercepts_))
    print("Number of Iterations for Which Estimator Ran : ", mlp.n_iter_)
    print("Name of Output Layer Activation Function : ", mlp.out_activation_)


    # accuracy = "Accuracy = " + str(mlp.score(X_test, y_test))

    # print(accuracy)

    return mlp, longestSmiles

def PrepSmilesString(longestSmiles):
    smilesString = input("Enter SMILES string: ")
    lb = LabelEncoder()

    temp = list(smilesString)
    temp = lb.fit_transform(temp)
    if len(temp) < longestSmiles:
        padSize = longestSmiles - temp.shape[0]
        temp = np.pad(temp, (0, padSize), 'constant')

        x = temp.reshape(1,-1)

    return x, smilesString

def Main():
    mlp, longestSmiles = TrainModel()
    prepedSmilesString, smilesString = PrepSmilesString(longestSmiles)
    prediction = mlp.predict(prepedSmilesString)
    LogP = Descriptors.MolLogP(MolFromSmiles(smilesString))
    print("Predicted PAMPA = " + str(prediction[0]) + "    LogP = " + str(LogP))


    fig, ax = plt.subplots(figsize=(6, 4))
    plt.axhline(y=10, color='r', linestyle='--')
    plt.axvline(x=1.5, color='r', linestyle='--')
    plt.axvline(x=4.5, color='r', linestyle='--')
    plt.scatter(LogP, prediction[0])
    ax.set_xlabel('LogP')
    ax.set_ylabel('Papp')

    plt.xlim([-5, 5])
    plt.ylim([-5, 40])
    plt.show()

Main()