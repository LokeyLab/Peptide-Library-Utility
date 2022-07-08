print("\nImporting Packages...")

from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors, rdFreeSASA
from rdkit.Chem import Crippen
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from ast import literal_eval
from tqdm import tqdm
from itertools import chain

from sklearn.preprocessing import StandardScaler, MinMaxScaler, RobustScaler

import tensorflow as tf
from tensorflow import keras
from keras.utils.vis_utils import plot_model
import visualkeras


print("\nImporting Packages Complete.")

def ReadTrainingData():
    print("\nReading Training Data File...")

    df = pd.read_csv("TrainingData/hex.csv")
    X = df.iloc[0:, 10]
    y = df.iloc[0:, 14]

    # print(df.to_string())

    print("\nReading Training Data File Complete.")

    return X, y


def PreprocessData(X, y):
    print("\nProcessing Training Data...\n")

    df = pd.DataFrame(columns=["Exact Mass", "ALogP", "TPSA", "SASA", "# Atoms", "# Rotatable Bonds",
                               "# Amide Bonds", "# Aromatic Rings","# HBA", "# HBD",  "Atoms", "PappE-6"])

    indexCount = X.size
    for x in tqdm(range(0, indexCount)):
        mol = Chem.MolFromSmiles(X[x])

        molH = Chem.AddHs(mol)

        AllChem.EmbedMolecule(molH)

        AllChem.MMFFOptimizeMolecule(molH)

        atomList = Chem.MolToMolBlock(molH).split("\n")

        newAtomList = []
        periodDict = {"H": 1, "He": 2, "Li": 3, "Be": 4, "B": 5, "C": 6, "N": 7, "O": 8, "F": 9, "Ne": 10, "Na": 11,
                      "Mg": 12, "Al": 13, "Si": 14, "P": 15, "S": 16, "Cl": 17, "Ar": 18, "K": 19, "Ca": 20, "Sc": 21,
                      "Ti": 22, "V": 23, "Cr": 24, "Mn": 25, "Fe": 26, "Co": 27, "Ni": 28, "Cu": 29, "Zn": 30,
                      "Ga": 31, "Ge": 32, "As": 33, "Se": 34, "Br": 35, "Kr": 36, "Rb": 37, "Sr": 38, "Y": 39,
                      "Zr": 40, "Nb": 41, "Mo": 42, "Tc": 43, "Ru": 44, "Rh": 45, "Pd": 46, "Ag": 47, "Cd": 48,
                      "In": 49, "Sn": 50, "Sb": 51, "Te": 52, "I": 53, "Xe": 54, "Cs": 55, "Ba": 56, "La": 57,
                      "Ce": 58, "Pr": 59, "Nd": 60, "Pm": 61, "Sm": 62, "Eu": 63, "Gd": 64, "Tb": 65, "Dy": 66,
                      "Ho": 67, "Er": 68, "Tm": 69, "Yb": 70, "Lu": 71, "Hf": 72, "Ta": 73, "W"	: 74, "Re": 75,
                      "Os": 76, "Ir": 77, "Pt": 78, "Au": 79, "Hg": 80, "Tl": 81, "Pb": 82, "Bi": 83, "Po": 84,
                      "At": 85, "Rn": 86, "Fr": 87, "Ra": 88, "Ac": 89, "Th": 90, "Pa": 91, "U"	: 94, "Np": 93,
                      "Pu": 94, "Am": 95, "Cm": 96, "Bk": 97, "Cf": 98, "Es": 99, "Fm": 100, "Md": 101, "No": 102,
                      "Lr": 103, "Rf": 104, "Db": 105, "Sg": 106, "Bh": 107, "Hs": 108, "Mt": 109, "Ds": 110,
                      "Rg": 111, "Cn": 112, "Nh": 113, "Fl": 114, "Mc": 115, "Lv": 116, "Ts": 117, "Og": 118 }

        for i in range(4, len(atomList)):
            atom = atomList[i].split(" ")
            tempAtom = []
    
            for j in atom:
                if j.strip():
                    tempAtom.append(j)

            atom = tempAtom[0:4]

            isAtom = atom[-1].isalpha()
    
            if isAtom:
                newAtomList.append(list(atom[0:3] + [periodDict[atom[3]]]))
            else:
                break

        exactMass = rdMolDescriptors.CalcExactMolWt(molH)
        aLogP = Crippen.MolLogP(molH)
        tPSA = rdMolDescriptors.CalcTPSA(molH)
        radii = rdFreeSASA.classifyAtoms(molH)
        sASA = rdFreeSASA.CalcSASA(molH, radii)
        numAtoms = rdMolDescriptors.CalcNumAtoms(molH)
        numRotatableBonds = rdMolDescriptors.CalcNumRotatableBonds(molH)
        numAmideBonds = rdMolDescriptors.CalcNumAmideBonds(molH)
        numAromaticRings = rdMolDescriptors.CalcNumAromaticRings(molH)
        numHBD = rdMolDescriptors.CalcNumHBA(molH)
        numHBA = rdMolDescriptors.CalcNumHBD(molH)

        data = {"Exact Mass" : exactMass,
                "ALogP" : aLogP,
                "TPSA": tPSA,
                "SASA": sASA,
                "# Atoms" : numAtoms,
                "# Rotatable Bonds": numRotatableBonds,
                "# Amide Bonds" : numAmideBonds,
                "# Aromatic Rings": numAromaticRings,
                "# HBA": numHBA,
                "# HBD": numHBD,
                "Atoms" : newAtomList,
                "PappE-6" : y[x]}

        df.loc[len(df.index)] = data

    print("\n" + df.to_string())

    df.to_csv("TrainingData/hexData.csv")

    print("\nProcessing Training Data Complete.")

    return df


def Train():

    # df = pd.read_csv(data)
    df = pd.read_csv("TrainingData/hexData.csv")

    X = np.asarray(df[["Exact Mass", "TPSA", "ALogP", "# Atoms", "# Rotatable Bonds",
                     "# Amide Bonds", "# Aromatic Rings"]])

    X2 = np.asarray(df["Atoms"].apply(literal_eval))

    for x in range(0, len(X2)):
        X2[x] = list(chain(*X2[x]))
        for i in range(0, len(X2[x])):
            X2[x][i] = np.float64(X2[x][i])
        X2[x] = np.stack(X2[x])

    largestAtomCount = len(max(X2, key=len))

    for j in range(0, len(X2)):
        X2[j] = np.pad(X2[j], (0, largestAtomCount - len(X2[j] + 50)), 'constant', constant_values= 0.0)

    X2 = np.stack(X2)

    y = np.asarray(df["PappE-6"])

    X = np.nan_to_num(X, copy=True, nan=0.0, posinf=None, neginf=None)
    X2 = np.nan_to_num(X2, copy=True, nan=0.0, posinf=None, neginf=None)
    y = np.nan_to_num(y, copy=True, nan=0.0, posinf=None, neginf=None)

    y = y.reshape(-1,1)

    print("X shape is " + str(X.shape))
    print(X.dtype)
    print("X2 shape is " + str(X2.shape))
    print(X2.dtype)
    print("y shape is " + str(y.shape))
    print(y.dtype)

    scaler = RobustScaler()

    X = scaler.fit_transform(X)
    X2 = scaler.fit_transform(X2)
    # y = scaler.fit_transform(y)

    inputA = keras.Input(shape=X.shape[1])
    inputB = keras.Input(shape=X2.shape[1])

    x = keras.layers.Dense(4, activation="relu")(inputA)

    x = keras.layers.RepeatVector(1)(x)

    x = keras.Model(inputs=inputA, outputs=x)

    atomsNodeCount = largestAtomCount

    x2 = keras.layers.RepeatVector(3)(inputB)
    x2 = keras.layers.Reshape((x2))
    x2 = keras.layers.Conv1D(atomsNodeCount, 3, activation="relu")(x2)
    x2 = keras.layers.BatchNormalization()(x2)
    x2 = keras.layers.MaxPooling1D(2)(x2)

    while atomsNodeCount > 4:
        atomsNodeCount = atomsNodeCount // 2
        x2 = keras.layers.Conv1D(atomsNodeCount, 3, activation="relu")(x2)
        x2 = keras.layers.BatchNormalization()(x2)
        x2 = keras.layers.MaxPooling1D(2)(x2)

    x2 = keras.Model(inputs=inputB, outputs=x2)

    combined = keras.layers.concatenate([x.output, x2.output])

    x3 = keras.layers.Dense(2, activation="relu")(combined)
    x3 = keras.layers.Dense(1, activation="linear")(x3)

    model = keras.Model(inputs=[x.input, x2.input], outputs=x3)

    print(model.summary())

    plot_model(model, "modela.png", show_shapes=True, show_layer_names=True)
    # visualkeras.layered_view(model, to_file='modelb.png', legend = True)

    model.compile(optimizer = keras.optimizers.Adam(
        learning_rate=0.001, beta_1=0.9, beta_2=0.999, epsilon=1e-7, amsgrad=False, name="Adam"),
        loss = keras.losses.MeanSquaredError(reduction="auto", name="mean_absolute_error"),
        metrics = ['accuracy'])


    print("\nTraining Model...\n")

    history = model.fit(
        x= [X,X2],
        y= y,
        batch_size = 32,
        epochs = 100,
        validation_split = 0.33,
    )

    score, acc = model.evaluate(
        x=[X, X2],
        y=y
    )

    print("\nTest score = ", score)
    print("Test accuracy = ", acc)

    fig, (ax1, ax2) = plt.subplots(2, 1)

    ax1.plot(history.history['accuracy'])
    ax1.plot(history.history['val_accuracy'])
    ax1.set_title('model accuracy')
    ax1.set_ylabel('accuracy')
    ax1.set_xlabel('epoch')
    ax1.legend(['train', 'test'], loc='upper left')

    ax2.plot(history.history['loss'])
    ax2.plot(history.history['val_loss'])
    ax2.set_title('model loss')
    ax2.set_ylabel('loss')
    ax2.set_xlabel('epoch')
    ax2.legend(['train', 'test'], loc='upper right')

    fig.tight_layout()
    plt.show()

def PampNetMain():
    X, y = ReadTrainingData()
    PreprocessData(X, y)
    # model = Train()


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
