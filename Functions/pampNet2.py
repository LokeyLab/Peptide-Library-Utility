print("\nImporting Packages...")
# Internal
from Functions import smilesToMatricies as stm

# External
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors, rdFreeSASA
from rdkit.Chem import Crippen
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from ast import literal_eval
from tqdm import tqdm
from itertools import chain

from sklearn.preprocessing import MinMaxScaler

import tensorflow as tf
from tensorflow import keras
from keras.utils.vis_utils import plot_model
import visualkeras

import pickle

print("\nImporting Packages Complete.")

def ReadTrainingData():
    print("\nReading Training Data File...")

    df = pd.read_csv("TrainingData/hex.csv")
    X = df.iloc[0:, 10]
    y = df.iloc[0:, 14]

    # print(X.to_string())

    print("\nReading Training Data File Complete.")

    return X, y

def PreprocessData(X,y):
    print("\nProcessing Training Data...\n")

    """
    df = pd.DataFrame(columns=["One-Hot", "PappE-6"])

    periodDict = {"H": 1, "He": 2, "Li": 3, "Be": 4, "B": 5, "C": 6, "N": 7, "O": 8, "F": 9, "Ne": 10, "Na": 11,
                  "Mg": 12, "Al": 13, "Si": 14, "P": 15, "S": 16, "Cl": 17, "Ar": 18, "K": 19, "Ca": 20, "Sc": 21,
                  "Ti": 22, "V": 23, "Cr": 24, "Mn": 25, "Fe": 26, "Co": 27, "Ni": 28, "Cu": 29, "Zn": 30,
                  "Ga": 31, "Ge": 32, "As": 33, "Se": 34, "Br": 35, "Kr": 36, "Rb": 37, "Sr": 38, "Y": 39,
                  "Zr": 40, "Nb": 41, "Mo": 42, "Tc": 43, "Ru": 44, "Rh": 45, "Pd": 46, "Ag": 47, "Cd": 48,
                  "In": 49, "Sn": 50, "Sb": 51, "Te": 52, "I": 53, "Xe": 54, "Cs": 55, "Ba": 56, "La": 57,
                  "Ce": 58, "Pr": 59, "Nd": 60, "Pm": 61, "Sm": 62, "Eu": 63, "Gd": 64, "Tb": 65, "Dy": 66,
                  "Ho": 67, "Er": 68, "Tm": 69, "Yb": 70, "Lu": 71, "Hf": 72, "Ta": 73, "W": 74, "Re": 75,
                  "Os": 76, "Ir": 77, "Pt": 78, "Au": 79, "Hg": 80, "Tl": 81, "Pb": 82, "Bi": 83, "Po": 84,
                  "At": 85, "Rn": 86, "Fr": 87, "Ra": 88, "Ac": 89, "Th": 90, "Pa": 91, "U": 94, "Np": 93,
                  "Pu": 94, "Am": 95, "Cm": 96, "Bk": 97, "Cf": 98, "Es": 99, "Fm": 100, "Md": 101, "No": 102,
                  "Lr": 103, "Rf": 104, "Db": 105, "Sg": 106, "Bh": 107, "Hs": 108, "Mt": 109, "Ds": 110,
                  "Rg": 111, "Cn": 112, "Nh": 113, "Fl": 114, "Mc": 115, "Lv": 116, "Ts": 117, "Og": 118}

    atomList = []

    SMILES_CHARS = [' ',
                    '#', '%', '(', ')', '+', '-', '.', '/',
                    '0', '1', '2', '3', '4', '5', '6', '7', '8', '9',
                    '=', '@',
                    'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q',
                    'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z',
                    '[', '\\', ']',
                    'a', 'b', 'c', 'e', 'g', 'i', 'l', 'n', 'o', 'p', 'r', 's',
                    't', 'u', 'v', 'w', 'x', 'y', 'z']

    smi2index = dict((c, i) for i, c in enumerate(SMILES_CHARS))

    index2smi = dict((i, c) for i, c in enumerate(SMILES_CHARS))

    def smiles_encoder(smiles, maxlen=200):
        smiles = Chem.MolToSmiles(Chem.MolFromSmiles(smiles))
        X = np.zeros((maxlen, len(SMILES_CHARS)))
        for i, c in enumerate(smiles):
            X[i, smi2index[c]] = 1
        return X

    for x in tqdm(range(0, X.size)):

        mol = Chem.MolFromSmiles(X[x])
        # print(X[x])
        smiles = Chem.MolToSmiles(mol)

        oneHot = smiles_encoder(smiles)

        # print(smiles)

        
        # molH = Chem.AddHs(mol)

        # AllChem.EmbedMolecule(molH)

        # AllChem.MMFFOptimizeMolecule(molH)

        # morganFingerprint = AllChem.GetMorganFingerprintAsBitVect(molH, radius=2, bitInfo={})

        # molBlock = AllChem.MolToXYZBlock(molH).split("\n")[2:]

        # print(morganFingerprint)


        aL = []
        xL = []
        yL = []
        zL = []

        # print(molBlock)

        for i in range(len(molBlock)-1):
            tempList = molBlock[i].split()
            a = periodDict[tempList[0]]
            aL.append(a)
            x = float(tempList[1])
            xL.append(x)
            y = float(tempList[2])
            yL.append(y)
            z = float(tempList[3])
            zL.append(z)

        # atomList.append()
        

        df.loc[len(df.index)] = {"One-Hot": oneHot, "PappE-6": y[x]}

    print("\n" + df.to_string())

    df.to_csv("TrainingData/hexData.csv")
    """

    pickleList = []

    for i in tqdm(range(0, len(X))):
        matricies = stm.smilesToMatricies(X[i])
        pickleList.append(matricies)

    fileName = 'TrainingData/pickledPeptides/cyclicPeptideData'
    outfile = open(fileName, 'wb')
    pickle.dump(pickleList, outfile)
    outfile.close()


    print("\nProcessing Training Data Complete.")

def Train():
    fileNameX = 'TrainingData/pickledPeptides/cyclicPeptideData'
    fileNameY = 'TrainingData/pickledPeptides/cyclicPeptideYData'

    infile = open(fileNameX, 'rb')
    matricies = pickle.load(infile)
    infile.close()

    scaler = MinMaxScaler()

    X = []
    maxLen = 0
    for i in range(0,len(matricies)):
        tempX = []
        tempX.append(scaler.fit_transform(np.array(matricies[i].adjacencyMatrix)))
        tempX.append(scaler.fit_transform(np.array(matricies[i].atomicNumberMatrix)))
        tempX.append(scaler.fit_transform(np.array(matricies[i].aromaticityMatrix)))
        tempX.append(scaler.fit_transform(np.array(matricies[i].chiralityMatrix)))
        tempX.append(scaler.fit_transform(np.array(matricies[i].distanceMatrix)))

        if len(matricies[i].adjacencyMatrix) > maxLen:
            maxLen = len(matricies[i].adjacencyMatrix)

        X.append(tempX)

    # print(maxLen)

    for i in range(0, len(matricies)):
        for j in range(0, 5):
            X[i][j] = np.pad(X[i][j], pad_width=(0, maxLen - len(X[i][j])), constant_values=0)

    """
    fig, axs = plt.subplots(2, 3)


    axs[0, 1].set_title("Adjacency")
    axs[0, 1].matshow(X[0][0])

    axs[0, 2].set_title("Atomic Number")
    axs[0, 2].matshow(X[0][1])

    # axs[1,0].set_title("Formal Charge")
    # axs[1,0].matshow(formalChargeMatrix)

    axs[1, 0].set_title("Aromaticity")
    axs[1, 0].matshow(X[0][2])

    axs[1, 1].set_title("Chirality")
    axs[1, 1].matshow(X[0][3])

    axs[1, 2].set_title("3D Distance")
    axs[1, 2].matshow(X[0][4])

    plt.tight_layout()
    plt.show()
    """

    df = pd.read_csv("TrainingData/hex.csv")
    y = scaler.fit_transform(np.asarray(df["PappE-6"]).reshape(-1, 1))

    X = np.nan_to_num(X, copy=True, nan=0.0, posinf=None, neginf=None)
    y = np.nan_to_num(y, copy=True, nan=0.0, posinf=None, neginf=None)

    print("X shape is " + str(X.shape))
    print(X.dtype)
    print("y shape is " + str(y.shape))
    print(y.dtype)


    X = np.stack(X)
    X = tf.transpose(X, [0, 2, 3, 1])
    print(X.shape)

    input = keras.Input(shape=X.shape[1:])

    x = keras.layers.Conv2D(filters = 32, kernel_size = 3, input_shape = X.shape[1:])(input)
    x = keras.layers.Activation('relu')(x)
    x = keras.layers.BatchNormalization(axis=-1)(x)
    x = keras.layers.MaxPooling2D(pool_size = 2)(x)

    x = keras.layers.Conv2D(filters=64, kernel_size=3, input_shape=X.shape[1:])(x)
    x = keras.layers.Activation('relu')(x)
    x = keras.layers.BatchNormalization(axis=-1)(x)
    x = keras.layers.MaxPooling2D(pool_size=2)(x)

    x = keras.layers.Flatten()(x)

    x = keras.layers.Dense(64)(x)
    x = keras.layers.Activation('relu')(x)
    x = keras.layers.BatchNormalization(axis=-1)(x)
    x = keras.layers.Dropout(0.9)(x)

    x = keras.layers.Dense(16)(x)
    x = keras.layers.Activation('relu')(x)
    x = keras.layers.BatchNormalization(axis=-1)(x)
    x = keras.layers.Dropout(0.75)(x)

    x = keras.layers.Dense(1)(x)
    x = keras.layers.Activation('linear')(x)

    model = keras.Model(inputs=input, outputs=x)

    print(model.summary())

    plot_model(model, "modela.png", show_shapes=True, show_layer_names=True)
    visualkeras.layered_view(model, to_file='modelb.png', legend=True)

    model.compile(optimizer = keras.optimizers.Adam(learning_rate = 5e-3, decay = 1e-4, name="Adam"),
                  loss = "mse",
                  metrics = 'accuracy')

    print("\nTraining Model...\n")

    history = model.fit(
        x=X,
        y=y,
        batch_size=32,
        epochs=20,
        shuffle=True,
        validation_split=0.2,
        verbose=1
    )


    score, acc = model.evaluate(
        x=X,
        y=y
    )

    print("\nTest score = ", score)
    print("Test accuracy = ", acc)

    fig, [ax1, ax2] = plt.subplots(2)

    ax1.plot(history.history['accuracy'])
    ax1.plot(history.history['val_accuracy'])
    ax1.set_title('Accuracy')
    ax1.set_xlabel('epoch')
    ax1.legend(['train', 'test'], loc='upper right')

    ax2.plot(history.history['loss'])
    ax2.plot(history.history['val_loss'])
    ax2.set_title('Loss')
    ax2.set_xlabel('epoch')
    ax2.legend(['train', 'test'], loc='upper right')

    fig.tight_layout()
    plt.savefig("Accuracy and Loss.png")
    plt.show()

def PampNetMain():
    # X, y = ReadTrainingData()
    # PreprocessData(X, y)
    model = Train()

