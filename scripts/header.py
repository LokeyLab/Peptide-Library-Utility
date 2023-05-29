#!/usr/bin/python
# -*- coding: utf-8 -*-

# Standard Libraries
import datetime
import itertools
import os.path
import PIL
from PIL import ImageTk, Image
import random
import sys

import tkinter as tk
from tkinter import ttk
from tkinter import messagebox
from tkinter import filedialog
from tkinter import font

import webbrowser

# External Libraries

import pathlib2

import platform

import matplotlib
platform = platform.system().lower()
if platform == 'darwin':
    matplotlib.use("TkAgg")

import matplotlib.pyplot as plt

from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

import numpy as np

import pandas as pd


import rdkit
from rdkit.Chem import Descriptors, Draw, rdMolDescriptors, AllChem, RemoveHs, AddHs, Crippen, \
    MolFromSmiles, PandasTools, rdFreeSASA

import tqdm

from scipy import stats

# Internal Libraries
from scripts.a_log_p import a_log_p_86
from scripts import bonding
from scripts import cheminformatics
from scripts import classes
from scripts import combinatronics
from scripts import gui
from scripts.gui import gui1, gui_root, gui_functions
from scripts.mcgowan_volume import mcgowan_volume
from scripts import setup
from scripts import subunit_builder
from scripts import tbui
from scripts import utilities

# Global Variables
global SUBUNIT_LIBRARY_POPULATED
SUBUNIT_LIBRARY_POPULATED = False

global SUBUNIT_LIBRARY
SUBUNIT_LIBRARY = {}

global SUBUNIT_STRING_DATA
SUBUNIT_STRING_DATA = []

global SUBUNIT_LIST
SUBUNIT_LIST = []

global SUBUNIT_COUNT
SUBUNIT_COUNT = 0

global INPUT_FILE
INPUT_FILE = ''

global OUTPUT_DIR
OUTPUT_DIR = ''

global OUTPUT_TYPES
OUTPUT_TYPES = []

global OUTPUT_IMG_TYPES
OUTPUT_IMG_TYPES = []

global CHEMINFORMATICS_OPTIONS
CHEMINFORMATICS_OPTIONS = []

global QUOTE_LIST
QUOTE_LIST = [
    "A house divided against itself cannot stand.",
    "Publish or Perish.",
    "Starting today, I decided to stop simping.",
    "Single life for me because I chose to be a cheater.",
    "The food of my people!!",
    "Like a ship passing in the night.",
    "DTF calculations.",
    "Baaabcockkkk.",
    "...Is that Comic Sans?",
    "Eh, I'ma go get high and watch 300.",
    "You wanna live free... Or die hard?",
    "...A lot of dudes here.",
    "We can't take Pan Pan in there.",
    "City Boys up!!",
    "Hot Girl Summer!!",
    "That's graphing 101.",
    "Alex doesn't look big.",
    "If you aren't grinding, what the fuck are you doing?",
    "We chiefin'??",
    "I understand him... But I don't respect him.",
    "We won.",
    "Patent a drug or fuck off.",
    "Real eyes realize real lies.",
    "Lab is church and I am your god.",
    "Weak men don't have a choice.",
    "King of the Brobot.",
    "I hate America.",
    "Chastity is the flowering of a man.",
    "John Timstone",
    "Measure once, cut twice.",
    "KRONK.",
    "Enzymes... They don't give a shit!",
    "What if I upload a picture of me and my glock?",
    "Smooth seas never made a smooth sailor.",
    "HE LIES STILL.",
    "FLEET FOCKETS.",
    "I wish I was a cat.",
    "Six months in the lab can save you one hour in the library.",
    "You have some chemistry background.",
    "Google is shit.",
    "It's the year of the DOG.",
    "I'm a DRC TA.",
    "Try my thigh.",
    "Bro I get so sad sometimes.",
    "Forced vaccinations are fascist.",
    "LIFE HAS NO LIMITS.",
    "Jesus is the son of God.",
    "Moon landing = Fake.",
    "Don't fear MAN, only fear God.",
    "Please keep the flies away from mommy.",
    "God has a plan.",
    "The age of colorblindness.",
    "your mental health doesn't matter.",
    "Thank God for ME.",
    "Life's a bitch and then you die.",
    "Without emotion, ANYTHING IS POSSIBLE.",
    "It's time to weed out the weak.",
    "CO Habitation.",
    "Oklahoma ostrich farm.",
    "If I cared, you would know.",
    "Sometimes you gotta get on your knees to stand on your feet.",
    "I'd rater die on my feet than lvie on my knees.",
    "Fedz watchin'.",
    "The nose always knows.",
    "Fun's over.",
    "OUR FUTURE.",
    "I want to see some children die.",
    "You're either IN or you're OUT.",
    "Cognac Boyz",
    "God's REAL.",
    "Bye bye Beef Boy.",
    "God Guns America.",
    "Demitrius is Mary-Ann now.",
    "Those knees look pretty BENDABLE.",
    "It was never me.",
    "It's hard to stay mad when you're winning.",
    "My word is my pride.",
    "Where's the lamb??",
    "RIP",
    "Let's crack some 40's and GET IT ON.",
    "End to end WE'RE ALL THE SAME",
    "Summer's over and it's gonna be a LONG WINTER.",
    "She wants a lot of things.",
    "I'm saving myself for MY WIFE.",
    "Just because you felt something doesn't mean it's real.",
    "They don't think it be like it is, but it do."
             ]

global CAP_AMINE
CAP_AMINE = False

global CAP_C_TERMINUS
CAP_C_TERMINUS = False