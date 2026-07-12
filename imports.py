#Imported Modules
import numpy as np
import os
import sys
from scipy import stats as stats
import astropy.units as units
from matplotlib import pyplot as plt
import importlib
import KapteynClustering as kc
import vaex, vaex.ml, vaex.ml.cluster
import cmasher as cml

os.chdir('/cosma/home/durham/dc-coll7/summerproject') #Forces notebook to operate in it's own directory so it can load other local .py files
