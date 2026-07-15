#Imported Modules
import numpy as np
import os
import sys
from scipy import stats as stats
import astropy.units as units
from matplotlib import pyplot as plt
import importlib
import vaex, vaex.ml, vaex.ml.cluster
import cmasher as cml
import yaml
import json
from matplotlib import colors
import pickle
from numba import jit
import seaborn as sns
import itertools

from scipy.cluster.hierarchy import linkage
from scipy.cluster.hierarchy import dendrogram
from scipy.cluster.hierarchy import fcluster

import KapteynClustering as kc
from KapteynClustering.Main import main as kc_main
import KapteynClustering.dynamics_funcs as dynf
import KapteynClustering.data_funcs as dataf
import KapteynClustering.param_funcs as paramf
import KapteynClustering.cluster_funcs as clusterf
from KapteynClustering.plot_funcs import plot_simple_scatter
from KapteynClustering.plotting_utils import plotting_scripts as ps 
from KapteynClustering import grouping_funcs as groupf
from KapteynClustering.plotting_utils.dendo_plotting import dendo_plot
from KapteynClustering.plot_funcs import group_plot


os.chdir('/cosma/home/durham/dc-coll7/summerproject') #Forces notebook to operate in it's own directory so it can load other local .py files
