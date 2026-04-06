
import numpy as np 
import pandas as pd 
import matplotlib.pyplot as plt
import seaborn as sns

from matplotlib.patches import Patch, Circle, Rectangle, Arrow, FancyArrowPatch
from matplotlib.colors import LinearSegmentedColormap, Normalize

from scipy.stats import pearsonr, spearmanr
from scipy.spatial.distance import pdist, squareform


class MantelPlot:
    '''Plot Mantel test results.'''

    def __init__(self):
        pass 

    def _color_config(self, p_value):
        '''Configure colors for plotting.'''

        cmap = LinearSegmentedColormap.from_list(
            "custom", ["#1E90FF", "#D3D3D3", "#FF6347"], N=256
        )
        norm = Normalize(vmin=-1, vmax=1)

        # p-value color mapping
        if p_value <= 0.001:
            p_color = "red"
            p_lable = "***"
        elif p_value <= 0.01:
            p_color = "orange"
            p_lable = "**"
        elif p_value <= 0.05:
            p_color = "yellow"
            p_lable = "*"
        else:
            p_color = "lightgrey"
            p_lable = "ns"

        return cmap, norm, p_color, p_lable
    
    def _line_config(self, r_value):
        pass

    def plot(
            self, 
            mantel_results_df: pd.DataFrame,
            figsize: tuple = (8, 6)
        ) -> plt.Figure:
        '''Plot Mantel test results.'''
        
        fig, ax = plt.subplots(constrained_layout=True, figsize=figsize)
