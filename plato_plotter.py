#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 20 18:59:06 2026

@author: Jonatan Rudrasingam
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.path import Path
from scipy.spatial import ConvexHull
import pandas as pd
from astropy.table import Table
from astropy.coordinates import SkyCoord
import astropy.units as u
import warnings
warnings.filterwarnings('ignore', category = RuntimeWarning)

plt.rcParams["savefig.dpi"] = 300

plt.close("all")

def wrap_coordinates(c):
    return (-c.l.wrap_at(180 * u.degree).value, c.b.value)

def rasterize_all_lines(fig=False):
    if not fig:
        figs = list(map(plt.figure, plt.get_fignums()))
    else:
        figs=[fig]
    
    for fig in figs:
        for axi in fig.axes:
            axi.set_rasterization_zorder(1)
            for linei in axi.lines:
                linei.set_zorder(linei.get_zorder()-100)

# Choose which FOV to plot
options = "south"
options = "north"

# Load the PLATO LOP (FOV)
lop = Table.read('lops2_healpix9_footprint.fits').to_pandas()
c = SkyCoord(lop['l'].values*u.deg, lop['b'].values*u.deg, frame = 'galactic')
c = wrap_coordinates(c)
lop['lwrap'] = -c[0]
lop['bwrap'] = c[1]

def plot(stars, label_stars = None, cluster_name = None, LOP = "south", 
         save_stars = None):

    if LOP == "north":
        l_diff, b_diff = 255.9375 - 81.56250, -24.62432 - 24.62432
        xlims = 48.49050572194119, 114.64548309571188
        ylims = -2.72107829269467, 52.3959750205521
        LOP_name = "LOPN1"
    else:
        l_diff, b_diff = 0, 0
        xlims = 222.8655057219412, 289.0204830957119
        ylims = -51.96971829269467, 3.1473350205520987
        LOP_name = "LOPS2"
    
    ncams = [6, 12, 18, 24]
    colors = ["lightsteelblue", "lightskyblue", "cornflowerblue", "royalblue"]
    
    # Plotting
    fs = 13
    
    # For the fast camera
    fc_sq = 1037
        
    plt.figure()
        
    for i in np.arange(len(ncams)):
        m = lop['ncam'] == ncams[i] #6, 12, 18, 24  
        points = lop[m][["l", "b"]].to_numpy()
        plt.plot(points[:,0] - l_diff, points[:,1] - b_diff, "o", markersize = 1, color = colors[i])
    
    points_shifted = lop[["l", "b"]].to_numpy()
    points_shifted[:, 0] -= l_diff
    points_shifted[:, 1] -= b_diff
    
    hull = ConvexHull(points_shifted)
    footprint_path = Path(points_shifted[hull.vertices])
    
    # Plot the FOV of the fast cameras (Note that not the entire FOV is used)
    theta = np.linspace(0, 2*np.pi, 100)
    x = np.sqrt(fc_sq/np.pi)*np.cos(theta)
    y = np.sqrt(fc_sq/np.pi)*np.sin(theta)
    plt.plot(x + (255.9375) - l_diff, y + (-24.62432) - b_diff, color = "purple")
        
    # Convert to Galactic coordinates (or alternativly supply list of stars in l and b)
    stars_l, stars_b = stars.galactic.l.degree, stars.galactic.b.degree
    stars_l = np.atleast_1d(stars_l)
    stars_b = np.atleast_1d(stars_b)
        
    star_coords = np.column_stack((stars_l, stars_b))
    is_inside = footprint_path.contains_points(star_coords)
    
    # Stars inside the LOP
    plt.plot(stars_l[is_inside], stars_b[is_inside], "*", color = "gold", 
             markersize = 15, markeredgecolor = "black")

    # Stars outside the LOP
    plt.plot(stars_l[~is_inside], stars_b[~is_inside], "*", color = "lightgray", 
             markersize = 5, alpha =0.5)
    
    # Now to insert the the name of the star(s).
    if label_stars is not None:
        labels = np.atleast_1d(label_stars)
        
        if cluster_name != None:
            if np.any(is_inside):
                avg_l = np.mean(stars_l[is_inside])
                avg_b = np.mean(stars_b[is_inside])
                plt.text(avg_l - 3.5, avg_b + 1.75, cluster_name, 
                         size = 15)
            else:
                print(f"Cluster {labels[0]} is entirely outside {LOP_name}")
        else:
            for i in range(len(stars_l)):
                if is_inside[i]:
                    plt.text(stars_l[i] - 3.5, stars_b[i] + 1.75, labels[i], size = 15)
                else:
                    label_name = labels[i] if i < len(labels) else labels[0]
                    print(f"{label_name} is not within {LOP_name}")
        
    plt.xlabel("$l \ (^\circ)$", fontsize = fs)
    plt.ylabel("$b \ (^\circ)$", fontsize = fs)
        
    plt.xlim(xlims[0], xlims[1])
    plt.ylim(ylims[0], ylims[1])
        
    plt.tick_params(axis = "both", which = "minor", labelsize = 12)
    
    rasterize_all_lines()
    
    if label_stars is not None and len(label_stars) == len(stars_l):
        names_save = label_stars
    else:
        names_save = np.arange(len(stars_l))
    
    output_table = Table()
    output_table["Name"] = names_save
    output_table["l"] = stars_l*u.deg
    output_table["b"] = stars_b*u.deg
    output_table["is_in_lop"] = is_inside
    
    # Print the table to console for a quick check
    print(output_table)
    
    if save_stars != None:
        output_table.write(f"{save_stars}", overwrite=True)
    