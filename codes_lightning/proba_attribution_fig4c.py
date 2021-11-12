from gooey import Gooey, GooeyParser
import pickle
import os
import copy
import numpy as np
import matplotlib.cm as cm
import alienlab.segment
import imageio
import mayavi.mlab as mlab
import json
from tkinter.filedialog import askdirectory
import matplotlib.pyplot as plt
import pandas as pd
from analyse_results_amplitude import select_bacteria
from mpl_toolkits.axes_grid1 import make_axes_locatable

from mpl_toolkits.axes_grid1 import ImageGrid


import glob

from segmentation_and_fits_update import segmentation

from matplotlib import rcParams

rcParams['mathtext.bf'] = 'Arial:bold'

control_stats = pickle.load( open( "statistics_extracted.pkl", "rb" ) )

colors = cm.tab20(np.linspace(0, 1, 30))

control = ["1","3","6", "11","13", "14", "15", "17", "20"] 

control_tapis = np.load('controle_tapis.npy')


# modelisation of the distribution with the hard threshold data
mean_distribution = np.stack([np.array([control_stats[type_tapis]["tau1low_mean"],
                                #control_stats[type_tapis]["tau2low_mean"],
                                control_stats[type_tapis]["tau1high_mean"],
                                control_stats[type_tapis]["tau2high_mean"]]) for type_tapis in control])

cov_distribution = np.stack([np.linalg.pinv(control_stats[type_tapis]["covariance"]) for type_tapis in control])


all_bact = 0
bact_kept = 0
attribution = {}
label_correspondance = {}
for c in control:
    attribution[c] = 0

remaining_images = {}

folders = glob.glob("G:/LIGHTNING_C/review/tapis_new/*")

for f in folders: 
    folder = glob.glob(f + "/*")[0]
  
    sp = os.path.split

    protein_name = sp(sp(folder)[0])[1]
    name = protein_name + '_' + sp(folder)[1]
    if True: #if the images have already been segmented with segmentation_and_fits_update.py
        segmented = np.load("proba_matrix/%s_segmented.npy"%name)
        all_tau_file = glob.glob(folder + "/*tau_val.pkl")[0]
        with open(all_tau_file, 'rb') as f:
           all_tau = pickle.load(f)
    else:
        all_tau, numbered_segmented, segmented, im_rgb, save_path = segmentation(folder, False, False, 5, 2, 2, 3, 1, False, False)
        np.save("proba_matrix/%s_segmented.npy"%name, segmented)


    im = copy.copy(segmented * 0) #image to fill
    ims = {}
    attrs = {}
    for type_tapis in control:
        ims[type_tapis] = copy.copy(im)
        attrs[type_tapis] = np.stack([copy.copy(im)]*3)

    allkey = all_tau['high'].keys() & all_tau['low'].keys()
    allkey_list = copy.copy(list(allkey))

    tau1high_list = []
    tau2high_list = []
    tau1low_list = []
    tau2low_list = []

    color_prot = {"1": (31/255,119/255,1),
                    "3": (0, 1, 1),
                    "6": (1, 0, 0),
                    "11": (1,0,2/3),
                    "13": (0, 1, 0),
                    "14": (255/250,102/255,0),
                    "15": (127/250,127/250,1),
                    "17": (1, 1, 0),
                    "20": (1,119/250,1)}


    dataframe = pd.DataFrame(index = control)
    tauframe = pd.DataFrame(index = ['tau1low','tau1high','tau2high'])

    for k in allkey_list: #k: bacterie '1.0', '2.0' ...
        all_bact += 1
        #apply selection threshold
        allkey, taus_1_low, taus_2_low, taus_1_high, taus_2_high = select_bacteria(k, all_tau, allkey, 
                                                                                        1, 2)
        #if bacteria has not been removed
        if k in allkey: 
            bact_kept += 1
            bacteria_mix = np.log10(np.stack([taus_1_low.mean(),# taus_2_low.mean(),
                                                taus_1_high.mean(), taus_2_high.mean()], axis = 0))
            tauframe[k]=bacteria_mix
            
            #compute probability
            ini = np.expand_dims(bacteria_mix, 0)
            a = ini - mean_distribution
            a = np.expand_dims(a, 1)
            b = np.transpose(a, (0, 2, 1))
            p = np.exp(-0.5*(a@cov_distribution@b))

            #probablity map for each bacteria
            for i, type_tapis in enumerate(control):
                ims[type_tapis][segmented == int(float(k))] = float(p[i])

            #class attribution
            dataframe[k] = p[:,0,0]
            ind_protein = np.argmax(p, axis = 0)[0][0]
            attribution[control[ind_protein]] += 1
            #coloring of the bacteria with label color and 0-1 intensity for 0-1 probability
            for channel in range(3):
                attrs[control[ind_protein]][channel][segmented == int(float(k))] = float(p[ind_protein])*color_prot[control[ind_protein]][channel]


            tau1high_list.append(taus_1_high.mean())
            tau2high_list.append(taus_2_high.mean())
            tau1low_list.append(taus_1_low.mean())
            tau2low_list.append(taus_2_low.mean())        

        else: 
            dataframe[k] = [0]*len(control)
            bacteria_mix = np.log10(np.stack([taus_1_low.mean(),# taus_2_low.mean(),
                                                taus_1_high.mean(), taus_2_high.mean()], axis = 0))
            tauframe[k]=bacteria_mix



    col = 3

    # Build the figure with the images
    f = plt.figure(figsize=(12, 10))

    grid = ImageGrid(f, 111, 
                    nrows_ncols=(3,3),
                    axes_pad=0.15,
                    share_all=True,
                    cbar_location="right",
                    cbar_mode="single",
                    cbar_size="4%",
                    cbar_pad=0.1,
                    )   
                    
    for i, type_tapis in enumerate(control):
            axi = grid[i]
            imageprob = ims[type_tapis]
            imcb = axi.imshow(imageprob, cmap='gray', vmin=0, vmax=1)
            axi.text(5, 25, r"$P_{%d,%d}$"%(int(type_tapis), int(protein_name)),color='white', fontsize=17)
            axi.axis("off")
    
    # Colorbar
    axi.cax.colorbar(imcb)
    axi.cax.toggle_label(True)



    dataframe.to_csv('proba_matrix/%s_proba.csv'%name)
    tauframe.to_csv('proba_matrix/%s_tau.csv'%name)
    plt.savefig("proba_matrix/%s_proba.pdf"%name)

    with open('proba_matrix/%s_attrib_soft.pkl'%name, 'wb') as output:
        pickle.dump(attrs, output)

    f, ax = plt.subplots(3,col, figsize = (11,10) ) 
    for i, type_tapis in enumerate(control):
        axi = ax[i//col][i%col]
        axi.imshow(attrs[type_tapis].transpose(1,2,0))
        axi.text(5, 25, r"$P_{%d,%d}$"%(int(type_tapis), int(protein_name)),color='white', fontsize=17)
        axi.axis("off")
    f.tight_layout()

    remaining_images[name] = attrs



    dataframe.to_csv('proba_matrix/%s_proba.csv'%name)
    tauframe.to_csv('proba_matrix/%s_tau.csv'%name)
    plt.savefig("proba_matrix/%s_attrib.pdf"%name)


with open('proba_matrix/remaining_images_soft.pkl', 'wb') as output:
    pickle.dump(remaining_images, output)
    
f, ax = plt.subplots(3,col, figsize = (11,10) )
for i, tapis in enumerate(control):
    ta = remaining_images[control_tapis[tapis]]
    im_col = np.stack([copy.copy(im)]*3)
    for k in ta.keys():
        im_col += ta[k]
    axi = ax[i//col][i%col]
    axi.imshow(im_col.transpose(1,2,0))
    axi.text(5, 25, r"$\bf{" + list(ta.keys())[i] + "}$",color='white', fontsize=25)
    axi.axis("off")
plt.tight_layout()
plt.savefig("proba_matrix/%s_criteres.pdf"%name)

