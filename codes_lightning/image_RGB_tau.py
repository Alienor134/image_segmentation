import os
import copy
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import alienlab.segment
import mayavi.mlab as mlab
import alienlab.plot
g = alienlab.plot.ShowFigure()

from analyse_results_amplitude import select_bacteria
from segmentation_and_fits_update import segmentation


def RGB_tau(folder):
    colors = cm.tab20(np.linspace(0, 1, 100))

    all_bact = 0
    bact_kept = 0
    #Perform segmentatioon
    all_tau, number_segmented, segmented, im_rgb, save_path = segmentation(folder, False, False, 5, 2, 2, 3, 1, False, False)

    R = copy.copy(segmented * 0)
    G = copy.copy(segmented * 0)
    B = copy.copy(segmented * 0)

    allkey = all_tau['high'].keys() & all_tau['low'].keys()
    allkey_list = copy.copy(list(allkey))

    tau1high_list = []
    tau2high_list = []
    tau1low_list = []
    tau2low_list = []

    for k in allkey_list: #k: bacterie '1.0', '2.0' ...
        all_bact += 1
        allkey, taus_1_low, taus_2_low, taus_1_high, taus_2_high = select_bacteria(k, all_tau, allkey, 
                                                                                        1, 2)
        #threshold parameters: 
        # std < 1
        # size > 2

        if k in allkey:
            bact_kept += 1

            """
            tau bourndaries for image normalization
                1low	1high	2high
            min -2,8	min -3,2	min -2,2
            max -0,5	max -1,3	max -0,5
            """
            #update tau image
            G[segmented == int(float(k))] = (np.log10(taus_1_high.mean()) + 3.2)/(3.2-1.3)
            R[segmented == int(float(k))] = (np.log10(taus_2_high.mean()) + 2.2)/(2.2-0.7)
            B[segmented == int(float(k))] = (np.log10(taus_1_low.mean()) + 2.8)/(2.8-0.5)

            tau1high_list.append(taus_1_high.mean())
            tau2high_list.append(taus_2_high.mean())
            tau1low_list.append(taus_1_low.mean())
            tau2low_list.append(taus_2_low.mean())    


       
    print("all_bact", all_bact)
    print("batc_kept", bact_kept, 100 * bact_kept/all_bact)

    plt.close('all')

    #RGB image from tau images
    im = np.array([R,G,B])
    im = np.transpose(im, [1, 2, 0])
 #   im_rgb = 2.5 * im_rgb
 #   im_rgb[im_rgb > 1 ] = 1
 #   im = im * im_rgb[:im.shape[0],:im.shape[1]]

    g.save_folder = folder + "/RGB"
    g.cmap = 'gray'
    g.save_name = "RGB"
    g.extension = ".tiff"
    g.col_num = 2
    g.title_list = ["R", "G", "B", "RGB"]
    fig = g.multi([R, G, B, im])
    g.saving(fig)

    plt.close('all')
    fig = plt.figure()
    plt.imshow(R, cmap='Reds')
    plt.savefig(folder + "/RGB/R.tiff")
    fig = plt.figure()
    plt.imshow(G, cmap='Greens')
    plt.savefig(folder + "/RGB/G.tiff")
    fig = plt.figure()
    plt.imshow(B, cmap='Blues')
    plt.savefig(folder + "/RGB/B.tiff")

    fig = plt.figure()
    plt.imshow(im)
    plt.savefig(folder + "/RGB/RGB.tiff")

    return im, im_rgb



