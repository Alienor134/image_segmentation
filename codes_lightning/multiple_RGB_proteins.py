import numpy as np
import matplotlib.pyplot as plt 
import os
import matplotlib
from image_RGB_tau import RGB_tau
import alienlab.plot

folder = "F:\Manips_Raja_Test/Manips"
liste_type_tapis = os.listdir(folder)

g = alienlab.plot.ShowFigure()
g.title_list = []

ims = []
ims_rgb = []

for type_tapis in ["3", "6", "11", "13", "17", "18", "19", "22"] :
        print(type_tapis)
        dir_path = folder + "/" +  type_tapis+ "/"
        liste_tapis = os.listdir(dir_path)
    
        for tapis in liste_tapis:
    
            print(tapis)
            im, im_rgb = RGB_tau(dir_path + "/" +  tapis)
            ims.append(im)
            ims_rgb.append(im_rgb)
            g.title_list.append(type_tapis)
    

ctrl = []
k = 2.5
for i , image in enumerate(ims_rgb):
     bact = ims[i]
     bis = np.copy(image)*k 
     bis[bis > 1] = 1 #thresholding the amplitude of the overlapped image of bacteria
     ctrl.append(bact * bis[:bact.shape[0],:bact.shape[1]]) #color the background image with tau false colors

g.col_num = 4
fig = g.multi(ctrl)
g.save_folder = "results_RGB"
g.save_name = "images_RGB"
g.saving(fig)

