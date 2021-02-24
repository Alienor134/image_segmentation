import subprocess
import numpy as np
import matplotlib.pyplot as plt 
import os
import pickle
import matplotlib
from gooey import Gooey, GooeyParser

from segmentation_and_fits_update import segmentation


@Gooey(optional_cols = 5, dump_build_config = True, default_size=(900, 600))
def main():
    parser = GooeyParser()


    parser.add_argument("folder", help = "Name of the folder you want to process", widget="DirChooser")
    #parser.add_argument("control", action="store_true", default = True,
    #                    help = "Run segmentation and fit on control")
    #parser.add_argument("mix", "--mixture", action="store_true", default = True,
    #                    help = "Run segmentation and fit on mix")

    parser.add_argument("-al", "--autolevel", type = int, default = 5, #2
                        help = "Size of structure element for autolevel")
    parser.add_argument("-c", "--contrast", type = int, default = 2, #10
                        help = "Size of structure element for contrast")
    parser.add_argument("-dd", "--disk_size", type = int, default = 2,
                        help = "Size of structure element for dilation of the mask")

    parser.add_argument("-mc", "--max_contrast", type = int, default = 3,
                        help = "Size of structure element for selection of max contrast")

    parser.add_argument("-sh", "--soft_hard", type = float, default = 1,
                        help = "Threshold")
                    

    parser.add_argument("-dm", "--dist_max", action="store_true",
                        help = "Use distance map for local maxima computation")
    parser.add_argument("-ds", "--dist_seg", action="store_true",
                        help = "Use distance map for watershed computation")      
                        
    parser.add_argument("-showit", "--showit", action="store_true", default = True,
                        help = "Plot intermediary steps to debug if checked")

    parser.add_argument("-ROI", "--select_roi", action="store_true", default = False,
                            help = "Manually select the ROI if checked")


    args = parser.parse_args()


    liste_type_tapis = os.listdir(args.folder)#[x[0] for x in os.walk(path)]

    for type_tapis in liste_type_tapis:
            print(type_tapis)
            dir_path = args.folder + "/" +  type_tapis+ "/"
            liste_tapis = os.listdir(dir_path)#[x[0] for x in os.walk(dir_path )]
        
            for tapis in liste_tapis:
        
                print(tapis)
                segmentation(dir_path + "/" +  tapis, args.dist_max, args.dist_seg, args.autolevel, args.contrast, args.disk_size, args.max_contrast, args.soft_hard, args.showit, args.select_roi)
        


if __name__ == '__main__':
    main()