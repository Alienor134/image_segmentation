import subprocess
import numpy as np
import matplotlib.pyplot as plt 
import os
import pickle
import matplotlib
from gooey import Gooey, GooeyParser
from alienlab import utils
import alienlab
import pandas as pd
from segmentation_and_fits_update import segmentation, reference_ellipses, fit_exponential, exp_decay

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

import mayavi.mlab as mlab
import numpy as np
import matplotlib.cm as cm
from matplotlib import rcParams
matplotlib.rcParams['font.sans-serif'] = "Arial"
matplotlib.rcParams['font.family'] = "sans-serif"


import scipy.stats

import json
import glob
from copy import copy

rcParams['mathtext.fontset'] = 'custom'
rcParams['mathtext.it'] = 'Arial:italic'
rcParams['mathtext.rm'] = 'Arial'
rcParams['mathtext.bf'] = 'Arial:bold'




@Gooey(required_col = 3, optional_cols = 5, default_size=(900, 600))

def select_bacteria(k, all_tau, allkey, threshold, bacteria_size):
    taus_1_low = np.array(all_tau['low'][k]['0']) #list tau_1_low of bacteria k
    taus_2_low = np.array(all_tau['low'][k]['1']) #list of tau_2_low of bacteria k
    taus_1_high = np.array(all_tau['high'][k]['0']) #list of tau_1_high of bacteria k
    taus_2_high = np.array(all_tau['high'][k]['1']) #liste of tau_2_high of bacteria k

    if taus_1_low.shape[0]>bacteria_size and taus_1_high.shape[0]>bacteria_size and  taus_2_high.shape[0]>bacteria_size and taus_2_low.shape[0]>bacteria_size:
# measure used for thresholding, the range of the values reached by tau reduced by the mean: 
        v_1_low  = -(np.log10(taus_1_low.max())-np.log10(taus_1_low.min()))/np.log10(taus_1_low.mean())
        v_1_high = -(np.log10(taus_1_high.max())-np.log10(taus_1_high.min()))/np.log10(taus_1_high.mean())  
        v_2_high = -(np.log10(taus_2_high.max())-np.log10(taus_2_high.min()))/np.log10(taus_2_high.mean())  
        # thresholding
        if np.abs(v_1_high) < threshold and np.abs(v_2_high) < threshold and np.abs(v_1_low ) < threshold: 

    
            if np.log10(taus_1_low.mean()) > 0:
                allkey.remove(k)
                print("POSITIVE: removing bacteria", k)
                
            elif np.log10(taus_1_high.mean()) >  0:
                allkey.remove(k)
                print("POSITIVE: removing bacteria", k)

            elif np.log10(taus_2_high.mean()) >  0:
                allkey.remove(k)
                print("POSITIVE: removing bacteria", k)
            
            elif np.log10(taus_2_low.mean()) > 4:
                allkey.remove(k)
                print("POSITIVE: removing bacteria", k)

            else:
                pass
        
        else:
            allkey.remove(k)    
            print("THRESHOLD: removing bacteria", k)
        
    else:
        allkey.remove(k)
        print("TOO SMALL: removing bacteria", k)                        
    
    return(allkey, taus_1_low, taus_2_low, taus_1_high, taus_2_high)        


def main():
    parser = GooeyParser()

    #parser.add_argument("folder", help = "Name of the folder you want to process", widget="DirChooser")

    parser.add_argument("-c", "--control", action="store_true", default = True,
                        help = "Results for the control")
    parser.add_argument("-m", "--mix", action="store_true", default = True,
                        help = "Results for the mix")

    parser.add_argument("-note", "--annotate", action="store_true", default = False,
                            help = "Manually select the ROI if checked")

    parser.add_argument("-t", "--threshold", type = float, default = 0.3,
                            help = "Threshold for data selection")
    parser.add_argument("-s", "--bacteria_size", type = int, default = 3,
                            help = "Minimal size of bacteria")

    args = parser.parse_args()

    result_name =  "_threshold_%.02f_size_%d"%(args.threshold, args.bacteria_size)
    #path = args.folder

    path = "F:\Manips_Raja_Test\Manips"
    liste_type_tapis = os.listdir(path)

    print(liste_type_tapis)
    dic = {}
    mean_val = {}
    colors = cm.tab20(np.linspace(0, 1, 30))
    for type_tapis in liste_type_tapis:

            dir_path = os.path.join(path, type_tapis) + "/"
            liste_tapis = os.listdir(dir_path)
            dic[type_tapis] = {}
            mean_val[type_tapis] = {}
            for tapis in liste_tapis:
                print(type_tapis, tapis)
                f = glob.glob(dir_path + tapis + "/*_tau_val.pkl")
                f = open(f[0],'rb')
                dict_item = pickle.load(f)
                dic[type_tapis][tapis] = dict_item
                mean_val[type_tapis][tapis]={}
                data_path_high = glob.glob(dir_path + tapis + "/*high_bin*.csv")
                for filep in data_path_high:
                    if "coords" not in filep:                    
                        mean_val[type_tapis][tapis]["high"] = pd.read_csv(filep)
                data_path_low = glob.glob(dir_path + tapis + "/*low_bin*.csv")
                for filep in data_path_low:
                    if "coords" not in filep:                    
                        mean_val[type_tapis][tapis]["low"] = pd.read_csv(filep)

    plt.close("all")

    statistics = {}

    c = 0
    for type_tapis in dic.keys(): #protein type, or mix 


            print(type_tapis)
            p = alienlab.plot.PlotFigure() #tool to plot graphs and multiple curves
            p.figsize = (7,7)
            p.fontsize = 16
            p.save_folder = "results" + result_name + "/"
            p.date = False
            p.extension = ".pdf"

            if type_tapis == "1":

                fig = plt.figure(figsize = (3,3))
                p.save_name = "test"
                p.saving(fig)
            framerate = {}
            limit = {}
            framerate['high'] = 1000
            framerate['low'] = 1000
            limit['high'] = [0, 200, 400, 600]
            limit['low'] = [0, 800, 1600, 2400]
            all_tau_fit = {}
                
            dic_t = dic[type_tapis]
            statistics[type_tapis] = {}
            statistics[type_tapis]["tau1low_list"] = []
            statistics[type_tapis]["tau2low_list"] = []
            statistics[type_tapis]["tau1high_list"] = []
            statistics[type_tapis]["tau2high_list"] = []
            statistics[type_tapis]["amplitude1low_list"] = []
            statistics[type_tapis]["amplitude2low_list"] = []
            statistics[type_tapis]["amplitude1high_list"] = []
            statistics[type_tapis]["amplitude2high_list"] = []
            statistics[type_tapis]["taumean1low_list"] = []
            statistics[type_tapis]["taumean2low_list"] = []
            statistics[type_tapis]["taumean1high_list"] = []
            statistics[type_tapis]["taumean2high_list"] = []



            bact_tot = 0
            bact_kept = 0
            for ta in dic_t.keys(): #ta : ta1, ta2, ta3 ...

                print(ta)
                all_tau = dic_t[ta]
                allkey = all_tau['high'].keys() & all_tau['low'].keys()
                allkey_list = copy(list(allkey))
                for k in allkey_list: #k: bacterie '1.0', '2.0' ...
                    
                    bact_tot += 1
                    allkey, taus_1_low, taus_2_low, taus_1_high, taus_2_high = select_bacteria(k, all_tau, allkey, 
                                                                                    args.threshold, args.bacteria_size)

                    if k in allkey:
                            bact_kept += 1
                            statistics[type_tapis]["tau1high_list"].append(taus_1_high.mean())
                            statistics[type_tapis]["tau2high_list"].append(taus_2_high.mean())
                            statistics[type_tapis]["tau1low_list"].append(taus_1_low.mean())
                            statistics[type_tapis]["tau2low_list"].append(taus_2_low.mean())


                try:
                    allkey.remove('0.0') #remove background
                except:
                    pass

                p.save_folder = "results"+ result_name +"/mean_fits/" + type_tapis + "/" + ta
                for film in ["high", "low"]:
                    all_tau_fit = fit_exponential(mean_val[type_tapis][ta][film], limit[film], framerate[film], 
                                                    all_tau_fit, film, allkey, p, showit = False, function_to_fit = exp_decay, amplitude = True)


                    statistics[type_tapis]["amplitude1" + film +"_list"] += [all_tau_fit[film][k][str(0)][0][0] for k in allkey]
                    statistics[type_tapis]["amplitude2" + film +"_list"] += [all_tau_fit[film][k][str(1)][0][0] for k in allkey] 
                    statistics[type_tapis]["taumean1" + film +"_list"] += [all_tau_fit[film][k][str(0)][0][1] for k in allkey]
                    statistics[type_tapis]["taumean2" + film +"_list"] += [all_tau_fit[film][k][str(1)][0][1] for k in allkey]

            protein = statistics[type_tapis]
            list_keys = list(protein.keys())
            hist_range = {'tau1low_list':(-2.5, -0),
                            'tau2low_list':(-3,3),
                            'tau1high_list':(-2.7, -1.3),
                            'tau2high_list':(-2, -1),
                            'amplitude1low_list':(1, 5),
                           'amplitude2low_list':(1, 5),
                            'amplitude1high_list':(1, 5),
                            'amplitude2high_list':(1, 5),
                            'taumean1low_list':(-2.5, -0),
                            'taumean2low_list':(-3,3),
                            'taumean1high_list':(-2.7, -1.3),
                            'taumean2high_list': (-2, -1)
                    }
            list_key = copy(list(protein.keys()))
            
            correspondance = {"tau1low_list": "1", "tau2low_list": "2", "tau1high_list": "3", "tau2high_list": "4"}
            correspondance_bis = {"tau1low_list": "1Ilow", "tau2low_list": "2IIlow", "tau1high_list": "3Ihigh", "tau2high_list": "4IIhigh"}
            
            for tau_i in ["tau1low_list", "tau2low_list", "tau1high_list", "tau2high_list"]:
                tau_list = np.nan_to_num(np.log10(np.array(protein[tau_i])))
                print(type_tapis, tau_i, tau_list)
                fitfunc  = lambda p, x: p[0]*exp(-0.5*((x-p[1])/p[2])**2)+p[3]
                errfunc  = lambda p, x, y: (y - fitfunc(p, x))

                init  = [1.0, 0.5, 0.5, 0.5]

                fig = plt.figure(figsize = (10,10))
                ax = plt.gca()
                p.save_folder = "results"+ result_name +"/distribution/Bacteria_%02d"%int(type_tapis) + "/Bacteria_%02d"%int(type_tapis) + "_" + correspondance_bis[tau_i] + "/"
                p.save_name = "tau_" + correspondance[tau_i] + "_bacteria_%02d"%int(type_tapis) 
                plt.xlabel(r"$\mathit{\ell_{" + correspondance[tau_i] + "\mathbf{" + type_tapis +"}}}$", fontsize = 36)
                plt.ylabel("Bacteria number", fontsize = 36)
                # We change the fontsize of minor ticks label 
                ax.tick_params(axis='both', which='major', labelsize=24, direction = 'in', top = True, right = True, length = 10 )
                _, bins, _ = plt.hist(tau_list, 40, density= False, range = hist_range[tau_i], alpha=1, facecolor = "white", edgecolor = "black")
                
                
                #mu, sigma = scipy.stats.norm.fit(tau_list)
                #best_fit_line = scipy.stats.norm.pdf(bins, mu, sigma)
                #print(type_tapis, len(tau_list))
                #D, pval = scipy.stats.kstest(tau_list, "norm")
                #plt.plot(bins, best_fit_line)
                #plt.title("mu = %f, sig = %f, D = %f, p = %f"%(mu, sigma, D, pval))

                p.saving(fig)
                plt.close('all')
                #protein[tau_i[:-5] + "_distribution"] = [mu, sigma, D, pval] 
                protein[tau_i[:-5] + "_mean"] = np.mean(tau_list)
                protein[tau_i[:-5] + "_std"] = np.std(tau_list)
            """
            for pos in ['1high', '2high', '1low', '2low']:
                plt.close('all')
                fig = plt.figure()
                p.xval = protein["taumean" + pos + '_list']
                p.yval = protein["amplitude" + pos + "_list"]
                print(type_tapis, len(p.xval), len(p.yval))
                p.xlabel = "tau" + pos + "_list"
                p.ylabel = "amplitude" + pos + "_list"
                p.xlim(hist_range["taumean" + pos + '_list'])
                p.ylim(hist_range["amplitude" + pos + '_list'])
                plt.scatter(np.log10(p.xval), np.log10(p.yval))
                p.save_name = "correlation_" + pos
                p.saving(fig)
            """
            statistics[type_tapis]["bacteria_kept"] = bact_kept
            statistics[type_tapis]["bacteria_tot"] = bact_tot
            statistics[type_tapis]["bacteria_ratio"] = bact_kept/bact_tot * 100
            
            mlab.points3d(np.log10(np.array(protein["tau1low_list"])),
                         np.log10(np.array(protein["tau1high_list"])),
                         np.log10(np.array(protein["tau2high_list"])), color = tuple(colors[c][:-1]), scale_factor=.02, name = type_tapis) 

            c+=1
            multivariate_distribution = np.array([np.log10(protein[it_key+'_list']) for it_key in ['tau1low', 'tau2low', 'tau1high', 'tau2high']])
            #statistics[type_tapis]["covariance"] = np.cov(multivariate_distribution)

    df = pd.DataFrame()
    for protein in statistics.keys():
        for stat_val in statistics[protein].keys():
            if "list" not in stat_val and "distribution" not in stat_val and "covariance" not in stat_val:
                df.loc[protein, stat_val] = statistics[protein][stat_val]
            
    df.to_csv("results" + result_name  + "/results.csv", sep = ",")
    mlab.savefig("3D_point_cloud.wrl")
      

    mlab.show()
    
    pickle.dump(statistics, open('statistics_extracted.pkl', 'wb'))
    with open("results" + result_name + "/data_3D.json", 'w') as fp:
        json.dump(statistics, fp)

    return statistics

if __name__ == '__main__':
    statistics = main()

