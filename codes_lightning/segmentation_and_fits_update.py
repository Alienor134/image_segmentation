
import alienlab.plot
from alienlab.improcessing import normalize, grey_to_rgb, make_binary
import alienlab.segment
from alienlab.fo import FramesOperator
import alienlab.utils
import matplotlib.pyplot as plt
import matplotlib

from tkinter.filedialog import askopenfilename

import copy
import time
import os
import numpy as np
import skimage.io
import imageio
import subprocess
import cv2
import pandas as pd
from scipy import optimize
import pickle
import glob
import matplotlib.cm as cm


def exp_decay(parameters, xdata):
    '''
    Calculate an exponential decay of the form:
    S = a * exp(-xdata/b)
    '''
    A = parameters[0]
    tau = parameters[1]
    y0 = parameters[2]
    return A * np.exp(-xdata/tau) + y0


def residuals(parameters,x_data,y_observed,func):
    '''
    Compute residuals of y_predicted - y_observed
    where:
    y_predicted = func(parameters,x_data)
    '''
    return func(parameters,x_data) - y_observed

def fit_exponential(items_dict_binned, lim, frate, all_tau, film, keys, p, showit, function_to_fit, amplitude = False):
    '''
    fit exponential curve to high and low regimes of light pulses. It's a two step fit, the second fit is performed on a decay range of 5*tau with tau predicted by the first fit. 
    '''

    all_tau[film] = {}
    key_list = copy.copy(keys)
    for k in key_list: #k corresponds to bacteria number k
        all_tau[film][k] = {}
        p.figsize = (12,12)
        for j in range(len(lim)-2): #j corresponds to slope number. j = 0: tau_1, j = 1: tau_2
            all_tau[film][k][str(j)] = []
            if showit:
                plt.close('all')
                fig = plt.figure()
            p.save_name = k + '_' + film + '%d'%j
            p.figsize = (10, 10)
            p.fontsize = 36
            p.fonttick = 24
            p.title = ""
            p.xlabel = "t(s)"
            p.ylabel = "Fluorescence intensity (a.u.)"
            p.extension = ".pdf"
            p.label_list = []
            p.xval = []
            p.yval = []
            p.color_list = []
            #collect data to fit
            try: 
                data =  items_dict_binned[k]['pixel_values']
            except:
                data =  items_dict_binned[k]
                data = np.expand_dims(data, 1)

            for i in range(data.shape[1]): #i corresponds to a pixel in the bacteria
                if j == 1:# and film == "low":
                    start = lim[j]
                    stop = lim[j]+200
                #elif j == 1 and film == "high":
                #    start = lim[j] + np.argmax(data[lim[j]: lim[j] + (lim[j+1] - lim[j])//2,i]) #the rise of tau_2_high is not accessible so we focus on the 
                    #following exponential decrease, by identifying the local maximum
                #   stop = lim[j+1]
                else:
                    start = lim[j]
                    stop = lim[j+1]

                x = np.linspace(0, stop - start, stop - start)

                x0 = [1e5, (stop - start)/8, 1e2] #initial values for the fit (A, tau, k)
                data_sample = data[start: start + (stop - start)//2,i] #select slope of interest depending on start and stop values

                #first fit
                OptimizeResult  = optimize.least_squares(residuals,  x0, bounds = (-1e9,1e9),
                                                    args = (x[0:(stop - start)//2], data_sample, function_to_fit))
                parameters_estimated = OptimizeResult.x
                tau = parameters_estimated[1]
                #conditions on tau too low or too high for the second, more accurate, fit, because we will fit on a signal that lasts 5*tau
                if tau >  (stop - start)//10: #if too high
                    tau =  (stop - start)//10
                if tau < 3: #if too low, increase it
                    tau = 5
                x0 = parameters_estimated #initial guess: parameters from previous fit
                #second fit
                OptimizeResult  = optimize.least_squares(residuals,  x0, bounds = (-1e9,1e9),
                                                    args = (x[0:int(tau*5)], data[start: start + int(tau*5),i], function_to_fit))
                parameters_estimated = OptimizeResult.x
                J =  OptimizeResult.jac
                cov = np.linalg.pinv(J.T.dot(J))
                sig = np.sqrt(np.abs((np.diagonal(cov))))

                y = function_to_fit(parameters_estimated, x[0: int(tau * 5)])
 
                # condition on fit variance to discard or not the signal of the pixel: 
                #if 3*sig[0] < np.abs(parameters_estimated[0]) and 3*sig[1] < parameters_estimated[1]  and tau > 0:
                if True:
                    if amplitude == True: 
                        all_tau[film][k][str(j)].append((parameters_estimated[0], parameters_estimated[1]/frate))

                    else: 
                         all_tau[film][k][str(j)].append(parameters_estimated[1]/frate)

                    if showit:
                        p.yval.append(data[start:stop - 5, i])
                        p.yval.append(y)

                        p.xval.append(x[0:(stop-start - 5)]/frate)
                        p.xval.append(x[0:int(tau*5)]/frate)

                        p.label_list.append("Pixel " + str(i))
                        p.label_list.append("Fit " + str(i) + ", tau = %0.3f"%np.log10(parameters_estimated[1]/frate) )
                        #plt.plot(x[0:int(tau*5)], data[start:start + int(tau*5), i], label = np.log10(parameters_estimated[1]/frate))

                        #plt.plot(x[0:int(tau*5)], y)
                else:
                    try: 
                        keys.remove(k)
                    except:
                        print("bacteria already removed")
            if showit: 
                p.legend = False
                p.color_list = cm.tab20(np.tile(np.linspace(0, 1, 20), 10))[::-1]
            

                fig, df, ax = p.plotting(p.xval, p.yval)
                p.xlabel = "t (s)"
                p.ylabel = "Fluorescence (a.u)"
                p.minor_ticks = False
                p.saving((fig, df))
    return all_tau

def reference_ellipses(n=1):
    plt.close('all')
    fig, (ax1, ax2) = plt.subplots(2)

    ax1.set_xlabel('tau_low')
    ax1.set_ylabel('tau_high')
    ax2.set_title('deviation')
    ax2.set_xlabel('scd_low')
    ax2.set_ylabel('scd_high')
    ax1.set_xlim([-3, 0])
    ax1.set_ylim([-4, -1])
    #ax1.set_xscale('log')
    #ax1.set_yscale('log')
    ax1.set_title('TAU')
    return fig, (ax1,ax2)



def segmentation(folder, dist_max, dist_seg, autolevel, contrast, disk_size, max_contrast, soft_hard, showit, select_roi):

    direc, tapis = os.path.split(folder)
    direc, type_tapis = os.path.split(direc)

    '''Import frames'''

    # Import video file in HQ and select ROI

    #file_path = folder + "/" + type_tapis[0:2] + "_" + tapis + "_low_hq.tif" # You can simply input the file path in cmd args

    file_list = glob.glob(folder + "/*.tif" )

    file_path = None
    for f in file_list: 
        if "low_hq" in f:
            file_path = f
    if file_path == None:
        print("Error in the file names")

    direc = os.path.split(file_path)[0]

    # Initialize plotting tools
    g = alienlab.plot.ShowFigure() #tool to show images and multiple images
    g.figsize = (15,7) #size of the image (increase it to increase resolution)
    g.save_folder = direc #where the images are saved 
    g.date = False #not include the date on the image's filename
    p = alienlab.plot.PlotFigure() #tool to plot graphs and multiple curves
    p.figsize = (15,7)
    p.save_folder = direc
    p.date = False
    p.extension = ".jpg"


    #read the stacked frame. dim = NxHxW (N images in the video, Heigt, Width)

    frames_full = skimage.io.imread(file_path) 
    FO = FramesOperator(frames_full) #This class allows to compute operations on the images rapidly. It computes and stores properties of the video
    im = normalize(FO.frames[0], 0, 1) #Selects the first frame and normalize it 
    im = grey_to_rgb(im)*255  #Transfer to rgb to allow pop-up window to open it

    #y, x = alienlab.io.select_roi(np.uint8(im))
    #FO.x, FO.y = x, y

    if select_roi == True:
        r = cv2.selectROI(np.uint8(im)) #Select the minimal region encompassing all bacteria. Main purpose: speedup computation to ROI of interest and easier visualisation commenter pour le reste
    else: 
        r = (130, 0, 250, 250) 
    print("ROI COLLECTED THANKS", r)
    FO.x, FO.y = (r[1], r[1] + r[3]), (r[0], r[0] + r[2])

    FO.crop() #crops all the frames to the ROI

    start_time = time.time()
    FO.compute_stats() #computes mean, min, max, std along the frame dimension and the pixel dimension
    FO.normalize(0, 1) #normalises all the frames in 0-1 by the same affine function (conservation of global maximum and minimum). Necessary for image processing libraries
    print("--- Computed frames statistics in %04f seconds ---" % (time.time() - start_time))

    #FO.global_stats: each array has size N, number of frames and represents the stats of each frame
    #FO.frames_stats: each array has size FO.x, FO.y and is an image representing the N frames stats overlayed

    if showit:
        p.title = 'statistics'
        p.xlabel = 'frame number'
        p.ylabel = 'amplitude'
        p.label_list = ['max', 'min', 'mean', 'std']
        fig = p.plotting(np.asarray(FO.inds), [FO.global_stats['max'], 
                            FO.global_stats['min'], 
                            FO.global_stats['mean']])
        p.save_name = 'frames_stats'
        p.saving(fig)

    ''' IMAGE SEGMENTATION '''

    # First we generate references images that will be used for image segmentation: one reference image for the binary mask and one for the watershed segmentation
    # To compute these images we will stack several frames to increase the contrast
    # for more info on watershed: https://scikit-image.org/docs/dev/auto_examples/segmentation/plot_watershed.html and Wikipedia


        
    start_time = time.time()
    FO.selected_inds = FO.select_frames(FO.global_stats['max'], FO.global_stats['max'].max()*0.8)
    # Select only images with high intensity to increase contrast and lower computation time
    #for more info on the local operations applied: https://scikit-image.org/docs/dev/api/skimage.filters.rank.html

    frames_contrast = FO.apply(skimage.filters.rank.enhance_contrast,  selem = skimage.morphology.disk(contrast)) # Apply contrast filter to selected frames
    frames_autolevel = FO.apply(skimage.filters.rank.autolevel, selem = skimage.morphology.disk(autolevel)) #Apply autolevel filter to selected frames 
    frame_contrast = np.sum(frames_contrast, axis = 0) # sum selected processed images
    frame_autolevel = np.sum(frames_autolevel, axis = 0)

    mask_contrast = make_binary(frame_contrast, soft_hard = soft_hard) * make_binary(frame_autolevel) # dilate the items in the binary mask to have larger items
    mask_contrast = skimage.morphology.binary_dilation(mask_contrast, selem = skimage.morphology.disk(disk_size))
    auto_contrast = normalize(mask_contrast * frame_autolevel)# restrain the autolevel processed image to the binary mask (avoid outliers where there are no bacteria)

    print("--- Computed binary mask in %04f seconds ---" % (time.time() - start_time))

    if showit:
        g.title_list = 'contrast', 'autolevel', 'reference image'
        fig = g.multi([frame_contrast, frame_autolevel,  auto_contrast])
        g.save_name = 'Segmentation reference'
        g.saving(fig)

    start_time = time.time()
    ref = auto_contrast
    mask = mask_contrast
    local_maxi = alienlab.segment.local_maxima(auto_contrast, max_contrast, g,
                                                    ref_distance = dist_max, mask = mask, show = showit) # Find the local maxima in the reference image.
                                                    # Each local maxima is the source of 1 bacteria 
                                                    # in the watershed transformation: there will be as many bacteria as maxima.
                                                    # Possible to play on the second arg (int) to increase or decrease the number of maxima identified
                                                    # as it is the minimal distance between 2 local maxima


    watershed_im_mask = alienlab.segment.watershed(ref, mask, local_maxi,
                                                        g, ref_distance = dist_seg, show = showit) #Watershed segmentation flowing from identified local maxima restrained to binary mask
    segmented = watershed_im_mask
    print("--- Computed segmentation in %04f seconds ---" % (time.time() - start_time))

    if showit:
        alienlab.segment.show_segmentation(FO, segmented, g)

    ''' CORRECT SEGMENTATION WITH LABELME - which we don't use'''

    start_time = time.time()
    labelme_im = normalize(frame_autolevel)
    im_std = grey_to_rgb(np.tile(normalize(FO.frames_stats['mean'], 0, 1), 2))

    im_rgb = grey_to_rgb(np.tile(labelme_im, 2)) #duplicate the image as it will be easier to visualise the labels in LabelMe

    save_path = os.path.join(os.getcwd(), p.save_folder, "im_ref_" + 
                                os.path.split(os.path.split(file_path)[0])[1] + ".png")

    imageio.imwrite(save_path, im_rgb) #Save the image to open it with LabelMe. 
    json_path = alienlab.utils.replace_extension(save_path, ".json")
    alienlab.segment.segmented_to_json(segmented, save_path, im_rgb) # Extract the polygons from the watershed output and save it as '.json' labels for LabelMe

    #Actually we don't use this part because it's manual, time consuming and the segmentation code performs ok
    #print("--- Ready to open with LabelMe in %04f seconds ---" % (time.time() - start_time))

    #subprocess.run(['labelme', save_path, '-O', json_path ]) # Run LabelMe to allow correction of the labels
    #open labelme


    numbered_segmented, segmented = alienlab.segment.json_to_segmented(json_path, labelme_im) # Transform the corrected labels in a corrected image segmentation output

    items = np.unique(segmented) #returns the set of values in items, corresponds to the values of the markers of local_maxima

    if showit:
        alienlab.segment.show_segmentation(FO, numbered_segmented, g)

    framerate = {}
    limit = {}
    framerate['high'] = 1000
    framerate['low'] = 1000
    limit['high'] = [0, 200, 400, 600]
    limit['low'] = [0, 800, 1600, 2400]

    all_tau = {}

    for film in ['low', 'high']:
        # Binned frames
        #file_path = askopenfilename(initialdir = direc, title = 'Select the binned video %s now'%film) # pops up a window to select your file
        #file_path = folder + "/" + type_tapis[0:2] + "_" + tapis +"_" +  film + "_bin.tif"
        file_path = None
        for f in file_list: 

            if film in f and "bin" in f:
                file_path = f
        if file_path == None:
            print("Error in the file names")

        frames_binned = skimage.io.imread(file_path)
        # Select area of interest based on HQ frames cropped
        factor = 4
        FB = FramesOperator(frames_binned)

        midx = (FB.x[1] - FB.x[0])//2
        midy = (FB.y[1] - FB.y[0])//2
        dx = (FO.x[1] - FO.x[0])//factor
        dy = (FO.y[1] - FO.y[0])//factor

        # Rough binned ROI coordinates
        x_ref = [midx-dx//2, midx+dx//2] # The binned ROI corresponding to the HQ ROI is in the middle of the binned video and has heigth and width a fourth of HQ's ROI
        y_ref = [midy-dy//2, midy+dy//2]

        FB.compute_stats()
        if showit:
            p.xval
            p.title = 'statistics'
            p.xlabel = 'frame number'
            p.ylabel = 'amplitude'
            p.label_list = ['max', 'min', 'mean', 'std']
            fig = p.plotting(np.asarray(FB.inds), [FB.global_stats['max'], 
                                FB.global_stats['min'], 
                                FB.global_stats['mean'], 
                                FB.global_stats['std']])
            p.save_name = 'frames_stats_' + film
            p.saving(fig)

        start_time = time.time()
        # Find shift between images downscaling HQ image and computing cross-correlation between the 2 ROI. Precisely the image reference for HQ 
        # and binned is an image build with the max of the value taken by this pixel along the frames
        compute_binning = skimage.transform.downscale_local_mean(FO.frames_stats['max'], (4, 4))[:x_ref[1]-x_ref[0], :y_ref[1]-y_ref[0]] # downscaling of HQ ROI. 1 binned pixel value is the mean of factor HQ pixels
        im_binned_max = FB.frames_stats['max'][x_ref[0]:x_ref[1], y_ref[0]:y_ref[1]] # reference binned image
        shift, error, diffphase = skimage.feature.register_translation(im_binned_max, compute_binning, 100) # find the x and y shift between the 2 images with cross-correlation 
        # subpixel shift is computed by up-scaling the Fourier transform to increase the resolution of the images (factor 100 here)
        # more info here https://scikit-image.org/docs/dev/auto_examples/transform/plot_register_translation.html 
        dshift = np.copy(shift)
        shift[0] += x_ref[0] # take the shift into account for the ROI
        shift[1] += y_ref[0]
        print("--- Aligned binned and HQ videos in %04f seconds ---" % (time.time() - start_time))

        if showit:

            g.title_list = 'Binned image', 'Downscaled image'
            fig = g.multi([im_binned_max, compute_binning])
            g.save_name = 'shift_%s'%film
            g.saving(fig)

        SHIFT = np.zeros((2,1))
        SHIFT[0] = shift[0] 
        SHIFT[1] = shift[1]
        SHIFT *= factor
        # Shift matrix will be used to compute binned coordinates

        ''' COMPUTE BACTERIA TIME TRAJECTORIES '''
        start_time = time.time()

        # Collect item labels

        # Item time trajectories with overlaps
        # create a dictionnary with one entry for each item:
        '''
        { '1.0': {'x_coords': np array, x coordinates in HQ}
                    'y_coords': np array,  y coordinates in HQ
                    'binned_coords': set, couples of (x,y) coordinates in binned video
                    'surface': number of pixels in the item in HQ
                    'pixel_values': array, size: (N, s) where N is number of frames and s surface
                    'mean': array, size N, mean value of the item intensity for each frame
                    'std':  array, size N, std value of the item intensity for each frame
                    'remains' : True, the item is present in this segmentation step
                    }
        '2.0': {'x_coords'...
                    ...
                        }
            }

        '''
        items_dict = {}
        for k in items:
            key = str(k)
            items_dict[key] = {}
            x_coords, y_coords = np.nonzero(segmented == k)
            items_dict[key]['x_coords'] = x_coords
            items_dict[key]['y_coords'] = y_coords
            binned_coords = (np.stack((x_coords, y_coords))+ SHIFT)//factor
            items_dict[key]['binned_coords'] = set(zip(*binned_coords))
            pixel_values = FO.frames[:,x_coords, y_coords]
            items_dict[key]['pixel_values'] = pixel_values
            items_dict[key]['surface'] = pixel_values.shape[1]
            items_dict[key]['mean'] = np.mean(pixel_values, axis = 1)
            items_dict[key]['std'] = np.std(pixel_values, axis = 1)
            items_dict[key]['remains'] = True


        # list all the item coordinates in the binned image
        items_coordinates = []
        for k in items: 
            if k != 0:
                key = str(k)
                items_coordinates += list(items_dict[key]['binned_coords'])

        # Pixel positions of bacterias    
        couples, counts = np.unique(np.asarray(items_coordinates), axis = 0, return_counts = True) #extract the set of couple (x, y) and the number of occurence of each couple
        # Remove pixels belonging to several classes
        outliers = couples[counts > 1] # remove if occurence is over 1 (overlap)
        outliers_set = set(zip(*(outliers.T))) #convert to pyhton set

        # Bacteria time trajectories without overlaps
        items_dict_binned = {} # for all infos
        panda_dict = {}
        panda_dict['frame'] = FB.inds # only for mean values, transfer to other software
        panda_xy_binned = {}

        # create a dictionnary with one entry for each item:
        '''
        { '1.0': {'x_coords': np array, x coordinates in binned}
                    'y_coords': np array,  y coordinates in binned
                    'binned_coords': set, couples of (x,y) coordinates in binned video
                    'surface': number of pixels in the item in binned
                    'pixel_values': array, size: (N, s) where N is number of frames and s surface
                    'mean': array, size N, mean value of the item intensity for each frame
                    'std':  array, size N, std value of the item intensity for each frame
                    'remains' : True, the item is present in this segmentation step, False if the item was completely discarded due to overlaps
                    }
        '2.0': {'x_coords'...
                    ...
                        }
            }

        '''
        for k in items:
            key = str(k)
            items_dict_binned[key] = {}
            items_dict_binned[key]['binned_coords'] = items_dict[key]['binned_coords'] - outliers_set
            try:
                x_coords_binned, y_coords_binned = zip(*items_dict_binned[key]['binned_coords'])
                x_coords_binned = np.array(x_coords_binned).astype(int)
                y_coords_binned = np.array(y_coords_binned).astype(int)
                items_dict_binned[key]['x_coords'] = x_coords_binned
                items_dict_binned[key]['y_coords'] = y_coords_binned
                items_dict_binned[key]['surface'] = x_coords_binned.shape[0] 
                pixel_values_binned = FB.frames[:,x_coords_binned, y_coords_binned]
                items_dict_binned[key]['pixel_values'] = pixel_values_binned
                items_dict_binned[key]['mean'] = np.mean(pixel_values_binned, axis = 1)
                items_dict_binned[key]['std'] = np.std(pixel_values_binned, axis = 1)
                items_dict_binned[key]['remains'] = True
                panda_dict[key] = np.mean(pixel_values_binned, axis = 1)
                panda_xy_binned[key] = (x_coords_binned.mean(), y_coords_binned.mean(), x_coords_binned.shape[0])
            except:
                items_dict_binned[key]['remains'] = False  # If there were only overlaps fr the binned pixels  
        print("--- Computed time trajectories in %04f seconds ---" % (time.time() - start_time))

        if showit: 
            plot_items_coordinates = []
            no_overlap_coordinates = []
            labels_binned = []
            for k in items: 
                key = str(k)
                
                if k != 0 and items_dict_binned[key]['remains'] == True:
                    plot_items_coordinates += list(items_dict[key]['binned_coords'])
                    reduced_set = list(items_dict_binned[key]['binned_coords'])
                    no_overlap_coordinates += list(reduced_set)
                    labels_binned += [k/255]*len(reduced_set)
                    
            couples, counts = np.unique(np.asarray(plot_items_coordinates), axis = 0, return_counts = True)
            no_overlap_coordinates = np.asarray(no_overlap_coordinates)

            im_max = FB.frames_stats['max']
            overlap = normalize(np.copy(im_max), 0, 1)
            overlap = grey_to_rgb(overlap)
            overlap[:,:,:] = 0
            overlap[couples[:,0].astype(int), couples[:,1].astype(int), 0] = 1


            no_overlap = normalize(np.copy(im_max), 0, 0.5)
            no_overlap = grey_to_rgb(no_overlap)

            no_overlap[:,:,0] = 0
            no_overlap[no_overlap_coordinates[:,0].astype(int), 
                    no_overlap_coordinates[:,1].astype(int), 0] = 1#labels_binned

            g.figsize = (20, 20)
            g.col_num = 2
            g.title_list = ['image HQ',  'no overlaps', 'max_binned', 'overlaps']
            x_0 = max(int(x_ref[0] + dshift[0]), 0)
            x_1 = min(int(x_ref[1] + dshift[0]), im_max.shape[1])
            y_0 = max(int(y_ref[0] + dshift[1]), 0)
            y_1 = min(int(y_ref[1] + dshift[1]), im_max.shape[1])
            fig = g.multi([im_rgb, 
                    no_overlap[x_0: x_1, y_0: y_1],
                    im_max[x_0: x_1, y_0: y_1], 
                    overlap[x_0: x_1, y_0: y_1]])
            g.save_name = 'binned_segmentation_%s'%film
            g.saving(fig)

        data = pd.DataFrame.from_dict(panda_dict, orient = 'index').T    # Convert to table
        data2 = pd.DataFrame.from_dict(panda_xy_binned, orient = 'index').T    # Convert to table

        front = os.path.split(direc)[1]
        end =  alienlab.utils.replace_extension(os.path.split(file_path)[1], '.csv')
        save_data = os.path.join(direc, front + '_' + end)
        save_data2 = os.path.join(direc, front + '_coords_' + end)

        data.to_csv(save_data, index = False) # save as comma seperated values format
        data2.to_csv(save_data2)

        keys = list(panda_dict.keys())
        keys.pop(0)
        keys.pop(0)

        all_tau = fit_exponential(items_dict_binned, limit[film], framerate[film], all_tau, film, keys, p, True, exp_decay)

    fig, (ax1,ax2) = reference_ellipses()
    allkey = all_tau['high'].keys() & all_tau['low'].keys()

    colors_j = {'0': 'r', '1': 'g', '2': 'b'}
    for k in allkey:
        for j in all_tau['high'][k].keys():
            taus_low = np.array(all_tau['low'][k][j])
            taus_high = np.array(all_tau['high'][k][j])


            ax1.scatter(np.log10(taus_low.mean()), np.log10(taus_high.mean()), color = colors_j[j])
            ax1.annotate(str(k), (np.log10(taus_low.mean()), np.log10(taus_high.mean())))



            ax2.scatter(taus_low.std()/taus_low.mean(), taus_high.std()/taus_high.mean(), color = colors_j[j])
            ax2.annotate(k,( taus_low.std()/taus_low.mean(), taus_high.std()/taus_high.mean()) )
            ax2.set_title('deviation')
            ax2.set_xlabel('scd_low')
            ax2.set_ylabel('scd_high')

    g.save_name = "Distribution"
    g.saving(fig)
    #plt.show() #Faire apparaître image tapis
    a_file = open(direc + "/" + tapis + "_tau_val.pkl", "wb")

    pickle.dump(all_tau, a_file)

    a_file.close()

    return(all_tau, numbered_segmented, segmented, im_std, save_path )


