import pygbif
import numpy as np
import matplotlib.pyplot as plt
import osgeo.gdal as gdal
from sklearn.neighbors import KernelDensity as kde

# # Pagination for exceeding the 300 parsing limit; 
def get_coord(sp_name, limit=999999, **kwargs):
    """
    Retrieve occurrence coordinate-s via streaming, using the search function from pygbif.occurrences module.
    Accepts all key-value pair arguments of the original function, such as continent and country.
    
    Params:
        sp_name (str): name of species.
        limit (int): the number of data intended to retrieve. Defaults to 999999.
        for more params, please refer to https://pygbif.readthedocs.io/en/latest/modules/occurrence.html
        
    Return:
        (list) A Nx2 nested list where each row is an observation.
        first item (index: 0) is the latidude and the latter is the longitude. [[latitude, longitude],[]...]
    """
    EOR = False
    sp_key = pygbif.species.name_backbone(sp_name)['speciesKey']
    filter_keys = ['decimalLatitude', 'decimalLongitude']
    off_count = 0
    new_list = []
    while EOR is False:
        search_output = pygbif.occurrences.search(taxonKey=sp_key, limit=limit, offset=off_count, hasCoordinate=True, **kwargs)
        occurrences_results = search_output['results']
        for item in occurrences_results:
            new_list.append([item[key] for key in filter_keys])
            new_list.pop() if len(new_list[-1]) == 0 else None
        off_count += 300
        EOR = search_output['endOfRecords']
        if 0 < limit <= 300:
            break
        limit -= 300
    return new_list

# # get WorldClim GeoTIFF
def bioclimap(coord_list, bcv_list, reso=10, model='ACCESS-CM2', scenario=585, pad=1, future=False):
    """
    Retrieve the bioclimate rasters from WorldClim (https://www.worldclim.org/data/index.html)
    and crop it to the specified padding from the outermost observations
    
    Params:
        coord_list (list[float]): A Nx2 list of species records.
        bcv_list (list[int]): A list of numbers [1:19] of the bioclimate variables listed in WorldClim.
        reso (int): Geospatial resolution. Options including 10, 5, 2.5, 30 (sec). Defaults to 10.
        model (str): Climate prediction models. Please see WorldClim for all the available options. Defaults to "ACCESS-CM2".
        scenario (int): Shared Socioeconomic Pathways (SSPs) climate scenarios. Options including 126, 245, 370, 585. Defaults to 585.
        pad (int, float): Padding from the outermost observations. Defaults to 1.
        future (bool): Assign "True" to retrieve future climate raster. Defaults to "False".
        
    Return:
        None. This function does not return anything but define a few global variables for inter-functional usage.
    """
    global array_drawer_84320700, geo_repo_84320700, extent_84320700
    array_drawer_84320700 = []
    ggt_item = gdal.Open(fr"/vsizip/vsicurl/https://geodata.ucdavis.edu/climate/worldclim/2_1/base/wc2.1_{reso}{'m' if reso != 30 else 's'}"
                    fr"_bio.zip/wc2.1_{reso}{'m' if reso != 30 else 's'}_bio_1.tif").GetGeoTransform()
    lat = [i[0] for i in coord_list]
    long = [i[1] for i in coord_list]
    long_min, long_max = min(long), max(long)
    lat_min, lat_max = min(lat), max(lat)
    x_BgId = int(np.floor(((long_min-pad)-ggt_item[0])/ggt_item[1]))
    y_BgId = int(np.floor(((lat_max+pad)-ggt_item[3])/ggt_item[5]))
    x_ExId = int(np.ceil(((long_max+pad)-(long_min-pad))/ggt_item[1]))+x_BgId
    y_ExId = int(np.ceil(((lat_min-pad)-(lat_max+pad))/ggt_item[5]))+y_BgId
    long_anc = ggt_item[0]+(x_BgId)*ggt_item[1]
    lat_anc = ggt_item[3]+(y_BgId)*ggt_item[5]
    ras_width, ras_height = x_ExId-x_BgId, y_BgId-y_ExId
    array_drawer_84320700.clear()
    for bcvar in bcv_list:
        if not future:
            gdal_box = gdal.Open(fr"/vsizip/vsicurl/https://geodata.ucdavis.edu/climate/worldclim/2_1/base/wc2.1_"
                            fr"{reso}{'m' if reso != 30 else 's'}_bio.zip/wc2.1_{reso}{'m' if reso != 30 else 's'}_bio_{bcvar}.tif")
            gdal_band = gdal_box.GetRasterBand(1)
        else:
            gdal_box = gdal.Open(fr"/vsicurl/https://geodata.ucdavis.edu/cmip6/{reso}{'m' if reso != 30 else 's'}/{model}/ssp{scenario}/wc2.1_"
                        fr"{reso}{'m' if reso != 30 else 's'}_bioc_{model}_ssp{scenario}_2021-2040.tif")
            gdal_band = gdal_box.GetRasterBand(bcvar)
        gdal_nan = gdal_band.GetNoDataValue()
        gdal_array = gdal_band.ReadAsArray()
        gdal_array[gdal_array <= gdal_nan] = np.nan
        gdal_array[np.isinf(gdal_array)] = np.nan
        sliced_array = gdal_array[y_BgId: y_ExId, x_BgId: x_ExId]
        array_drawer_84320700.append(sliced_array)
    geo_repo_84320700 = [ggt_item[0], ggt_item[1], ggt_item[2], ggt_item[3], ggt_item[4], ggt_item[5],
                long_min, long_max, lat_max, lat_min, x_BgId, x_ExId, y_BgId, y_ExId, lat_anc, long_anc, ras_width, ras_height]
    extent_84320700=(geo_repo_84320700[-3], geo_repo_84320700[-3] + geo_repo_84320700[-2] * geo_repo_84320700[1], 
                   geo_repo_84320700[-4] - geo_repo_84320700[-1] * geo_repo_84320700[5], geo_repo_84320700[-4])

# # clusters and outliers
def radius_cluster(coord_list, radius=1):
    """
    Calculate the cluster-ness or density of an observation.
    
    Params:
        coord_list (list[float]): A Nx2 list of species records.
        radius (int, float): The radius to be grouped into a cluster, in the unit of degree. Defaults to 1.
        
    Return:
        (list) A list of the observation-count in a cluster, following the sequence of the input coordinate list.
    """
    coord_array = np.asarray(coord_list)
    c_list = []
    for coorows in coord_array:
        pythag = np.sqrt(np.sum((coorows - coord_array)**2, axis=1))
        c_list.append(np.count_nonzero(pythag <= radius))
    return c_list

# # normalise to 0-1
def get_norm_const(coord_list, bins=30, cv_check=False):
    """
    Obtain the optimum (habitat suitability of 1) from the training set to normalise the testing (i.e., prediction) set.
    A histogram with kernel density estimation will be plot and users will be prompted for customised optimum input.
    Defaults to the mode of the dataset as the function assumes a normal distribution of species with respect to the respective bioclimatic variable.
    
    Params:
        coord_list (list[float]): A Nx2 list of species records.
        bins (int): Number of histogram bins.
    
    Return:
        None. This function does not return anything but define a few global variables for inter-functional usage.
    """
    # try:
    global norm_rng_84320700
    norm_rng_84320700 = [] # reset every run
    coord_arr = np.asarray(coord_list)
    indices_arr = np.floor((coord_arr-[geo_repo_84320700[-4], geo_repo_84320700[-3]])/[geo_repo_84320700[5], geo_repo_84320700[1]]).astype(int)
    row_i = indices_arr[:, 0]
    col_i = indices_arr[:, 1]
    out_check = True
    for arr in array_drawer_84320700:
        enviro_vlist = arr[row_i, col_i]
        enviro_vlist = enviro_vlist[~np.isnan(enviro_vlist)]
        reshp_env = enviro_vlist.reshape(-1, 1)
        while out_check:
            the_min = min(enviro_vlist)
            the_max = max(enviro_vlist)
            the_mean = np.mean(enviro_vlist)
            the_median = np.median(enviro_vlist)
            # # get mode bin
            bin_width = (the_max-the_min)/bins
            binedge_list = np.arange(the_min, the_max+bin_width, bin_width)
            bin_count, _ = np.histogram(enviro_vlist, bins=binedge_list.tolist())
            bin_count = bin_count.astype(int)
            mode_index = np.argmax(bin_count)
            optimum = np.median([v for v in enviro_vlist if binedge_list[mode_index] <= v < binedge_list[mode_index+1]])
            scaling_range = [the_min, the_max, optimum]
            # # plot the graph
            if not cv_check:
                spl_x = np.linspace(np.nanmin(arr)-1, np.nanmax(arr), 400).reshape(-1, 1)
                kde_func = kde().fit(reshp_env)
                kde_y = kde_func.score_samples(spl_x)
                plt.plot(spl_x, np.exp(kde_y), color='purple', label=f'Kernel Density')
                plt.axvline(optimum, linestyle='dashed', color='red', label=f'mode: {optimum}')
                plt.axvline(the_median, linestyle='dashed', color='blue', label=f'median: {the_median}')
                plt.axvline(the_mean, linestyle='dashed', color='green', label=f'mean: {the_mean}')
                plt.hist(enviro_vlist, bins=binedge_list, density=True)
                plt.legend()
                plt.show()
                out_remover = input('To remove outliers, input to-be-retained [inclusive range] in '
                                    '--> start:end <-- seperated by colon. Leave blank to retain all values: ')
                if not out_remover: break
                start_HS, end_HS = map(float, out_remover.split(':'))
                enviro_vlist = [x for x in enviro_vlist if start_HS <= x <= end_HS]
            else: break
        if not cv_check:
            input_values = (input('min (int or float only. Leave blank to default): '),
                input('max (int or float only. Leave blank to default): '),
                input('optimum (int or float only. Leave blank to default): '))
        else: input_values = ((),(),(the_median)) #!
        scaling_range = [float(val) if val else var for var, val in zip(scaling_range, input_values)]
        norm_rng_84320700.append(scaling_range)
    # except:
    #     print("No BioClim raster was found. Please call the bioclimap function before proceed.")

def normalise(thres=0.5, extrm=True, cv_check=False):
    """
    Normalise the array values.
    Both extreme ends of the bioclim-variables were deemed unsuitable and the two-tailed 5% were remove.
    
    Params:
        extrm (bool): Remove the two-tailed extreme 5% if set to True. Defaults to true.
        
    Return:
        (numpy.array): An array of normalised habitat suitability.
    """
    # try:
    array_lab = array_drawer_84320700 # resetting
    for arr, arr_value in zip(array_lab, norm_rng_84320700):
        right_mask = arr >= arr_value[2]; left_mask = arr < arr_value[2]
        np.putmask(arr, right_mask, (arr-arr_value[1])/(arr_value[2]-arr_value[1]))
        np.putmask(arr, left_mask, (arr-arr_value[0])/(arr_value[2]-arr_value[0]))
        if extrm:
            np.putmask(arr, arr<0.05, 0)
    masking_arr = np.minimum.reduce(array_lab)
    np.putmask(masking_arr, masking_arr>0, 1)
    final_arr = masking_arr*(np.mean(array_lab, axis=0))
    if not cv_check:
        return final_arr
    else:
        bi_class = np.ones(arr.shape)
        np.putmask(bi_class, final_arr<thres, -1) #-1 & 1 for model grading
        return bi_class
    # except:
    #     print("Scaling values for normalisation are not defined. Please call the related function before proceed.")

def model_grade(cv_norm_arr, valid_set, target_list):
        presence_arr = np.ones(cv_norm_arr.shape)
        test_arr = np.asarray(valid_set)
        indicies_arr = np.floor((test_arr-[geo_repo_84320700[-4], geo_repo_84320700[-3]])/[geo_repo_84320700[5], geo_repo_84320700[1]]).astype(int)
        id_list = indicies_arr.tolist()
        compare = None
        for i in id_list:
            if i == compare:
                id_list.remove(i)
                compare = i
        id_arr = np.asarray(id_list)
        row_i = id_arr[:, 0]
        col_i = id_arr[:, 1]
        presence_arr[row_i, col_i] = 2
        mul_matrix = cv_norm_arr*presence_arr
        # # performance-quantifying-metrics matrix
        # # # Only presence-only metric (i.e., sensitivity) a.t.m.
        true_posi = np.count_nonzero(mul_matrix == 2)
        false_nega = np.count_nonzero(mul_matrix == -2)
        sensitivity = true_posi/(true_posi+false_nega)
        target_list.append(sensitivity)

# # cross validation function
def crossvalid(coord_list, kfold=None, test_perc=None, N=10, thres=0.5):
    sensi_list = []
    shuff_copy = coord_list
    if kfold and test_perc:
        print("Method conflict. Please assign only value to kfold or test_perc.")
        return
    elif not kfold and not test_perc:
        print("Nethod Error: No cross-validate method selected.")
        return
    
    if kfold:
        if not isinstance(kfold, int):
            print("Please ensure k-fold is integer.")
            return 
        test_slicer = round(len(shuff_copy)/kfold)
        np.random.shuffle(shuff_copy)
        train_form = []
        k_bgid = 0; k_exid = test_slicer
        for k in range(kfold):
            fold_test = shuff_copy[k_bgid:k_exid]
            train_latt = shuff_copy[k_exid:] if k is not kfold-1 else []
            fold_train = train_form + train_latt
            get_norm_const(fold_train, cv_check=True)
            cv_arr = normalise(cv_check=True)
            model_grade(cv_arr, fold_test, sensi_list)
            # loopers 
            train_form.extend(fold_test)
            k_bgid += test_slicer
            k_exid = k_exid + test_slicer if k != kfold-1 else None
        avrg_sensi = np.mean(sensi_list)
        print(f"Average sensitivity: {avrg_sensi}")
    
    if test_perc:        
        if test_perc <1 or test_perc >50:
            print("percentage cant be 0, negative, and has a max set at 50%")
            return
        deno = 100/test_perc
        shuff_uniq = set()
        while len(shuff_uniq) != N:
            np.random.shuffle(shuff_copy)
            test_set = shuff_copy[:len(coord_list)//deno]
            train_set = shuff_copy[len(coord_list)//deno:]
            if test_set not in shuff_uniq:
                shuff_uniq.add(tuple(sorted(test_set)))
                get_norm_const(train_set, cv_check=True)
                cv_arr = normalise(thres=thres, cv_check=True)
                # # model grading
                model_grade(cv_arr, test_set, sensi_list)
        avrg_sensi = np.mean(sensi_list)
        print(f"Average sensitivity: {avrg_sensi}")

# example demonstration
if __name__ == '__main__':
    
    bioclims = [1, 12]
    # get some occurrences data
    exmpl_coord = get_coord("Bombus affinis", continent='North America', year='1970,2020')
    latidudes = [i[0] for i in exmpl_coord]
    longitudes = [i[1] for i in exmpl_coord]
    bioclimap(exmpl_coord, bioclims, reso=10)
    # cross validate
    crossvalid(exmpl_coord, kfold=10)

    # Once satisfy with the performance, use the model to predict future distribution
    future_coord = get_coord("Bombus affinis", continent='North America', year='2021,2024')
    fut_lat = [i[0] for i in future_coord]
    fut_long = [i[1] for i in future_coord]
    bioclimap(future_coord, bioclims, reso=10, pad=2, future=True)
    get_norm_const(future_coord)
    final_future = normalise()
    # calculate clustering density
    dens_colour = radius_cluster(future_coord, radius=2)
    
    # map plotting
    scatter_cmap = plt.get_cmap('cool', 7)
    scatt_plot = plt.scatter(fut_long, fut_lat, c=dens_colour, cmap=scatter_cmap)
    legd_elements = scatt_plot.legend_elements()
    handles = legd_elements[0][::2]
    labels = [f"Up to {x}" for x in legd_elements[1][1::2]]
    plt.legend(handles, labels, title='Occurrences density')
    plt.imshow(final_future, cmap='viridis',
            extent=extent_84320700)
    plt.colorbar(label='Habitat Suitability Index')
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')
    plt.title('Rusty Patched Bumblebee BioClim')
    plt.show()