import numpy as np
import osgeo.gdal as gdal

def ras_reso(coords, div=1000, pad=5): # in coordinate degree
    global GT_0109
    coray = np.asarray(coords)
    long = coray[:,0]; lat = coray[:,1]
    dstnc = [max(long)-min(long), max(lat)-min(lat)]
    ratio1 = float(max(dstnc)/min(dstnc))
    if not ratio1.is_integer:
        if round(ratio1)-ratio1 > 0: # round up
            i = dstnc.index(max(dstnc))
            x = max(dstnc)+(min(dstnc)*round(ratio1)-max(dstnc))
        else:
            i = dstnc.index(min(dstnc))
            x = min(dstnc)+(max(dstnc)/round(ratio1)-min(dstnc))
        dstnc[i] = x
        ratio1 = round(ratio1)
    new_width = dstnc[1]/div #shape need class!
    long_pin = min(long)-pad*new_width
    lat_pin = min(lat)-pad*new_width #not to be confused with GT's negative y value
    GT_0109 = [long_pin, new_width, 0, lat_pin, 0, -new_width]
    shape_x_0109 = [div+2*pad if dstnc[0] <= dstnc[1] else div*ratio1+2*pad]
    shape_y_0109 = [div+2*pad if dstnc[1] <= dstnc[0] else div*ratio1+2*pad]
    GT_0109.append(shape_x_0109[0]); GT_0109.append(shape_y_0109[0])

def rescale(coords, arr, GT, nan=None):
    # GT = geo-transformed set of tiff
    # arr = numpy array to be resized
    # width = size of cell of new grid (in lat/long)
    if len([coords][0]) != 2:
        print("list of coordinates must comprised pairs of first-longitude, then-latitude, respectively")
    elif not isinstance(arr, np.ndarray):
        print("please parse in the target tif in the form of NumPy's array")
    elif arr.ndim != 2:
        print("only 2D array is allowed")
    elif len(GT) != 6:
        print("ensure gdal's geotransform info is included")
    elif GT[1] != -GT[5]:
        print("why in the Earth do you have a rectangular grid?")
    elif "GT_0109" not in globals():
        ras_reso(coords)
    elif GT[0] > GT_0109[0] or GT[3] < GT_0109[3]:
        print("input raster was too small for rescaling") #!!!

    # parallels = np.arange(GT_0109[0], GT_0109[0]+GT_0109[6]*GT_0109[1]+GT_0109[1]*0.1, GT_0109[1])
    # ^^^stop num was set so it's always smaller than width
    arr = arr[int((GT[3]-GT_0109[3])/-GT[5]): int((GT[3]-(GT_0109[3]+GT_0109[5]*GT_0109[7]))/-GT[5])+1, int((GT_0109[0]-GT[0])/GT[1]): int((GT_0109[0]+GT_0109[1]*GT_0109[6])-GT[0])/GT[1]+1]
    arr[arr == nan] = np.nan
    GT = list(GT); GT.append(arr.shape[1]); GT.append(arr.shape[0])
    # generate a grid of centres for finer grid
    fine_GT = min(GT, GT_0109, key=lambda x:x[1]); coarse_GT = max(GT, GT_0109, key=lambda x:x[1])
    mid_longs = np.array(np.arange(fine_GT[0]+fine_GT[1]/2, fine_GT[0]+fine_GT[1]/2+fine_GT[6]*fine_GT[1], fine_GT[1]))
    #cut out OOR rows
    arr = arr[:, 1:] if mid_longs[0] < coarse_GT[0] else arr
    arr = arr[:, :-1] if mid_longs[-1] > coarse_GT[0]+coarse_GT[1]*coarse_GT[-2] else arr
    fine_GT[-2] = arr.shape[1]
    
    mid_lats = np.array(np.arange(fine_GT[3]+fine_GT[5]/2, fine_GT[3]+fine_GT[5]/2+fine_GT[5]*fine_GT[7], fine_GT[5]))
    arr = arr[1:, :] if mid_lats[0] < coarse_GT[3] else arr
    arr = arr[:-1, :] if mid_lats[-1] > coarse_GT[3]+coarse_GT[5]*coarse_GT[-1] else arr
    fine_GT[-1] = arr.shape[0]
    
    find_col = (mid_longs-coarse_GT[0])/coarse_GT[1]
    col1 = round(find_col); col2 = round(find_col)+1 if np.float64(round(find_col)) == find_col else round(find_col)
    find_row = (mid_lats-coarse_GT[3])/coarse_GT[5]
    row1 = round(find_row); row2 = round(find_row)+1 if np.float64(round(find_row)) == find_row else round(find_row)
    c1r1_r, c1r1_c = np.meshgrid(row1, col1); c2r1_r, c2r1_c = np.meshgrid(row1, col2)
    c1r2_r, c1r2_c = np.meshgrid(row2, col1); c2r2_r, c2r2_c = np.meshgrid(row2, col2)
    
    if coarse_GT[1] == GT[1]: #upscale
        c1r1 = arr[c1r1_r, c1r1_c]; c2r1 = arr[c2r1_r, c2r1_c]; c1r2 = arr[c1r2_r, c1r2_c]; c2r2 = arr[c2r2_r, c2r2_c]
        new_arr = np.nansum([c1r1, c2r1, c1r2, c2r2], axis=0)/4
        return new_arr
    
    elif coarse_GT[1] == GT_0109[1]: #downscale
        arr = np.float64(arr) # to avoid mismatch of value type
        # count of points included in greater grid
        c1r1_n, c2r1_n, c1r2_n, c2r2_n = [np.zeros(((coarse_GT[-2], coarse_GT[-1])))]*4
        [exec("np.add.at(c{i}r{j}_n, (c{i}r{j}_r, c{i}r{j}_c), 1)") for i in (1,2) for j in (1,2)]
        # nan-bit classifier
        nan_mask = np.where(np.isnan(arr), 0, 1)
        c1r1_b, c2r1_b, c1r2_b, c2r2_b = [np.zeros((coarse_GT[-2], coarse_GT[-1]))]*4
        [exec("np.add.at(c{i}r{j}_b, (c{i}r{j}_r, c{i}r{j}_c), nan_mask)") for i in (1,2) for j in (1,2)]
        # environmental values
        c1r1, c2r1, c1r2, c2r2 = [np.zeros((coarse_GT[-2], coarse_GT[-1]))]*4
        numerised_arr = np.where(np.isnan(arr), 0, arr) # preparation for add.at()
        [exec("np.add.at(c{i}r{j}, (c{i}r{j}_r, c{i}r{j}_c), numerised_arr)") for i in (1,2) for j in (1,2)]
        [exec("c{i}r{j}_b /= 4*c{i}r{j}_n") for i in (1,2) for j in (1,2)]; nan_mask = sum([c1r1_b, c2r1_b, c1r2_b, c2r2_b]) 
        [exec("c{i}r{j} /= 4*c{i}r{j}_n") for i in (1,2) for j in (1,2)]; new_arr = sum([c1r1, c2r1, c1r2, c2r2])
        new_arr = np.where(nan_mask >= 0.5, new_arr, np.nan)
        return new_arr