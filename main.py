import pygbif
import numpy as np
import matplotlib.pyplot as plt
import osgeo.gdal as gdal
from sklearn.neighbors import KernelDensity as kde
from bioclim import get_coord, bioclimap, radius_cluster, get_norm_const, normalise, model_grade, crossvalid

def main():
    
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
    ##!
    from bioclim import extent_84320700
    plt.legend(handles, labels, title='Occurrences density')
    plt.imshow(final_future, cmap='viridis',
            extent=extent_84320700)
    plt.colorbar(label='Habitat Suitability Index')
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')
    plt.title('Rusty Patched Bumblebee BioClim')
    plt.show()
    
if __name__ == '__main__':
    print("example demonstration")
    main()