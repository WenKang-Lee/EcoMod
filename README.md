# EcoMod
Greetings, fellow ecologists. This is a simple species distribution modelling package, with currently only having the BIOCLIM function.

## installation
Use the package manager pip to install the module.
Due to the issue of osgeo package, I've included the whl of GDAL for Win py3.11 in this package.
You might need to look for other relevant source if you're running on a different system, such as Mac.
```bash
pip install YOUR_PYTHON_FOLDER/GDAL-3.4.3-cp311-cp311-win_amd64.whl
pip install EcoMod
```

## Usage
```python
import EcoMod

# get occurrences data & WorldClim GeoTiff
bioclims = [1, 12]
# # |: repeat for testing set
coordinates = EcoMod.get_coord("sp. name")
latitudes = [i[0] for i in coordinates]
longitudes = [i[1] for i in coordinates]
EcoMod.bioclimap(coordinates, bioclims, reso=10)
# # # [1] call to get the optimum from the training set.
EcoMod.get_norm_const(coordinates)
# # :| [2] normalise the prediction (i.e., testing, present) map
clusters = EcoMod.radius_cluster(coordinates, 2) # optional
normed_array = EcoMod.normalise()

# view the prediction map using matplotlib
# # assign 'extent_84320700' (variable) to the 'extent' param
# # to crop out the Region of Focus
...
plt.imshow(normed_array, extent=extent_84320700)
plt.show()
```

## Roadmap
I'll next be working on a more advance and statistically comprehensive feature like MaxEnt and its related target group background.
The aesthetic part of the map (i.e, visualisation) is also part of the roadmap.
This project is also following the guidelines and suggestions in Soley-Guardia et al.'s (2024) publication on the top ten hazards that were seen when modelling species distribution.

## Contributing
Contributions are welcome. Feel free to open an issue for any suggestion or submit a pull request.

## License
This project is licensed under the conventional MIT License - see the LICENSE file for details.

## Acknowledgement
This eco-modelling package was inspired by Timothée Poisot's "BioClim from scratch" tutorial. Special thanks to Professor Poisot (@tpoisot) for creating such an informative guide on modelling the bioclimatic distribution of species.

## References
1. Mateo, R.G., Croat, T.B., Felicísimo, A.M. & Muñoz, J. (2010) Profile or group discriminative techniques? Generating reliable species distribution models using pseudo-absences and target-group absences from natural history collections._Diversity and Distribution_, 16(1), 84-94. https://doi.org/10.1111/j.1472-4642.2009.00617.x
2. Poisot, T. (2021) Building the BioClim model from scratch_.[\[video\]]YouTube.(https://www.youtube.com/watch?v=81MblyNWz3Q)
3. Soley-Guardia, M., Alvarado-Serrano, D.F. & Anderson, R.P. (2024) _Top ten hazards to avoid when modeling species distributions: a didactic guide of assumptions, problems, and recommendations._Ecography_, e06852. https://doi.org/10.1111/ecog.06852