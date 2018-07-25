# RAP-Aviation
RAP plots for aviation nowcasting

July 25, 2018 (5:11 PM) - Ellrod Index

Ellrod Index is now working. Convergence and deformation are computing using MetPy, whereas vertical wind shear is computed using an approximated difference
method. Right now, plots are only available at isobaric levels. Eventually, I would like to be able to interpolate these to flight levels (e.g. FL300), but
the isobaric levels are a start. Also note that the data is smoothed using scipy.ndimage.filters.gaussian_filter in order to remove some of the noisiness
that can be a problem in high-resolution model output. The amount of smoothing can be controlled via the "sigma" variable.

July 25, 2018 - Initial Build

Right now, the only thing this code does is plot wind speeds and geopotential heights at each pressure level using RAP GRIB data from NCEP. Long-term, I hope to use this code to have some real-time and short-term forecast RAP products for clear air turbulence, mountain turbulence, mountain obscuration, and icing products for use in aviation forecasting. The to-do list is below. Feel free to contact me via my Twitter (@jgodwinWX) if you have any questions, comments, concerns, or requests.

To-Do List:
1. Clear Air Turbulence Products: Ellrod Index, Divergence Tendency, Richardson Number, Turbulent Kinetic Energy, Graphical Turbulence Guidance, and Ellrod-Knox Index.
2. Interpolate information to flight levels instead of pressure levels.
3. Mountain Wave Turbulence
4. Mountain Obscuration
5. Icing
6. Forecast
7. Overlay PIREPs

Completed To-Do List Items:
1. Clear Air Turbulence Products: Ellrod Index (2018/07/25).

Product Descriptions:
-Ellrod Index: https://aviationweather.gov/exp/ellrod/info.php

Sample plots can be found at: http://www.jasonsweathercenter.com/aviation/
