# RAP-Aviation
RAP plots for aviation nowcasting

July 26, 2018 - Divergence Tendency and Richardson Number

Divergence tendency seems to be working. The computation should be fairly straightforward since we already have the convergence values at the analysis time
and the forecast time (DIVERGENCE TREND = -CONVERGENCE_T2 + CONVERGENCE_T1). The values seems pretty small with not many contours being plotted, but I would
also estimate that these values are generally pretty small in the summer with not much in the way of ageostrophic wind tendencies compared to in the cold
season. The addition of this product required some tweaking of how variables are imported into pygrib since we now need two grib files. The same basic 
algorithms are used, but in the spirit of DRY ("don't repeat yourself"), they were written as functions with everything being passed. The assert statements
were also set up as functions. The addition of the functions should allow an easier transition into forecast products.

Richardson Number was a bit tricky, but was accompolished by first computing the Brunt-Vaisala Frequency, then dividing it by the square of the vertical wind
shear. This product seems to be the most useful so far, and in just a couple days of casually looking, appears to line up the best with turbulence PIREPs.
This is probably because Richardson Number is the only product being computed that considers both dynamic effects (wind shear) and thermodynamic effects
(Brunt-Vaisala Frequency is related to static stability) that lead to the development of turbulent eddies.

July 25, 2018 (5:11 PM) - Ellrod Index

Ellrod Index is now working. Convergence and deformation are computing using MetPy, whereas vertical wind shear is computed using an approximated difference
method. Right now, plots are only available at isobaric levels. Eventually, I would like to be able to interpolate these to flight levels (e.g. FL300), but
the isobaric levels are a start. Also note that the data is smoothed using scipy.ndimage.filters.gaussian_filter in order to remove some of the noisiness
that can be a problem in high-resolution model output. The amount of smoothing can be controlled via the "sigma" variable.

July 25, 2018 - Initial Build

Right now, the only thing this code does is plot wind speeds and geopotential heights at each pressure level using RAP GRIB data from NCEP. Long-term, I hope to use this code to have some real-time and short-term forecast RAP products for clear air turbulence, mountain turbulence, mountain obscuration, and icing products for use in aviation forecasting. The to-do list is below. Feel free to contact me via my Twitter (@jgodwinWX) if you have any questions, comments, concerns, or requests.

To-Do List:
1. Clear Air Turbulence Products: Ellrod Index*, Divergence Tendency*, Richardson Number*, Turbulent Kinetic Energy, Graphical Turbulence Guidance, and Ellrod-Knox Index.
2. Interpolate information to flight levels instead of pressure levels.
3. Mountain Wave Turbulence
4. Mountain Obscuration
5. Icing
6. Forecast
7. Overlay PIREPs

*Completed To-Do List Items:
1. Clear Air Turbulence Products: Ellrod Index (2018/07/25), Richardson Number (2018/07/25), Divergence Tendency (2018/07/26).

Product Descriptions:
-Ellrod Index: https://aviationweather.gov/exp/ellrod/info.php

Sample plots can be found at: http://www.jasonsweathercenter.com/aviation/
