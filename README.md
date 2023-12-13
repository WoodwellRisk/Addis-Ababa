# distribution_plots

Code in this repository is used to make distribution plots of climate variables.

distribution_plot_single_city: Use this plot when making a single distribution of a single variable for a single location using an ensemble of climate model data. The script is currently set up to run optimally using bias-adjusted and downscaled 0.5-deg CMIP6 data, but can be altered to work with a variety of climate model data. In its current form, the script assumes that all climate model input are on the same grid, with data from the grid cell nearest to the user-specified lat and lon will be used in the distribution.
