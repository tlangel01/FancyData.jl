# FancyData.jl

Some functions that I use regularly when working with DataFrames (DataFrames.jl) and uncertainties (Measurement.jl)

## Install via

- `] add https://github.com/tlangel01/FancyData.jl.git`

## Features

- `mes()`: Converts a Measurement to a String with format `value(uncertainty)`
- `mean_std()`: Calculates the mean with standard deviation of a given array
- `wmean()`: Calculates the weighted mean of a given array of Measurements
- `tableDF()`: Prints a DataFrame to a TeX-formatted table
- `readDF()`, `writeDF()`: Reads a delimited file to a DataFrame and vice versa
- `readfits()`: Fitlistreader for HDTV-based XML-files
