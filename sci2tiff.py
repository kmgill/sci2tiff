
import sys
import os
import numpy as np
from libtiff import TIFFimage
import argparse
from scipy.misc import imresize
from astropy.io import fits
import vicar
import re

def load_fits_matrix(input_file, band=0):
    hdu_list = fits.open(input_file)

    if band + 1 > len(hdu_list):
        raise Exception("Band '%d' not found in data"%band)

    pixel_matrix = hdu_list[band].data

    pixel_matrix = np.flipud(pixel_matrix)
    return pixel_matrix


def load_vic_matrix(input_file):
    pixel_matrix, value_pairs = vicar.load_vic(input_file)
    return pixel_matrix

def load_image_matrix(input_file, band=0):

    if input_file.lower().find(".fits") >= 0:
        pixel_matrix = load_fits_matrix(input_file, band=band)
    else:  # Otherwise assume it to be vicar. TODO: Don't make this assumption
        pixel_matrix = load_vic_matrix(input_file)

    return pixel_matrix

"""
Detect rows of null values that are a result of Cassini non-lossy compression
"""
def detect_null_stripes(pixel_matrix):

    stripes = []

    height = pixel_matrix.shape[0]
    width = pixel_matrix.shape[1]

    for y in range(0, height-1, 1):
        for x in range(width - 2, -1, -1):
            pixel_value_prev_row = pixel_matrix[y-1][x]
            pixel_value = pixel_matrix[y][x]
            pixel_value_next_row = pixel_matrix[y+1][x]
            if not np.isnan(pixel_value) and pixel_value != 0.0 and pixel_value_prev_row != 0.0 and pixel_value_next_row != 0.0:
                if x < width - 2:
                    stripes.append((y, x+1))
                break

    return stripes

def fill_stripe(pixel_matrix, row, start_x):
    for x in range(start_x, pixel_matrix.shape[1]):
        prev_row_value = pixel_matrix[row-1][x]
        next_row_value = pixel_matrix[row+1][x]
        fill_value = np.mean([prev_row_value, next_row_value])
        pixel_matrix[row][x] = fill_value


def fill_stripes(pixel_matrix, stripes):
    for stripe in stripes:
        fill_stripe(pixel_matrix, stripe[0], stripe[1])


def get_data_min_max(input_file,
            fill_null_stripes=False,
            fillsat=False,
            dohisteq=False,
            minpercent=None,
            maxpercent=None,
            band=0,
            trim=0):

    pixel_matrix = load_image_matrix(input_file, band=band)

    pixel_matrix = process_data(pixel_matrix,
                                None,
                                None,
                                fill_null_stripes,
                                fillsat,
                                dohisteq,
                                minpercent,
                                maxpercent,
                                None,
                                trim=trim,
                                stretch=False,
                                outputformat=None)

    pixel_min = np.nanmin(pixel_matrix)
    pixel_max = np.nanmax(pixel_matrix)
    return pixel_min, pixel_max


def histeq(pixel_matrix, nbr_bins=65536):
    im = pixel_matrix
    imhist, bins = np.histogram(im.flatten(), nbr_bins, normed=True)
    cdf = imhist.cumsum()  # cumulative distribution function
    cdf = 65535 * cdf / cdf[-1]  # normalize

    # use linear interpolation of cdf to find new pixel values
    im2 = np.interp(im.flatten(), bins[:-1], cdf)
    im2 = im2.reshape(im.shape)

    return np.array(im2, np.uint16)

def build_output_filename(input_file):
    return input_file[:input_file.lower().rindex(".")]


def process_data(pixel_matrix,
            force_input_min=None,
            force_input_max=None,
            fill_null_stripes=False,
            fillsat=False,
            dohisteq=False,
            minpercent=None,
            maxpercent=None,
            std_mult=None,
            resize=None,
            trim=0,
            stretch=True,
            outputformat="uint16"):

    if type(trim) == int and trim > 0:
        pixel_matrix = pixel_matrix[trim:-trim, trim:-trim]


    if force_input_min is not None:
        pixel_min = force_input_min
    else:
        pixel_min = np.nanmin(pixel_matrix)

    if force_input_max is not None:
        pixel_max = force_input_max
    else:
        pixel_max = np.nanmax(pixel_matrix)


    if fill_null_stripes is True:
        stripes = detect_null_stripes(pixel_matrix)
        fill_stripes(pixel_matrix, stripes)

    # The min/max percent stuff isn't correct. TODO: Make it correct.
    if minpercent is not None:
        diff = pixel_min + ((pixel_max - pixel_min) * (minpercent / 100.0))
        pixel_matrix[pixel_matrix < diff] = diff
        pixel_min = diff

    if maxpercent is not None:
        diff = pixel_min + ((pixel_max - pixel_min) * (maxpercent / 100.0))
        pixel_matrix[pixel_matrix > diff] = diff
        pixel_max = diff

    if fillsat is True:
        inds = np.where(np.isnan(pixel_matrix))
        pixel_matrix[inds] = np.nanmax(pixel_matrix)


    if std_mult is not None:
        std = pixel_matrix.std()
        pixel_min = pixel_matrix.mean() - std * std_mult
        pixel_max = pixel_matrix.mean() + std * std_mult

    if stretch is True:
        pixel_matrix = pixel_matrix - pixel_min
        pixel_matrix = pixel_matrix / (pixel_max - pixel_min)
        pixel_matrix[pixel_matrix < 0] = 0

    if resize is not None:
        if re.match("\A\d+x\d+\Z", resize) is not None:
            new_size = map(int, re.split("x", resize))
        elif re.match("\A\d+[\.]*[\d]*%\Z", resize) is not None:
            new_size = (np.array((pixel_matrix.shape[1], pixel_matrix.shape[0])) * (float(resize[:resize.index("%")]) / 100.0)).astype(np.int32)
        else:
            print "Invalid resize specifier:", resize
            sys.exit(0)

        pixel_matrix = imresize(pixel_matrix, size=new_size, interp='bicubic', mode='F')

    if outputformat == "byte" or outputformat == "uint8":
        pixel_matrix = pixel_matrix * 256.0
        pixel_matrix = pixel_matrix.astype(np.uint8)
    elif outputformat == "uint16":
        pixel_matrix = pixel_matrix * 65535.0
        pixel_matrix = pixel_matrix.astype(np.uint16)
    elif outputformat == "uint32":
        pixel_matrix = pixel_matrix * 4294967295.0
        pixel_matrix = pixel_matrix.astype(np.uint32)
    elif outputformat == "float16":
        pixel_matrix = pixel_matrix.astype(np.float16)
    elif outputformat == "float32":
        pixel_matrix = pixel_matrix.astype(np.float32)

    if dohisteq is True:
        pixel_matrix = histeq(pixel_matrix, nbr_bins=65535)

    return pixel_matrix


def sci2tiff(input_file,
            force_input_min=None,
            force_input_max=None,
            fill_null_stripes=False,
            fillsat=False,
            dohisteq=False,
            minpercent=None,
            maxpercent=None,
            std_mult=None,
            resize=None,
            band=0,
            trim=0,
            outputformat="uint16"):

    pixel_matrix = load_image_matrix(input_file, band=band)
    pixel_matrix = process_data(pixel_matrix,
                                force_input_min,
                                force_input_max,
                                fill_null_stripes,
                                fillsat,
                                dohisteq,
                                minpercent,
                                maxpercent,
                                std_mult,
                                resize,
                                trim,
                                outputformat=outputformat)

    output_filename = build_output_filename(input_file)
    print "Writing", output_filename

    # Create output tiff
    tiff = TIFFimage(pixel_matrix, description='')
    tiff.write_file(output_filename, compression='none', verbose=False)
    return True


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--data", help="Input file(s)", required=True, type=str, nargs='+')
    parser.add_argument("-i", "--intensity", help="Match intensities when stretching", required=False, action="store_true")
    parser.add_argument("-f", "--fill", help="Fill null stripes", required=False, action="store_true")
    parser.add_argument("-s", "--fillsat", help="Fill saturated pixels to match max value", required=False, action="store_true")
    parser.add_argument("-e", "--histeq", help="Apply histogram equalization", required=False, action="store_true")
    parser.add_argument("-x", "--maxpercent", help="Clamp values to maximum percent (0-100)", type=float, default=None)
    parser.add_argument("-n", "--minpercent", help="Clamp values to minimum percent (0-100)", type=float, default=None)
    parser.add_argument("-r", "--resize", help="Resize image to WidthxHeight", type=str, default=None)
    parser.add_argument("-b", "--band", help="Data band", type=int, default=1)
    parser.add_argument("-t", "--trim", help="Trim borders", type=int, default=0)
    parser.add_argument("-o", "--outputformat", help="Output format", type=str, default="uint16")
    parser.add_argument("-X", "--maxvalue", help="Clamp values to maximum value", type=float, default=None)
    parser.add_argument("-N", "--minvalue", help="Clamp values to minimum value", type=float, default=None)
    parser.add_argument("-S", "--std", help="Clamp values to within multiple of standard deviations of mean", type=float, default=None)
    args = parser.parse_args()

    input_files = args.data
    match_intensities = args.intensity
    dohisteq = args.histeq
    fill = args.fill
    fillsat = args.fillsat
    maxpercent = args.maxpercent
    minpercent = args.minpercent
    resize = args.resize
    band = args.band
    trim = args.trim
    outputformat = args.outputformat
    std_mult = args.std

    maxvalue = args.maxvalue
    minvalue = args.minvalue

    input_min = None
    input_max = None

    # TODO: This doesn't account for changes made due to clamping and histogram. Fix that.
    if match_intensities is True:
        data_limits = []
        for input_file in input_files:
            data_min, data_max = get_data_min_max(input_file, band=band, trim=trim)
            data_limits += [data_min, data_max]

        input_min = float(np.array(data_limits).min())
        input_max = float(np.array(data_limits).max())
    else:
        input_min = minvalue
        input_max = maxvalue

    for input_file in input_files:
        sci2tiff(input_file,
                force_input_min=input_min,
                force_input_max=input_max,
                fill_null_stripes=fill,
                fillsat=fillsat,
                dohisteq=dohisteq,
                minpercent=minpercent,
                maxpercent=maxpercent,
                std_mult=std_mult,
                resize=resize,
                band=band,
                trim=trim,
                outputformat=outputformat
                )


