# sci2tiff

Provides conversion and some data handling capabilities for scientific image formats. Mostly this is a script to satisfy requirements of the developer (me) as part of my image processing workflow. 

It has currently been tested on calibrated Cassini VICAR files that are output from Cisscal and with Hubble-derived fits data.

```
usage: sci2tiff.py [-h] -d DATA [DATA ...] [-i] [-f] [-s] [-e] [-x MAXPERCENT]
                   [-n MINPERCENT] [-r RESIZE] [-b BAND] [-t TRIM]

optional arguments:
  -h, --help            show this help message and exit
  -d DATA [DATA ...], --data DATA [DATA ...]
                        Input file(s)
  -i, --intensity       Match intensities when stretching
  -f, --fill            Fill null stripes
  -s, --fillsat         Fill saturated pixels to match max value
  -e, --histeq          Apply histogram equalization
  -x MAXPERCENT, --maxpercent MAXPERCENT
                        Clamp values to maximum percent (0-100)
  -n MINPERCENT, --minpercent MINPERCENT
                        Clamp values to minimum percent (0-100)
  -r RESIZE, --resize RESIZE
                        Resize image to WidthxHeight
  -b BAND, --band BAND  Data band
  -t TRIM, --trim TRIM  Trim borders
```
