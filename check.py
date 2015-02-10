#!/usr/bin/env python

try:
    import pdb
# import ast
    import os
    import sys
    import datetime
    import dateutil
    import dateutil.tz
    import re

    import dicom
    import uuid
    import math
    import numpy  # ; numpy.test()
    import argparse

    from dicom.sequence import Sequence
    from dicom.dataset import Dataset

    import scipy  #  ; scipy.test()
    if scipy.__version__[2] == 7:
        scipy.pkgload('signal')
        scipy.pkgload('ndimage')
        scipy.pkgload('fftpack')
    else:
        from scipy.fftpack import fftn, ifftn, fftshift, ifftshift
        from scipy import ndimage
        from scipy import signal


except ImportError:
    raise ImportError("Import failed.")

print "Python imports successful."
