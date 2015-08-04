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
    import math
except ImportError:
    raise ImportError("Import failed for standard modules.")

try:
    import dicom
    import uuid
    import argparse

    from dicom.sequence import Sequence
    from dicom.dataset import Dataset
except ImportError:
    raise ImportError(
        "Import failed for Pip imports: dicom, uuid, or argparse.")

try:
    import numpy  # ; numpy.test()
    import scipy  # ; scipy.test()
    if scipy.__version__[2] == 7:
        scipy.pkgload('signal')
        scipy.pkgload('ndimage')
        scipy.pkgload('fftpack')
    else:
        from scipy.fftpack import fftn, ifftn, fftshift, ifftshift
        from scipy import ndimage
        from scipy import signal


except ImportError:
    raise ImportError("Import failed for numpy or scipy.")

print "Python imports successful."
