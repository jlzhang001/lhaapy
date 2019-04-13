import matplotlib.pyplot as plt
import numpy as np
import astropy.units as u
from astropy.coordinates import SkyCoord
from lhaapy.irf import EffectiveAreaTable2D
from lhaapy.maps import WcsGeom, MapAxis, WcsNDMap
from lhaapy.spectrum.models import PowerLaw
from lhaapy.image.models import SkyGaussian
from lhaapy.utils.random import get_random_state
from lhaapy.cube import make_map_exposure_true_energy, MapFit, MapEvaluator
from lhaapy.cube.models import SkyModel

filename = "$LHAASO_DATA/wcda/data/1.fits"

aeff = 
