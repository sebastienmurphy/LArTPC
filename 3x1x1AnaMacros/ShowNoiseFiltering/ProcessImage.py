import os
from glob import glob
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec
import ROOT
from datetime import datetime
from dateutil import tz

data=np.load("image.npy")
imgplot=plt.imshow(data)
plt.show()

#lum_img = data[:, :, 0]
#plt.imshow(lum_img)
#Print data
