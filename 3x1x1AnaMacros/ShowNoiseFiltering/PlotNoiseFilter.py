import ROOT
import matplotlib.pyplot as plt
from ROOT import gROOT, TCanvas, TF1, gApplication
from ROOT import TFile, TColor
from ROOT import TH1F,TH2F,TH3F
from ROOT import kBlack, kBlue, kRed, kCyan, kMagenta
import numpy as np
import sys
import math
import matplotlib.pyplot as plt
from ROOT import  TPolyMarker3D
import getopt
import os  
from cStringIO import StringIO # Python 2
import argparse
from ROOT import TEveManager, TRandom
from ROOT import TEveViewer,TEveBox,TGLViewer,TEveLine,TEveBox,TEveShape,TEveElement
from ROOT import TPolyMarker3D, TPolyLine3D
from ROOT import TEvePointSet,TEveLine,TEveBoxSet
from ROOT import TGLVector3,TGLOrthoCamera,TEveGeoTopNode,TGeoNode
from ROOT import TEveRGBAPalette, TEveRGBAPaletteOverlay,TLatex
from ROOT import TGeoManager as gGeo

from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle

def smooth(y, box_pts):
    box = np.ones(box_pts)/box_pts
    y_smooth = np.convolve(y, box, mode='same')
    return y_smooth


NpAmp=np.load("NpyFigs/imageFulReco_ev2.npy")
NpAmpRot= np.rot90(NpAmp[:320])
NpAmpRotFlip=np.flipud(NpAmpRot)

NpAmpNoisy=np.load("NpyFigs/imageNoisy_ev2.npy")
NpAmpNoisyRot=np.rot90(NpAmpNoisy[:320])
NpAmpNoisyRotFlip=np.flipud(NpAmpNoisyRot)
#ax=ax.imshow(np.rot90(NpAmpNoisy[:320]),cmap=plt.cm.bone_r,aspect=320.0/1667.0,clim=[0,20],origin="upper")
#ax.invert_yaxis() 
#plt.show()


#NpAmpNoisy=np.array(AmpNoisy)
#


## t = np.arange(0.0, 2.0, 0.01)

## s1 = np.sin(2 * np.pi * t)
## s2 = np.exp(-t)
## s3 = s1 * s2
FontSize=18

fig, axs = plt.subplots(1, 3, sharey=False, facecolor="white", figsize=(15,5))

# Plot each graph, and manually set the y tick values
axs[0].imshow(NpAmpNoisyRotFlip,cmap=plt.cm.bone_r,aspect=315.0/(1667.0*.4),clim=[0,20], origin='upper', interpolation='nearest',extent=[15,320,0,1667*.4])
axs[0].set_xlabel('channel numnber', fontsize=FontSize)
axs[0].set_ylim(5,600)
axs[0].set_xlim(25,300)
axs[0].set_ylabel('drift time (us)', fontsize=FontSize)
#plt.xlim(0,320)
#Axis=plt.gca()
#XYaxis.set_x_ticks()

#axs[0].xlim(0,320)
axs[1].imshow(NpAmpRotFlip,cmap=plt.cm.Reds,aspect=315.0/(1667.0*.4),clim=[0,20],extent=[15,320,0,1667*.4])
axs[1].set_xlabel('channel numnber', fontsize=FontSize)
axs[1].set_ylim(5,600)
axs[1].set_xlim(25,300)
Axis=plt.gca()
#Axis
#axs[2].imshow(np.rot90(NpAmp[:320]),cmap=plt.cm.bone_r,aspect=320.0/1667.0,clim=[0,20])

color_list=plt.cm.BuGn(np.linspace(0,1,100))
n_waveform=1
for i in range(190, 190+n_waveform):axs[2].plot(np.arange(0,1667,1)*.4,smooth(NpAmpNoisy[i],10),color="gray",linewidth=2.0)
for i in range(177, 177+n_waveform):axs[2].plot(np.arange(0,1667,1)*.4,smooth(NpAmp[i],10),color="red")
axs[2].set_xlabel('drift time (us)', fontsize=FontSize)
axs[2].set_ylabel('amplitude (ADC)',fontsize=FontSize)
plt.ylim(-5,28)
plt.xlim(5,600)

fig.add_subplot()

plt.subplots_adjust(bottom=0.15, wspace=0.2)
plt.show()
