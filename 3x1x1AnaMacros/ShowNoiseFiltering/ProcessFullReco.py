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


ADC2CHARGE = 67.

run="/Users/sebastienmurphy/311Data/RecoFull/840-1-RecoFull-Parser.root"

#   
f = ROOT.TFile(run)
tree = f.Get("analysistree/anatree")
NEntries= tree.GetEntries()
print NEntries
driftvelocity=0.16*1000
#tree.Print()

first_event=2
last_event=3
iev_itr=0
for iev in range(first_event,last_event) : #loop on events
    tree.GetEntry(iev)
    NTracks=tree.NumberOfTracks
    iev_itr+=1
    print("-------- Event:",iev,"number of tracks:",NTracks, "--------------")
    boxes=[]
    boxes_view1=[]
    index_hit_inside_track = 0
    for itrk in range(NTracks) :
        Nhits_per_view_0=tree.Track_NumberOfHitsPerView[itrk*2]
        Nhits_per_view_1=tree.Track_NumberOfHitsPerView[itrk*2+1]
        Nhits=tree.Track_NumberOfHits[itrk]
        start_x=tree.Track_StartPoint_X[itrk]
        start_y=tree.Track_StartPoint_Y[itrk]
        start_z= tree.Track_StartPoint_Z[itrk]
        end_x=tree.Track_EndPoint_X[itrk]
        end_y=tree.Track_EndPoint_Y[itrk]
        end_z=tree.Track_EndPoint_Z[itrk]
        trk_length=tree.Track_Length_StraightLine[itrk]
        print "track number:", itrk," has: ",Nhits," hits", " and length:",trk_length, "hit vector between:", index_hit_inside_track, "and:", Nhits+index_hit_inside_track
#        if trk_length<90 :
#            index_hit_inside_track+=Nhits
#            continue
        print("------------------", "total:",Nhits, "View0:",Nhits_per_view_0, "View1:",Nhits_per_view_1, iev, index_hit_inside_track, Nhits+index_hit_inside_track,"---------")
        for ihit in range(index_hit_inside_track,Nhits+index_hit_inside_track) :
#        for ihit in range(index_hit_inside_track,10+index_hit_inside_track) :
            hit_x= tree.Track_Hit_X[ihit]
            hit_y= tree.Track_Hit_Y[ihit]
            hit_z= tree.Track_Hit_Z[ihit]
            iview=tree.Track_Hit_View[ihit]
            hit_integral=tree.Hit_ChargeIntegral[ihit]/ADC2CHARGE
            hit_dx=tree.Track_Hit_ds_3DPosition[ihit]
            hit_dx_local=tree.Track_Hit_ds_LocalTrackDirection[ihit]
            hit_dqdx=hit_integral/hit_dx_local
            drift=-(hit_x/driftvelocity)*1000
            hit_start_time= tree.Hit_StartTime[ihit]
            hit_end_time= tree.Hit_EndTime[ihit]
            hit_peak_time= tree.Hit_PeakTime[ihit]
            hit_width= tree.Hit_Width[ihit]
            hit_amplitude= tree.Hit_Amplitude[ihit]
            hit_channel= tree.Hit_Channel[ihit]
            if (math.isnan(hit_integral)): continue
            x1=hit_channel
            x2=hit_channel+1
            y1=hit_start_time
            y2=hit_end_time
            print ihit, iview, x1, x2
            if iview==0:
                rect = Rectangle((x1,y1),1,y2-y1,color="green")
                boxes.append(rect)
#                print "+++ Append View 0:", x1, x2
            if iview==1:
                rect_view1 = Rectangle((x1,y1),1,y2-y1,color="yellow")
                boxes_view1.append(rect_view1)
#                print "---Append View1:", x1, x2
                
        index_hit_inside_track+=Nhits

    Amp=[]
#    print tree.RecoWaveforms_NumberOfChannels
#
#    print 
    tick_itr=0
    ch_itr=0
    for ich in np.array(tree.RecoWaveform_Channel):
        ADC=[]
#        print ich, tree.RecoWaveform_Channel[ich]
        for it in range (tick_itr, tick_itr+tree.RecoWaveform_NTicks[ch_itr]):
            ADC.append(tree.RecoWaveform_ADC[it])
            tick_itr+=1
        #ADC.reverse()
        ch_itr+=1
        Amp.append(ADC)
#      print "append to Amp"
#         print tree.RecoWaveform_ADC[it]
#print ADC[0]
#print Amp[0]

#Amp=np.load("imageFulReco_ev2.npy")
NpAmp=np.array(Amp)

#print boxes
#print NpAmp
#print type(NpAmp)

np.save("imageFulReco_ev2.npy",Amp)

#imgplot=plt.imshow(np.rot90(NpAmp[:320]),cmap=plt.cm.BuPu_r,aspect=320.0/1667.0)
#plt.show()
#pc = PatchCollection(boxes, facecolor="r", alpha=0.5,
#                         edgecolor="None")
pc = PatchCollection(boxes, color="green", alpha=0.9)
pc_v1 = PatchCollection(boxes_view1, color="red", alpha=0.9)
#pc = PatchCollection(boxes)


ax=plt.subplot(311)
plt.imshow(np.rot90(NpAmp[:320]),cmap=plt.cm.bone_r,aspect=320.0/1667.0,clim=[0,20])
#ax.invert_yaxis()
#plt.imshow(NpAmp[:320],cmap=plt.cm.BuPu_r,aspect=1667.0/320.0,clim=[-10,20])
#ax.add_collection(pc)

ax_v1=plt.subplot(312)
plt.imshow(np.rot90(NpAmp[320:]),cmap=plt.cm.BuPu_r,aspect=320.0/1667.0,clim=[-10,20])


# color_list=plt.cm.jet(np.linspace(0,1,320))
n_waveform=10
#color_list=plt.cm.Blues(np.linspace(0,1,100))
color_list=plt.cm.BuGn(np.linspace(0,1,100))
 
ax_w=plt.subplot(313)
for i in range(50, 50+n_waveform):
    plt.plot(np.arange(0,1667,1),smooth(Amp[i],30),color=color_list[i])
#ax_v1.add_collection(pc_v1)
#ax.add_patch(Rectangle((0,0), 100, 100, color="green"))




plt.subplots_adjust(bottom=0.1, right=0.8, top=0.9)
cax = plt.axes([0.85, 0.1, 0.075, 0.8])
plt.colorbar(cax=cax)
plt.show()

#    print tree.RecoWaveforms_NumberOfChannels
#    print tree.RecoWaveform_NTicks[0]
#    print tree.RecoWaveform_Channel[0]
