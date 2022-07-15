#!/usr/local/bin/python3

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors, ticker,cm
import sys

#--------------------------------------------------------------------------------------------------
def plot_column(col):

    zet  = np.zeros((nt,nk), dtype='d')
    qit  = np.zeros((nt,nk), dtype='d')
    Frim = np.zeros((nt,nk), dtype='d')
    rhoi = np.zeros((nt,nk), dtype='d')
    Di   = np.zeros((nt,nk), dtype='d')
    Dhmax= np.zeros((nt,nk), dtype='d')
    muiB = np.zeros((nt,nk), dtype='d')   # column 13  (interporated from lookupTable)
    qc   = np.zeros((nt,nk), dtype='d')
    qr   = np.zeros((nt,nk), dtype='d')
    Wup  = np.zeros((nt,nk), dtype='d')
    Rliq = np.zeros((nt,nk), dtype='d')   #2D array for convenience only (only value at k=nk is relevant)
    Rsol = np.zeros((nt,nk), dtype='d')   #2D array for convenience only (only value at k=nk is relevant)

# extract from 'data'.  Last index is the column, correponding to order in 'out_p3.dat'
    for t in range(nt):
        Wup[t,:]  = data[col][(t)*nk:(t+1)*nk, 1]
        zet[t,:]  = data[col][(t)*nk:(t+1)*nk, 2]
        qc[t,:]   = data[col][(t)*nk:(t+1)*nk, 3] *1.e+3
        qr[t,:]   = data[col][(t)*nk:(t+1)*nk, 4] *1.e+3
        qit[t,:]  = data[col][(t)*nk:(t+1)*nk, 5] *1.e+3
        muiB[t,:] = data[col][(t)*nk:(t+1)*nk, 6]
        rhoi[t,:] = data[col][(t)*nk:(t+1)*nk, 7]
        Di[t,:]   = data[col][(t)*nk:(t+1)*nk, 8] *1.e+3   #convert to [mm]
        Dhmax[t,:]= data[col][(t)*nk:(t+1)*nk, 9] *1.e+3   #convert to [mm]
        Rliq[t,:] = data[col][(t)*nk:(t+1)*nk, 10]
        Rsol[t,:] = data[col][(t)*nk:(t+1)*nk, 11]
#        Frim[t,:] = data[col][(t)*nk:(t+1)*nk, 7]

    # apply mask:
    mask = qit > 0.001
    Wup  = np.where(Wup>0., Wup, -1.)
    Rliq = np.where(Rliq>0., Rliq, -1.)
    Rsol = np.where(Rsol>0., Rsol, -1.)
    muiB = np.where(muiB>0., muiB, 0.001)
    muiB = np.where(mask, muiB, -99.)
    #zet  = np.where(mask, zet,  -99.)
    Frim = np.where(mask, Frim, -99.)
    Di   = np.where(mask, Di,   -99.)
    Dhmax= np.where(mask, Dhmax,-99.)

    cb_pad = 0.03
    lnwid  = 1.
    lncol  = 'black'

    im0 = ax[0].contourf(range(nt), z_agl, np.transpose(Wup), levels=[0.,1.,2.,3.,4.,5.], vmin=0., vmax=8., cmap='Greys')
    fig.colorbar(im0, ax=ax[0], shrink=0.8, pad=cb_pad)
    ax_b = ax[0].twinx()
    ax_b.plot(range(nt), Rliq, color='red')
    ax_b.plot(range(nt), Rsol, color='blue')
    ax_b.set_ylim([0.,80.])
    ax_b.tick_params(axis="y",direction="in", pad=-22)

    im1 = ax[1].contourf(range(nt), z_agl, np.transpose(qc), levels=levs_Qc,   cmap='Greens')
    fig.colorbar(im1, ax=ax[1], shrink=0.8, pad=cb_pad)
    ax[1].contour( range(nt), z_agl, np.transpose(qc), levels=[0.01,0.01001], linewidths=lnwid, linestyles='dotted', colors='Green')
    ax[1].contour(range(nt), z_agl, np.transpose(qr), levels=levs_Qc,   cmap='Reds')
    ax[1].contour( range(nt), z_agl, np.transpose(qr), levels=[0.01,0.01001], linewidths=lnwid, linestyles='dashed', colors='Red')

    im2 = ax[2].contourf(range(nt), z_agl, np.transpose(qit), levels=levs_Q,   cmap='Blues') #
    fig.colorbar(im2, ax=ax[2], shrink=0.8, pad=cb_pad)
    ax[2].contour( range(nt), z_agl, np.transpose(qit), levels=[0.001,0.001001], linewidths=lnwid, linestyles='dashed', colors=lncol)

    im3 = ax[3].contourf(range(nt), z_agl, np.transpose(Di), levels=levs_Di, colors=ccol_zet)
    fig.colorbar(im3, ax=ax[3], shrink=0.8, pad=cb_pad)
    ax[3].contour( range(nt), z_agl, np.transpose(qit), levels=[0.001,0.001001], linewidths=lnwid, linestyles='dashed', colors=lncol)

    im4 = ax[4].contourf(range(nt), z_agl, np.transpose(zet), levels=levs_Z, colors=ccol_zet)
    fig.colorbar(im4, ax=ax[4], shrink=0.8, pad=cb_pad)


    xloc = 93
    yloc = 12.5
    labsize = 12
    ax[0].text(xloc,yloc,'m s$^{-1}$', size=labsize)
    ax[1].text(xloc,yloc,'g kg$^{-1}$', size=labsize)
    ax[2].text(xloc,yloc,'g kg$^{-1}$', size=labsize)
    ax[3].text(xloc,yloc,'mm',     size=labsize)
    ax[4].text(xloc,yloc,'dBZ',    size=labsize)

    xloc = 5
    yloc = 10.5
    labsize = 18
    ax[0].text(xloc,yloc,'$w$', size=labsize)
    ax[0].text(67,yloc,'$R_{liq}$', size=labsize, color='red')
    ax[0].text(67,yloc-1.8,'$R_{sol}$', size=labsize, color='blue')
    ax[0].text(64,yloc-3.6,'(mm h$^{-1}$)', size=12)
    ax[1].text(xloc,yloc,'$Q_c,Q_r$', size=labsize)
    ax[2].text(xloc,yloc,'$Q_{i,tot}$', size=labsize)
    ax[3].text(xloc,yloc,'$D_m$',       size=labsize)
    ax[4].text(xloc,yloc,'$Z_e$',       size=labsize)


    xaxislabel = 'Time (min)'
    yaxislabel = 'Height (km)'
    if col==0:
        ax[0].yaxis.set_label_text(yaxislabel, size=14)
        ax[1].yaxis.set_label_text(yaxislabel, size=14)
        ax[2].yaxis.set_label_text(yaxislabel, size=14)
        ax[3].yaxis.set_label_text(yaxislabel, size=14)
        ax[4].yaxis.set_label_text(yaxislabel, size=14)

    bottom_ax = ax.shape[0]-1
    ax[bottom_ax].xaxis.set_label_text(xaxislabel, size=14)
    ax[0].set_title(plotTitle[col], size=16, fontweight='bold')

#--------------------------------------------------------------------------------------------------


#-------- read in data from cld1d; store as individual field arrays:
outputDir = './'
path1 = '../'

if len(sys.argv)==1:
    plotTitle0 = 'TITLE'
else:
    plotTitle0 = sys.argv[1]

data0 = np.loadtxt(path1 + 'out_p3.dat')
data = [data0]
plotTitle = [plotTitle0]


nk = 41                         # number of vertical levels
nt = int(data[0].shape[0]/nk)   # number of output times
#z_agl = data[0][0:nk,0]        # array of elevations AGL [m]  (column 0)
z_agl = data[0][0:nk,0]/1000.   # array of elevations AGL [km]  (column 0)


#--------  plotting intervals:
#levs_Z    = [-20., 0., 10., 20., 30., 40., 50., 60., 70., 80.]
levs_Z = [-30,-10,0,5,10,15,20,25,30,35,40,45,50,55,60,65,70]
#levs_Q    = [0.001, 0.01, 0.1, 0.2, 0.3, 0.4, 0.5, 1., 2., 3., 4., 5.,7.,10.]
levs_Q    = [0.1, 0.2, 0.3, 0.4, 0.5, 1., 2., 3., 4., 5.,6.,7.,8.,9.,10.]
levs_Qc   = [0.1, 0.2, 0.3, 0.4, 0.5, 1., 2., 3., 4., 5.]
levs_Fr   = [0., 0.2, 0.4, 0.6, 0.8, 1.0]
levs_Vi   = [0.001, 0.2, 0.5, 1., 2., 3., 4., 5., 7., 10., 15., 20., 25.]
#levs_muB  = [0., 0.5, 1., 2., 3., 4., 5., 10., 15., 20., 25.]
levs_Di   = [1.e-3,1.,2.,3.,4.,5.,10.,15.,20.,25.,30.,35.,40.,45.,50.,60.,70.,80.,90.]
levs_Dh   = levs_Di
levs_rhoi = [0.001, 100., 200., 300., 400., 500., 600., 700., 800., 900., 1000.]
levs_lami = [0.000, 25., 50., 75., 100., 200., 300., 400., 500.]
#levs_mui = [0., 2., 4., 6., 8., 10., 12., 14., 16., 18., 20.]

ccol_zet  = ['lightgrey','darkgrey','grey','lightskyblue','dodgerblue','blue','limegreen','forestgreen',
	     'darkgreen','yellow','orange','salmon','red','darkred','darkviolet','Indigo','darkblue','black']
levs_mui = [0.,1.,2.,3.,4.,5.,6.,7.,8.,9.,10.,11.,12.,13.,14.,15.,16.,17.,99.]

fig, ax = plt.subplots(nrows=5, ncols=1, figsize=(9,14), sharex=True)

plot_column(0)

# adjust subplots
#plt.subplots_adjust(hspace=0.10)
#plt.subplots_adjust(wspace=0.01)

plt.tight_layout()

#--------------------------
plt.savefig('./fig.png', format='png')
plt.show()
