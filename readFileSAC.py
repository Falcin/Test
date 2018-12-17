from obspy import *
import datetime
import matplotlib.pyplot as plt
import numpy as np
import os
import glob
import scipy.signal as sg



files=[]
files.extend(glob.glob("*HHZ*"))
for i in range(len(files)):
    print("Tapez ",i," pour selectionner le fichier : ",files[i])

index=int(input())

print(files[index])

s=read(files[index])

signal=s[0].data
stats=s[0].stats
fs=stats['sampling_rate']
starttime=stats['starttime'].datetime
endtime=stats['endtime'].datetime
npts=stats['npts']
time_vect=[starttime+datetime.timedelta(microseconds=i*(1/fs*1000000)) for i in range(len(signal))]

fmax=20
ratioDecimate=int(fs/(2*fmax))
w_size=250
w=sg.kaiser(w_size,18)
w=w*w_size/ sum([pow(j,2)for j in w])
freq,time,spec=sg.spectrogram(sg.decimate(signal,ratioDecimate,zero_phase=False),
                             fs=fs/ratioDecimate,
                             nperseg=w_size,
                             noverlap=0.96*w_size,
                             nfft=1.5*w_size,
                             window=w,
                             scaling="density")
specDB=10*np.log10(spec)

f=plt.figure(figsize=(15,8))
f.subplots_adjust(hspace=0.25, top=0.95)
ax1 = f.add_subplot(311)
ax2 = f.add_subplot(312)
ax3 = f.add_subplot(313)
ax1.plot(time_vect,signal)
ax1.set_xlabel("Time(s)")
ax1.set_ylabel("Amplitude")
ax2.quadmesh=ax2.pcolormesh(time,freq,(specDB),shading="gouraud", vmax=(np.amax(specDB)))
ax2.set_xlabel("Time (s)")
ax2.set_ylabel('Frequency (Hz)')
#colorbar(ax2, orientation=horizontal)
#ax2.set_title("Spectrogram of an "+ligneCatalogue["class"]+" event")
#quadmesh.set_clim(0,np.amax(specDB)/1)
ax3.plot(freq,np.sum(abs(specDB),axis=1))
ax3.set_xlabel("Frequency (Hz)")
ax3.set_ylabel("Energy")

plt.show()
#plt.figure(figsize=(15,7))
#plt.plot(time_vect,signal,linewidth=0.1)
#plt.xticks(rotation=30)
#plt.show()
