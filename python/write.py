import numpy as np
import struct
import os
import shutil
import tkinter as tk
from tkinter import filedialog

def select_folder():
    root = tk.Tk()
    root.withdraw()
    folder_path = filedialog.askdirectory(
        title='Please select the save location for the dat file.',
    )
    return folder_path if folder_path else None

if __name__ == '__main__':

    name0 = select_folder()
    name='/output'
    #num='5'
    type1='.dat'
    type2='.segy'
    folder=name0+'/dat'
    if not os.path.exists(folder):
        os.makedirs(folder)

    for num in ['1','2','3','4']:
        for d in ['x','y','z']:
            f1=name0+name+num+'_'+d+type1
            f2=name0+name+num+'_'+d+type2
            fm=name0+'/dat/'+name+num+'_'+d+type1
    # Read seismic data
            f = open(f1)
            data = f.readlines()
            nTrace=int(data[0].split(' ')[0])       #Number of seismic traces
            nSample=int(data[0].split(' ')[1])      # Number of sampling points per trace

            S = np.zeros((nTrace, nSample))   # Seismic data
            row_z = 0
            for i in range(nTrace):
                S[row_z, :] = data[i+1].split(' ')[0:nSample]
                row_z += 1
            f.close()
            shutil.move(f1,fm)


        # Write seismic data and add volume header, trace header, and other information.

            fSegy = open("segy11", "rb")    # Open a SEGY format data file containing complete file header information.
            fSegy_size = os.path.getsize("segy11")
            fileHandle = open(f2, 'wb')
            cardHeader = fSegy.read(3200)
            fileHandle.write(cardHeader)
            reelHeader = fSegy.read(400)
            traceHeader = fSegy.read(240)
            traceHeader0 = bytearray(traceHeader)

            reelHeader0 = bytearray(reelHeader)
            reelHeader0[12:14] = struct.pack('<i', nTrace)[0:2]
            reelHeader0[16:18] = struct.pack('<i', 100)[0:2]
            reelHeader0[20:22] = struct.pack('<i', nSample)[0:2]

            tempValue = reelHeader0[20:22]
            decValue = int.from_bytes(tempValue , 'little')
            print("Number of sampling points per trace：")
            print(decValue)

            fileHandle.write(bytes(reelHeader0))


            for i in range(nTrace):
                traceHeader0[0:4] = struct.pack('<i', i + 1)
                traceHeader0[8:12] = struct.pack('<i', i + 1)
                traceHeader0[12:16] = struct.pack('<i', i + 1)
                traceHeader0[114:116] = struct.pack('<i', nSample)[0:2]
                traceHeader0[116:118] = struct.pack('<i', 100)[0:2]
                fileHandle.write(bytes(traceHeader0))
                for j in range(nSample):
                    fileHandle.write(struct.pack('<f', S[i, j]))
            fSegy.close()
            fileHandle.close()

            f = open(f2, "rb")
            f_size = os.path.getsize(f2)
            print("Number of seismic traces：")
            print((f_size - 3600) / (nSample * 4 + 240))
            f.close()