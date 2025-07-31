import numpy as np
import matplotlib.pyplot as plt
import time
import struct
import os
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
    name = 'output'
    type1 = '.dat'
    t0=40
    dt=20
    nt=200
    for num in range(t0, nt+1, dt):
        name1 = name0 + '/' + str(num) + '/'
        #for d in ['x','y','z']:
        for d in ['z']:
            f1 = name1 + name + '_' + d + type1
            # Read seismic data
            f = open(f1)
            data = f.readlines()
            x = int(data[0].split(' ')[0])  # Number of seismic traces
            y = int(data[0].split(' ')[1])  # Number of sampling points per trace

            S = np.zeros((x, y))  # Seismic data
            row_z = 0
            for i in range(x):
                S[row_z, :] = data[i + 1].split(' ')[0:y]
                row_z += 1
            f.close()

            max = np.max(abs(S)) * 0.8
            #max=0.05
            #if num >=800: max=0.02


            plt.contourf(S.T, levels=np.linspace(-max, max, 101), extend='both', cmap="gray")
            plt.colorbar()
            plt.title('t='+str(num))
            plt.show()
            time.sleep(0.5)
            plt.close()
