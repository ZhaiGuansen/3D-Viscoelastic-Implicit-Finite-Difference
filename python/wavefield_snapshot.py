import numpy as np
import matplotlib.pyplot as plt
import struct
import os

if __name__ == '__main__':
    name0 = './all1/timeshot/2,62/'
    name = 'output'
    type1 = '.dat'
    t0=100
    dt=50
    nt=1000
    for num in range(t0, nt+1, dt):
        name1 = name0 + str(num) + '/'
        #for d in ['x','y','z']:
        for d in ['x']:
            f1 = name1 + name + '_' + d + type1
            # 读取地震数据
            f = open(f1)
            data = f.readlines()
            x = int(data[0].split(' ')[0])  # 地震道数
            y = int(data[0].split(' ')[1])  # 每道采样点数

            S = np.zeros((x, y))  # 地震数据
            row_z = 0
            for i in range(x):
                S[row_z, :] = data[i + 1].split(' ')[0:y]
                row_z += 1
            f.close()

            max = np.max(abs(S)) * 0.8
            max=0.01
            #if num >=800: max=0.02


            plt.contourf(S, levels=np.linspace(-max, max, 101), extend='both', cmap="RdBu_r")
            plt.colorbar()
            plt.title('t='+str(num))
            plt.show()
