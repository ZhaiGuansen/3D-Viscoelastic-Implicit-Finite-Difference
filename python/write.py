import numpy as np
import struct
import os

if __name__ == '__main__':
    name0 = './all3/'
    name='output'
    #num='5'
    type1='.dat'
    type2='.segy'
    for num in ['1','2','3','4','5']:
        for d in ['x','y','z']:
            f1=name0+name+num+'_'+d+type1
            f2=name0+name+num+'_'+d+type2
    # 读取地震数据
            f = open(f1)
            data = f.readlines()
            nTrace=int(data[0].split(' ')[0])       #地震道数
            nSample=int(data[0].split(' ')[1])      # 每道采样点数

        #nTrace=data[0].split('')[0]       #地震道数
        #nSample=data[1].split('')[0]   # 每道采样点数
            S = np.zeros((nTrace, nSample))   # 地震数据
            row_z = 0
            for i in range(nTrace):
                S[row_z, :] = data[i+1].split(' ')[0:nSample]
                row_z += 1
            f.close()

        # 写入地震数据，并添加卷头、道头等信息
        # nTrace = 2301
        # nSample = 751
            fSegy = open("segy11", "rb")    # 打开一份包含完整文件头信息的segy格式数据
            fSegy_size = os.path.getsize("segy11")    # 计算文件所占字节数
        # print((fSegy_size - 3600) / (400 * 4 + 240))    # 计算地震道数
            fileHandle = open(f2, 'wb')
            cardHeader = fSegy.read(3200)    # 读取3200字EBCDIC文件头
            fileHandle.write(cardHeader)     # 写入EBCDIC文件头
            reelHeader = fSegy.read(400)     # 读取400字节二进制文件头
            traceHeader = fSegy.read(240)
            traceHeader0 = bytearray(traceHeader)

        # 修改相关参数
            reelHeader0 = bytearray(reelHeader)
            reelHeader0[12:14] = struct.pack('<i', nTrace)[0:2]    # 修改地震道数
            reelHeader0[16:18] = struct.pack('<i', 100)[0:2]      # 修改采样时间间隔 单位为微秒(us)
            reelHeader0[20:22] = struct.pack('<i', nSample)[0:2]   # 修改地震道道采样点数

        # 测试修改是否正确
            tempValue = reelHeader0[20:22]
            decValue = int.from_bytes(tempValue , 'little')
            print("每道采样点数：")
            print(decValue)

            fileHandle.write(bytes(reelHeader0))    # 写入二进制文件头

        # 写入数据体
            for i in range(nTrace):
                # fSegy.seek(240 * i + 3600 + nSample * 4 * i, 0)    # 每一道都跳过240字节道头
                traceHeader0[0:4] = struct.pack('<i', i + 1)
                traceHeader0[8:12] = struct.pack('<i', i + 1)
                traceHeader0[12:16] = struct.pack('<i', i + 1)
                traceHeader0[114:116] = struct.pack('<i', nSample)[0:2]  # 数据道采样点数
                traceHeader0[116:118] = struct.pack('<i', 100)[0:2]     # 微秒（us）形式的采样间隔
                fileHandle.write(bytes(traceHeader0))
                for j in range(nSample):
                    fileHandle.write(struct.pack('<f', S[i, j]))
            fSegy.close()
            fileHandle.close()

            f = open(f2, "rb")
            f_size = os.path.getsize(f2)    # 计算文件所占字节数
            print("地震道数：")
            print((f_size - 3600) / (nSample * 4 + 240))    # 计算地震道数
            f.close()