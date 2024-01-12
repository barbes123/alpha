import json
import sys
from os.path import exists
import matplotlib.pyplot as plt
import math
import os.path
from datetime import datetime
time_format = "%Y-%m-%d %H:%M:%S"

import numpy as np
sys.path.append(os.path.abspath('../lib'))
from libLineEq import  *


# import pandas as pd
def SortX(entry):
    return entry['x']

def OpenJsonFromFile(fname):
    try:
        with open('{}'.format(fname),'r') as f:
            js = json.load(f)
            print('File {} is okay'.format(fname))
    except:
        print('no carbon CS file')
    return js

def PrepeareJsonToPlot(js, energy_units = 1e3, norm = 1.0):
    dic_lit = []
    for entry in js:
        datapoint1 = {'x': 0, 'y': 0, 'yerr':0}
        datapoint1['x'] = entry['energy']*energy_units
        datapoint1['y'] = entry['cs']*norm
        datapoint1['y_err'] = entry['cs_err']
        dic_lit.append(datapoint1)
    return dic_lit

def PrepeareJsonToPlotSeries(js, series = 0, norm = 1):
    dic_lit = []
    for entry in js:
        if entry['plotTag'] == series:
            datapoint2 = {'x': 0, 'y': 0, 'err':0, 'y_sub':0, 'y_minus':0, 'y_plus':0, 'label':''}
            datapoint2['x'] = entry['energy_calib']
            datapoint2['y'] = entry['cs']*norm
            datapoint2['err'] = datapoint2['y'] * entry['cs_err']
            datapoint2['y_minus'] = datapoint2['y'] - datapoint2['err']
            datapoint2['y_plus'] = datapoint2['y'] + datapoint2['err']
            datapoint2['label'] = entry['runNumber']
            dic_lit.append(datapoint2)
    print(dic_lit.sort(key=SortX))
    return dic_lit

def GetCS(js, energy):
    cs = 0
    cs_prev = 0
    for carbon in js:
        if energy <= carbon['energy']:
            # cs = carbon['cs']
            cs = abs(carbon['cs'] + cs_prev)/2
            # print(energy, ' ', cs_prev, ' ', cs)
            break
        else:
            cs_prev = carbon['cs']
    return cs

def JsonCSToAscii(js, fname):
    try:
        with open('{}'.format(fname),'w') as fout:
            for element in js:
                fout.write(element['energy'], element['cs'])
    except:
        pass

def ListCSToAscii(list, fname):
    try:
        with open('{}'.format(fname), 'w') as fout:
            for element in list:
                fout.writelines('{} {} \n'.format(element['x'], element['y']))

    except:
        pass

def GetCSfromLit(fname, energy_units): #energy_units = 1e6
    dic_lit = []
    with open('{}'.format(fname), 'r') as fjson:
        js_lit = json.load(fjson)

    dic_lit = []
    for entry in js_lit:
        # if entry['energy'] * 1e6 > 4.6e6:
        #     continue
        datapoint1 = {'x': 0, 'y': 0}
        datapoint1['x'] = entry['energy'] * energy_units
        datapoint1['y'] = entry['cs']
        try:
            datapoint1['y_err'] = entry['cs_err']
            datapoint1['x_err'] = entry['energy_err']
        except:
            datapoint1['y_err'] = 0
            datapoint1['x_err'] = 0
            pass
        dic_lit.append(datapoint1)
    print(dic_lit.sort(key=SortX))
    return dic_lit

def LoadRingRatios(r2r = 12, infile = '../CS/ring/ring_sim.dat'):
    index = 0
    if r2r == 12:
        index = 1
    elif r2r == 13:
        index = 2
    elif r2r == 23:
        index = 3
    else:
        print('wrong ring-2-ring combination')
        return

    list_dic = []

    if not exists(infile):
        print('no ring_sim.dat')
        sys.exit()


    with open('{}'.format(infile), 'r') as fring:
        for line in fring:
            if line[0] == '#' or line.strip() == "":
                continue
            data_line = line.split()
            list = [data_line[0], data_line[index]]
            list_dic.append(list)
        #
        # for el in list_dic:
        #     np.array()

        return list_dic

def get_closest_element(list, colum_nbr, element):
    col = list[:, colum_nbr]
    difference_array = np.absolute(col - element)
    pos = difference_array.argmin()
    return list[pos][0], col[pos], pos

def RingRing2Energy(list, r2r = 12, ratio = 1):
    # print(list)
    index = 0
    if r2r == 12:
        index = 1
    elif r2r == 13:
        index = 2
    elif r2r == 23:
        index = 3
    else:
        print('wrong ring-2-ring combination')
        return
    # print(list)
    #Find First min
    col = list[:, index]
    energy, value, pos = get_closest_element(list,index,ratio)
    p1 = (energy, value)

    #Find second min
    col = np.delete(col, pos)
    energy, value, pos = get_closest_element(list, index, ratio)
    p2 = (energy, value)
    points = [p1, p2]
    pol = LineEquation(points)
    # print ((ratio - pol[0])/pol[1])
    return (ratio - pol[0])/pol[1]

def convert_json_to_np_array(js, keys):
    list_dic = []
    for datum in js:
        list = []
        for key in keys:
            list.append(datum[key])
        list_dic.append(list)
    return np.array(list_dic, dtype='float32')

def get_contaminant_cs(list, ener):
    energy, value, pos = get_closest_element(list, 0, ener)
    # print(energy, list[pos][1], pos)
    p1 = (energy, list[pos][1])

    list = np.delete(list, pos, 0)
    energy, value, pos = get_closest_element(list, 0, ener)
    p2 = (energy, list[pos][1])
    # print(energy, list[pos][1], pos)


    pol = LineEquation([p1, p2])
    ## print(pol)
    ## print (pol[0] + pol[1]*ener)
    return (pol[0] + pol[1]*ener)


def SetUpData_2022_w39(run):
    if run['source'] == 'alpha' and run['status'] == 'data' and run['target'] == 'N1-CaF':
       return 1, 'series1'
       # return 1, 'N1-1'
    elif run['source'] == 'alpha' and run['status'] == 'data' and run['target'] == 'N2-CaF' and run['runNumber'] >= 525 and run['runNumber'] <= 840:
    # elif run['status'] == 'data' and run['target'] == 'N2-CaF' and run['energy'] < 5000 and run['energy'] > 2250:
        return 2, 'series2'
    elif run['source'] == 'alpha' and run['status'] == 'data' and run['target'] == 'N2-CaF' and run['energy'] < 2250 and run['energy'] > 1700:
        return 3, 'N2-3'
    elif run['source'] == 'alpha' and run['status'] == 'data' and run['target'] == 'N2-CaF' and run['energy'] < 1700:
        return 4, 'N2-4'
    # elif run['status'] == 'data' and run['target'] == 'N2-CaF' and run['energy'] > 5000:
    elif run['source'] == 'alpha' and run['status'] == 'data' and run['target'] == 'N2-CaF' and run['runNumber'] >= 979 and run['runNumber'] <= 1273:
        return 5, 'series5'
    elif run['source'] == 'alpha' and run['status'] == 'data' and run['target'] == 'N8-Cu12mg' and run['runNumber'] >= 1274 and run['runNumber'] <= 1353:
        return 6, 'series6'
    elif run['source'] == 'alpha' and run['status'] == 'data' and run['target'] == 'N5-C13' and run[
        'runNumber'] >= 105 and run['runNumber'] <= 173:
        return 7, 'series7'
    elif run['source'] == 'alpha' and run['status'] == 'data' and run['target'] == 'N7-C13' and run[
        'runNumber'] >= 174 and run['runNumber'] <= 256:
        return 8, 'series8'
    elif run['source'] == 'alpha' and run['status'] == 'data' and run['target'] == 'N5-C13' and run[
        'runNumber'] >= 64 and run['runNumber'] <= 104:
        return 9, 'series9'
    elif run['source'] == 'alpha' and run['status'] == 'data' and run['target'] == 'N3-CaF':
        return 99, 'none'
    else:
         return 0, 'None'

def SetUpData_2023_w43(run):

    if run['status'] == 'data' and run['target'] == 'N3-13C' and run['runNumber'] <= 111:
        return 1, 'series1'
    elif run['status'] == 'data' and run['target'] == 'N3-13C' and run['runNumber'] >= 112 and run['runNumber'] <= 209:
        return 2, 'series2'
    elif run['runNumber'] >= 210 and run['runNumber'] <= 222:
        return 3, 'series3'
    elif run['runNumber'] >= 223 and run['runNumber'] <= 249:
        return 4, 'series4'
    elif run['runNumber'] >= 250 and run['runNumber'] <= 341:
        return 5, 'series5'
    elif run['runNumber'] >= 342 and run['runNumber'] <= 445:
        return 6, 'series6'
    elif run['runNumber'] >= 446 and run['runNumber'] <= 472:
        return 7, 'series7'
    elif run['status'] == 'data' and run['runNumber'] >= 473 and run['runNumber'] <= 672:
        return 8, 'series8'
    elif run['status'] == 'data' and run['runNumber'] >= 673 and run['runNumber'] <= 1000:
        return 9, 'series9'
    else:
         return 0, 'None'

# def SetUpData_2023_w17(run):
#     # print('i am in here')
#     # print(run['runNumber'])
#     if run['runNumber'] == 113 or run['runNumber'] == 114 or run['runNumber'] == 115 or run[
#         'runNumber'] == 116:
#         # run['status'] = 3
#         return 3, '3'
#     elif run['runNumber'] == 228 or run['runNumber'] == 229 or run['runNumber'] == 230 or run[
#         'runNumber'] == 231:
#         # run['status'] = 5
#         return 5, '5'
#     elif run['runNumber'] >= 47 and run['runNumber'] <= 54:
#         # run['status'] = 1
#         return 1, '1'
#     elif run['runNumber'] >= 55 and run['runNumber'] <= 141:
#         # run['status'] = 2
#         return 2, '2'
#     elif run['runNumber'] >= 144 and run['runNumber'] <= 226:
#         # run['status'] = 4
#         return 4, '4'
#     elif run['runNumber'] > 226 and run['runNumber'] < 300 and run['runNumber'] != 242:
#         # run['status'] = 5
#         return 5, '5'
#     else:
#         # run['status'] = -1
#         return -1, '-1'

    #
    # if run['status'] == 1:
    #     return 1, '1'
    # elif run['status'] == 2:
    #     return 2, '2'
    # elif run['status'] == 3:
    #     return 3, '3'
    # elif run['status'] == 4:
    #     return 4, '4'
    # elif run['status'] == 5:
    #     return 5, '6'
    # elif run['status'] == 6:
    #     return 6, '6'
    # else:
    #      return 0, 'None'
#
# R2R = np.loadtxt('../CS/ring/ring_sim.dat')
#
# print(RingRing2Energy(R2R, 12, 2.97154))
# print(RingRing2Energy(R2R, 13,  2.91841))
# print(RingRing2Energy(R2R, 23, 0.982121))


# R12 = LoadRingRatios(12)
# print(R12)
# RinrRing2Energy(LoadRingRatios(12))

def SetUpData_2023_w29(run):
    if run['runNumber'] >= 316 and run['runNumber'] <= 502 and run['source'] == 'alpha' and run['target'] == '13C':
        return 4, 'series4'
    elif run['runNumber'] >= 512 and run['runNumber'] <= 736 and run['source'] == 'alpha' and run['target'] == 'N1-27Al':
        return 1, 'series1'
    elif run['runNumber'] >= 737 and run['runNumber'] <= 880 and run['source'] == 'alpha' and run['target'] == 'N1-27Al' and run['energy'] > 4450:
        return 2, 'series2' # REMOVE ENERGY CONDITION IT WAS FOR IAEA
    elif run['runNumber'] >= 882 and run['runNumber'] <= 964 and run['source'] == 'alpha' and run['target'] == 'N1-27Al':
        return 3, 'series3'
    # elif run['status'] == 'data1' and run['target'] == 'N1-27Al':
    #      return 1, 'series1'
    # elif run['status'] == 'data2' and run['target'] == 'N1-27Al':
    #     return 2, 'series2'
    # elif run['status'] == 'data' and run['target'] == '13C':
    #     return 3, 'series3'
    else:
        return 0, 'None'


def SetUpData_2023_w17(run):
    if run['runNumber'] >= 241 and run['runNumber'] <= 300 and run['source'] == 'alpha' and run['target'] == 'N7-27Al' and run['runNumber'] != 299 and run['runNumber'] != 242 and run['energy']>4500:
        return 6, 'series6'
    # elif run['runNumber'] >= 512 and run['runNumber'] <= 736 and run['source'] == 'alpha' and run['target'] == 'N1-27Al':
    #     return 1, 'series1'
    # elif run['runNumber'] >= 737 and run['runNumber'] <= 880 and run['source'] == 'alpha' and run['target'] == 'N1-27Al':
    #     return 2, 'series2'
    # elif run['runNumber'] >= 882 and run['runNumber'] <= 964 and run['source'] == 'alpha' and run['target'] == 'N1-27Al':
    #     return 3, 'series3'
    # elif run['status'] == 'data1' and run['target'] == 'N1-27Al':
    #      return 1, 'series1'
    # elif run['status'] == 'data2' and run['target'] == 'N1-27Al':
    #     return 2, 'series2'
    # elif run['status'] == 'data' and run['target'] == '13C':
    #     return 3, 'series3'
    else:
        return 0, 'None'

def GetTargetDegradation(js, run1, run2, nm1, nm2):
    bl1 = False
    bl2 = False
    for el in js:
        if el['runNumber'] == run1:
            time1 = datetime.strptime(el['tstart'], time_format)
            bl1 = True
        elif el['runNumber'] == run2:
            time2 = datetime.strptime(el['tstop'], time_format)
            bl2 = True
    if not bl1 or not bl2:
        print('time from run {} or run {} is not available'.format(run1, run2))
        return 0.0

    dtime = (time2 - time1).total_seconds()
    p1 = (0, nm1)
    p2 = (dtime, nm2)
    pol = LineEquation([p1, p2])
    # print(pol)
    # print (pol[0] + pol[1]*ener)
    # return (pol[0] + pol[1] * ener)
    # print (pol[0] + pol[1] * ener)
    return time1, pol
    # dtime = (time2 - time1).total_seconds()
    # nm = (nm2-nm1)/dtime
    # print(nm)

def GetCurrentDegradationFactor(timeZero, pol, time_str):
    time = datetime.strptime(time_str, time_format)
    # time0 = datetime.strptime(timeZero, time_format)
    current_time = (time - timeZero).total_seconds()
    # print(pol[0]+pol[1]*current_time)
    return pol[0]+pol[1]*current_time
    # print((time2 - time1).total_seconds())

    return nm #normalization gradient per second

def SubtractionCS(js_data, js_bg):
    list_data = js_data.dumps()


    def get_energy_col(js):
        list_ene = []
        for element in js:
            list_ene.append(element['energy'])
        return list_ene

    list_ene_data = []
    list_ene_bg = []
    list_ene_data = get_energy_col(list_data)
    list_ene_bg = get_energy_col(list_bg)

    for energy in list_ene_data:
        closest_value = min(list_ene_bg, key=lambda x: abs(energy - x))
        # print('{} : {}'.format(closest_value, energy))

def GetTaCS(js_ta_cs, ener):
    list_ene = []
    for element in js_ta_cs:
        list_ene.append(element['energy'])
    closest_energy = min(list_ene, key=lambda x: abs(ener - x))

    NeutronsFromTa = 0
    for element in js_ta_cs:
        if element['energy'] == closest_energy:
            NeutronsFromTa = element['cs'] * element['thickness']
            break

    return  NeutronsFromTa
    # return closest_value