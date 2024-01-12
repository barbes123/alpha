import json
import sys
from os.path import exists
from os.path import abspath
import matplotlib.pyplot as plt
import math
import subprocess
import shutil

sys.path.append(abspath('../lib/'))

from libCS import *
from libCS import ListCSToAscii
# import pandas as pd
# import joblib
# import pickle
# import crosshair
import matplotlib.backend_bases
import mouse
import numpy as np
import matplotlib.patches as patches

pol = [0, 1]
# pol = [-80.902859308, 1.02559837796973]
pol_LE = [-75.48320, 1.024358] #not so bad for energy 4180 for #Balakrishnan
# pol = [-22.046391, 1.01365979] #Balakrishnan only 3 peaks
# pol = [-67.270729, 1.0229244439] #Peters, 3 peaks, different target, high energy
pol_HE = [9.79489, 1.0087671053] #Peters, 4 peaks, different target, high perfect energy
# pol = [-69.0368242657, 1.022670384004826] #Peters + Balakrishnan all points#not perfect

# pol_ser2=[-294.926363831903, 1.07482926894846]
pol_ser2=[-90.971638, 1.0279958891]
# pol_ser2 = [0, 1]
pol_ser7_ser8 = [-95.955139488, 1.0214478885781]#13C
# pol = [0, 1]

blSubtractContaminants  = False
blErrorsCS              = True
blCrossHair = False
blPlotHighEnergy = False
blPlotOxygenCS = False
blPlotCarbonCS = True
blPlotRunLables = True


carbon_fraction = 0.03
oxygen_fraction = 0.00

errNeutronEff = 0.07
errTarget = 0.05
neutron_eff = 0.37
barn = 1e24 #cm2

normalization = 1#8
ax = plt
if blCrossHair:
    fig, ax = plt.subplots()
    ax = fig.add_axes([0, 0, 1, 1])

fileJson = 'json/new_my_log.json'

CSfolder = '../CS'

fileCS_13C = '../CS/CS13C/Harispolus.json'
fileCS_18O_kunz = '../CS/CS18O/kunz.json'
fileCS_18O_bair1973 = '../CS/CS18O/bair1973.json'
fileCS_18O_bair1962 = '../CS/CS18O/bair1962.json'
fileCS_18O_hansen = '../CS/CS18O/hansen.json'

js_carbon = OpenJsonFromFile(fileCS_13C)
js_oxygen_kunz = OpenJsonFromFile(fileCS_18O_kunz)
js_oxygen_bair1973 = OpenJsonFromFile(fileCS_18O_bair1973)
js_oxygen_bair1962 = OpenJsonFromFile(fileCS_18O_bair1962)
js_oxygen_hansen = OpenJsonFromFile(fileCS_18O_hansen)

class BlittedCursor:
    """
    A cross-hair cursor using blitting for faster redraw.
    """
    def __init__(self, ax):
        self.ax = ax
        self.background = None
        self.horizontal_line = ax.axhline(color='k', lw=0.8, ls='--')
        self.vertical_line = ax.axvline(color='k', lw=0.8, ls='--')
        # text location in axes coordinates
        self.text = ax.text(0.02, 0.55, '', transform=ax.transAxes)
        self._creating_background = False
        ax.figure.canvas.mpl_connect('draw_event', self.on_draw)

    def on_draw(self, event):
        self.create_new_background()


    def set_cross_hair_visible(self, visible):
        need_redraw = self.horizontal_line.get_visible() != visible
        self.horizontal_line.set_visible(visible)
        self.vertical_line.set_visible(visible)
        self.text.set_visible(visible)
        return need_redraw

    def create_new_background(self):
        if self._creating_background:
            # discard calls triggered from within this function
            return
        self._creating_background = True
        self.set_cross_hair_visible(False)
        self.ax.figure.canvas.draw()
        self.background = self.ax.figure.canvas.copy_from_bbox(self.ax.bbox)
        self.set_cross_hair_visible(True)
        self._creating_background = False

    def on_mouse_move(self, event):
        if self.background is None:
            self.create_new_background()
        if not event.inaxes:
            need_redraw = self.set_cross_hair_visible(False)
            if need_redraw:
                self.ax.figure.canvas.restore_region(self.background)
                self.ax.figure.canvas.blit(self.ax.bbox)
        else:
            self.set_cross_hair_visible(True)
            # update the line positions
            x, y = event.xdata, event.ydata
            self.horizontal_line.set_ydata([y])
            self.vertical_line.set_xdata([x])
            self.text.set_text('x=%1.2f, y=%1.2f' % (x, y))

            self.ax.figure.canvas.restore_region(self.background)
            self.ax.draw_artist(self.horizontal_line)
            self.ax.draw_artist(self.vertical_line)
            self.ax.draw_artist(self.text)
            self.ax.figure.canvas.blit(self.ax.bbox)


def AddLabels(list):
    for i in list:
        plt.text(i['x'], i['y'], i['x'])
def GetOxygenCS(energy):
    if energy <= 1970:
        return GetCS(js_oxygen_kunz, energy)
    elif energy >1970 and energy<= 2504:
        return GetCS(js_oxygen_bair1973, energy)
    elif energy > 2504 and energy <= 5140:
        return GetCS(js_oxygen_bair1962, energy)
    elif energy > 5140 and energy <= 12500:
        return GetCS(js_oxygen_hansen, energy)
    else:
        print('Energy {} is not found in carbon CS'.format(energy))
        return -999
# def SortX(entry):
#     return entry['x']

if not exists(fileJson):
    print('no fileJson')
    sys.exit()

with open('{}'.format(fileJson),'r') as fjson:
    js_data =json.load(fjson)




def PlotDataLiterature(datum, strColor = 'r', strStyle = 'None', blErr = True, strLabel = 'label', ptSize=20, ptLinewidth=3 ):
    ax.scatter([i['x'] for i in datum], [i['y'] for i in datum], s=ptSize, color=strColor, label = strLabel)
    ax.plot([i['x'] for i in datum], [i['y'] for i in datum], color=strColor, linestyle=strStyle, linewidth=ptLinewidth)
    if blErrorsCS:
        plt.errorbar([i['x'] for i in datum], [i['y'] for i in datum], [i['y_err'] for i in datum], color=strColor, linestyle=strStyle)
def PlotData(plot_here, data, strColor = 'r', strStyle = 'None', blErr = True, strLabel = 'label', ptSize=20, ptLinewidth=3, blLabel = False):
    # ln, = plot_here.plot
    plot_here.plot([i['x'] for i in data], [i['y_sub'] for i in data], color=strColor, linestyle=strStyle)
    plot_here.scatter([i['x'] for i in data], [i['y'] for i in data], s=ptSize, color=strColor, label=strLabel)
    plot_here.plot([i['x'] for i in data], [i['y'] for i in data], color=strColor)
    # plot_here.text([i['x'] for i in data], [i['y'] for i in data], 'mytext')
    # plot_here.text([i['x'] for i in data)]  .text(i, j, f'({i}, {j})') for (i, j) in zip(x, y)]
    # [plot_here.text(i['x'],i['y'], 'mytext') for (i['x']) in data]
    if blLabel:
        for i in data:
            plot_here.text(i['x'],i['y'], i['label'])

    if blErr:
        # uncomment line below to draw error bars
        # plot_here.errorbar([i['x'] for i in data], [i['y'] for i in data], [i['err'] for i in data], color=strColor)
        plot_here.fill_between([i['x'] for i in data], [i['y_minus'] for i in data],
                         [i['y_plus'] for i in data], color=strColor, alpha=0.2)
    return plot_here

########################################
time0_ser1, pol_ser1 =  GetTargetDegradation(js_data, 519, 409, 7, 7) #time0_series1 - 'zero' time of run 409 Ecalib=5235
# time0_ser1, pol_ser1 =  GetTargetDegradation(js_data, 266, 519, 7, 7) #time0_series1 - 'zero' time of run 266
time0_2, pol_2 =  GetTargetDegradation(js_data, 525, 840, 6, 6) #time0_2 - 'zero' time of run 525
# print(time0_2)



# sys.exit()


#Calculations of CS
for run in js_data:
    if run['status'] != 'data':
        continue
    # print(run['status'])
    # sys.exit()

    # if run['energy'] > 5000:

    # if run['target'] == 'N1-CaF' or (run['target'] == 'N2-CaF' and run['energy'] > 5000):
    #     pol = pol_HE
    # else:
    #      pol = pol_LE

    # pol = pol_LE
    run['energy_calib'] = round(pol[0]+run['energy']* pol[1], 2)

    # if run['target'] == 'N2-CaF':
    #     run['cs_norm'] = 1.0

    blNonZero = True
    if run['beamParticles'] == 0:
        print(run['runNumber'], 'beamParticles is zero')
        blNonZero = False
    if run['beam'] == 0:
        print(run['beam'], 'beam is zero')
        blNonZero = False
    if run['thickness'] == 0:
        print(run['runNumber'], 'thickness is zero')
        blNonZero = False
    if run['charge'] == 0:
        print(run['runNumber'], 'charge is zero')
        blNonZero = False

    plot_index = -1;
    data_tag = ''
    plot_index, data_tag = SetUpData_2022_w39(run)
    run['plotTag'] = data_tag

    if run['plotTag'] == 'series1':
        run['cs_norm'] = 10
        # if run['runNumber'] <= 409:
        #     run['cs_norm'] = 7
        # else:
        #     run['cs_norm'] = GetCurrentDegradationFactor(time0_ser1, pol_ser1, run['tstart'])
    elif run['plotTag'] == 'series2':
        # run['cs_norm'] = 4
        run['cs_norm'] = GetCurrentDegradationFactor(time0_2, pol_2, run['tstart'])
        # run['cs_norm'] = 6
    elif run['plotTag'] == 'N2-3':
        pass
    elif run['plotTag'] == 'N2-4':
        pass
    elif run['plotTag'] == 'series5':
        run['cs_norm'] = 6.2
    elif run['plotTag'] == 'series6':
        run['cs_norm'] = 1
    elif run['plotTag'] == 'series7' or run['plotTag'] == 'series8':
        run['cs_norm'] = 5.2
    elif run['plotTag'] == 'N3-X':
        pass
    else:
        pass

    if run['plotTag'] == 'series7' or run['plotTag'] == 'series8':
        pol = pol_ser7_ser8
    elif run['plotTag'] == 'series2':
        pol = pol_ser2
    elif run['target'] == 'N1-CaF' or (run['target'] == 'N2-CaF' and run['energy'] > 5000):
        pol = pol_HE
    else:
         pol = pol_LE



    if blNonZero:

        NeutronsFromCarbon = 0
        NeutronsFromOxygen = 0

        if blSubtractContaminants:
            carbon_thickness = run['thickness'] * carbon_fraction
            oxygen_thickness = run['thickness'] * oxygen_fraction
            NeutronsFromCarbon = GetCS(js_carbon, run['energy_calib']) * 1e-24 * run['beamParticles'] * carbon_thickness
            NeutronsFromOxygen = GetOxygenCS(run['energy_calib']) * 1e-24 * run['beamParticles'] * oxygen_thickness
            # run['cs_sub'] = run['cs_norm'] * (run['neutrons'] - NeutronsFromCarbon - NeutronsFromOxygen) / run[
            #     'beamParticles'] / run['thickness'] / neutron_eff * barn

            run['cs_sub'] = run['cs_norm'] * (run['neutrons'] - NeutronsFromCarbon - NeutronsFromOxygen)  / run['beamParticles'] / run['thickness'] / run[
                'charge'] / neutron_eff * barn

            # print(run['energy'], GetOxygenCS(run['energy']))
            # print('Nn {}, Ncarbon {}, Noxygen {}, Sub {} '.format(run['neutrons'], NeutronsFromCarbon, NeutronsFromOxygen, run['neutrons'] - NeutronsFromCarbon - NeutronsFromOxygen))

        run['cs'] = run['cs_norm']*run['neutrons']/run['beamParticles']/run['thickness']/run['charge']/neutron_eff*barn

        if blErrorsCS:
            run['cs_err'] = math.sqrt( 1/run['neutrons'] + 1/run['beamParticles'] ) + errTarget + errNeutronEff
        # run['cs_err_abs'] = run['cs_err'] *  run['cs']
        # print(run['runNumber'], ' ', run['energy'], ' cs ', run['cs'], run['cs_err'], ' Nn ', run['neutrons'], ' Beam ', run['beam'], 'Norm ', run['cs_norm'])

x_axis = []
y_axis = []

# plot_data = {'x':0, 'y':0}
plot_data = []
plot_data1 = []
plot_data2 = []
plot_data3 = []
plot_data4 = []
plot_data5 = []
plot_data6 = []

#Prepearing Plot
for run in js_data:
    if run['status'] != 'data' or run['source'] != 'alpha':
        continue

    if run['cs'] == 0:
        continue
    # print(run['status'])
    datapoint = {'x':0, 'y':0}
    # datapoint['x'] = run['energy']*1000
    # datapoint['x'] = (run['energy'] * 1.02227 - 63.208)*1000
    # datapoint['x'] = run['energy']
    datapoint['x'] = run['energy_calib']

    datapoint['y'] = run['cs']
    datapoint['y_sub'] = run['cs_sub']
    datapoint['err'] = run['cs']*run['cs_err']
    datapoint['y_minus'] = run['cs'] - datapoint['err']
    datapoint['y_plus'] = run['cs'] + datapoint['err']
    datapoint['label'] = run['runNumber']
    if datapoint['err'] < 0:
        datapoint['err']= abs(datapoint['err'])
        print('after neutron number correction negative value')

    if run['plotTag'] == 'series1':
        plot_data1.append(datapoint)
    elif run['plotTag'] == 'series2':
        plot_data2.append(datapoint)
    elif run['plotTag'] == 'N2-3':
        plot_data3.append(datapoint)
    elif run['plotTag'] == 'N2-4':
        plot_data4.append(datapoint)
    elif run['plotTag'] == 'series5':
        plot_data5.append(datapoint)
    elif run['plotTag'] == 'series6':
        plot_data6.append(datapoint)

print(plot_data.sort(key = SortX))
print(plot_data1.sort(key = SortX))
print(plot_data2.sort(key = SortX))
print(plot_data3.sort(key = SortX))
print(plot_data4.sort(key = SortX))
print(plot_data5.sort(key = SortX))
print(plot_data6.sort(key = SortX))

new_js_data = json.dumps(js_data, indent=3)

with open('{}'.format('result_2022_w39.json'),'w') as foutjson:
    foutjson.write(new_js_data)

# ListCSToAscii(plot_data, 'ascii/data_27Al.txt')
# ListCSToAscii(plot_data1, 'ascii/data1_27Al.txt')
# ListCSToAscii(plot_data2, 'ascii/data2_27Al.txt')
# ListCSToAscii(flynn_data2, 'ascii/flynn_data2.txt')
# ListCSToAscii(flynn_data, 'ascii/flynn_data.txt')


print('= Settings =')
print('blSubtractContaminants ', blSubtractContaminants)

if blSubtractContaminants:
    print('carbon_fraction: ', carbon_fraction)
    print('oxygen_fraction: ', oxygen_fraction)

print('blErrorsCS', blSubtractContaminants)
print('blCrossHair', blCrossHair)
print('blPlotHighEnergy', blPlotHighEnergy)
print('blPlotOxygenCS', blPlotOxygenCS)
print('blPlotCarbonCS', blPlotCarbonCS)



if blPlotCarbonCS:
    # carbon_plot = PrepeareJsonToPlot(js_carbon)
    PlotDataLiterature(PrepeareJsonToPlot(js_carbon, 1), 'whitesmoke', '--', False, '13C Harispolus', 12)

    # plt.scatter([i['x'] for i in carbon_plot], [i['y'] for i in carbon_plot], s=12, color='gray', label='13C Harispolus')
    # plt.plot([i['x'] for i in carbon_plot], [i['y'] for i in carbon_plot], color='gray', linestyle='--')

    # plt.scatter([i['x'] for i in plot_data3], [i['y'] for i in plot_data3], s=20, color='orange', label='data 2 ')
    # plt.plot([i['x'] for i in plot_data3], [i['y'] for i in plot_data3], color='orange')
    # plt.errorbar([i['x'] for i in plot_data3], [i['y'] for i in plot_data3], [i['err'] for i in plot_data3],
    #              color='orange')
    # plt.plot([i['x'] for i in plot_data3], [i['y_sub'] for i in plot_data3], color='orange', linestyle='--')
    # plt.fill_between([i['x'] for i in plot_data3], [i['y_minus'] for i in plot_data3],
    #                  [i['y_plus'] for i in plot_data3], color='orange', alpha=0.2)

blCS1 = True
if blCS1:
    PlotData(ax, plot_data1,'lime','--', blErrorsCS, 'N1-CaF-ser1',20, 3, False)
    PlotData(ax, plot_data2,'b','--', blErrorsCS, 'N2-CaF-ser2',20, 3, False)
    PlotData(ax, plot_data3,'steelblue','--', blErrorsCS, 'N2-3-CaF')
    PlotData(ax, plot_data4,'deepskyblue','--', blErrorsCS, 'N2-4-CaF')
    PlotData(ax, plot_data5,'cyan','--', blErrorsCS, 'N2-CaF-ser5',20, 3, False)
    PlotData(ax, plot_data6,'black','--', blErrorsCS, 'Cu:ser6:2022',20, 3, False)

########################## Literature data ##########################
    PlotDataLiterature(GetCSfromLit('{}/CS19F/Peters2016.json'.format(CSfolder),1 ), 'm','--',True, 'Peters2016',40, 3.5)
    PlotDataLiterature(GetCSfromLit('{}/CS19F/Wrean2000.json'.format(CSfolder), 1), 'hotpink','dotted',True, 'Wrean 2000',12)
    PlotDataLiterature(GetCSfromLit('{}/CS19F/Balakrishnan.json'.format(CSfolder), 1), 'y','--',False, 'Balakrishnan',12)
    PlotDataLiterature(GetCSfromLit('{}/CS19F/Peters2016B.json'.format(CSfolder), 1), 'g','None',True, 'Peters2016-B',12)

    if blPlotHighEnergy == True:
        PlotDataLiterature(GetCSfromLit('{}/CS19F/Ismail1993.json'.format(CSfolder), 1), 'peru','--',True, 'Ismail-1993',12)
        PlotDataLiterature(GetCSfromLit('{}/CS19F/Gladun.json'.format(CSfolder), 1), 'plum','--',True, 'Gladun',12)
        PlotDataLiterature(GetCSfromLit('{}/CS19F/Norman2015.json'.format(CSfolder), 1), 'coral','--',False, 'Norman',12)
########################## Literature data ##########################


    if blPlotOxygenCS:
        kunz_plot = PrepeareJsonToPlot(js_oxygen_kunz, 1)
        bair1973_plot = PrepeareJsonToPlot(js_oxygen_bair1973, 1)
        bair1962_plot = PrepeareJsonToPlot(js_oxygen_bair1962, 1)

        plt.scatter([i['x'] for i in kunz_plot], [i['y'] for i in kunz_plot], s=12, color='black', label='18O Hansen')
        plt.plot([i['x'] for i in kunz_plot], [i['y'] for i in kunz_plot], color='black', linestyle='--')

        plt.scatter([i['x'] for i in bair1973_plot], [i['y'] for i in bair1973_plot], s=12, color='black',  label='18O Bair 1973')
        plt.plot([i['x'] for i in bair1973_plot], [i['y'] for i in bair1973_plot], color='black', linestyle='--')

        plt.scatter([i['x'] for i in bair1962_plot], [i['y'] for i in bair1962_plot], s=12, color='black', label='18O Bair 1962')
        plt.plot([i['x'] for i in bair1962_plot], [i['y'] for i in bair1962_plot], color='black', linestyle='--')





    # AddLabels(williamson_data)
    # AddLabels(plot_data1)
    # print(plot_data2)


    # ax.legend(loc='upper left', fontsize = 16)
    # ax.legend(loc='lower right', fontsize = 18)
    # plt.title(r"$^{nat}$Cu($\alpha$, n)", fontsize = 18)
    plt.title(r"CS ($\alpha$, n)", fontsize=18)

    plt.xlabel('Incident $\u03B1$-particle Energy, MeV ', fontsize=15)
    plt.ylabel('cross-section, barn ', fontsize=15)
    plt.yscale('log')

    plt.yticks(fontsize=16)
    plt.xticks(fontsize=16)


    plt.xlim(2000, 8000)
    # plt.ylim(1e-5,1e-3)
    plt.ylim(0.001,0.9)

    #for series1
    # plt.xlim(4500, 6500)
    # plt.ylim(0.04,1.2)

    # plt.xlim(3600, 4500)
    subprocess.run(['python3', 'ring2energy.py'])
    shutil.copy('tmp_result_2022_w39.json', 'result_2022_w39.json')
    if (blCrossHair):
        blitted_cursor = BlittedCursor(ax)
        fig.canvas.mpl_connect('motion_notify_event', blitted_cursor.on_mouse_move)
        t = ax.transData
        matplotlib.backend_bases.MouseEvent("motion_notify_event", ax.figure.canvas, *t.transform((0.5, 0.5)))._process()
    plt.show()

    # matplotlib.pyplot.tra
