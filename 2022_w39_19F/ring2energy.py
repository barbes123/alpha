import json
import sys
from os.path import exists
import matplotlib.pyplot as plt
import numpy as np

from os.path import abspath
sys.path.append(abspath('../lib/'))
from libCS import *


ringring = [12, 13, 23]
fileJson = 'result_2022_w39.json'

if not exists(fileJson):
    print('no fileJson')
    sys.exit()

def ConvertRatioToEnergy(js):
    R2R = np.loadtxt('../CS/ring/ring_sim.dat')
    for datum in js:
        i = 0
        # print(ringring[0])
        ener_list = []
        for el in datum['ring2ring']:
            ener = round(RingRing2Energy(R2R,ringring[i],el),2)
            ener_list.append(ener)
            # print (datum['runNumber'],' ', el, ' ', ener)
            i+=1
        datum['ring_ener'] = ener_list
    new_js = json.dumps(js, indent=3)
    return new_js

with open('{}'.format(fileJson),'r') as fjson:
    js_data =json.load(fjson)


js_data_new = ConvertRatioToEnergy(js_data)


with open('{}'.format('tmp_result_2022_w39.json'), 'w') as foutjson:
    foutjson.write(js_data_new)