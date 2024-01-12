import json
lut_dic = []
unit_constant = 1e-3
infile = 'hensen_170'
energy_kev = 1e3

if infile == 'kunz':
    unit_constant = 1e-6
elif infile == 'bair1973' or infile == 'bair1962' or infile == 'Harispolus' or infile == 'hansen':
    unit_constant = 1e-3
    energy_kev = 1000


with open('{}.dat'.format(infile), 'r') as fin:
    print(fin)
    for line in fin:
        if line[0] == '#':
            continue
        lut_line = line.split()
        if not lut_line:
            continue
        par_list = {'energy': 0, 'cs': 0.0, 'cs_err':0.0}  # r1/r2, r1/r3, r2/r3

        par_list['energy'] = float(lut_line[0])*energy_kev
        par_list['cs'] = float(lut_line[2])*unit_constant
        try:
            par_list['cs_err'] = float(lut_line[3])*unit_constant
        except:
            pass

        lut_dic.append(par_list)
# print(lut_dic)
lut_json = json.dumps(lut_dic, indent=3)

with open('{}.json'.format(infile), 'w') as fout:
    fout.write(lut_json)


