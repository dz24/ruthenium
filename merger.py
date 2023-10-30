import os
import numpy as np
import MDAnalysis as mda
import shutil
from common import read_adresses, read_trajtxt, read_homos
from op_from_hmo import finder

inp = '../proc/adress.txt'
data = read_adresses(inp)
store = '../proc/trajs/'

for idx, (key, item) in enumerate(data.items()):
    # if item[:3] in ('000', '001', '002', '003', '004', '005', '006'):
    #     continue
    orig = f'../data/{item}'
    dest = store + key
    if not os.path.exists(orig):
        print(f'WARNING: {orig} does not exist! Ignoring this loop.')
        continue

    if not os.path.isdir(dest):
        os.makedirs(dest)
    if not os.path.isdir(dest + '/accepted'):
        os.makedirs(dest + '/accepted')

    ordertxt = np.loadtxt(orig +'/order.txt')
    shutil.copyfile(orig +'/order.txt', dest + '/order.txt')
    trajtxt, adressset = read_trajtxt(orig +'/traj.txt')
    homos = read_homos(trajtxt, orig)

    ## here we merge the xyz files.
    # load them up:
    mdas = {}
    natoms = []
    for adress in adressset:
        adress_file = f'{orig}/traj/{adress}'
        mdas[adress] = mda.Universe(adress_file)
        natoms.append(mdas[adress].atoms.n_atoms)
    # if homos we save those as well
    mda_homos1 = {}
    mda_homos2 = {}
    for homo in homos:
        adress_file = f'{orig}/traj/{homo}'
        if 'home1' in adress_file:
            mda_homos1[homo] = mda.Universe(adress_file)
        else:
            mda_homos2[homo] = mda.Universe(adress_file)

    # trajfile = dest + '/accepted/merged.xyz'
    homofile = dest + '/accepted/merged-homos.xyz'

    keys = list(trajtxt.keys())
    # with mda.Writer(trajfile, natoms[0]) as write:
    #     for i in range(len(trajtxt[keys[0]])):
    #         step, filename, index, vel = (trajtxt[j][i] for j in keys)
    #         pos = mdas[filename].trajectory[int(index)]
    #         write.write(mdas[filename].atoms)

    # write homos
    print('howdy', keys[0])
    # with mda.Writer(homofile, 716) as write:
    with mda.Writer(homofile, 194) as write:
        for i in range(len(trajtxt[keys[0]])):
            step, filename, index, vel = (trajtxt[j][i] for j in keys)
            if 'second' not in filename:
                hfilename1 = filename[:-4] + '-home1.xyz'
                hfilename2 = filename[:-4] + '-home2.xyz'
                if hfilename1 not in mda_homos1 or hfilename2 not in mda_homos2:
                    continue
                pos1 = mda_homos1[hfilename1].trajectory[int(index)]
                pos2 = mda_homos2[hfilename2].trajectory[int(index)]
                merg = mda.Merge(mda_homos1[hfilename1].atoms,
                                 mda_homos2[hfilename2].select_atoms('name X'))
                whada = finder(merg.atoms.positions)
                merg2 = mda.Merge(mda_homos1[hfilename1].atoms[:193],
                                 merg.select_atoms(f'bynum {whada[3]+1}'))
                whada2 = finder(merg2.atoms.positions)
                write.write(merg2.atoms)
            else:
                pos = mdas[filename].trajectory[int(index)]
                copy = mdas[filename].copy()
                copy.atoms[0].name = 'X'
                copy.atoms[0].element= 'X'
                agroup = mda.AtomGroup([copy.atoms[0]]*1)
                merg = mda.Merge(mdas[filename].atoms, agroup)
                write.write(merg2.atoms)

    # mda.Universe(

    # else:
    #     print(f'WARNING: {dest} already exists! Ignoring this loop.')
    #     continue

    # print(key, item, orig, os.path.exists(orig))
    # os.path.exists(orig)

    print('iteration:', idx)
    # if idx == 3:
    #     break
