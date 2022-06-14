"""Merge subpaths to one whole path using traj.txt."""
import os.path
from os import getcwd, mkdir
import numpy as np
import matplotlib.pyplot as plt
from pyretis.inout.formats.xyz import (read_xyz_file,
                                       convert_snapshot,
                                       write_xyz_trajectory)
from pyretis.pyvisa.common import read_traj_txt_file
from op_from_hmo import finder, oh_finder
from movieMaker import (createTitusMovie,
                        checkSaveFolders)
from read_write import reader


def combine_xyz(traj, name, out, plot=False):
    """Combine subpaths to make one whole path.xyz."""
    traj_dic = read_traj_txt_file(traj + '/traj.txt')
    traj_txt = np.loadtxt(traj + '/traj.txt', dtype='str')
    x_op_l, o_op_l = [], []
    if not os.path.exists(f'./0-dump/{name}'):
        for _, values in traj_dic.items():
            base = values[0][:-4]
            traj_l = list(read_xyz_file(f'./{traj}/traj/{values[0]}'))
            traj_s1 = list(read_xyz_file(f'./{traj}/traj/{base}-home1.xyz'))
            traj_s2 = list(read_xyz_file(f'./{traj}/traj/{base}-home2.xyz'))
            for j in [line for line in traj_txt if line[1] in values[0]]:
                box, xyz, vel, names = convert_snapshot(traj_l[int(j[2])])
                _, xyz_1, vel_1, names_1 = convert_snapshot(traj_s1[int(j[2])])
                _, xyz_2, vel_2, names_2 = convert_snapshot(traj_s2[int(j[2])])
                x = np.concatenate((xyz_1,
                                    xyz_2[names_2.index('X'):]))

                x_op, _, spin, loc, _ = finder(xyz_1, xyz_2)
                o_idx, _, o_op = oh_finder(xyz_1[:names_2.index('X')])
                names_1[o_idx] = 'N'
                x_op_l.append(x_op)
                o_op_l.append(o_op)

                if spin == 's1':
                    names_1[loc] = 'F'
                else:
                    names_2[loc] = 'F'

                names = names_1 + [nam for nam in names_2 if nam not in ('Ru', 'H', 'O')]

                vel = np.concatenate((vel, vel_1[names_1.index('X'):],
                                      vel_2[names_2.index('X'):]))
                write_xyz_trajectory(f'{out}/{name}', x,
                                     vel, names, box, j[0])
    with open(f'{out}/order_dz.txt', 'w', encoding='utf-8') as writer:
        for idx, (op_x, op_o) in enumerate(zip(x_op_l, o_op_l)):
            writer.write(f'{idx}\t{op_x:.4f}\t{op_o:.4f}\n')
    # WRITE o_op_l, x_op_l...!
    if plot:
        plotter(f'{out}/{name}', o_op_l, x_op_l, check=True)


def plotter(traj, o_op, x_op, check=False):
    """Plot the first ru atom position throughout the traj wrt origin."""
    frame, pos = [], []
    if check:
        for idx, snapshot in enumerate(read_xyz_file(traj)):
            _, xyz, _, _ = convert_snapshot(snapshot)
            frame.append(idx)
            pos.append(np.linalg.norm(xyz[0]))
        plt.plot(frame, [i/np.max(pos) -1 for i in pos])
    if o_op and x_op:
        plt.plot(frame, o_op, label='OH-oxygen')
        plt.plot(frame, x_op, label='OP')
    if check or (o_op and x_op):
        plt.show()


DATA_PATH = './1-data'
STRG_PATH = './2-movis'
ens_l = ['004']
for ens in ens_l:
    pathens = np.loadtxt(f'{DATA_PATH}/{ens}/pathensemble.txt', dtype='str')
    acc = [line for line in pathens if 'ACC' in line and 'cs' in line]
    # subfolder = [line[0] for line in acc]
    for idx, cycle in enumerate(acc):
        sub, A, B = cycle[0], cycle[3], cycle[5]
        traj_path = f'{DATA_PATH}/{ens}/traj/traj-acc/{str(sub)}/'
        folder_name = f'{STRG_PATH}/{ens}-{idx}-{A}{B}'
        path_name =  f'{ens}-{idx}-{A}{B}.xyz'
        exist = checkSaveFolders(folder_name)
        if exist:
            # combine_xyz(traj_path, path_name, folder_name, plot=True)
            createTitusMovie(path_name, folder_name, folder_name, 0, 2000)
            reader(folder_name +'/' + path_name, folder_name + f'/{path_name[:-4]}_noX.xyz')
            print('hi')                           
        exit("ape")
    # print(acc)
    # print(subfolder)

# for i in ensembles:
# pathensemble.txt reader..
# TRAJ_PATH = './1-data/007/traj/traj-acc/39/'
# combine_xyz(TRAJ_PATH, PATH plot=True)
# createTitusMovie('100739.xyz', './0-dump', './0-dump', 0, 2000)
# reader(inp, out) 
