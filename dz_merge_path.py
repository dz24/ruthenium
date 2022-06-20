"""Merge subpaths to one whole path using traj.txt."""
import os.path
import numpy as np
import matplotlib.pyplot as plt
from pyretis.inout.formats.xyz import (read_xyz_file,
                                       convert_snapshot,
                                       write_xyz_trajectory)
from pyretis.pyvisa.common import read_traj_txt_file
from op_from_hmo import finder, oh_finder
from movieMaker import (createTitusMovie,
                        checkSaveFolders)
from read_write import o_history


def combine_xyz(traj, name, out, plot=False):
    """Combine subpaths to make one whole path.xyz."""
    traj_dic = read_traj_txt_file(traj + '/traj.txt')
    traj_txt = np.loadtxt(traj + '/traj.txt', dtype='str')
    x_op_l, o_op_l, o_idx_l = [], [], []
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
            o_idx, _, o_op, _ = oh_finder(xyz_1[:names_2.index('X')])
            names_1[o_idx] = 'N'
            x_op_l.append(x_op)
            o_op_l.append(o_op)
            o_idx_l.append(o_idx)

            if spin == 's1':
                names_1[loc] = 'F'
            else:
                names_2[loc] = 'F'
            names = names_1 + [nam for nam in names_2 if nam not in ('Ru', 'H', 'O')]

            vel = np.concatenate((vel, vel_1[names_1.index('X'):],
                                  vel_2[names_2.index('X'):]))
            write_xyz_trajectory(f'{out}/{name}', x,
                                 vel, names, box, step=j[0])
    with open(f'{out}/order_{name[:-4]}.txt', 'w', encoding='utf-8') as writer:
        writer.write('# idx\top_x\top_o\to_idx\n')
        for idx, (op_x, op_o) in enumerate(zip(x_op_l, o_op_l)):
            writer.write(f'{idx}\t\t{op_x:.4f}\t{op_o:.4f}\t{o_idx_l[idx]:.0f}\n')


def plotter(traj, o_op, x_op, check=False):
    """Plot the first ru atom position throughout the traj wrt origin."""
    frame, pos = [], []
    if check:
        for idx, snapshot in enumerate(read_xyz_file(traj)):
            _, xyz, _, _ = convert_snapshot(snapshot)
            frame.append(idx)
            pos.append(np.linalg.norm(xyz[0]))
        plt.plot(frame, [i/np.max(pos) - 1 for i in pos])
    if o_op and x_op:
        plt.plot(frame, o_op, label='OH-oxygen')
        plt.plot(frame, x_op, label='OP')
    if check or (o_op and x_op):
        plt.show()


def make_movies(data_path, strg_path, ens_l):
    """Use all the diff functions to create movies."""
    for ens in ens_l:
        pathens = np.loadtxt(f'{data_path}/{ens}/pathensemble.txt', dtype='str')
        acc = [line for line in pathens if 'ACC' in line and 'cs' in line]
        for idx, cycle in enumerate(acc):
            print(f'ens:{ens}, cycle: {idx}')
            sub, A, B, length = cycle[0], cycle[3], cycle[5], cycle[6]
            traj_path = f'{data_path}/{ens}/traj/traj-acc/{str(sub)}/'
            folder_name = f'{strg_path}/{ens}-{idx}-{A}{B}-{length}'
            path_name = f'{ens}-{idx}-{A}{B}-{length}.xyz'
            out_path = os.path.join(folder_name, path_name)
            out_path_O = path_name[:-4] + '-Ohist.xyz'

            print('traj_path', traj_path)
            print('path_name', path_name)
            print('folder_name', folder_name)
            make_movie(traj_path, path_name, folder_name)
            exit('ape')

            if not os.path.exists(out_path):
                print('combining xyz.')
                combine_xyz(traj_path, path_name, folder_name, plot=False)
            if not os.path.exists(out_path[:-4] + '-Ohist.xyz'):
                print('making o-history xyz.')
                order_file = os.path.join(folder_name, 'order_' + path_name[:-4] + '.txt')
                o_history(out_path, order_file, out_path[:-4] + '-Ohist.xyz')
            if not os.path.exists(out_path[:-4] + '-Ohist_moviePBC.xyz'):
                print('making titus movie')
                createTitusMovie(out_path_O, folder_name, folder_name, 0, 2000)


def make_movie(in_path, out_path_name, out_folder, extra_o_idx_l=[]):
    """Make single movie."""
    out_path = os.path.join(out_folder, out_path_name)
    checkSaveFolders(out_folder)
    print('bling')
    if not os.path.exists(out_path):
        print('combining xyz.')
        combine_xyz(in_path, out_path_name, out_folder)
    if not os.path.exists(out_path[:-4] + '-Ohist.xyz'):
        print('making o-history xyz.')
        order_file = os.path.join(out_folder,
                                  'order_' + out_path_name[:-4] + '.txt')
        o_history(out_path, order_file, out_path[:-4] + '-Ohist.xyz',
                  extra_o_idx_l)
    if not os.path.exists(out_path[:-4] + '-Ohist_moviePBC.xyz'):
        print('making titus movie')
        out_path_O = out_path_name[:-4] + '-Ohist.xyz'
        createTitusMovie(out_path_O, out_folder, out_folder, 0, 2000)


DATA_PATH0 = './1-data'
STRG_PATH0 = './2-movis'
ENS_L0 = ['004']
# make_movies(DATA_PATH0, STRG_PATH0, ENS_L0)
make_movie(None, '006-15-LR-264.xyz', '2-movis/006-15-LR-264', [40])
