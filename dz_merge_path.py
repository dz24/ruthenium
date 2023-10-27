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
    # traj_dic = read_traj_txt_file(traj + '/traj.txt')
    # traj_txt = np.loadtxt(traj + '/traj.txt', dtype='str')
    x_op_l, o_op_l, o_idx_l, h_idx_l, x_idx_l = [], [], [], [], []
    o_h_dist_l = []
    o_h_l = [[] for i in range(66)]
    h2o_idxes = []
    print('whdafa')
    traj_n = traj + '/accepted/' + name
    print('a1', traj_n)
    traj_l = list(read_xyz_file(traj_n))
    # print('tiger', traj_l)
    for frame in traj_l:
        box, xyz, vel, names = convert_snapshot(frame)
        o_idx, _, o_op, at_ar, h2o_idxes_0, o_h_dist = oh_finder(xyz[:names.index('X')])
        x_op, _, spin, x_idx, _ = finder(xyz)
        x_op_l.append(x_op)
        o_op_l.append(o_op)
        o_idx_l.append(o_idx)
        h_idx_l.append(at_ar[int(o_idx)][0][0])
        x_idx_l.append(x_idx)
        o_h_dist_l.append(o_h_dist)
        h2o_idxes_0 += [h_tup[0] for h_tup in at_ar[o_idx]]
        h2o_idxes = list(set(h2o_idxes_0 + h2o_idxes))

    if 3 in [len(i) for i in o_h_l]:
        o_h_l_index = [len(i) for i in o_h_l].index(3)
        h2o_idxes = list(set(h2o_idxes + o_h_l[o_h_l_index]))
        if  o_h_l_index not in h2o_idxes:
            h2o_idxes.append(o_h_l_index)

    with open(f'{out}/order_{name[:-4]}.txt', 'w', encoding='utf-8') as writer:

        writer.write('# idx\top_x\top_o\top_oh\to_idx\tx_idx\th_idx\twater_idx\n')
        for idx, (op_x, op_o) in enumerate(zip(x_op_l, o_op_l)):
            writer.write(f'{idx}\t\t{op_x:.4f}\t{op_o:.4f}\t{o_h_dist_l[idx]:.4f}\t')
            writer.write(f'{o_idx_l[idx]:.0f}\t\t{x_idx_l[idx]:.0f}\t\t')
            writer.write(f'{h_idx_l[idx]:.0f}\t\t')
            for h2o_idx in h2o_idxes:
                writer.write(f'{h2o_idx:.0f}\t')
            writer.write('\n')


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

def make_movie(in_path, out_path_name, out_folder, extra_o_idx_l=[]):
    """Make single movie."""
    out_path = os.path.join(out_folder, out_path_name)
    checkSaveFolders(out_folder)
    order_file = os.path.join(out_folder,
                                'order_' + out_path_name[:-4] + '.txt')
    in_file = in_path + '/accepted/' + out_path_name
    print('making orderp')
    if not os.path.exists(out_path_name):
        combine_xyz(in_path, out_path_name, out_folder)

    if not os.path.exists(out_path[:-4] + '-Ohist.xyz'):
        print('making o-history xyz.')
        print('pumpkin 1', in_path)
        print('pumpkin 2', out_path)
        o_history(in_file, order_file, out_path[:-5] + '-Ohist.xyz')
    if not os.path.exists(out_path[:-4] + '-Ohist_moviePBC.xyz'):
        print('making titus movie')
        out_path_O = out_path_name[:-5] + '-Ohist.xyz'
        createTitusMovie(out_path_O, out_folder, out_folder, 0, 2000, order_file)


# DATA_PATH0 = './1-data/3-tis_sh_cs/2-cs/'
# DATA_PATH0 = './1-data/3-tis_sh_cs/2-cs/'
# STRG_PATH0 = './2-movis/3-tis_sh_cs/'
# ENS_L0 = ['002']
# make_movies(DATA_PATH0, STRG_PATH0, ENS_L0)
# make_movie('./1-data/006/traj/traj-acc/19', '006-9-LL-301.xyz', './0-dump/0-dump/')
# make_movie('/home/danielzh/dump/plot/4/005/65', '006-9-LL-301.xyz', './0-dump/0-dump/')
# make_movie('./1-data/4-sim_4/007/traj/traj-acc/215/', '007.xyz', './0-dump/0-dump/')
make_movie('dump', 'merged-homos.xyz', './dump_out/')
# make_movie('./1-data/006/traj/traj-acc/66', '006-15-LR-264.xyz', './2-movis/1-sim_1/006-15-LR-264/')
