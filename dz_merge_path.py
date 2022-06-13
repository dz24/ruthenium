"""Merge subpaths to one whole path using traj.txt."""
import os.path
import numpy as np
import matplotlib.pyplot as plt
from pyretis.inout.formats.xyz import (read_xyz_file,
                                       convert_snapshot,
                                       write_xyz_trajectory)
from pyretis.pyvisa.common import read_traj_txt_file
from op_from_hmo import finder, oh_finder
from movieMaker import createTitusMovie


def combine_xyz(traj, plot=False):
    """Combine subpaths to make one whole path.xyz."""
    traj_dic = read_traj_txt_file(traj + '/traj.txt')
    traj_txt = np.loadtxt(traj + '/traj.txt', dtype='str')
    name = ''.join(i for i in traj if i in '0123456789') + '.xyz'
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
                write_xyz_trajectory(f'./0-dump/{name}', x,
                                     vel, names, box, j[0])
    if plot:
        plotter(f'./0-dump/{name}', o_op_l, x_op_l)


def plotter(traj, o_op, x_op, check=False):
    """Plot the first ru atom position throughout the traj wrt origin."""
    frame, pos = [], []
    if check:
        for idx, snapshot in enumerate(read_xyz_file(traj)):
            _, xyz, _, _ = convert_snapshot(snapshot)
            frame.append(idx)
            pos.append(np.linalg.norm(xyz[0]))
        plt.plot(frame, [i/np.max(pos) for i in pos])
    if o_op and x_op:
        plt.plot(frame, o_op, label='OH-oxygen')
        plt.plot(frame, x_op, label='OP')
    if check or (o_op and x_op):
        plt.show()


# TRAJ_PATH = './1-data/007/traj/traj-acc/39/'
# combine_xyz(TRAJ_PATH, plot=True)
createTitusMovie('100739.xyz', './0-dump', './0-dump', 0, 2000)
