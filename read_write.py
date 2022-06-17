"""remove Xs in a periodicly generated traj for jmol."""
import numpy as np
from pyretis.inout.formats.xyz import (read_xyz_file,
                                       convert_snapshot,
                                       write_xyz_trajectory)


def reader(inp, out):
    """remove Xs in a periodicly generated traj for jmol."""
    # a = [[], [], [], []]
    arr = []
    cnt = 0
    next0 = False
    with open(inp, encoding='utf-8') as file:
        with open(out, 'w', encoding='utf-8') as writer:
            for line in file:
                split = line.rstrip().split()
                next0 = next0 if len(split) != 1 else not next0
                cnt += 1 if len(split) == 1 and next0 else 0
                if len(split) != 1:
                    arr.append(line)
                if len(split) == 1 and next0 and cnt > 1:
                    arr_nox = [kaka for kaka in arr if 'X' not in kaka]
                    writer.write(f'{len(arr_nox)}\n')
                    writer.write(f'step {cnt-2}\n')
                    for write_l in arr_nox:
                        writer.write(f'{write_l}')
                    arr = []

def o_history(inp, order_path, out):
    """remove Xs in a periodicly generated traj for jmol."""
    order = np.loadtxt(order_path)
    o_idxes = set([i[3] for i in order])
    traj = list(read_xyz_file(inp))
    if len(traj) != len(order):
        print(f'wrong length between {traj}')
        print(f'and  {order_path}')
        exit('ape6')
    # elif len(o_idxes) == 1:
    #     print('The OH- does not change. no need to make history based trajectory.')
    #     # exit('ape5')
    #     return
        
    for idx, (order_idx, frame) in enumerate(zip(order, traj)):
        box, xyz, vel, names = convert_snapshot(frame)
        for o_idx in o_idxes:
            if names[int(o_idx)] == 'O':
                names[int(o_idx)] = 'B'
        write_xyz_trajectory(out, xyz, vel, names, box, idx)


# INPUT_ = './100739_moviePBC.xyz'
# OUTPUT = './100739_moviePBC_noX.xyz'
# reader(INPUT_, OUTPUT)
