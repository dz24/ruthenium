"""remove Xs in a periodicly generated traj for jmol."""
import numpy as np
from pyretis.inout.formats.xyz import (read_xyz_file,
                                       convert_snapshot,
                                       write_xyz_trajectory)
from op_from_hmo import oh_finder

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
    """rearrange..."""
    order = np.loadtxt(order_path)
    o_idx_0, x_idx_0 = int(order[0][4]), int(order[0][5])
    traj = list(read_xyz_file(inp))
    if len(traj) != len(order):
        print(f'wrong length between {traj}')
        print(f'and  {order_path}')
        exit('ape6')
    
    print('pirate', inp, out)
    for idx, (order_idx, frame) in enumerate(zip(order, traj)):
        box, xyz, vel, names = convert_snapshot(frame)
        # o_idx, _, _, at_ar, _ = oh_finder(xyz[:names.index('X')])
        o_idx, x_idx, h_idx = order[idx][4:7]
        
        # Switch xyz between active OH, X for jmol color.
        if int(o_idx) != o_idx_0:
            xyz_temp = xyz[int(o_idx)].copy()
            xyz[int(o_idx)] = xyz[o_idx_0]
            xyz[o_idx_0] = xyz_temp
        if int(x_idx) != x_idx_0:
            xyz_temp = xyz[int(x_idx)].copy()
            xyz[int(x_idx)] = xyz[x_idx_0]
            xyz[x_idx_0] = xyz_temp

        write_xyz_trajectory(out, xyz, vel, names, box, idx)

# INPUT_ = './100739_moviePBC.xyz'
# OUTPUT = './100739_moviePBC_noX.xyz'

# reader(INPUT_, OUTPUT)
