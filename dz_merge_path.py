"""Merge subpaths to one whole path using traj.txt."""
import os.path
import numpy as np
from pyretis.inout.formats.xyz import (
    read_xyz_file,
    convert_snapshot,
    write_xyz_trajectory,
)
from op_from_hmo import finder, oh_finder
from movieMaker import createTitusMovie, checkSaveFolders


def o_history(inp, order_path, out):
    """rearrange..."""
    order = np.loadtxt(order_path)
    o_idx_0, x_idx_0 = int(order[0][4]), int(order[0][5])
    traj = list(read_xyz_file(inp))
    if len(traj) != len(order):
        print(f'wrong length between {traj}')
        print(f'and  {order_path}')
        exit('ape6')
    
    for idx, (order_idx, frame) in enumerate(zip(order, traj)):
        box, xyz, vel, names = convert_snapshot(frame)
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

def orderp_mov(traj_n, out, plot=False):
    """Combine subpaths to make one whole path.xyz."""
    x_op_l, o_op_l, o_idx_l, h_idx_l, x_idx_l = [], [], [], [], []
    o_h_dist_l = []
    o_h_l = [[] for i in range(66)]
    h2o_idxes = []
    traj_l = list(read_xyz_file(traj_n))
    for frame in traj_l:
        box, xyz, vel, names = convert_snapshot(frame)
        o_idx, o_op, at_ar, h2o_idxes_0, o_h_dist = oh_finder(xyz[:names.index("X")])
        x_op, x_idx,= finder(xyz)
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
        if o_h_l_index not in h2o_idxes:
            h2o_idxes.append(o_h_l_index)

    with open(f"{out}/order.txt", "w", encoding="utf-8") as writer:
        writer.write("# idx\top_x\top_o\top_oh\to_idx\tx_idx\th_idx\twater_idx\n")
        for idx, (op_x, op_o) in enumerate(zip(x_op_l, o_op_l)):
            writer.write(f"{idx}\t\t{op_x:.4f}\t{op_o:.4f}\t{o_h_dist_l[idx]:.4f}\t")
            writer.write(f"{o_idx_l[idx]:.0f}\t\t{x_idx_l[idx]:.0f}\t\t")
            writer.write(f"{h_idx_l[idx]:.0f}\t\t")
            for h2o_idx in h2o_idxes:
                writer.write(f"{h2o_idx:.0f}\t")
            writer.write("\n")


def make_movie(infile, out_folder, extra_o_idx_l=[]):
    """Make single movie."""
    outfile = os.path.join(out_folder, "merged-homos.xyz")
    checkSaveFolders(out_folder)
    order_file = os.path.join(out_folder, "order.txt")

    if not os.path.exists(order_file):
        print("making order.txt.")
        orderp_mov(infile, out_folder)

    if not os.path.exists(outfile[:-4] + "-Ohist.xyz"):
        print("making o-history xyz.")
        o_history(infile, order_file, outfile[:-5] + "-Ohist.xyz")

    if not os.path.exists(outfile[:-4] + "-Ohist_moviePBC.xyz"):
        print("making titus movie")
        outfile_O = outfile[:-5] + "-Ohist.xyz"
        createTitusMovie(outfile_O, out_folder, out_folder, 0, 2000, order_file)
    
    # delete old


make_movie("dump/accepted/merged-homos.xyz", "./dump_out/")
