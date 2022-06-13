"""The order parameter for the hysteresis example."""
import numpy as np
from MDAnalysis.analysis import distances

# Assign global variables:
BOX = [12.4138, 12.4138, 12.4138, 90.0, 90.0, 90.0]
AT, O, RU = 193, 64, 2
NAMES = ['RU']*RU + ['O']*O + ['H']*(AT-O-RU)


def finder(xyz1, xyz2):
    """Calculate something ."""
    at_ar = [0] * (RU+O)
    dic = {i: [] for i in range(RU+O)}
    bad = False

    # Generate atom_array and dict:
    at_ar, dic = finder_a(at_ar, dic, xyz1, 's1')
    at_ar, dic = finder_a(at_ar, dic, xyz2, 's2')

    # Find the atom a_idx and excess electron x_idx:
    if 6 in at_ar[0:2]:
        # Excess electron exists in one of the Ru atoms
        a_idx = 0 if at_ar[0] == 6 else 1
        x_idx = np.argmax([i['dist'] for i in dic[a_idx]])
    elif (at_ar[0] == at_ar[1] == 5):
        # Excess electron exists in one of the oxygen atoms.
        # Get list for oxygen with max count(x):
        o_l = [i for i, j in enumerate(at_ar) if j == np.max(at_ar)]
        a_idx0, x_idx0, dist = 0, 0, 0
        for a_idx0 in o_l:
            x_idx0 = np.argmax([i['dist'] for i in dic[a_idx0]])
            if dist < dic[a_idx0][x_idx0]['dist']:
                dist = dic[a_idx0][x_idx0]['dist']
                a_idx, x_idx = a_idx0, x_idx0
    else:
        print('It did not seem like we converged..')
        bad = True

    if not bad:
        loc = dic[a_idx][x_idx]['x_loc']
        spin = dic[a_idx][x_idx]['spin']
        traj = xyz1 if spin == 's1' else xyz2

        # Ru1-el, Ru2-el, Ru1-Ru2
        dist1 = distances.distance_array(traj[0], traj[loc], box=BOX)
        dist2 = distances.distance_array(traj[1], traj[loc], box=BOX)
        dist3 = distances.distance_array(traj[0], traj[1], box=BOX)
        OP = (dist1[0][0]-dist2[0][0])/dist3[0][0]
        return OP, dist1[0][0], spin, loc, bad
    else:
        return None, None, None, None, bad, None


def finder_a(atom_array, dic, xyz0, spin):
    """Calculate something ."""
    for i in range(len(xyz0[AT:][:, 0])):
        dist_arr = distances.distance_array(xyz0[AT+i],
                                            xyz0[:RU+O],
                                            box=BOX)[0]
        loc = np.argmin(dist_arr)
        atom_array[loc] += 1
        dic[loc].append({'x_loc': i+AT,
                         'dist': dist_arr[loc],
                         'spin': spin,
                         'a_loc': loc})
    return atom_array, dic


def oh_finder(xyz):
    at_ar = [[] for i in range(RU+O)]
    for i in range(len(xyz[RU+O:])):
        dist_arr = distances.distance_array(xyz[RU+O + i], xyz[:RU+O], box=BOX)[0]
        at_ar[np.argmin(dist_arr)].append((RU+O + i, np.min(dist_arr)))
    o_idx = [len(i) for i in at_ar].index(1)
    h_idx, o_OP = at_ar[o_idx][0]
    # ru_o = distances.distance_array(xyz[0], xyz[o_idx], box=BOX)[0]
    dist1 = distances.distance_array(xyz[0], xyz[o_idx], box=BOX)
    dist2 = distances.distance_array(xyz[1], xyz[o_idx], box=BOX)
    dist3 = distances.distance_array(xyz[0], xyz[1], box=BOX)
    ru_o = (dist1[0][0]-dist2[0][0])/dist3[0][0]
    return o_idx, o_OP, ru_o



