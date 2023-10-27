"""The order parameter for the hysteresis example."""
import numpy as np
from MDAnalysis.analysis import distances

# Assign global variables:
BOX = [12.4138, 12.4138, 12.4138, 90.0, 90.0, 90.0]
AT, O, RU = 193, 64, 2
NAMES = ["RU"] * RU + ["O"] * O + ["H"] * (AT - O - RU)


def finder(xyz):
    """Calculate something ."""
    at_ar = [0] * (RU + O)
    dic = {i: [] for i in range(RU + O)}
    bad = False

    # Generate atom_array and dict:
    at_ar, dic = finder_a(at_ar, dic, xyz, "s1")

    # Find the atom a_idx and excess electron x_idx:
    if 6 in at_ar[0:2]:
        # Excess electron exists in one of the Ru atoms
        a_idx = 0 if at_ar[0] == 6 else 1
        x_idx = np.argmax([i["dist"] for i in dic[a_idx]])
    elif at_ar[0] == at_ar[1] == 5:
        # Excess electron exists in one of the oxygen atoms.
        # Get list for oxygen with max count(x):
        o_l = [i for i, j in enumerate(at_ar) if j == np.max(at_ar)]
        a_idx0, x_idx0, dist = 0, 0, 0
        for a_idx0 in o_l:
            x_idx0 = np.argmax([i["dist"] for i in dic[a_idx0]])
            if dist < dic[a_idx0][x_idx0]["dist"]:
                dist = dic[a_idx0][x_idx0]["dist"]
                a_idx, x_idx = a_idx0, x_idx0
    else:
        if sum(at_ar) == 1:
            a_idx = at_ar.index(1)
            x_idx = np.argmax([i["dist"] for i in dic[a_idx]])
        else:
            print("It did not seem like we converged..")
            bad = True

    if not bad:
        loc = dic[a_idx][x_idx]["x_loc"]
        # Ru1-el, Ru2-el, Ru1-Ru2
        dist1 = distances.distance_array(xyz[0], xyz[loc], box=BOX)
        dist2 = distances.distance_array(xyz[1], xyz[loc], box=BOX)
        dist3 = distances.distance_array(xyz[0], xyz[1], box=BOX)
        OP = (dist1[0][0] - dist2[0][0]) / dist3[0][0]
        return OP, loc
    else:
        return None, None


def finder_a(atom_array, dic, xyz0, spin):
    """Calculate something ."""
    for i in range(len(xyz0[AT:][:, 0])):
        dist_arr = distances.distance_array(xyz0[AT + i], xyz0[: RU + O], box=BOX)[0]
        loc = np.argmin(dist_arr)
        atom_array[loc] += 1
        dic[loc].append(
            {"x_loc": i + AT, "dist": dist_arr[loc], "spin": spin, "a_loc": loc}
        )
    return atom_array, dic


def oh_finder(xyz):
    at_ar = [[] for i in range(RU + O)]
    for i in range(len(xyz[RU + O :])):
        dist_arr = distances.distance_array(xyz[RU + O + i], xyz[: RU + O], box=BOX)[0]
        at_ar[np.argmin(dist_arr)].append((RU + O + i, np.min(dist_arr)))
    o_idx = [len(i) for i in at_ar].index(1)
    h_idx, o_op = at_ar[o_idx][0]
    o_h_dist = distances.distance_array(xyz[o_idx], xyz[h_idx], box=BOX)[0]

    # find the surrounding ru h2o molecules
    h2o_indexes = []
    ru1_d = distances.distance_array(xyz[0], xyz[2:66], box=BOX)[0]
    ru2_d = distances.distance_array(xyz[1], xyz[2:66], box=BOX)[0]
    oo_idxes_1 = [
        j for j, bol in enumerate([0] * 2 + [i < 2.5 for i in ru1_d]) if bol == 1
    ]
    oo_idxes_2 = [
        j for j, bol in enumerate([0] * 2 + [i < 2.5 for i in ru2_d]) if bol == 1
    ]
    h2o_indexes += oo_idxes_1 + oo_idxes_2
    for i in oo_idxes_1 + oo_idxes_2:
        for h_idxes in at_ar[i]:
            h2o_indexes += [h_idxes[0]]
    h2o_indexes += [o_idx] + [h_tup[0] for h_tup in at_ar[o_idx]]

    return o_idx, o_op, at_ar, h2o_indexes, o_h_dist[0]
