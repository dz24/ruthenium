from os import getcwd, mkdir
from copy import *
from operator import itemgetter
import numpy as np
import time


def checkSaveFolders(saveDirectory):
    temp = ""
    exist = False
    for folders in saveDirectory.split("/"):
        try:
            mkdir(temp + folders)
            print(temp + folders, "not found, creating directory...")
            exist = False
        except:
            pass
            exist = True
        temp += folders + "/"
    return exist


def Create_jmolspt(L, cutoffs, saveFolder, fileName, order):
    OP_txt = np.loadtxt(order)
    h2o_idx = [0]+[1]+list(OP_txt[0][4:])
    o_idx_0, x_idx_0 = int(OP_txt[0][4]), int(OP_txt[0][5])
    o_idx_1 = set([int(i[4]) for i in OP_txt if int(i[4]) != o_idx_0])

    # Color the complex (that does not become OH:
    o_complex_idx = []
    if sum([i < 66 for i in OP_txt[0][7:]]) > 12 and OP_txt[0][-1] < 66:
        o_complex_idx = [OP_txt[0][-1]]

    f = open(fileName[:-4]+"_jmol.spt", "w")
    dir = getcwd()  # this give the working directory
    boxdimensions = str(L).replace("[", "").replace("]", "").replace(",", "")
    string = "load '"
    string += dir + "/"
    string += fileName[:-4]+"_moviePBC.xyz'\n"
    f.write(string)
    string = "set PerspectiveDepth false\n"
    string += "select all; set measure angstroms\n"
    string += "set autobond false;\n"
    string += "boundbox off\n"
    f.write(string)
    string = ""
    string = "connect 1.5 2.5 (_O) (_H) hbond\n"
    string += "hbonds 0.01\n"
    string += "color hbonds yellow\n"
    f.write(string)
    string = ""
    f.write(f"select not atomno=[{' '.join([str(int(i+1)) for i in h2o_idx])}]; hbonds off\n")
    string += "connect (_Ru) (_O) DELETE\n"
    string += "cpk 25%\n"
    string += f"select atomno={o_idx_0+1}; color blue\n"
    string += f"select atomno={x_idx_0+1}; color green\n"
    for o_idx in o_idx_1:
        string += f"select atomno={o_idx+1}; color peachpuff\n"
    if o_complex_idx:
        string += f"select atomno={o_complex_idx[0]+1}; color peachpuff\n"
    f.write(string)
    f.write(f"select not atomno=[{' '.join([str(int(i+1)) for i in h2o_idx])}]; wireframe only")
    f.close()



def DISTANCE(c1, c2, L=None):
    vector = c1 - c2
    if L is not None: vector -= L * np.around(vector / L)
    d = np.sqrt(sum(vector * vector))
    return d

def CONNECTED(at1, at2, cutoffs, L=False):
    c1, el1 = at1[0], at1[1]
    cutoff1 = cutoffs[el1]
    c2, el2 = at2[0], at2[1]
    cutoff2 = cutoffs[el2]
    d = DISTANCE(c1, c2, L)
    return d < cutoff1 + cutoff2

def EXTRACTNEIGHBORSFROMLIST(atom, leftover, cutoffs, L):
    indexleftover = 0
    extract = []
    while indexleftover < len(leftover):
        secatom = leftover[indexleftover]
        if CONNECTED(atom, secatom, cutoffs, L):
            extract += [secatom]
            del leftover[indexleftover]
        else:
            indexleftover += 1
    return extract, leftover


def MOLECLIST4(atomlist, L, cutoffs):
    moleclist = []
    tot_idx = []
    leftover = deepcopy(atomlist)
    xyz = np.array([i[0] for i in leftover])
    for idx, atom in enumerate(leftover):
        if idx not in tot_idx:
            names = np.array([cutoffs[i[1]] for i in leftover[idx+1:]])
            vector = xyz[idx] - xyz[idx+1:]
            vector -= L * np.around(vector / L)
            d = np.sqrt(np.sum(vector * vector, axis=1))
            neighbors = d < cutoffs[leftover[idx][1]] + names
            neighbors = (len(leftover) - len(neighbors))*[False] + list(neighbors)
            neighbors_idx = [idx] + [i for i in range(len(neighbors)) if neighbors[i] and i not in tot_idx]
            moleclist += [[leftover[i] for i in neighbors_idx]]
            tot_idx = list(set(neighbors_idx + list(tot_idx)))
    return moleclist


def MIRRORCOORDINATES(mol, L, N, Mirrors, cutoffs):
    firstat = mol[0]
    mirror = [firstat]
    del mol[0]
    imirror = 0
    listoftrans = []
    while len(mol) > 0:
        atom = mirror[imirror]
        neighbors, mol = EXTRACTNEIGHBORSFROMLIST(atom, mol, cutoffs, L)
        for ni in neighbors:
            cni = ni[0]
            cat = atom[0]
            vector = cni - cat
            trans = np.around(vector / L)
            cni -= trans * L
            ni[0] = cni
            mirror += [ni]
            if list(trans) != [0, 0, 0] and list(trans) not in [list(item) for item in
                                                                listoftrans]: listoftrans += [
                trans]
        imirror += 1
    if Mirrors:
        mol = deepcopy(mirror)
        for at in mol:
            for trans in listoftrans:
                newat = deepcopy(at)
                newat[0] += trans * L
                newat[2] += N
                mirror += [newat]
    return mirror


def WRITEMOLECLIST(g, moleclist, counter, commentline):
    g.write("********** counter=" + str(counter) + " ****************\n")
    g.write(commentline)
    for mol in moleclist:
        g.write("nat=" + str(len(mol)) + " ")
        for at in mol:
            g.write(at[1])
        for at in mol:
            g.write(" " + str(at[2]))
        g.write("\n")


def WRITEMIRROR2MOV(mirrorlist, framecount, f):
    reorderedmirror = sorted(mirrorlist, key=itemgetter(2))
    f.write(str(len(mirrorlist)) + "\n")
    f.write(str(framecount) + "\n")
    for at in reorderedmirror:
        f.write(f'{at[1]}\t{at[0][0]}\t{at[0][1]}\t{at[0][2]}\n')


def WRITEFRAME(N, f, g, framecount, L, strcoordinates, cutoffs, centerindex, Mirrors, commentline):
    shift = np.array([0., 0., 0.])

    for i in centerindex:
        x, y, z = strcoordinates[i - 1].split()[1:4]
        xyz = np.array([float(x), float(y), float(z)])
        shift += xyz

    shift /= max(1, len(centerindex))
    shift -= L / 2.
    atomlist = []

    for atindex in range(N):
        element, x, y, z = strcoordinates[atindex].split()[0:4]
        xyz = np.array([float(x), float(y), float(z)])
        xyz -= shift
        xyz -= L * np.floor(xyz / L)
        atom = [xyz, element, atindex + 1]
        atomlist += [atom]

    moleclist = MOLECLIST4(atomlist, L, cutoffs)
    WRITEMOLECLIST(g, moleclist, framecount, commentline)
    mirrorlist = []

    for mol in moleclist:
        mol2 = MIRRORCOORDINATES(mol, L, N, Mirrors, cutoffs)
        mirrorlist += mol2

    WRITEMIRROR2MOV(mirrorlist, framecount, f)


def createTitusMovie(fileName, loadDirectory, saveDirectory, frameStart, frameEnd, order):

    L = 12.4138
    selectframes = range(frameStart, frameEnd)
    centerindex = [1,2]
    Mirrors = True

    #Bondlengths
    dH = 0.91 / 2.
    dO = 1.81 / 2.
    dN = 1.81 / 2.
    dRu = 0
    dX = 0. / 2

    # ifile = loadDirectory+"/"+fileName
    f = open(fileName, "r")
    lines = f.readlines()
    N = int(lines[0])
    L = np.array([L, L, L])
    nframes = int(len(lines) / (N + 2))

    # cutoffs = {"H": dH, "O": dO,"Ru": dRu, "X": dX, "N": dN, "F": dX, "B": dO, "He": dH}
    cutoffs = {"H": dH, "O": dO,"Ru": dRu, "X": dX}

    try:
        mkdir(saveDirectory)
        print(saveDirectory, "not found, creating directory:", saveDirectory)
    except:
        pass

    f = open(fileName[:-4] + "_moviePBC.xyz", "w")
    g = open(fileName[:-4] + "_moleclist.txt", "w")
    lcount = 0
    framecount = 0
    for l in range(nframes):
        loadBar = "["
        if l in selectframes:
            framecoordinates = lines[lcount + 2:lcount + 2 + N]
            commentline = lines[lcount + 1]
            WRITEFRAME(N, f, g, framecount, L, framecoordinates, cutoffs, centerindex, Mirrors, commentline)
        for j in range(1, 31):
            if (j <= (framecount / nframes) * 30):
                loadBar = loadBar + "▋"
            else:
                loadBar = loadBar + " "
        loadBar = loadBar + "]"
        print("\rGenerating movie from '"+fileName+"':" + loadBar + "  " + str(
            round((framecount / nframes) * 100, 1)) + "% complete", end="")
        framecount += 1
        lcount = lcount + 2 + N
    f.close()
    g.close()

    print("")
    print("Molecular saved in " + saveDirectory + " as '"+fileName+"_moleclist.txt'")
    print("Moviefile saved to: " + saveDirectory + " as '"+fileName+"_moviePBC.xyz'")
    print("Trajectory files saved in:", loadDirectory)

    # new code: create jmol.spt from scratch
    Create_jmolspt(L, cutoffs, saveDirectory, fileName, order)
