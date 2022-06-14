# from trajManipulation import *
# from manipulateOrderParameter import *
from os import getcwd, mkdir
from copy import *
from operator import itemgetter
import numpy as np
from numba import jit
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


def Create_jmolspt(L, cutoffs, saveFolder, fileName):
    f = open(saveFolder+"/"+fileName[:-4]+"_jmol.spt", "w")
    dir = getcwd()  # this give the working directory
    boxdimensions = str(L).replace("[", "").replace("]", "").replace(",", "")
    string = "load '"
    string += dir + "/"
    string += saveFolder+"/"+fileName[:-4]+"_moviePBC.xyz' {1 1 1} unitcell { "
    string += boxdimensions
    string += "  90.00  90.00  90.00}\n"
    f.write(string)
    string = "set PerspectiveDepth false\n"
    string += "select all; set measure angstroms\n"
    string += "cpk off\n"
    string += "set autobond false;\n"
    string += "boundbox on\n"
    f.write(string)
    string = ""
    for aa in cutoffs.keys():
        for bb in cutoffs.keys():
            string += "connect " + str(cutoffs[aa] + cutoffs[bb]) + " (_" + aa + ") (_" + bb + ")\n "
    f.write(string)
    string = "connect 1.5 2.5 (_O) (_H) hbond\n"
    string += "wireframe 0.2\n"
    string += "hbonds 0.01\n"
    string += "color hbonds yellow\n"
    string += "cpk 25%\n"
    f.write(string)
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


def MOLECLIST(atomlist, L, cutoffs):
    moleclist = []
    leftover = deepcopy(atomlist)
    while len(leftover) > 0:
        mol = []
        mol += [leftover[0]]
        del leftover[0]
        iat = 0
        while iat < len(mol):
            atom = mol[iat]
            neighbors, leftover = EXTRACTNEIGHBORSFROMLIST(atom, leftover, cutoffs, L)
            mol += neighbors
            iat += 1
        moleclist += [mol]
    return moleclist


def MOLECLIST2(atomlist, L, cutoffs):
    moleclist = []
    leftover = deepcopy(atomlist)
    xyz = np.array([i[0] for i in leftover])
    for idx, atom in enumerate(leftover):
        names = np.array([cutoffs[i[1]] for i in leftover[idx+1:]])
        vector = xyz[idx] - xyz[idx+1:]
        vector -= L * np.around(vector / L)
        d = np.sqrt(np.sum(vector * vector, axis=1))
        neighbors = d < cutoffs[leftover[0][1]] + names
        print(idx, sum(neighbors))
    print(xyz)

    # names = [i[1] for i in leftover][1:]
    # vector = xyz[0] - xyz[1:]
    # vector -= L * np.around(vector / L)
    # d = np.sqrt(np.sum(vector * vector, axis=1))
    # neighbors = [d[idx] < cutoffs[leftover[0][1]] + cutoffs[name] for idx, name in enumerate(names)]

    print(sum(neighbors))

    #print(d < cutoffs
    # print(xyz[0] - xyz[1:] - L*np.around()
    exit('booger')

    while len(leftover) > 0:
        mol = []
        mol += [leftover[0]]
        del leftover[0]
        iat = 0
        while iat < len(mol):
            atom = mol[iat]
            extract, leftover0, exx = [], [], []
            for sec in leftover:
                vector = atom[0] - sec[0]
                vector -= L * np.around(vector / L)
                d = np.sqrt(sum(vector * vector))
                exx.append(d < cutoffs[atom[1]]+ cutoffs[sec[1]])
            # exx = [DISTANCE(atom[0], sec[0], L) < cutoffs[atom[1]]+ cutoffs[sec[1]] for sec in leftover]
            for i, j in zip(exx, leftover):
                if i:
                    extract.append(j)
                else:
                    leftover0.append(j)
            leftover = leftover0
            neighbors = extract
            # neighbors, leftover = EXTRACTNEIGHBORSFROMLIST(atom, leftover, cutoffs, L)
            mol += neighbors
            iat += 1
        moleclist += [mol]
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
    start = time.time()
    reorderedmirror = sorted(mirrorlist, key=itemgetter(2))
    print(time.time() - start, 'part a')
    f.write(str(len(mirrorlist)) + "\n")
    f.write(str(framecount) + "\n")
    for at in reorderedmirror:
        f.write(f'{at[1]}\t{at[0][0]}\t{at[0][1]}\t{at[0][2]}\n')


def WRITEFRAME(N, f, g, framecount, L, strcoordinates, cutoffs, centerindex, Mirrors, commentline):
    shift = np.array([0., 0., 0.])

    start = time.time()
    for i in centerindex:
        x, y, z = strcoordinates[i - 1].split()[1:4]
        xyz = np.array([float(x), float(y), float(z)])
        shift += xyz

    # print(time.time() - start, 'first loop')
    shift /= max(1, len(centerindex))
    shift -= L / 2.
    atomlist = []

    start = time.time()
    for atindex in range(N):
        element, x, y, z = strcoordinates[atindex].split()[0:4]
        xyz = np.array([float(x), float(y), float(z)])
        xyz -= shift
        xyz -= L * np.floor(xyz / L)
        atom = [xyz, element, atindex + 1]
        atomlist += [atom]
    # print(time.time() - start, 'second loop')

    moleclist = MOLECLIST2(atomlist, L, cutoffs)
    WRITEMOLECLIST(g, moleclist, framecount, commentline)
    mirrorlist = []

    # start = time.time()
    for mol in moleclist:
        mol2 = MIRRORCOORDINATES(mol, L, N, Mirrors, cutoffs)
        mirrorlist += mol2
    # print(time.time() - start, 'third loop')

    # start = time.time()
    WRITEMIRROR2MOV(mirrorlist, framecount, f)
    # print(time.time() - start, 'fourth loop')
    # exit('ape')


def createTitusMovie(fileName, loadDirectory, saveDirectory, frameStart, frameEnd):

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

    ifile = loadDirectory+"/"+fileName
    f = open(ifile, "r")
    lines = f.readlines()
    N = int(lines[0])
    L = np.array([L, L, L])
    nframes = int(len(lines) / (N + 2))

    cutoffs = {"H": dH, "O": dO,"Ru": dRu, "X": dX, "N": dN, "F": dX}

    try:
        mkdir(saveDirectory)
        print(saveDirectory, "not found, creating directory:", saveDirectory)
    except:
        pass

    f = open(saveDirectory + "/" + fileName[:-4] + "_moviePBC.xyz", "w")
    g = open(saveDirectory + "/" + fileName[:-4] + "_moleclist.txt", "w")
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
    Create_jmolspt(L, cutoffs, saveDirectory, fileName)


def makeMovie(fileName, loadDirectory=".", saveDirectory=".", ensemble="007", frameStart=0, frameEnd=2000, markOH=False, tracking=True ,**movieTypes):
    originalFileName = fileName
    checkSaveFolders(saveDirectory)
    ordered = False

    print("Merging HOMO-files...")
    mergePaths(loadDirectory, "trajB", ensemble)
    mergePaths(loadDirectory, "trajF", ensemble)
    loadFrom = loadDirectory + "/" + fileName + ".xyz"
    print("Generating full trajectory...")
    combinePath(loadDirectory, loadFrom, atomAmount=714)

    for movies in movieTypes:
        if(movieTypes[movies] == "Full"):
            fileName = originalFileName+"_Merged.xyz"
            loadFrom = loadDirectory+"/"+fileName
            if (markOH):
                print("Replacing OH...")
                replaceOH(loadFrom, fileName, saveDirectory, loadDirectory, ordered, 193, 714, atomAmount=714)
                ordered = True
                fileName = fileName[:-4]+"_Replaced_OH.xyz"

            createTitusMovie(fileName, loadDirectory, saveDirectory, frameStart, frameEnd, markOH, tracking)


        if(movieTypes[movies] == "Electronless"): #Make sure to do full anyways for now
            fileName = originalFileName+"_Merged.xyz"
            loadFrom = loadDirectory+"/"+fileName
            fileName = fileName[:-4]+"_Electronless.xyz"
            saveAs = loadDirectory+"/"+fileName
            print("Stripping electrons...")
            stripElectrons(loadFrom, saveAs, cutoff=193, atomAmount=714)
            loadFrom = saveAs
            if (markOH):
                print("Replacing OH...")
                replaceOH(loadFrom, fileName, saveDirectory, loadDirectory, ordered, 193, 193, atomAmount=193)
                ordered = True
                fileName = fileName[:-4]+"_Replaced_OH.xyz"

            createTitusMovie(fileName, loadDirectory, saveDirectory, frameStart, frameEnd, markOH, tracking)

        if(movieTypes[movies] == "OrderParameter"):
            fileName = originalFileName + "_Merged.xyz"
            loadFrom = loadDirectory + "/" + fileName
            fileName = fileName[:-4] + "_IsolatedOP.xyz"
            saveAs = loadDirectory + "/" + fileName
            print("Isolating order parameter...")
            isolateOrderParameter(loadFrom, saveAs, atomAmount=193, particleAmount=714)
            loadFrom = saveAs
            if (markOH):
                print("Replacing OH...")
                replaceOH(loadFrom, fileName, saveDirectory, loadDirectory, ordered, bondlistStart=193, bondlistEnd=194, atomAmount=194)
                ordered = True
                fileName = fileName[:-4] + "_Replaced_OH.xyz"

            createTitusMovie(fileName, loadDirectory, saveDirectory, frameStart, frameEnd, markOH, tracking)
