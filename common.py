import os

def pathens_reader(inp, dic, ensnum):
    print('birdie', ensnum)
    with open(inp, 'r') as read:
        for idx, line in enumerate(read):
            if 'ACC' not in line or 'ld' in line or 're' in line:
                continue
            rip = line.rstrip().split()
            # paths are identified by opmax-opmin-len
            plen, opmax, opmin = rip[6], rip[9], rip[10]
            lmr = rip[3] + rip[4] + rip[5]
            weight = float(rip[-1])
            key = f'{lmr}-{opmax}-{opmin}-{lmr}'
            adress = f'{ensnum}/{rip[0]}'
            if key not in dic:
                dic[key] = {'len': plen, 'opmax': opmax, 'opmin': opmin,
                            f'{ensnum}_w': weight, 'adress': [adress],
                            f'{ensnum}_f': 1}
            else:
                dic[key][f'{ensnum}_w'] = weight
                if f'{ensnum}_f' in dic[key]:
                    dic[key][f'{ensnum}_f'] += 1
                else:
                    dic[key][f'{ensnum}_f'] = 1
                dic[key]['adress'].append(adress)


def write_pathens(path_dic, ensnums):
    with open('../proc/cpens_data.txt', 'w') as write:
        write.write('\t#\txxx\tlen\tmax OP\t' + '\t'.join(ensnums) + '\n')
        dkeys = ['len', 'opmax', 'opmin'] + ensnums
        for idx, (key, item) in enumerate(path_dic.items()):
            plen, opmin, opmax = item['len'], item['opmax'], item['opmin']
            toprint = [f"\t{idx}\t{plen}\t{opmax}\t{opmin}\t\t"]
            for dkey in ensnums:
                freq = f'{dkey}_f'
                if freq in item:
                    toprint.append(str(item[freq]))
                else:
                    toprint.append('----')
            for dkey in ensnums:
                weig= f'{dkey}_w'
                if weig in item:
                    toprint.append(str(item[weig]))
                else:
                    toprint.append('----')
            print(idx, toprint)
            write.write('\t'.join(toprint) + '\n')

def write_adresses(path_dic, ensnums):
    with open('../proc/adress.txt', 'w') as write:
        write.write('\t#\txxx\tlen\tmax OP\t' + '\t'.join(ensnums) + '\n')
        dkeys = ['len', 'opmax', 'opmin'] + ensnums
        for idx, (key, item) in enumerate(path_dic.items()):
            plen, opmin, opmax = item['len'], item['opmax'], item['opmin']
            toprint = [f"\t{idx}\t{plen}\t{opmax}\t{opmin}\t\t"]
            toprint += item['adress']
            write.write('\t'.join(toprint) + '\n')

def read_adresses(inp):
    dic = {}
    with open(inp, 'r') as read:
        for line in read:
            if '#' in line:
                continue
            rip = line.rstrip().split()
            dic[rip[0]] = rip[-1]
    return dic 

def read_trajtxt(inp):
    keys = ['step', 'filename', 'index', 'vel']
    nums = list(range(len(keys)))
    traj = {i: [] for i in keys}
    adresses = []
    with open(inp, 'r') as read:
        for idx, line in enumerate(read):
            if '#' in line:
                continue
            rip = line.rstrip().split()
            for num, key in zip(nums, keys):
                traj[key].append(rip[num])
            adresses.append(rip[1])
    return traj, set(adresses)

def read_homos(trajtxt, orig):
    xyzs = list(set(trajtxt['filename']))
    files = [i for i in os.listdir(orig + '/traj') if 'home' in i]
    # for xyz in xyzs:
    #     print('api', xyz[:-4])
    # # print(xyzs)
    # print(orig)
    # print('snow', files)
    return files 
