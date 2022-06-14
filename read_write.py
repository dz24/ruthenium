"""remove Xs in a periodicly generated traj for jmol."""


def reader(inp, out):
    """Read text and write files."""
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


# INPUT_ = './100739_moviePBC.xyz'
# OUTPUT = './100739_moviePBC_noX.xyz'
# reader(INPUT_, OUTPUT)
