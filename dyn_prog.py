def dyn_prog(i, j, cts, protein, table):
    if table[i][j] != None:
        return table[i][j]
    table[i][j] = (dyn_prog(i-1, j-1, cts, protein, table)+1) if cts[j-1] == protein[i-1] else 0     return table[i][j]


#cts (complete translated sequence) est la sequence de l'ADN complete traduite en acide amine
def algo_q3(cts, protein):
    table = []
    for y in range(0, len(protein)+1):
        table.append([])
        for x in range(0, len(cts)+1):
            table[y].append(None)

    for x in range(0, len(cts)+1):
        table[0][x] = 0
    for y in range(0, len(protein)+1):
        table[y][0] = 0

    #Remplir la table
    for y in range(1, len(protein)+1):
        for x in range(1, len(cts)+1):
            dyn_prog(y, x, cts, protein, table)

    return table
