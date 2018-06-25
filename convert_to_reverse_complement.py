import numpy

# stop codant{'TAA', 'TAG', 'TGA'} sont presentes par *
dict_com = {
    'A':'T',
    'T':'A',
    'U':'A',
    'G':'C',
    'C':'G',
    }

def convert_to_reverse_complement(seq):
    complement = ''
    for i in range (0, len(seq)):
        complement = complement + dict_com[seq[i]]
    return complement[::-1]

def chevauchement_sequence(mot1, mot2, match=4, mismatch=-4, indel=-8):
    d = numpy.empty((len(mot1) + 1, len(mot2) + 1)).astype(int) #Stockage des valeurs de la grille
    fleches = numpy.empty((len(mot1) + 1, len(mot2) + 1, 3)).astype(int) #Stockage des fleches de la grille [gauche, diagonale, haut]
    for i in range(0, len(mot1) + 1):
        for j in range(0, len(mot2) + 1):
            d[i][j] = 0
    current_max = [0, 0, 0] #Coordonnees x et y du max, ainsi que sa valeur
    for i in range(1, len(mot1) + 1):
        for j in range(1, len(mot2) + 1):
            if (mot1[i - 1] == mot2[j - 1]):
                d[i][j] = d[i - 1][j - 1] + match
                fleches[i][j] = [0, 1, 0]
            else:
                maxi = max(d[i][j-1] + indel, d[i - 1][j - 1] + mismatch, d[i - 1][j] + indel)
                if (maxi == d[i][j-1] + indel):
                    d[i][j] = d[i][j-1] + indel
                    fleches[i][j] = [1, 0, 0]
                elif (maxi == d[i - 1][j - 1] + mismatch):
                    d[i][j] = d[i - 1][j - 1] + mismatch
                    fleches[i][j] = [0, 1, 0]
                elif (maxi == d[i - 1][j] + indel):
                    d[i][j] = d[i - 1][j] + indel
                    fleches[i][j] = [0, 0, 1]
    for i in range (1, len(mot1) + 1):
        if (d[i][len(mot2)] > current_max[2]):
            current_max[0] = i
            current_max[1] = len(mot2)
            current_max[2] = d[i][len(mot2)]
    for j in range (1, len(mot2) + 1):
        if (d[len(mot1)][j] > current_max[2]):
            current_max[0] = len(mot1)
            current_max[1] = j
            current_max[2] = d[len(mot1)][j]
    str1 = ''
    str2 = ''
    x = current_max[0] - 1
    y = current_max[1] - 1
    while (x >= 0 and y >= 0):
        if numpy.array_equal(fleches[x][y], [1, 0, 0]): #gauche
            str1 = str1 + mot1[x]
            str2 = str2 + '_'
            y = y - 1
        elif numpy.array_equal(fleches[x][y], [0, 1, 0]): #diagonale
            str1 = str1 + mot1[x]
            str2 = str2 + mot2[y]
            y = y - 1
            x = x - 1
        elif numpy.array_equal(fleches[x][y], [0, 0, 1]): #haut
            str1 = str1 + '_'
            str2 = str2 + mot2[y]
            x = x - 1
        else:
            str1 = str1 + mot1[x]
            str2 = str2 + mot2[y]
            y = y - 1
            x = x - 1
    return current_max

file = open('reads_with_reverse.fq', 'r')
read = file.read()

read = read.split()
seq = []
for i in range (1, len(read), 4):
    seq.append(read[i])

def assemblage_fragment(my_lst):
    l = len(my_lst)-1
    matrix = numpy.zeros(shape=(20,20))
    for i in range(0, len(matrix)):
        for j in range(0, len(matrix)):
            if (i != j):
                current_max = chevauchement_sequence(my_lst[i],my_lst[j])  #[prefixe = 0/suffixe = 1, score]
                if (current_max[0] == len(my_lst[i])):
                    matrix[i][j] = current_max[2]
                else:
                    matrix[i][j] = 0
    print matrix
    return matrix

def convert_edgelist(matrix):
    for i in range(0, len(matrix)):
        for j in range(0, len(matrix)):
            if (matrix[i][j] >= 80):
                print 'seq' + str(i) + '	' + 'seq' + str(j) + '	' + str(matrix[i][j])

matrix = assemblage_fragment(seq)
convert_edgelist(matrix)
#print convert_to_reverse_complement('ACTGCA')
