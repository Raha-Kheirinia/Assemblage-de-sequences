# Assemblage-de-se-quences


Assemblage de séquences
1. Quelle est la différence entre un alignement maximisant le
chevauchement et l'alignement global ?
L'alignement global est conçu pour comparer des séquences homologues (apparentées)
sur toute leur longueur. On cherche donc le plus long chemin entre (0,0) et
(n,m). L’alignement maximisant le chevauchement, en revanche, est conçu pour
rechercher dans la séquence A et la séquence B un alignement préfixe/suffixe. On
compare des segments de toutes longueurs et retient celui qui maximise le score de
similitude sur les segments.

2. Quelles doivent être les valeurs de la première ligne et de la
première colonne de la table de programmation dynamique
V ?
Les valeurs de la première colonne seront : ! !
Les valeurs de la première ligne seront : ! !
Ces valeurs permettront à l’alignement de commencer n’importe où sur une des
deux séquences.
3. Quelles sont les équations de récurrence à utiliser pour remplir
la table de programmation dynamique ?
∀j V(0, j) = 0
∀i V(i,0) = 0
Zohreh Kheirinia Tanguy Invernizzi
(match)
(dismatch)
V(i, j) = max


4. Comment peut-on retrouver l’alignement avec le meilleur
chevauchement a partir de la table de programmation dynamique
?
Pour retrouver le meilleur alignement à partir de la table de programmation dynamique,
on cherche tout d’abord la valeur maximum dans les deux bords ! et
! , de couleur rouge sur la table exemple. En effet, pour trouver un suffixe et un
préfixe, l’alignement doit forcément être construit a partir de la dernière ligne ou de la
dernière colonne de la table, car elles correspondent à la fin de l’une de deux
séquences.
Figure 1 : Table exemple de programmation dynamique du meilleur chevauchement
Après avoir trouvé le maximum des deux bords étudiés, il suffit de remonter les
pointeurs tant que ! et ! . Quand l’une des deux conditions n’est plus vraie,
cela signifie que l’on a atteint la fin de l’une des deux séquences. Dans le cas où ! ,
c’est ! qui sera notre préfixe et ! qui occupera la place du suffixe. Dans le cas où
! , c’est ! qui sera notre préfixe et ! qui occupera la place du suffixe.

5. Application de l’algorithme en Python : TP1_Partie1.py
Pré-requis techniques :
• Python 2.7.13 ou équivalent
• Numpy
Le programme TP1_Partie1.py supporte les fichiers .fq et .fasta. Les fichiers
doivent absolument respecter les conventions de ces formats, tel qu’elles ont été utilisées
dans les fichiers exemples du TP (geneX.fasta et reads.fq). Seules les deux pre-
(imax, j)
(i, jmax)
D E F A B
0 0 0 0 0 0
C 0 -4 -4 -4 -4 -4
B 0 -4 -8 -8 -8 0
A 0 -4 -8 -12 -4 -8
D 0 4 -4 -12 -12 -8
E 0 -4 8 0 -8 -16
F 0 -4 0 12 4 -4
x ≥ 0 y ≥ 0
x ≤ 0
x y
y ≤ 0 y x


mières séquences du fichier seront prises en compte. Pour étudier deux séquences
contenues dans un fichier, ouvrez une ligne de commande, placez-vous dans le dossier
contenant TP1_Partie1.py et votre fichier. Entrez la commande suivante, en remplaçant
<nomdufichier> par le nom du fichier a étudier :
python TP1_Partie1.py <nomdufichier>
Le programme TP1_Partie1.py supporte également l’étude directe de deux
séquences entrées dans la ligne de commande. Entrez la commande suivante, en remplaçant
<sequence1> et <sequence2> par les deux séquences a étudier :
python TP1_Partie1.py <sequence1> <sequence2>
Assemblage de fragments

1. Pour chaque paire de reads ! , calculer le score de
l’alignement correspondant au chevauchement maximal entre
! et ! .
Le but de l’algorithme étant d’obtenir des alignements pour entre différentes
paires de reads, il est inutile d’essayer d’aligner une séquence avec elle-même. La diagonale
de la matrice résultant de l’exécution de l’algorithme sera donc vide. De même,
aligner ! avec ! puis ! avec ! est inutile. On choisit donc une convention : la partie
supérieure droite obtient les scores dont les reads ! sont préfixes, et la partie inférieure
gauche obtient les scores dont les reads ! sont préfixes.
Figure 2 : Matrice obtenue après utilisation de l’algorithme calculant chaque paire de reads.
{Rx, Ry}
Rx Ry
Rx Ry Ry Rx
Rx
Ry


2. En déduire le graphe orienté de chevauchement de l’ensemble
des reads.

Figure 3 : Graphe orienté issu de l’algorithme de calcul de chevauchement de l’ensemble des
reads.
Le graphe orienté issu de la matrice, comme il est possible de le remarquer sur
la figure 3, possède des cycles. Cette caractéristique nous empêche de créer un enchainement
de séquences cohérent. Pour remédier a ce problème, il faut supprimer les
arêtes ayant un score trop faible.

(a) En déduire le graphe orienté de chevauchement de l’ensemble
des reads, en considérant un seuil minimum de 80.


Comme nous pouvons le voir sur la figure 4, le fait de filtrer l’ensemble des
scores pour ne garder que ceux supérieurs a 80 à pour effet de créer deux graphes orienté
acyclique différents.

(b) En déduire l’ensemble des reads reverse
On peux d’hors et déjà supposer que ces deux sous-ensembles correspondent
au reads reverse et forward. Si l’on considère que le read @READS_2 est définitivement
sur le brin forward, et que celui-ci est contenu dans le graphe du haut, on peut en déduire
que le graphe du bas correspond aux reads reverse.

3. Remplacer les reads reverse par leurs équivalents complémentaires
inverses

(a) En déduire le graphe orienté de chevauchement de l’ensemble
des reads, en considérant un seuil minimum de 80.
En remplaçant les séquences 5, 6, 8, 9, 11, 14 et 17 par leurs séquences
complémentaire inverse, nous obtenons le graphe suivant :
Figure 5 : Graphe orienté issu de l’algorithme de calcul de chevauchement de l’ensemble des
reads forwards avec les reads reverse complémentaires inverses, avec un seuil minimum de 80.
(b) Appliquer la réduction transitive sur le graphe obtenu.

(c) En déduire la séquence du fragment génomique séquencé
et sa longueur.
La séquence du fragment génomique sera de longueur
! où
! est le nombre de sommets du graphe réduit transitivement,
! sont les cotés du graphes et
! est la longueur de chevauchement du coté du graphe
considéré.
Pour retrouver la séquence du fragment, il suffit d’assembler un à un
chacune des séquences, dans l’ordre ou elles apparaissent dans le graphe
obtenu de la réduction transitive.
Recherche d’introns et Blast
1. Identifier la position de la protéine X
La premiere étape pour l’identification de la position de la protéine X est la traduction
de la sequence nucléotidique du fichier sequence.fasta. Pour cela, on va se
servir de la table de code génétique. Le cadre de lecture est une région de l’ADN située
entre un codon START (ou Methyonine) et un codon STOP (TAA, TAG, TGA).

traduction, on a remplacé les codons STOP par des caractères “ * ”.
Traduction :
Étape 1) On divise la séquence dans leurs codons
Étape 2) On traduit chaque codon en acide aminé
Étape 3) On assemble la chaîne des acides aminés

(a) Dans quel cadre de lecture se trouve le codon START de la séquence protéique ?
En traduisant la séquence protéique dans les 3 cas de lecture, on trouve le
codon START dans le deuxième cas de lecture. Il se trouve dans l’intervalle de [60: 115]
de la séquence protéine.
Voici la traduction de la séquence de protéine dans 3 cas différents :

(b) Algo- rithme
L’algorithme demandé ne considère pas les in- dels. On initialise donc la matrice avec des zéros sur la première ligne et la première colonne. Dans la matrice de score, on vérifie seule- ment la diagonale pour vérifier s’il y a un match. Si c’est le cas, on ajoute 1, sinon on passe au prochain caractère . Il est important de considérer uniquement les chemins des diagonales.

(c) Intervalle de position contenant les exons
Les intervalles exon de geneX et leur correspondance dans la séquence de protéine et la
séquence ARNm sont les suivants :
Pour trouver les intervalles correspondant dans la séquence ADN(ARNm), on a multiplié les intervalles
de séquence protéine par 3, auxquelles on as ajouté 1 ou 2, en fonction de son cas de
lecture. On a ensuite réalisé une concaténation entre les parties trouvées.
sequence_dna[66*3+1:121*3+1]
+ sequence_dna[287*3+2:324*3+2]
+ sequence_dna[480*3+1:488*3+1]
Enfin, on traduit cette sequence a ARNm en remplaçant les caractères 'T' par des ‘U' :
______________________________________ARN___________________________________
AUGUGCCAGCGUUGUGGUUUAAAACUAAUAGUAAUAAUAUGCUUCUUUGUUCAGUUGGCUAGAGAUUUACUACAUCCGUCCUUGGAAGAGGAAAAGAAAAAACAUAAAAAGAAACGCCUAGUACA
AAGUCCAAAUUCUUACUUUAUGGAUGUAAAAUGUCCAGGUGCUACAAGAUCACCACGGUUUUCAGCCAUGCUCAGACAGUGGUUCUUUGUGUAGGUUGUUCAACAGUGUUGUGCCAGCCUACAGG
AGGAAAGGCCAGACUCACAGAAGGAUUGUUCAUUUAGAAGAAAGCAACAC
____________________________________________________________________________
2. Identifier le nom de la protéine X
Nom : Ribosomal protein S27 [ Homo sapiens (human) ]
Description : Le gène RPS27 comprend 1,39 kb et se compose de 4 exons. Ce
gène est un membre du CCDS humain: CCDS1059.
