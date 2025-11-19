import sys
import numpy as np
from numpy.matlib import empty


def extraction_reads_fastq(fichier_fastq):
    """
    Extrait les reads d'un fichier FASTQ et vérifie qu'ils ont tous la même longueur.

    Paramètres:
        fichier_fastq : chemin vers le fichier FASTQ

    Retourne:
        list: liste des reads

    Erreur si:
        ValueError: si un read a une taille différente des autres
        FileNotFoundError: si le fichier n'existe pas
    """
    Reads = []
    longueur_ref = 0
    compteur = 0

    try:
        with open(fichier_fastq, 'r') as fichier:
            while True:
                # Lire les 4 lignes d'un read FASTQ
                ligne = fichier.readline()
                if not ligne:  # Fin du fichier
                    break

                # Ignorer header (@)
                sequence = fichier.readline().strip()  # Enlever le \n
                fichier.readline()  # Ignorer ligne "+"
                fichier.readline()  # Ignorer qualité

                # Vérifier la longueur des reads
                if compteur == 0:
                    longueur_ref = len(sequence)
                else:
                    if len(sequence) != longueur_ref:
                        raise ValueError(f"ERREUR : Read de taille différente détecté. "
                                         f"Attendu: {longueur_ref}, Trouvé: {len(sequence)}")

                Reads.append(sequence)
                compteur += 1

        print(f"Extraction réussie : {compteur} reads extraits")
        print(f"Longueur des reads : {longueur_ref} bp")

    except FileNotFoundError:
        print(f"ERREUR : Le fichier '{fichier_fastq}' n'existe pas")
        raise
    # Convertir la liste en numpy array
    return np.array(Reads), longueur_ref




def overlap(A, B):
    """
    Calcule la longueur du plus long suffixe de A qui est un préfixe de B.

    Paramètres:
        A (str): première séquence
        B (str): deuxième séquence

    Retourne:
        int: longueur du chevauchement maximal
    """
    o = 0
    lenA = len(A)
    lenB = len(B)
    max_possible = min(lenA, lenB)

    for k in range(1, max_possible + 1):
        # Compare le suffixe de longueur k de A avec le préfixe de longueur k de B
        if A[lenA - k:lenA] == B[0:k]:
            o = k

    return o


def matrice_adjacence(Reads):
    """
    Construit la matrice d'adjacence représentant les chevauchements entre reads.

    Paramètres:
        F : tableau de chaînes de caractères (reads)

    Retourne:
        list: matrice d'entiers à deux dimensions (liste de listes)
    """
    n = len(Reads)
    # Initialiser la matrice n x n
    M = np.zeros((n, n), dtype=int)

    for i in range(n):
        for j in range(n):
            if i == j:
                M[i, j] = -1
            else:
                M[i, j] = overlap(Reads[i], Reads[j])

    return M


def glouton_layout_matrice(M, len_read):
    """
    Trouve un chemin hamiltonien approximatif en utilisant un algorithme glouton.

    Paramètres:
        M (numpy.ndarray): matrice d'adjacence des chevauchements

    Retourne:
        numpy.ndarray: tableau de triplets (i, j, poids) représentant le chemin
    """
    n = M.shape[0]
    chemin = []
    m = 0

    while m < n :
        i_max = 0
        j_max = 0
        max_val = -1

        # Chercher le poids maximal dans la matrice
        for i in range(n):
            for j in range(n):
                if i != j and M[i, j] > max_val:
                    i_max = i
                    j_max = j
                    max_val = M[i, j]


        # Enregistrer l'arc avec son poids
        if chemin==[]:
            chemin.append([i_max, j_max, max_val])
        elif max_val!= 0 and (chemin[-1][0]!=j_max or chemin[-1][1]!=i_max) and max_val!=len_read: #evite les chevauchement nul et les boucles à 2 reads
            #print(f'on ajoute {max_val}')
            chemin.append([i_max, j_max, max_val])

        # Supprimer les arcs liés à i_max et j_max
        for k in range(n):
            M[i_max, k] = -1
            M[k, j_max] = -1

        m += 1

    # Convertir la liste en numpy array
    return np.array(chemin)


def glouton_layout_matrice_optimise(M, len_read):
    """
    Trouve un chemin hamiltonien approximatif en utilisant un algorithme glouton optimisé.

    Améliorations:
    1. Suivi des degrés entrants/sortants pour éviter les cycles
    2. Détection de cycles pour éviter les petites boucles
    3. Construction de chaînes linéaires (chaque read max 1 prédécesseur et 1 successeur)

    Paramètres:
        M (numpy.ndarray): matrice d'adjacence des chevauchements
        len_read (int): longueur des reads

    Retourne:
        numpy.ndarray: tableau de triplets (i, j, poids) représentant le chemin
    """
    n = M.shape[0]
    chemin = []

    # Suivi des degrés pour construire un chemin (pas de cycles)
    degre_sortant = np.zeros(n, dtype=int)  # Nombre de successeurs
    degre_entrant = np.zeros(n, dtype=int)  # Nombre de prédécesseurs

    # Copie de la matrice pour pouvoir la modifier
    M_copy = M.copy()

    arcs_rejetes_degre = 0
    arcs_rejetes_cycle = 0
    m = 0

    while m < n:
        i_max = 0
        j_max = 0
        max_val = -1

        # Chercher le poids maximal dans la matrice
        for i in range(n):
            for j in range(n):
                if i != j and M_copy[i, j] > max_val:
                    i_max = i
                    j_max = j
                    max_val = M_copy[i, j]

        # Vérifier les conditions d'ajout
        if max_val != 0 and max_val != len_read:
            # Vérifier les contraintes de degré
            if degre_sortant[i_max] >= 1 or degre_entrant[j_max] >= 1:
                arcs_rejetes_degre += 1
            # Vérifier si on crée un cycle
            elif creer_cycle(chemin, i_max, j_max):
                arcs_rejetes_cycle += 1
            else:
                # Ajouter l'arc au chemin
                chemin.append([i_max, j_max, max_val])
                degre_sortant[i_max] += 1
                degre_entrant[j_max] += 1

                if len(chemin) % 100 == 0:
                    print(f"  Progression: {len(chemin)} arcs ajoutés...")

        # Supprimer les arcs liés à i_max et j_max
        for k in range(n):
            M_copy[i_max, k] = -1
            M_copy[k, j_max] = -1

        m += 1

    print(f"\nStatistiques:")
    print(f"  Arcs rejetés (degré): {arcs_rejetes_degre}")
    print(f"  Arcs rejetés (cycle): {arcs_rejetes_cycle}")
    print(f"  Arcs acceptés: {len(chemin)}")

    return np.array(chemin)


def creer_cycle(chemin, i, j):
    """
    Vérifie si ajouter l'arc (i, j) créerait un cycle dans le chemin actuel.

    Principe:
    1. Vérifier les cycles directs (A->B déjà présent, on veut ajouter B->A)
    2. Remonter depuis j pour voir si on peut atteindre i
    Si oui, ajouter (i, j) créerait un cycle.
    """
    if not chemin:
        return False

    # 1. Vérifier les cycles directs (boucles à 2 sommets)
    for arc in chemin:
        src, dst, _ = arc
        # Si l'arc inverse existe déjà (j -> i existe et on veut ajouter i -> j)
        if src == j and dst == i:
            return True

    # 2. Construire un dictionnaire des prédécesseurs
    predecesseurs = {}
    for arc in chemin:
        src, dst, _ = arc
        predecesseurs[dst] = src

    # 3. Remonter depuis j
    courant = j
    visite = set()

    while courant in predecesseurs:
        if courant in visite:  # Cycle détecté dans le chemin existant
            break
        visite.add(courant)
        courant = predecesseurs[courant]

        if courant == i:  # On a atteint i en remontant depuis j
            return True

    return False

'''
seq1 = "ABCDEF"
seq2 = "DEFGHI"
print(f"overlap('{seq1}', '{seq2}') = {overlap(seq1, seq2)}")  # Devrait retourner 3 (DEF)
'''


if __name__ == "__main__":
    # Vérifier qu'un fichier a été fourni en argument
    if len(sys.argv) < 2:
        print("Usage: python extraction_reads.py <fichier_fastq>")
        print("Exemple: python extraction_reads.py sequences.fastq")
        sys.exit(1)

    # Récupérer le nom du fichier depuis les arguments
    nom_fichier = sys.argv[1]

    try:
        # Extraire les reads
        Reads, longueur_read = extraction_reads_fastq(nom_fichier)

        # Construire la matrice de chevauchement
        print(f"\nConstruction de la matrice de chevauchement...")
        matrice = matrice_adjacence(Reads)
        print(f"Matrice {matrice.shape[0]}x{matrice.shape[1]} créée")

        # Afficher un échantillon de la matrice
        print(f"\nÉchantillon de la matrice (3x3 premiers éléments):")
        taille = min(10, matrice.shape[0])
        print(matrice[:taille, :taille])

        # Trouver le chemin hamiltonien
        print(f"\nRecherche du chemin hamiltonien avec algorithme glouton...")
        chemin = glouton_layout_matrice_optimise(matrice, longueur_read)
        print(f"Chemin trouvé avec {len(chemin)} arcs")

        # Afficher les premiers arcs du chemin
        print(f"\nPremiers arcs du chemin (i -> j, poids):")
        for k in range(min(10, len(chemin))):
            i, j, poids = chemin[k]
            print(f"  {i} -> {j} (chevauchement: {poids})")

        if len(chemin) > 10:
            print(f"  ... et {len(chemin) - 10} autres arcs")

    except Exception as e:
        print(f"Erreur lors de l'extraction : {e}")
        sys.exit(1)
