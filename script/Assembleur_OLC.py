import sys
import os
import numpy as np
from numpy.matlib import empty
from python_tsp.heuristics import solve_tsp_local_search
import time




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
                #else:
                #    if len(sequence) != longueur_ref:
                #        raise ValueError(f"ERREUR : Read de taille différente détecté. "
                #                         f"Attendu: {longueur_ref}, Trouvé: {len(sequence)}")

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
    2. Vérifier si on peut atteindre i depuis j en suivant les successeurs
       (si j mène à i, alors ajouter i->j créerait un cycle)
    """
    if not chemin:
        return False

    # 1. Vérifier les cycles directs (boucles à 2 sommets)
    for arc in chemin:
        src, dst, _ = arc
        # Si l'arc inverse existe déjà (j -> i existe et on veut ajouter i -> j)
        if src == j and dst == i:
            return True

    # 2. Construire un dictionnaire des successeurs
    successeurs = {}
    for arc in chemin:
        src, dst, _ = arc
        successeurs[src] = dst

    # 3. Suivre les successeurs depuis j pour voir si on peut atteindre i
    # Si oui, ajouter i->j créerait un cycle
    courant = j
    visite = set()

    while courant in successeurs:
        if courant in visite:  # Cycle détecté dans le chemin existant
            break
        visite.add(courant)
        courant = successeurs[courant]

        if courant == i:  # On a atteint i en partant de j
            return True  # Ajouter i->j créerait un cycle

    return False

def reorganiser_chemin(chemin):
    """
    Réorganise le chemin pour obtenir une séquence ordonnée de reads.

    Le chemin glouton peut contenir des arcs non connectés ou dans le désordre.
    Cette fonction reconstruit des chaînes linéaires ordonnées.

    Algorithme:
    1. Construire un graphe successeurs[i] = (j, poids)
    2. Trouver les points de départ (reads sans prédécesseur)
    3. Pour chaque point de départ, suivre la chaîne jusqu'au bout
    4. Retourner la chaîne la plus longue

    Paramètres:
        chemin: array de triplets (i, j, poids)

    Retourne:
        list: chemin réorganisé [(i, j, poids), ...]
    """
    if len(chemin) == 0:
        return []

    # Construire le graphe des successeurs et prédécesseurs
    successeurs = {}  # successeurs[i] = (j, poids)
    predecesseurs = set()  # ensemble des reads qui ont un prédécesseur
    tous_reads = set()

    for arc in chemin:
        i, j, poids = arc
        successeurs[i] = (j, poids)
        predecesseurs.add(j)
        tous_reads.add(i)
        tous_reads.add(j)

    # Trouver les points de départ (reads sans prédécesseur)
    points_depart = tous_reads - predecesseurs

    print(f"\nRéorganisation du chemin:")
    print(f"  Nombre d'arcs: {len(chemin)}")
    print(f"  Nombre de reads impliqués: {len(tous_reads)}")
    print(f"  Points de départ trouvés: {len(points_depart)}")

    # Construire toutes les chaînes possibles
    chaines = []
    for depart in points_depart:
        chaine = []
        courant = depart
        visite = set()

        # Suivre la chaîne tant qu'il y a un successeur
        while courant in successeurs and courant not in visite:
            visite.add(courant)
            suivant, poids = successeurs[courant]
            chaine.append([courant, suivant, poids])
            courant = suivant

        if chaine:
            chaines.append(chaine)
            print(f"  Chaîne trouvée: {len(chaine)} arcs, départ={depart}")

    # Retourner la chaîne la plus longue
    if not chaines:
        print("  ATTENTION: Aucune chaîne valide trouvée!")
        return []
    print(f'Il y a {len(chaines)} chaines différentes')
    #chaine_max = max(chaines, key=len)
    #print(f"  Chaîne sélectionnée: {len(chaine_max)} arcs")

    return chaines


def consensus(Reads, chemin_brut):
    """
    Construit les séquences consensus à partir des reads et du chemin brut.

    Paramètres:
        Reads: array de reads (séquences)
        chemin_brut: array de triplets (i, j, poids) sortant de l'algo glouton

    Retourne:
        list: liste des séquences consensus assemblées (str)
    """
    if len(chemin_brut) == 0:
        print("ERREUR: Chemin vide, impossible de construire un consensus")
        return []

    # Étape 1: Réorganiser le chemin (On le fait UNE FOIS pour tout le graphe)
    # Cela renvoie une liste de chaînes (liste de listes)
    toutes_les_chaines = reorganiser_chemin(chemin_brut)

    if len(toutes_les_chaines) == 0:
        print("ERREUR: Impossible de réorganiser le chemin")
        return []

    seqs = []
    print(f"\nAssemblages des consensus ({len(toutes_les_chaines)} chaînes trouvées):")

    # Étape 2: Boucler sur chaque chaîne trouvée pour créer sa séquence
    for index, chaine_ordonnee in enumerate(toutes_les_chaines):

        # Premier assemblage pour cette chaîne
        i0, j0, p0 = chaine_ordonnee[0]
        seq = Reads[i0] + Reads[j0][p0:]

        # Assemblages suivants
        for k in range(1, len(chaine_ordonnee)):
            i, j, p = chaine_ordonnee[k]
            # On ajoute la partie non-chevauchante du read suivant
            seq = seq + Reads[j][p:]

        print(f"  - Consensus {index + 1}: {len(seq)} bp (formé de {len(chaine_ordonnee) + 1} reads)")
        seqs.append(seq)

    print("Fini")
    return seqs


def tsp_layout(M):
    """
    Utilise un solveur TSP avec l'heuristique 2-opt pour trouver un ordre optimal des reads,
    en maximisant les chevauchements.
    """
    n = M.shape[0]

    # On cherche à MINIMISER un coût.
    # Donc on convertit overlap → coût (inverse)
    max_overlap = np.max(M[M >= 0])   # ignore -1
    cost_matrix = max_overlap - M

    # Important : les diagonales doivent être très grandes (interdit de rester sur soi-même)
    for i in range(n):
        cost_matrix[i, i] = 999999999

    # Résolution TSP avec 2-opt (local search avec perturbation_scheme="two_opt")
    print("Résolution TSP avec l'heuristique 2-opt...")
    chemin, cout_total = solve_tsp_local_search(
        cost_matrix,
        perturbation_scheme="two_opt",
        max_processing_time=None  # Pas de limite de temps
    )

    print("Ordre TSP obtenu (premiers 20 reads) :", chemin[:20], "...")
    print("Coût total :", cout_total)

    # On convertit l'ordre TSP → arcs (i -> j, overlap)
    arcs = []
    for k in range(len(chemin) - 1):
        i = chemin[k]
        j = chemin[k + 1]
        arcs.append([i, j, M[i, j]])

    return np.array(arcs)

'''
seq1 = "ABCDEF"
seq2 = "DEFGHI"
print(f"overlap('{seq1}', '{seq2}') = {overlap(seq1, seq2)}")  # Devrait retourner 3 (DEF)
'''

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

        # Construire le consensus
        print(f"\n{'=' * 60}")
        print("CONSTRUCTION DU CONSENSUS")
        print('=' * 60)
        sequence_consensus = consensus(Reads, chemin)

        # Définition du chemin et du nom de fichier
        # Note: Le chemin remonte d'un dossier (..) puis va dans results/OCL_result
        dossier_sortie = "../results/OCL_result/"
        nom_fichier_sortie = "OLC_Result.fasta"
        chemin_complet = os.path.join(dossier_sortie, nom_fichier_sortie)

        try:
            # Créer le dossier s'il n'existe pas (exist_ok=True évite une erreur s'il existe déjà)
            os.makedirs(dossier_sortie, exist_ok=True)

            # Écriture du fichier
            with open(chemin_complet, "w") as f_out:
                for i, seq in enumerate(sequence_consensus):
                    # Création de l'entête demandé (j'ajoute un index pour différencier les seqs)
                    header = f">Seq{i + 1}_lenght{len(seq)}"
                    f_out.write(f"{header}\n")
                    f_out.write(f"{seq}\n")

            print(f"\nSUCCÈS : Fichier sauvegardé sous :")
            print(f"  {os.path.abspath(chemin_complet)}")

        except OSError as e:
            print(f"\nERREUR : Impossible de créer le fichier ou le dossier.")
            print(f"Détails : {e}")

        #print(sequence_consensus)
        if sequence_consensus:
            print(f"\n{'=' * 60}")
            print("RÉSULTAT FINAL")
            print('=' * 60)
            print(f"Longueur de la séquence consensus: {len(sequence_consensus)} bp")
            print(f"\nDébut de la séquence (100 premiers caractères):")
            print(sequence_consensus[:100])
            if len(sequence_consensus) > 200:
                print(f"\nFin de la séquence (100 derniers caractères):")
                print(sequence_consensus[-100:])



    except Exception as e:
        print(f"Erreur lors de l'extraction : {e}")
        sys.exit(1)
'''


def format_temps(secondes):
    """
    Convertit un temps en secondes en format lisible (heures:minutes:secondes).

    Paramètres:
        secondes (float): temps en secondes

    Retourne:
        str: temps formaté
    """
    heures = int(secondes // 3600)
    minutes = int((secondes % 3600) // 60)
    secs = secondes % 60

    if heures > 0:
        return f"{heures}h {minutes}min {secs:.2f}s"
    elif minutes > 0:
        return f"{minutes}min {secs:.2f}s"
    else:
        return f"{secs:.2f}s"

if __name__ == "__main__":

    # Démarrer le chronomètre
    temps_debut = time.time()

    # --- 1. Gestion des arguments et de l'option --heuristique ---

    # Vérifie si le flag est présent
    use_tsp = "--heuristique" in sys.argv

    # Récupère les arguments en ignorant le flag pour trouver le fichier
    # args_propres contient [nom_script, nom_fichier_fastq]
    args_propres = [arg for arg in sys.argv if arg != "--heuristique"]

    if len(args_propres) < 2:
        print("Usage: python extraction_reads.py <fichier_fastq> [--heuristique]")
        print("Exemple: python extraction_reads.py sequences.fastq --heuristique")
        sys.exit(1)

    nom_fichier = args_propres[1]

    try:
        # Extraire les reads
        Reads, longueur_read = extraction_reads_fastq(nom_fichier)

        # Construire la matrice de chevauchement
        print(f"\nConstruction de la matrice de chevauchement...")
        matrice = matrice_adjacence(Reads)
        print(f"Matrice {matrice.shape[0]}x{matrice.shape[1]} créée")

        # --- 2. Choix de l'algorithme (Glouton vs TSP) ---

        if use_tsp:
            print(f"\n{'=' * 60}")
            print("RECHERCHE DU CHEMIN : HEURISTIQUE TSP")
            print('=' * 60)
            # Appel de la fonction TSP
            chemin = tsp_layout(matrice)

            # Nom de fichier spécifique pour l'heuristique
            nom_fichier_sortie = "OLC_heurist_result.fasta"

        else:
            print(f"\n{'=' * 60}")
            print("RECHERCHE DU CHEMIN : GLOUTON OPTIMISÉ")
            print('=' * 60)
            # Appel de la fonction Glouton par défaut
            chemin = glouton_layout_matrice_optimise(matrice, longueur_read)

            # Nom de fichier par défaut
            nom_fichier_sortie = "OLC_Result.fasta"

        print(f"Chemin trouvé avec {len(chemin)} arcs")

        # Afficher les premiers arcs du chemin
        print(f"\nPremiers arcs du chemin (i -> j, poids):")
        for k in range(min(10, len(chemin))):
            i, j, poids = chemin[k]
            print(f"  {i} -> {j} (chevauchement: {poids})")

        # Construire le consensus
        print(f"\n{'=' * 60}")
        print("CONSTRUCTION DU CONSENSUS")
        print('=' * 60)

        sequence_consensus = consensus(Reads, chemin)

        # --- 3. Sauvegarde du résultat ---

        dossier_sortie = "../results/OCL_result/"
        chemin_complet = os.path.join(dossier_sortie, nom_fichier_sortie)

        try:
            # Créer le dossier s'il n'existe pas
            os.makedirs(dossier_sortie, exist_ok=True)

            # Écriture du fichier
            with open(chemin_complet, "w") as f_out:
                for i, seq in enumerate(sequence_consensus):
                    # Entête dynamique selon la séquence
                    header = f">Seq{i + 1}_lenght{len(seq)}"
                    f_out.write(f"{header}\n")
                    f_out.write(f"{seq}\n")

            print(f"\nSUCCÈS : Fichier sauvegardé sous :")
            print(f"  {os.path.abspath(chemin_complet)}")

        except OSError as e:
            print(f"\nERREUR : Impossible de créer le fichier ou le dossier.")
            print(f"Détails : {e}")

        # Affichage des stats finales
        if sequence_consensus:
            print(f"\n{'=' * 60}")
            print("RÉSULTAT FINAL")
            print('=' * 60)
            print(f"Nombre de contigs générés : {len(sequence_consensus)}")

            # Affichage détaillé
            for i, seq in enumerate(sequence_consensus):
                print(f"\n--- Séquence {i + 1} ({len(seq)} bp) ---")
                print(f"Début : {seq[:100]}")
                if len(seq) > 200:
                    print(f"Fin   : {seq[-100:]}")

    except Exception as e:
        print(f"Erreur critique lors de l'exécution : {e}")
        sys.exit(1)

    finally:
        # Calculer et afficher le temps d'exécution total
        temps_fin = time.time()
        temps_total = temps_fin - temps_debut

        print(f"\n{'=' * 60}")
        print("TEMPS D'EXÉCUTION")
        print('=' * 60)
        print(f"Temps total : {format_temps(temps_total)}")
        print(f"Temps brut  : {temps_total:.3f} secondes")
        print("=" * 60)
