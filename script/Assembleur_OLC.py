import sys


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
    reads = []
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

                reads.append(sequence)
                compteur += 1

        print(f"Extraction réussie : {compteur} reads extraits")
        print(f"Longueur des reads : {longueur_ref} bp")

    except FileNotFoundError:
        print(f"ERREUR : Le fichier '{fichier_fastq}' n'existe pas")
        raise

    return reads




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
        reads = extraction_reads_fastq(nom_fichier)



    except Exception as e:
        print(f"Erreur lors de l'extraction : {e}")
        sys.exit(1)
