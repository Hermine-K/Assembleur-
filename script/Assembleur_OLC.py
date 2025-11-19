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


seq1 = "ABCDEF"
seq2 = "DEFGHI"
print(f"overlap('{seq1}', '{seq2}') = {overlap(seq1, seq2)}")  # Devrait retourner 3 (DEF)

