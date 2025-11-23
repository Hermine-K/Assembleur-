# Cours HAU902I : Bioinformatique avancé  
### Projet : Algorithme de Reconciliation | Fait par : Conceptia Dagba Allade ; Hermine Kiossou ; Homero Sanchez


# Algorithme d’Assemblage

L’algorithme d’assemblage repose sur le principe de reconstruire une séquence d’ADN ou d’ARN à partir de fragments courts issus du séquençage. L’objectif est d’exploiter les chevauchements entre ces reads pour retrouver, autant que possible, la séquence originale. Cette étape, réalisée entièrement in silico, suit immédiatement le séquençage d’un organisme, d’une population clonale (comme une culture bactérienne) ou d’un mélange complexe.

Dans ce projet, nous avons réimplémenté un assemblage de type OLC (Overlap–Layout–Consensus). Cette approche nous a permis de suivre en détail chaque phase de l’algorithme (du calcul des chevauchements à la construction du consensus) et d’évaluer son fonctionnement sur des jeux de données réels.

## 1. Overlap-Layout-Consensus

Dans cette méthode, on exploite le calcul des chevauchements pour construire un graphe permettant d’identifier les lectures qui se chevauchent (overlap), de les organiser en une séquence continue (layout) puis de corriger les erreurs afin d’obtenir une séquence consensus.

* La phase **Overlap** consiste à calculer les chevauchements optimaux entre les reads par **alignement semi-global**, à construire progressivement un graphe orienté reliant les reads par des arêtes pondérées selon la qualité du chevauchement, et à exclure les reads entièrement contenus dans d’autres. On peut coder le graph de chevauchement sous la forme de **matrice d’adjacence** où la case **M[i][j]** stocke le poids du chevauchement du read i avec le read j.

* La phase **Layout** vise à trouver un ordre optimal des reads en résolvant un problème de type TSP, que l’on peut aborder soit par des méthodes exactes (comme la programmation dynamique ou le branch and bound), soit par des heuristiques plus rapides (glouton, plus proche voisin, k-opt, Lin Kernighan, etc.), en adaptant le problème asymétrique en un TSP symétrique.

* La phase **consensus** consiste à effectuer un alignement multiple des séquences pour en obtenir une version optimisée en score, mais comme l’alignement exact devient rapidement impraticable au-delà d’une dizaine de séquences, on utilise généralement des méthodes heuristiques.

## 2. Avantages et Inconvénients 

### a. Avantages

* Très adapté aux longues lectures : Cette approche fonctionne particulièrement bien avec les longues lectures.

* Assemblage généralement plus continus : Les longues lectures peuvent traverser des zones répétées du génome sans se casser en morceaux, ce qui permet de reconstruire des séquences plus longues et avec moins de coupures.

* Modèle explicite de l’assemblage : Le graphe OLC offre une vision claire des chevauchements entre les lectures et facilite les étapes de correction ou de vérification.

* Résistant aux erreurs systématiques : Les erreurs aléatoires des lectures longues sont largement corrigées lors de la phase de consensus final.

---

### b. Inconvénients

* Coût computationnel élevé : Comparer toutes les lectures entre elles demande beaucoup de temps et de mémoire, ce qui nécessite des optimisations.

* Peu adapté aux très grands jeux de données à courtes lectures : Avec de très nombreuses lectures courtes, l’OLC devient trop lourd et les graphes de Bruijn sont bien plus efficaces.

* Sensibilité aux régions très répétées : Certaines répétitions complexes peuvent créer des ambiguïtés dans le graphe, même avec des lectures longues.

* Pipeline plus complexe : Le processus OLC comporte plusieurs étapes distinctes qu’il peut être plus difficile de paramétrer et d’optimiser que dans un pipeline de Bruijn.


## 3. Pipeline de l'outil d'assemblage OLC

![Pipeline OLC](images/pipeline_olc.png)

## 4. Choix du langage de programation 
















⚠ :

      ALGORITHME: reconciliation 
      ENTREE: Arbre_gene, Arbre_espece
      SORTIE: Affiche Arbre reconcilié
      
   	\\Associe chaque feuille de l'arbre de gènes à sa fauille correspondante dans l'arbre d'espece
   	FONCTION initialize_mapping(Arbre_gene, Arbre_espece)
   		POUR Feuille DANS Arbre_gene FAIRE
   			Nom_gene <- Feuille.name  
   			Nom_espece <- Nom_gene[0:3] \\recupere les noms des espece dans le nom des genes
   			Feuille_espece <- Arbre_espece.recherche_feuille(Nom_espece) //fonction qui trouve la feuille avec le nom donné
   			Feuille.Ajouter_feature("M", Feuille_espece) //Associe le M(g) de chaque Feuille
   		FIN POUR
   
   	\\Calcule du dernier ancêtre commun (LCA) de deux noeuds
   	FONCTION compute_lca(Noeud_1, Noeud_2)
   		SI Noeud_1 = Noeud_2 ALORS
   			RETOURNER Noeud_1
   		FIN SI
   		\\ Construire le chemin depuis Noeud_1 jusqu'à la racine
   		Chemin_1 <- []
   		Noeud_courant <- Noeud_1
   		TANT QUE Noeud_courrant != None ALORS
   			Chemin_1.Ajouter(Noeud_courant)
   		FIN TANT QUE
   		Chemin_1.reverse() \\interverti l'ordre de la liste
   
   		\\ Construire le chemin depuis Noeud_2 jusqu'à la racine
                   Chemin_2 <- []
                   Noeud_courant <- Noeud_2
                   TANT QUE Noeud_courrant != None ALORS
                           Chemin_2.Ajouter(Noeud_courant)
   		FIN TANT QUE
                   Chemin_2.reverse() \\interverti l'ordre de la liste
   
   		\\Trouver le LCA
   		LCA <- None
   		POUR i ALLANT de 0 À Min(|Chemin_1|, |Chemin_2|) FAIRE
   			SI Chemin_1[i] = Chemin_2[i] ALORS
   				LCA = Chemin_1[i]
   			SINON FAIRE
   				ARRET
   			FIN SI
   		FIN POUR
   		RETOURNER LCA
   
   	\\Calcule le M(g) pour tout les noeuds internes et classifie les evenement en duplication ou spéciation
   	FONCTION compute_mapping_and_classify(Arbre_gene)
   		POUR Noeud DANS Arbre_gene \\la recherche se fait des feuille vers la racine de telle sorte qu’un nœud est traité lorsque tous ses fils sont traités
   			SI NON Noeud.est_feuille() ALORS
   				Enfant_1, Enfant_2 <- Noeud.recuperer_enfants()
   				\\Calcule M(g) de Noeud
   				M_1, M_2 <- recuperer_feature("M") \\Recupere le M(g) des 2 enfants
   				LCA_Noeud <- compute_lca(M_1,M_2) \\Calculer le LCA dans l'arbre d'espèces
   				Noeud.Ajouter_feature("M", LCA_Noeud) //Associe le M(g) de chaque Noeud
   				\\Classifie entre duplication et spéciation
   				SI LCA_Noeud = M_1 OU LCA_Noeud = M_2 ALORS
   					 Noeud.Ajouter_feature("Type", "DUPLICATION") //Associe le type DUPLICATION au noeud
   				SINON
   					 Noeud.Ajouter_feature("Type", "SPECIATION") //Associe le type SPECIATION au noeud
   	\\Appel des fonction et Affichage finale		
   	initialize_mapping(Arbre_gene, Arbre_espece)
   	compute_mapping_and_classify(Arbre_gene)
   	Afficher(Arbre_gene)



## Références

[1]  Réconciliation phylogénétique. (2023, février 12). Wikipédia, l'encyclopédie libre. Page consultée le 20:21, février 12, 2023 à partir de http://fr.wikipedia.org/w/index.php?title=R%C3%A9conciliation_phylog%C3%A9n%C3%A9tique&oldid=201335036.

[2]  ETE Toolkit Documentation: https://etetoolkit.org/docs/latest/tutorial/tutorial_trees.html.
