# Cours HAU902I : Bioinformatique avancé  
### Projet : Algorithme de Reconciliation | Fait par : Conceptia Dagba Allade ; Hermine Kiossou ; Homero Sanchez


# Algorithme de Reconciliation 

La réconciliation phylogénétique consiste à faire correspondre l’arbre des espèces et celui des gènes.
Certains gènes suivent une évolution différente du reste du génome, à cause d’événements comme la duplication, la perte ou la spéciation, ce qui peut compliquer la construction de l’arbre des espèces [1].

Cette méthode modélise ces événements pour mieux comprendre l’histoire du génome.
Notre algorithme réalise cette réconciliation avec la librairie ETE3 [2], qui exploite le format Newick, un standard bioinformatique pour représenter la structure, les distances et les nœuds d’un arbre phylogénétique.


## 1. Explication du Code  


### a. Parse Arguments et Load  Tree

La première étape de notre algorithme de réconciliation consiste à analyser les arguments passés au script.
Cette fonction permet d’indiquer, au moment de l’exécution, les fichiers contenant les arbres d’espèces et les arbres de gènes, ou bien de les fournir directement via la ligne de commande.

Elle s’appuie sur les librairies os (`import os`), sys (`import sys`), argparse (`import argparse`).

La fonction gère plusieurs options :

* `--loss` : active la prise en compte des pertes de gènes

* `--verif` : lance une vérification des arbres fournis

* `-h` : affiche l’aide automatique générée par argparse

En complément, la fonction **load Tree** permet de charger les arbres de gènes et d’espèces à partir d’un fichier ou d’une chaîne au format Newick, format standard pour représenter la structure hiérarchique d’un arbre phylogénétique.

---

### b. Initialize Mapping

La fonction initialize_mapping() constitue la première étape du processus de réconciliation. Elle permet d’établir une correspondance initiale entre les feuilles de l’arbre des gènes et celles de l’arbre des espèces. Cette fonction ajoute une étiquette `M(g)` à chaque feuille de l'arbre de gène. 

---

### c. Compute Lca

La fonction `compute_lca` permet de trouver le dernier ancêtre commun (LCA) de deux nœuds dans l’arbre des espèces.
Elle fonctionne en plusieurs étapes :

* Pour chaque nœud, elle construit la liste de tous ses ancêtres jusqu’à la racine.

* Elle inverse ces listes pour aller de la racine vers les nœuds.

* Elle compare les listes élément par élément pour identifier le dernier nœud partagé, qui correspond au LCA.

Cette fonction est utilisée dans notre algorithme pour déterminer le point d’origine commun de deux gènes dans l’arbre des espèce.

---

### d. Compute Mappings and Classify

La fonction `compute_mappings_and_classify()` correspond à la deuxième étape de l’algorithme de réconciliation, appelée phase montante. Elle sert à déterminer la correspondance M(g) pour les nœuds internes de l’arbre de gènes et à classer chaque nœud comme un événement de duplication ou de spéciation.

L’arbre de gènes est parcouru en post-ordre (des feuilles vers la racine). Pour chaque nœud interne, la fonction récupère les correspondances de ses enfants dans l’arbre des espèces, détermine leur dernier ancêtre commun (LCA) et assigne ce LCA au nœud courant. Si le LCA correspond à l’un des enfants, l’événement est une duplication ; sinon, une spéciation.

---

### e. Display Tree ASCII et Reconciliation

La fonction `display_tree_ascii()` affiche dans le terminal l’arbre de gènes réconcilié sous forme ASCII, accompagné des annotations associées à chaque nœud.
Pour chaque nœud, elle affiche soit :

* le nom du gène et sa correspondance d’espèce pour les feuilles,

* soit le type d’événement (duplication ou spéciation) et la correspondance M(g) pour les nœuds internes.

Elle produit également un résumé final du nombre de duplications et de spéciations identifiées.


La fonction principale `reconciliation()` appelle l’ensemble des fonctions citées plus haut. 

---

### f. Option Verif and Loss

La fonction `option_verif_et_loss()` gère les options facultatives `--verif` et `--loss` du programme.
Elle utilise la méthode `reconcile()` de la librairie ETE3 pour effectuer une réconciliation automatique entre l’arbre de gènes et l’arbre d’espèces, puis affiche les événements de spéciation et de duplication identifiés.
Selon les options activées, elle peut aussi afficher graphiquement l’arbre de gènes original (`--verif`) ou l’arbre réconcilié avec pertes (`--loss`). Ces options servent à vérifier notre implémentation. 

La fonction `main()` constitue le point d’entrée du programme. Elle analyse les arguments, vérifie les options, puis lance soit `option_verif_et_loss()` pour les cas spécifiques, soit la réconciliation standard via `reconciliation()`.

## 2. Exécution du script  

Pour exécuter le programme, il faut utiliser **Python 3.12.4** et s’assurer que le module **ete3** est installé :  
 `pip install ete3`  

Les bibliothèques **os**, **sys** et **argparse** sont incluses par défaut dans Python sinon il faut les installer.  

Le script se lance depuis le terminal en indiquant les arbres de gènes et d’espèces, soit sous forme de chaînes **Newick**, soit à partir de fichiers.  

Les options `--verif` et `--loss` permettent respectivement d’afficher l’arbre de gènes original et l’arbre réconcilié intégrant les pertes.  

#### Exemples d’utilisation :

* Avec des chaînes Newick directement 
`python3 reconciliation.py "(((AAA1,BBB1)1,CCC1)2,((CCC2,DDD1)3,DDD2)4)5;" "(((AAA,BBB)6,CCC)7,DDD)8;"`  

* Avec des fichiers
`python3 reconciliation.py gene_tree.nwk species_tree.nwk`

* Avec l’option verif
` python3 reconciliation.py gene_tree.nwk species_tree.nwk --verif` 

* Avec l’option loss
` python3 reconciliation.py gene_tree.nwk species_tree.nwk --loss`

* Avec les deux options
`python3 reconciliation.py gene_tree.nwk species_tree.nwk --verif --loss`


⚠ : Les arbres Newick doivent toujours avoir trois lettres pour les noms d'espèces et des noeuds internes numérotés, comme dans le premier exemple fourni sinon les options loss et vérif ne pourront pas être utilisées. 


## 3. Algorithme de Reconciliation


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
