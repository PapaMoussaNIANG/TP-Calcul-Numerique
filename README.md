# Équation de Poisson 1D : solveurs directs et itératifs en C

Travail pratique de calcul numérique, Master 1 Calcul Haute Performance et Simulation, Université Paris-Saclay, 2025-2026. Projet individuel. Le squelette du TP, c'est-à-dire le Makefile, les en-têtes et les fonctions d'écriture des résultats, a été fourni par l'enseignant.

## Le problème

Le code résout l'équation de la chaleur stationnaire en une dimension, avec des températures imposées aux deux bords. Une discrétisation par différences finies sur une grille uniforme donne un système linéaire tridiagonal. La matrice est stockée au format bande général de LAPACK, en priorité colonne.

La solution analytique est connue. Chaque solveur est donc jugé sur son erreur relative avant par rapport à cette solution.

## Solveurs directs

Programme `tpPoisson1D_direct`. L'argument choisit la méthode :

| Argument | Méthode |
|---|---|
| `0` | factorisation LU de LAPACK `dgbtrf`, puis résolution `dgbtrs` |
| `1` | factorisation LU tridiagonale écrite à la main, `dgbtrftridiag`, en O(n), puis `dgbtrs` |
| `2` | `dgbsv`, qui enchaîne factorisation et résolution |

## Solveurs itératifs

Programme `tpPoisson1D_iter`. L'argument choisit la méthode :

| Argument | Méthode |
|---|---|
| `0` | Richardson avec le pas optimal α = 2 / (λmin + λmax), les valeurs propres de la matrice étant calculées par leur formule analytique |
| `1` | Richardson préconditionné par Jacobi |
| `2` | Richardson préconditionné par Gauss-Seidel |

Les produits matrice-vecteur passent par BLAS. L'arrêt se fait quand la norme du résidu passe sous 10⁻³, ou après 1000 itérations. L'historique du résidu est écrit dans `RESVEC.dat`.

## Compiler et lancer

Il faut gcc, make, BLAS et LAPACK. Sous Ubuntu :

```sh
sudo apt-get install gcc make libblas-dev liblapacke-dev
```

Le Makefile lit les options de compilation dans un fichier `<nom-de-machine>.mk`. Le fichier `shams.mk` est fourni ; on peut le copier pour une autre machine ou le désigner directement :

```sh
mkdir -p bin
make HOSTNAME=shams all

./bin/tpPoisson1D_direct 1   # LU tridiagonale
./bin/tpPoisson1D_iter 0     # Richardson, pas optimal
```

La taille du problème se règle avec la variable `nbpoints` dans le `main` de chaque programme. Chaque exécution écrit le second membre, la solution exacte, la grille et la solution calculée dans des fichiers `.dat`.

Un Dockerfile est aussi fourni :

```sh
docker build -f docker/Dockerfile -t poisson1d .
```

## Rapport

`rapport.pdf` présente l'analyse et les résultats obtenus.

## Organisation

```
src/        sources C : bibliothèque Poisson 1D, méthodes de Richardson, programmes principaux
include/    en-têtes
docker/     Dockerfile
Makefile    cibles de compilation et d'exécution
shams.mk    options de compilation pour une machine
```
