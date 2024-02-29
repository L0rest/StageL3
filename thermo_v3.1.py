import tkinter as tkt
from tkinter import *
from tkinter import ttk

from matplotlib.backends._backend_tk import NavigationToolbar2Tk
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg)

import sys

sys.path.extend(['../'])
from sympy import exp, nsolve, Symbol
import matplotlib.pyplot as plt

HeaderString = """ Pour un cristal avec une espèce, on entre le fichier """

# =====================================================================
# % Valeurs des énergies de liaison    %%
# % Les énergies sont en eV            %%
# % Les concentrations sont en at ppm  %%
# =====================================================================
# Constantes utiles:
kB = 8.6173303e-5  # en eV/K
energie_reference = {'H': -6.76824288 / 2., 'C': -18.19943489 / 2., 'O': -9.87320241 / 2.,
                     'N': -16.64877294 / 2., 'B': -80.45763600 / 12.}

# Plusieurs choix dans le calcul de la concentration:
concentration = {1: [10 ** (t / 10) for t in range(-60, 0)], 2: [1e-8],
                 3: [1e-8, 5e-8, 1e-7, 5e-7, 1e-6, 5e-6, 1e-5, 5e-5, 1e-4, 3e-4, 5e-4, 7e-4, 9e-4, 1e-3, 3e-3, 5e-3,
                     7e-3, 9e-3, 1e-2, 3e-2, 5e-2, 7e-2, 9e-2, 1e-1],
                 4: [1e-13, 5e-12, 1e-12, 5e-11, 1e-11, 5e-10, 1e-9, 5e-9, 1e-8, 5e-8, 1e-7, 5e-7, 1e-6, 5e-6, 1e-5,
                     5e-5, 1e-4, 5e-4, 1e-3, 5e-3, 1e-2, 5e-2, 1e-1]}

""" Lecture des données d'entrée du code:
    nbY    : Nombre d'atome interstitiels dans les différents amas
    nbV    : Nombre de lacunes dans les différents amas
    E      : énergies de liaison ou DFT
    deg    : Contient le nom de l'amas = la multiplicité des configurations, dégénérescence des amas
    deftype: Contient le nom de l'amas
    oxyde  : Permet de lui dire si c'est un oxyde, et donc on doit lire les charges
    chge   : Contient la charge de l'amas
    """


def lecture(nomfile, debug=False, oxyde=False):  # Lecture des données d'entrée du code
    nbY = []  # Nombre d'atomes interstitiels dans les différents amas
    nbV = []  # Nombre de lacunes
    E = []  # énergies de liaison ou DFT
    deg = []  # Contient le nombre de cas identiques = la multiplicité des configurations, dégénérescence des amas
    deftype = []  # Contient le nom de l'amas
    charge = []  # contient la charge de l'amas

    with open(nomfile, 'r', encoding='utf8') as inp:
        for line in inp:
            data = line.split()
            if debug:  # Attention dans les données initiales nbV et nbY sont interverties dans le fichier data_test !!!!
                nbY.append(int(data[0]))
                nbV.append(int(data[1]))
            else:
                nbV.append(int(data[0]))
                nbY.append(int(data[1]))
            deg.append(int(data[2]))
            E.append(float(data[3]))
            if oxyde:
                charge.append(float(data[4]))
                deftype.append(str(data[5]))
            else:
                deftype.append(str(data[4]))
    print("Il y a", len(E), "défauts")
    return nbY, nbV, deg, E, deftype, charge


# =====================================================================

""" Conversion des données DFT en énergies de liaison """


def convNRJ2bind(nbY, nbV, deg, E, deftype, nat, atom="C", DFT=False, nrjREF=0.0):
    if DFT:  # &Kevin-debut& dans cette version il faut avoir les énergies du bulk, de la mono-lacune, et du soluté seul dans le site le plus stable!! &Kevin-fin&
        bulk = False
        lac = False
        ## &Kevin-debut&
        monosolute = False
        ## &Kevin-fin&
        for i in range(len(E)):
            if (nbY[i] == 0) and (nbV[i] == 0):
                Ebulk = E[i]
                print(f"     nrj du bulk: {Ebulk} eV")
                bulk = True
            ## &Kevin-debut&
            if (nbY[i] == 1) and (nbV[i] == 0) and (monosolute is False):
                indice_monosol = 0
                for j in range(len(E)):
                    if (nbY[j] == 1) and (nbV[j] == 0) and (j != i):
                        if E[j] < E[i]:
                            ## Le monosoluté j est plus stable que le monosoluté i !
                            E_monosol = E[j]
                            indice_monosol = j
                        else:
                            ## Le monosoluté i est plus stble que la monosoluté j !
                            E_monosol = E[i]
                            indice_monosol = i
                print(
                    f"Le monosoluté le plus stable a été trouvé: en position {indice_monosol - 1}, nom : {deftype[indice_monosol]}")
                print(f"     nrj du monosoluté le plus stable: {E_monosol} eV")
                monosolute = True
            ## &Kevin-fin&
            if (nbY[i] == 0) and (nbV[i] == 1):
                print(f"La monolacune a été trouvée: en position {i + 1}")
                Elac = E[i]
                print(f"     nrj de la monolacune: {Elac} eV")
                EFv = E[i] - Ebulk * (nat - nbV[i]) / nat
                print(f"NRJ de formation de la monolacune V{int(nbV[i])}: {EFv:.4f} eV")
                lac = True
        if not bulk:
            print(" attention pas l'énergie du bulk")
        if not lac:
            print(" attention pas l'énergie de la lacune")
        if not monosolute:
            print(" attention pas l'énergie de la lacune")
        print("=" * 40)

        num = []  # numero du défaut dans la liste lue utile pour après
        Etemp = []
        NRJgc = []
        ## &Kevin-debut&
        for i in range(len(E)):  # Seul le bulk ne doit pas être inclus dans la liste des défauts
            ## &Kevin-fin&
            ## &Kevin-debut&
            Eref = energie_reference[atom]
            if not ((nbY[i] == 0) and (nbV[i] == 0)):
                ## &Kevin-debut&
                NRJb = E[i] + (float(nbV[i]) + float(nbY[i]) - 1.0) * Ebulk - nbV[i] * Elac - nbY[i] * E_monosol
                ## &Kevin-fin&
                ## &Kevin-fin&
                Etemp.append(NRJb)
                Egc_t = E[i] - (Ebulk) - (Eref)
                NRJgc.append(Egc_t)
                num.append(i)
        Eb = [Etemp[i] for i in range(len(Etemp))]
        nY = [nbY[num[i]] for i in range(len(Etemp))]
        nV = [nbV[num[i]] for i in range(len(Etemp))]
        g = [deg[num[i]] for i in range(len(Etemp))]
        dtype = [deftype[num[i]] for i in range(len(Etemp))]
        Egc = [NRJgc[i] for i in range(len(Etemp))]
        Edft = [E[num[i]] for i in range(len(Etemp))]
        for i in range(len(Etemp)):
            print(f"nbY = {nY[i]} nbV = {nV[i]} g = {g[i]}, numero {i + 1} ",
                  f"Edft = {Edft[i]}, Eb = {Eb[i]}, Egc = {Egc[i]} ",
                  f"Défaut: {dtype[i]}")
    else:
        print('Ici on lit les énergies de liaison')
        nY = nbY.copy()
        nV = nbV.copy()
        g = deg.copy()
        ## &Kevin-debut&
        EFv = 2.058  # Il faut définir l'énergie de la lacune
        ## &Kevin-fin&
        Eb = E.copy()
        dtype = deftype.copy()
    print("nY, nV à la fin de convNRJtoBIND:", nY, nV)
    return nY, nV, g, Eb, dtype, EFv


###============================================================

""" Calcul des énergies de formation des amas """


def calculpotentielv1(nY, nV, Eb, dtype, EFv, fileout, muMIN=-12, muMAX=-3, atom="H"):
    # EN CANONIQUE : Les concentrations nominales sont fixées
    # fileout = 'Ef_amas_fct_T_fct_muC_2.dat'
    muY = muMIN
    with open(fileout, 'w', encoding='utf8') as output:
        output.write("###################################\n")
        output.write(f'# valeur des énergies de formation des amas en fonction de mu{atom} \n')
        output.write("###################################\n")
        output.write(f"# T     {atom}_o     {atom}_M     {atom}_t     V  ",
                     "   V2(1)     V2(2)     V2(3)     C2(1)       C2(2) ",
                     "      C2(3)       C2(4)       C2(5)       C3(1)     ,"
                     "  C3(2)       C3(3)       C3(4)       C4(1)       C4(2)",
                     "       C4(3)       C4(4)       C5(1)       C5(2)       C5(3)",
                     "       C5(4)       C6(1)       C6(2)       VC(1)       VC(2) ",
                     "      VC(3)       VC(4)       VC3(1)       VC3(2)       VC4(1)",
                     "       VC4(2)       VC5(1)       VC6(1)       V2C(1)       V2C(2) ",
                     "V2C(3)       V2C(4)       V2C(5)       V2C(6)\n")
        while muY < muMAX:
            clustersEf = [-Eb[i] + nV[i] * EFv + nY[i] * (-8.502 - muY) for i in range(len(nY))]
            output.writelines(
                ("{:15.4e} " * (len(clustersEf) + 1)).format(muY, *[clustersEf[i] for i in range(len(clustersEf))]))
            output.writelines("\n")
            muY += 0.01


###============================================================

""" nY: nombre d'atomes interstitiels
    nV: nombre de lacunes
    Ti: température initiale
    Tf: température finale
    dT: le pas de calcul

    EN CANONIQUE: Les concentrations nominales sont fixées
    Cvtot: concentration nominale en lacunes dans le système
    Ytot:  concentration nominale en oxygène dans le système
    Ns: nombre de site substitutionnels dans le système
    Cv = exp(-EFv / (kB * T))
"""


## &Kevin-debtu&
## Il y a maintenant 2 fonctions similaires "calcul_grand_canonique", une qui permet de calculer en isotherme (calcul_grand_canonique_isoth), l'autre en isoconcentration (calcul_grand_canonique)
def calcul_grand_canonique(nY, nV, g, Eb, dtype, charge, Ti, Tf, dT, conc, EFv, fileout, ensemble, atom, oxyde):
    global figs
    figs = []  # On réinitialise la liste des graphes, au cas où on lance les tracés plusieurs fois de suite!
    Yc = Symbol('Yc')
    Ytot = Symbol('Ytot')
    symvar = [Yc]
    if ensemble == "GCvy":
        Cvtot = Symbol('Cvtot')
        symvar.append(Cvtot)
    if oxyde:
        Ccy = Symbol('Ccy')
        symvar.append(Ccy)
        if ensemble == "GCvy":
            Ccv = Symbol('Ccv')
            symvar.append(Ccv)

    with open(fileout, 'w', encoding='utf8') as output:
        i = 0  # Compteur
        solvables = 0
        for Ytot in conc:  # Boucle sur les concentrations
            T = Ti  # température initiale en kelvin
            print('\n# FLAG = {} \n'.format(i))
            print(f"# Etude en fonction de ({atom}) de la concentration des amas V{atom}n C{atom}tot = {Ytot} \n")
            output.write('\n# FLAG = {} \n'.format(i))
            output.write('# CCtot = {} \n'.format(Ytot))
            output.write("### T(K)")
            output.write((" {} " * (len(dtype))).format(*[dtype[i] for i in range(len(dtype))]))
            output.write("\n")
            Temp = []
            amas = []
            amasnum = []
            while T < Tf:  # Boucle sur les T
                Cv = exp(- EFv / (kB * T))
                ## &Kevin-debut&
                clusters = [g[i] * exp(- Eb[i] / (kB * T)) * (Yc ** nY[i]) * (Cv ** nV[i]) for i in range(len(nY))]
                ##clusters.append(Cv)
                ## &Kevin-fin&

                # équations des concentrations nominales en fonction des concentrations des différents amas: eq1 pour oxy et eq2 pour lac
                # en canonique (concentrations fixées), bizarre dans le cas de la lacune!
                eqY = sum(nY[i] * clusters[i] for i in range(len(nY))) - Ytot
                # print("clusters:", clusters)
                print("nY:", nY, ", Yc:", Yc, ", Ytot:", Ytot)
                print("eqY:", eqY)
                eqs, v1, v2 = eqY, Yc, Ytot  # "GCy": # on ne controle que la concentration en interstitiel: Y (canonique)
                if ensemble == "GCvy":  # "GCvy": on contrôle la concentration en interstitiel et en lacunaire: Y et V (grand-canonique)
                    eqV = sum(nV[i] * clusters[i] for i in range(len(nV))) - Cvtot
                    print("eqV:", eqV, ", Cvtot", Cvtot, ", Cv", Cv)
                    eqs, v1, v2 = [eqY, eqV], [Yc, Cvtot], [0.99 * Ytot, 0.9 * Cv]
                if oxyde:
                    eqCY = sum(nY[i] * charge[i] for i in range(len(nY)))
                    print("eqCY:", eqCY)
                    if ensemble == "GCvy":
                        eqCV = sum(nV[i] * charge[i] for i in range(len(nY)))
                        print("eqCV:", eqCV)
                        eqs, v1, v2 = [eqY, eqV, eqCY, eqCV], [Yc, Cvtot, Ccy, Ccv], [0.99 * Ytot, 0.9 * Cv, 0, 0]
                    else:
                        eqs, v1, v2 = [eqY, eqCY], [Yc, Ccy], [0.99 * Ytot, 0]

                try:
                    values = nsolve(eqs, v1, v2, modules=['mpmath'], prec=9, verify=True)
                    if ensemble == "GCy" and not oxyde:  # Pour résoudre un conflit de typage
                        values = [values]

                    print("values: ", values)
                    # Stockage des concentrations des différents amas contenant de l'oxygène et des lacunes
                    clusternum = clusters
                    for x in range(len(values)):
                        clusternum = [item.subs(symvar[x], values[x]) for item in clusternum]
                    output.writelines(
                        ("{:15.5e} " * (len(clusters) + 1)).format(T, *[clusternum[i] for i in range(len(clusters))]))
                    # output.writelines(("{:15.4e} " * (len(clusters) + 1)).format(T, *[clusternum[i] for i in range(len(clusters))]))
                    output.writelines("\n")
                    Temp.append(T)  # Pour tracer des variables diverses
                    amas.append(clusters)
                    amasnum.append(clusternum)
                    solvables += 1
                except Exception:
                    print("Oooooops... I did it again!")
                T += dT
            i += 1
            print(solvables)
            output.writelines("\n\n")
            if solvables > 0:
                traceCONC(Temp, Ytot, atom, amas, amasnum, dtype)
            else:
                print(
                    "Impossible de tracer le graphe pour cette valeur de concentration, toutes les valeurs de températures ont provoqué un échec du solveur!")


def calcul_grand_canonique_isoth(nY, nV, g, Eb, dtype, charge, Ti, Tf, dT, EFv, fileout, ensemble, atom, oxyde):
    ## Cette fonction fait l'inverse de la première : elle boucle d'abord sur les températures, puis sur les concentrations
    global figs
    figs = []  # On réinitialise la liste des graphes, au cas où on lance les tracés plusieurs fois de suite!
    Yc = Symbol('Yc')
    Ytot = Symbol('Ytot')
    symvar = [Yc]
    if ensemble == "GCvy":
        Cvtot = Symbol('Cvtot')
        symvar.append(Cvtot)
    if oxyde:
        Ccy = Symbol('Ccy')
        symvar.append(Ccy)
        if ensemble == "GCvy":
            Ccv = Symbol('Ccv')
            symvar.append(Ccv)

    with open(fileout, 'w', encoding='utf8') as output:
        i = 0  # Compteur
        solvables = 0
        Temp = list(range(Ti, Tf, dT))
        for T in Temp:  # Boucle sur les températures
            print('\n# FLAG = {} \n'.format(i))
            print(f"# Etude en fonction de ({atom}) de la concentration des amas V{atom}n C{atom}tot = {Ytot} \n")
            output.write('\n# FLAG = {} \n'.format(i))
            output.write('# T = {} K \n'.format(T))
            output.write("### CCtot")
            output.write((" {} " * (len(dtype))).format(*[dtype[i] for i in range(len(dtype))]))
            output.write("\n")
            amas = []
            amasnum = []
            for Ytot in concentration[4]:  # Boucle sur les concentrations (valeurs linéaires)
                Cv = exp(- EFv / (kB * T))
                ## &Kevin-debut&
                clusters = [g[i] * exp(- Eb[i] / (kB * T)) * (Yc ** nY[i]) * (Cv ** nV[i]) for i in range(len(nY))]
                ##clusters.append(Cv)
                ## &Kevin-fin&

                # équations des concentrations nominales en fonction des concentrations des différents amas: eq1 pour oxy et eq2 pour lac
                # en canonique (concentrations fixées), bizarre dans le cas de la lacune!
                eqY = sum(nY[i] * clusters[i] for i in range(len(nY))) - Ytot
                # print("clusters:", clusters)
                print("nY:", nY, ", Yc:", Yc, ", Ytot:", Ytot)
                print("eqY:", eqY)
                eqs, v1, v2 = eqY, Yc, Ytot  # "GCy": # on ne controle que la concentration en interstitiel: Y (canonique)
                if ensemble == "GCvy":  # "GCvy": on contrôle la concentration en interstitiel et en lacunaire: Y et V (grand-canonique)
                    eqV = sum(nV[i] * clusters[i] for i in range(len(nV))) - Cvtot
                    print("eqV:", eqV, ", Cvtot", Cvtot, ", Cv", Cv)
                    eqs, v1, v2 = [eqY, eqV], [Yc, Cvtot], [0.99 * Ytot, 0.9 * Cv]
                if oxyde:
                    eqCY = sum(nY[i] * charge[i] for i in range(len(nY)))
                    print("eqCY:", eqCY)
                    if ensemble == "GCvy":
                        eqCV = sum(nV[i] * charge[i] for i in range(len(nY)))
                        print("eqCV:", eqCV)
                        eqs, v1, v2 = [eqY, eqV, eqCY, eqCV], [Yc, Cvtot, Ccy, Ccv], [0.99 * Ytot, 0.9 * Cv, 0, 0]
                    else:
                        eqs, v1, v2 = [eqY, eqCY], [Yc, Ccy], [0.99 * Ytot, 0]

                try:
                    values = nsolve(eqs, v1, v2, modules=['mpmath'], prec=9, verify=True)
                    if ensemble == "GCy" and not oxyde:  # Pour résoudre un conflit de typage
                        values = [values]

                    print("values: ", values)
                    # Stockage des concentrations des différents amas contenant de l'oxygène et des lacunes
                    clusternum = clusters
                    for x in range(len(values)):
                        clusternum = [item.subs(symvar[x], values[x]) for item in clusternum]
                    output.writelines(("{:15.5e} " * (len(clusters) + 1)).format(Ytot, *[clusternum[i] for i in
                                                                                         range(len(clusters))]))
                    # output.writelines(("{:15.4e} " * (len(clusters) + 1)).format(T, *[clusternum[i] for i in range(len(clusters))]))
                    output.writelines("\n")
                    amas.append(clusters)
                    amasnum.append(clusternum)
                    solvables += 1
                except Exception:
                    print("Oooooops... I did it again!")
                i += 1
            print(solvables)
            output.writelines("\n\n")
            if solvables > 0:
                traceCONC_isoth(T, concentration[4], atom, amas, amasnum, dtype)
            else:
                print(
                    "Impossible de tracer le graphe pour cette valeur de concentration, toutes les valeurs de températures ont provoqué un échec du solveur!")


###============================================================

""" Tracé des concentrations en fonction de la temperature: """


def traceCONC(T, Ytot, atom, clusters, clusternum, dtype):
    global figs
    fig1, ax1 = plt.subplots()
    lacune = ["LAC", "lac", "VAC", "vac", "Vac", "Lac", "1NN", "1nn", "2NN", "2nn"]

    for i in range(len(clusters[0])):
        data = [clusternum[j][i] for j in range(len(T))]
        if dtype[i] in lacune:
            ax1.plot(T, data, '+--', lw=1, label=dtype[i])
        elif atom + "o" in dtype[i]:
            ax1.plot(T, data, 'h--', lw=2, label=dtype[i])
        elif atom + "t" in dtype[i]:
            ax1.plot(T, data, 'v--', lw=2, label=dtype[i])
        elif "V2" in dtype[i]:
            ax1.plot(T, data, 's-', lw=2, label=dtype[i])
        ## &Kevin-debut&
        elif "V3" in dtype[i]:
            ax1.plot(T, data, '^-', lw=2, label=dtype[i])
        ## &Kevin-fin&
        elif "V" + atom in dtype[i]:
            ax1.plot(T, data, 'd-', lw=2, label=dtype[i])
        else:
            ax1.plot(T, data, 'x-', lw=2, label=dtype[i])

    plt.yscale('log')
    plt.gca().legend(loc='lower left', bbox_to_anchor=(0, 0.5), ncol=3, fontsize=8)
    plt.title("Total [" + atom + "] = " + str(Ytot), loc='left')
    ax1.set_xlabel('Temperature (in K)')
    ax1.set_ylabel('Concentrations')
    ## &Kevin-debut&
    plt.ylim([1e-13, 1])
    plt.xlim([T[0], T[-1]])
    ## &Kevin-fin&
    figs.append([fig1, atom, Ytot])  # Une astuce pour faire passer atom et Ytot dans les fonctions trace() et record()!


# ========================================================================
## &Kevin-debut&
def traceCONC_isoth(T, Ytot, atom, clusters, clusternum, dtype):
    global figs
    fig1, ax1 = plt.subplots()

    lacune = ["LAC", "lac", "VAC", "vac", "Vac", "Lac", "1NN", "1nn", "2NN", "2nn"]

    for i in range(len(clusters[0])):
        data = [clusternum[j][i] for j in range(len(Ytot))]
        if dtype[i] in lacune:
            ax1.plot(Ytot, data, '+--', lw=1, label=dtype[i])
        elif atom + "o" in dtype[i]:
            ax1.plot(Ytot, data, 'h--', lw=2, label=dtype[i])
        elif atom + "t" in dtype[i]:
            ax1.plot(Ytot, data, 'v--', lw=2, label=dtype[i])
        elif "V2" in dtype[i]:
            ax1.plot(Ytot, data, 's-', lw=2, label=dtype[i])
        ## &Kevin-debut&
        elif "V3" in dtype[i]:
            ax1.plot(Ytot, data, '^-', lw=2, label=dtype[i])
        ## &Kevin-fin&
        elif "V" + atom in dtype[i]:
            ax1.plot(Ytot, data, 'd-', lw=2, label=dtype[i])
        else:
            ax1.plot(Ytot, data, 'x-', lw=2, label=dtype[i])
    plt.xscale('log')
    plt.yscale('log')
    plt.gca().legend(loc='lower left', bbox_to_anchor=(0, 0.5), ncol=3, fontsize=8)
    plt.title("T = " + str(T) + " K", loc='left')
    ax1.set_xlabel('Concentration totale [' + atom + ']')
    ax1.set_ylabel('Concentrations des amas')
    plt.xlim([1e-13, 1e-1])
    plt.ylim([1e-13, 1e-1])
    plt.grid(linestyle='--')
    figs.append([fig1, atom, T])  # Une astuce pour faire passer atom et T dans les fonctions trace() et record()!


## &Kevin-fin&

""" Ti: température initiale
    Tf: température finale
    dT: le pas de calcul
    nx, ny, nz: taille de la boîte, en nombre d'atomes
    nat: nombre d'atomes dans la boîte de simulation
"""


# ========================================================================

def trace():  # Affichage des fonctions
    ## &Kevin-debut&
    ## Changement de tous les incréments de data, du fait de l'ajout des conditions v7 et vg
    global select
    data = recuperation()  # Une ruse pour transmettre plusieurs données de fonction en fonction sans provoquer de bug!
    mini = data[0]
    maxi = data[1] + 1  # Pour ne pas rater le dernier point sur le graphe!
    pas = data[2]
    x = data[3]
    y = data[4]
    z = data[5]
    nb_atome_maille = data[6]
    atome = data[7]
    canonique = data[8]
    cristal = data[9]
    cas = data[10]
    oxyde_flag = data[11]
    DFTnrj = data[12]
    isoT_flag = data[13]
    nbY = data[14]
    nbV = data[15]
    deg = data[16]
    E = data[17]
    deftype = data[18]
    charge = data[19]
    debug = data[20]
    print("nbY:", nbY)
    print("nbV:", nbV)
    print("degré:", nbY)
    print("E:", E)
    print("deftype:", deftype)
    ## &Kevin-fin&
    ## &Kevin-debut&
    nat = x * y * z * nb_atome_maille  # Maintenant, c'est adaptatif au nombre d'atomes par maille que l'on entre
    ## &Kevin-fin&
    # Attention, il ne faut pas oublier d'inverser la lecture des nB et nV
    if debug:
        DFTnrj = False
        namefile = "data_test"
        nameoutput = "C_amas_diff_T_C_canonique_C_canonique_V_gc_avc_V2.dat"
        EoREF = 0.
    else:
        namefile = "data_" + cristal + "_" + atome
        nameoutput = atome + "_test.dat"
        EoREF = energie_reference[atome]  # énergie de l'atome donné
    conc = concentration[cas]  # 1 = final, 2 = test, 3 = fixe

    nbY, nbV, deg, E, deftype, charge = lecture(namefile, debug, oxyde=False)
    nY, nV, g, Eb, dtype, EFv = convNRJ2bind(nbY, nbV, deg, E, deftype, nat, atome, DFTnrj, nrjREF=EoREF)
    print(namefile)
    ## &Kevin-debut&
    if isoT_flag:
        calcul_grand_canonique_isoth(nY, nV, g, Eb, dtype, charge, mini, maxi, pas, EFv, nameoutput, canonique, atome,
                                     oxyde_flag)
        tracer_isoth()
    else:
        calcul_grand_canonique(nY, nV, g, Eb, dtype, charge, mini, maxi, pas, conc, EFv, nameoutput, canonique, atome,
                               oxyde_flag)
        tracer()
    ## &Kevin-fin&


def tracer():  # Pour afficher le graphe sélectionné
    global select
    # Ouvrir une fenêtre pour afficher le graphe
    fenetreGraphe = tkt.Toplevel()
    fenetreGraphe.title("Graphe")
    fenetreGraphe.geometry("800x600")
    canvas = FigureCanvasTkAgg(figs[select % len(figs)][0], master=fenetreGraphe)
    canvas.get_tk_widget().pack()
    canvas.draw()
    toolbar = NavigationToolbar2Tk(canvas, fenetreGraphe)
    toolbar.update()
    canvas.get_tk_widget().pack()
    fenetreGraphe.mainloop()


## &Kevin-debut&
def tracer_isoth():  # Pour afficher le graphe sélectionné
    global select
    # Ouvrir une fenêtre pour afficher le graphe
    fenetreGraphe = tkt.Toplevel()
    fenetreGraphe.title("Graphe Isotherme")
    fenetreGraphe.geometry("800x600")
    canvas = FigureCanvasTkAgg(figs[select % len(figs)][0], master=fenetreGraphe)
    canvas.get_tk_widget().pack()
    canvas.draw()
    toolbar = NavigationToolbar2Tk(canvas, fenetreGraphe)
    toolbar.update()
    canvas.get_tk_widget().pack()
    fenetreGraphe.mainloop()


## &Kevin-fin&

def shift_left():
    global select
    select -= 1
    tracer()


def shift_right():
    global select
    select += 1
    tracer()


def record():  # Pour enregistrer le graphe dans le dossier courant
    global select
    figs[select][0].savefig("fig_" + figs[select][1] + "_" + str(figs[select][2]) + ".pdf", dpi=600)


def record_all():
    global select
    stock = select
    for i in range(len(figs)):
        select = i
        record()
    select = stock


def recuperation():  # Fonction récupérant les valeurs à partir des données saisies
    v1 = int(temperature_initiale.get())
    v2 = int(temperature_finale.get())
    v3 = int(temperature_pas.get())
    v4 = int(dim_x.get())
    v5 = int(dim_y.get())
    v6 = int(dim_z.get())
    ## &Kevin-debut&
    v7 = int(nb_atom_maille.get())
    ## &Kevin-fin&
    va = selection_y.get()
    vb = selection_c.get()
    vc = selection_r.get()
    vd = selection_o.get()
    ve = selection_x.get()
    vf = selection_d.get()
    ##&Kevin-debut&
    vg = selection_type_g.get()
    ##&Kevin-fin&
    debug = False

    if v1 > v2:
        v2, v1 = v1, v2
        tkt.messagebox.showwarning("Températures",
                                   "La température initiale doit être plus élevée que la température initiale !")
        tkt.messagebox.showinfo("Températures", "Les valeurs ont été automatiquement inversées.")
    if v3 <= 0:
        tkt.messagebox.showerror("Températures", "Le pas de calcul doit être strictement positif !")
    ## &Kevin-debut&
    if v7 <= 0:
        tkt.messagebox.showerror("Nombre atomes par maille",
                                 "Le nombre d'atomes par maille doit être strictement positif !")

    if va == "H: Dihydrogène":
        va = "H"
    elif va == "C: Carbone allotrope diamant":
        va = "C"
    elif va == "N: Diazote":
        va = "N"
    elif va == "O: Dioxygène":
        va = "O"
    elif va == "Dihydrogène, version de débugage":
        va = "H"
        debug = True
    else:
        tkt.messagebox.showerror("Atomes", "L'atome sélectionné est incorrect!")

    if vb == "GCy: lacunaire uniquement":
        vb = "GCy"
    elif vb == "GCvy: lacunaire et interstitiel":
        vb = "GCvy"
    else:
        tkt.messagebox.showwarning("Canonique",
                                   "L'ensemble canonique sélectionné est incorrect; le mode lacunaire uniquement est sélectionné par défaut.")

    if vc == "Be: Béryllium":
        vc = "Be"
    ## &Kevin-debut&
    elif vc == "Ti: Titane":
        vc = "Ti"
    ## &Kevin-fin&
    else:
        tkt.messagebox.showwarning("Atomes",
                                   "Le cristal sélectionné est incorrect; le Béryllium est sélectionné par défaut.")

    if vd == "Exponentiel borné":
        vd = 1
    elif vd == "Test unique":
        vd = 2
    elif vd == "Valeurs linéaires":
        vd = 3
    else:
        tkt.messagebox.showerror("Concentrations",
                                 "Le mode de concentration sélectionné est incorrect; les valeurs linéaires sont sélectionnées par défaut.")

    if ve == "Oxydes":
        ve = True
    elif ve == "Non chargés":
        ve = False
    else:
        ve = False
        tkt.messagebox.showwarning("Charges",
                                   "L'option sélectionnée est incorrecte; les charges sont par défaut désactivées.")

    if vf == "DFT":
        vf = True
    elif vf == "Energy Binding":
        vf = False
    else:
        vf = False
        tkt.messagebox.showwarning("énergies de liaison",
                                   "L'option sélectionnée est incorrecte; les énergies de liaison sont par défaut désactivées.")
    ## &Kevin-debut&
    if vg == "Isotherme":
        vg = True
    elif vg == "Iso-concentration":
        vg = False
    else:
        vg = True
        tkt.messagebox.showwarning("Mode de graphiques",
                                   "L'option sélectionnée est incorrecte; les graphes sont par défaut isothermes.")
    ## &Kevin-fin&
    ## &Kevin-debut&
    nbY, nbV, deg, E, deftype, charge = [], [], [], [], [], []
    if va in ["H", "C", "N", "O"] and vb in ["GCy", "GCvy"] and vc == ["Ti"] and vd in [1, 2, 3] and v3 > 0:
        nbY, nbV, deg, E, deftype, charge = lecture("data_Ti_" + va, False, False)

    # return [v1, v2, v3, v4, v5, v6, va, vb, vc, vd, ve, vf, nbY, nbV, deg, E, deftype, charge, debug] ##Ce qu'il y avait au début
    return [v1, v2, v3, v4, v5, v6, v7, va, vb, vc, vd, ve, vf, vg, nbY, nbV, deg, E, deftype, charge, debug]
    ## &Kevin-fin&


def lire_fichier():
    # Ouvrir une fenêtre pour afficher les données du fichier
    fenetreFichier = tkt.Toplevel()
    fenetreFichier.title("Fichier")
    fenetreFichier.geometry("800x600")
    fichier = open("data_Ti_N", "r")
    texte = fichier.read()
    fichier.close()
    texte = texte.replace("\n", "\n\n")

    textZone = tkt.Label(fenetreFichier, font=("Arial", 12), padx=20, pady=20)
    textZone.pack()
    textZone.config(text=texte)


figs = []  # Pour stocker notre série de graphes
select = 0  # Cette variable nous era utile pour faire défiler les graphes!
fenetre = Tk()
fenetre.title("Calcul des énergies")
fenetre.geometry("1300x750")
fenetre.columnconfigure(0, weight=1)
fenetre.rowconfigure(0, weight=1)
fenetre.rowconfigure(1, weight=1)
fenetre.rowconfigure(2, weight=1)

"""y_scrollbar = ttk.Scrollbar(cadre0, orient = VERTICAL, command = my_canvas.yview)
y_scrollbar.pack(side = RIGHT, fill = Y)
my_canvas.configure(yscrollcommand = y_scrollbar.set)
my_canvas.bind("<Configure>", lambda e: my_canvas.config(scrollregion= my_canvas.bbox(ALL)))"""

cadre1 = tkt.LabelFrame(fenetre, text="Paramétrage des variations de température, en Kelvin", padx=20, pady=20, font=("Arial", 10), bd=5, relief="groove")
cadre1.rowconfigure(0, weight=1)
cadre1.rowconfigure(1, weight=1)
cadre1.rowconfigure(2, weight=1)
cadre1.rowconfigure(3, weight=1)
cadre1.columnconfigure(0, weight=1)
cadre1.columnconfigure(1, weight=1)
cadre1.columnconfigure(2, weight=1)
cadre1.columnconfigure(3, weight=1)
cadre1.columnconfigure(4, weight=1)
cadre1.grid(row=0, column=0, sticky=NSEW)

temperature_initiale = tkt.Entry(cadre1, width=5, font=("Arial", 12), justify="center")
temperature_initiale.insert(0, "500")
titre_initiale = tkt.Label(cadre1, text="Entrez la température initiale :", font=("Helvetica", 11))
temperature_finale = tkt.Entry(cadre1, width=5, font=("Arial", 12), justify="center")
temperature_finale.insert(0, "1000")
titre_finale = tkt.Label(cadre1, text="Entrez la température maximale :", font=("Arial", 11))
temperature_pas = tkt.Entry(cadre1, width=5, font=("Arial", 12), justify="center")
temperature_pas.insert(0, "100")
titre_pas = tkt.Label(cadre1, text="Entrez le pas de calcul :", font=("Arial", 11))
stock_dim = tkt.Canvas(cadre1)
dim_x = tkt.Entry(stock_dim, width=2, font=("Arial", 12), justify="center")
dim_y = tkt.Entry(stock_dim, width=2, font=("Arial", 12), justify="center")
dim_z = tkt.Entry(stock_dim, width=2, font=("Arial", 12), justify="center")
dim_x.insert(0, "1")
dim_y.insert(0, "1")
dim_z.insert(0, "1")
titre_dim = tkt.Label(cadre1, text="Entrez les dimensions x, y et z du cristal :", font=("Arial", 11))
## &Kevin-debut&
nb_atom_maille = tkt.Entry(cadre1, width=2, font=("Arial", 12), justify="center")
titre_atom = tkt.Label(cadre1, text="Entrez le nombre d'atomes par maille :", font=("Arial", 11))
nb_atom_maille.insert(0, "1")
## &Kevin-fin&

cadre2 = tkt.LabelFrame(fenetre, text="Sélection des atomes", padx=20, pady=20, font=("Arial", 10), bd=5, relief="groove")
cadre2.grid(row=1, column=0, sticky="nsew")

selection_atome = tkt.Label(cadre2, text="Sélectionnez la molécule à faire pénétrer dans le cristal :")
atomes_y = ['H: Dihydrogène', 'C: Carbone allotrope diamant', 'N: Diazote', 'O: Dioxygène',
            'Dihydrogène, version de débugage']
selection_y = ttk.Combobox(cadre2, values=atomes_y)
selection_y.current(1)

selection_cristal = tkt.Label(cadre2, text="Sélectionnez le cristal :")
## &Kevin-debut&
atomes_c = ['Be: Béryllium', 'Ti: Titane']
## &Kevin-fin&
selection_r = ttk.Combobox(cadre2, values=atomes_c)
selection_r.current(1)

selection_canonique = tkt.Label(cadre2, text="Sélectionnez l'ensemble canonique :")
canoniques = ['GCy: lacunaire uniquement', 'GCvy: lacunaire et interstitiel']
selection_c = ttk.Combobox(cadre2, values=canoniques)
selection_c.current(1)

selection_oxyde = tkt.Label(cadre2, text="Sélectionnez la charge des atomes :")
oxydes = ['Oxydes', 'Non chargés']
selection_x = ttk.Combobox(cadre2, values=oxydes)
selection_x.current(1)

selection_conc = tkt.Label(cadre2, text="Sélectionnez l'intervalle des valeurs de concentration :")
concentrations = ['Exponentiel borné', 'Test unique', 'Valeurs linéaires']
selection_o = ttk.Combobox(cadre2, values=concentrations)
selection_o.current(1)

selection_dft = tkt.Label(cadre2, text="Sélectionnez le mode de traitement des énergies de liaison :")
debugs = ['DFT', 'Energy Binding']
selection_d = ttk.Combobox(cadre2, values=debugs)
selection_d.current(1)

## &Kevin-debut&
selection_type_graphe = tkt.Label(cadre2, text="Sélectionnez le mode de représentation des calculs :")
types_g = ['Isotherme', 'Iso-concentration']
selection_type_g = ttk.Combobox(cadre2, values=types_g)
selection_type_g.current(1)
## &Kevin-fin&


etape1 = tkt.Button(cadre2, text="Vérifier les données saisies", command=lire_fichier)
etape2 = tkt.Button(cadre2, text="Tracer le graphe!", command=trace)
arrow1 = tkt.Button(cadre2, text="Graphe précédent", command=shift_left)
arrow2 = tkt.Button(cadre2, text="Graphe suivant", command=shift_right)
caught = tkt.Button(cadre2, text="Enregistrer le graphe dans le dossier courant", command=record)
cauall = tkt.Button(cadre2, text="Enregistrer tous les graphes", command=record_all)

tot_label = tkt.Label(cadre2, text=' ')

titre_initiale.grid(row=0, column=0, sticky=S)
temperature_initiale.grid(row=1, column=0, sticky=N)
titre_finale.grid(row=0, column=2, sticky=S)
temperature_finale.grid(row=1, column=2, sticky=N)
titre_pas.grid(row=0, column=4, sticky=S)
temperature_pas.grid(row=1, column=4, sticky=N)
titre_dim.grid(row=2, column=1, sticky=S)
dim_x.grid(row=0, column=0, sticky=NSEW)
dim_y.grid(row=0, column=1, sticky=NSEW)
dim_z.grid(row=0, column=2, sticky=NSEW)
stock_dim.grid(row=3, column=1, sticky=N)
## &Kevin-debut&
titre_atom.grid(row=2, column=3, sticky=S)
nb_atom_maille.grid(row=3, column=3, sticky=N)
## &Kevin-fin&

selection_atome.grid(row=1, column=1)
selection_y.grid(row=1, column=2)
selection_oxyde.grid(row=1, column=3)
selection_x.grid(row=1, column=4)
selection_cristal.grid(row=2, column=1)
selection_r.grid(row=2, column=2)
selection_canonique.grid(row=2, column=3)
selection_c.grid(row=2, column=4)
selection_conc.grid(row=3, column=1)
selection_o.grid(row=3, column=2)
selection_dft.grid(row=3, column=3)
selection_d.grid(row=3, column=4)
etape1.grid(row=4, column=1)
etape2.grid(row=4, column=2)
arrow1.grid(row=5, column=1)
arrow2.grid(row=5, column=2)
## &Kevin-debut&
selection_type_graphe.grid(row=5, column=3)
selection_type_g.grid(row=5, column=4)
## &Kevin-fin&
caught.grid(row=6, column=1)
cauall.grid(row=6, column=2)
## &Kevin-debut&
tot_label.grid(row=6, column=3)
## &Kevin-fin&

texte0 = tkt.Label(fenetre,
                   text="Version: 0.3.1 du 14/12/2023 par Kevin Gautier \n Versions précédentes: 0.2 du 06/05/2023 par Gabriel Faraut \n 0.1 du 25/06/2021 par Damien Connétable \n E_mail: damien.connetable@ensiacet.fr \n CNRS CIRIMAT")
texte0.grid(row=2, column=0, sticky="nsew")

fenetre.resizable(width=True, height=True)
fenetre.mainloop()  # Tout
