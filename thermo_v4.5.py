import sys
import tkinter as tkt
from tkinter import ttk, messagebox, NSEW, N, S, NS, END, filedialog

import matplotlib.pyplot as plt
from matplotlib import ticker
from matplotlib.backends._backend_tk import NavigationToolbar2Tk
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg)
from sympy import exp, nsolve, Symbol

sys.path.extend(['../'])

HeaderString = """ Pour un cristal avec une espèce, on entre le fichier """

# =====================================================================
# % Valeurs des énergies de liaison    %%
# % Les énergies sont en eV            %%
# % Les concentrations sont en at ppm  %%
# =====================================================================
# Constantes utiles:
kB = 8.6173303e-5  # en eV/K
energie_reference = {'H': -6.76824288 / 2., 'C': -18.19943489 / 2., 'O': -9.87320241 / 2., 'N': -16.64877294 / 2.,
                     'B': -80.45763600 / 12.}

# Plusieurs choix dans le calcul de la concentration:
concentration = {1: [10 ** (t / 10) for t in range(-60, 0)], 2: [1e-8],
                 3: [1e-8, 5e-8, 1e-7, 5e-7, 1e-6, 5e-6, 1e-5, 5e-5, 1e-4, 3e-4, 5e-4, 7e-4, 9e-4, 1e-3, 3e-3, 5e-3,
                     7e-3, 9e-3, 1e-2, 3e-2, 5e-2, 7e-2, 9e-2, 1e-1],
                 4: [1e-13, 5e-12, 1e-12, 5e-11, 1e-11, 5e-10, 1e-9, 5e-9, 1e-8, 5e-8, 1e-7, 5e-7, 1e-6, 5e-6, 1e-5,
                     5e-5, 1e-4, 5e-4, 1e-3, 5e-3, 1e-2, 5e-2, 1e-1]}

"""
Lecture des données d'entrée du code:
nbY    : Nombre d'atome interstitiels dans les différents amas
nbV    : Nombre de lacunes dans les différents amas
E      : énergies de liaison ou DFT
deg    : Contient le nom de l'amas = la multiplicité des
         configurations, dégénérescence des amas
deftype: Contient le nom de l'amas
oxyde  : Permet de lui dire si c'est un oxyde, et donc on doit lire
         les charges
chge   : Contient la charge de l'amas
"""


def lecture(nomfile, debug=False, oxyde=False):  # Lecture des données d'entrée du code
    nbY = []  # Nombre d'atomes interstitiels dans les différents amas
    nbV = []  # Nombre de lacunes
    E = []  # Énergies de liaison ou DFT
    deg = []  # Le nombre de cas identiques = multiplicité des configurations, dégénérescence des amas
    deftype = []  # Contient le nom de l'amas
    charge = []  # contient la charge de l'amas

    try:
        with open(nomfile, 'r', encoding='utf8') as inp:
            for line in inp:
                if not (line.replace(' ', '') == '\n' or 'flag' in line):
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
    except Exception:
        messagebox.showerror(
            "Impossible de tracer le(s) graphe(s) demandé(s)!, modifiez les paramètres d'entrée et/ou le fichier sélectionné")


# =====================================================================

""" Conversion des données DFT en énergies de liaison """


def convNRJ2bind(nbY, nbV, deg, E, deftype, nat, atom="C", DFT=False, nrjREF=0.0):
    if DFT:  # &Kevin-debut& dans cette version il faut avoir les énergies du bulk, de la mono-lacune, et du soluté seul dans le site le plus stable!! &Kevin-fin&
        bulk = False
        lac = False
        monosolute = False
        EFv = 0.0

        for i in range(len(E)):
            if (nbY[i] == 0) and (nbV[i] == 0):
                Ebulk = E[i]
                print(f"     nrj du bulk: {Ebulk} eV")
                bulk = True

            if (nbY[i] == 1) and (nbV[i] == 0) and (monosolute is False):
                indice_monosol = 0
                for j in range(len(E)):
                    if (nbY[j] == 1) and (nbV[j] == 0) and (j != i):
                        if E[j] < E[i]:
                            # Le monosoluté j est plus stable que le monosoluté i !
                            E_monosol = E[j]
                            indice_monosol = j
                        else:
                            # Le monosoluté i est plus stble que la monosoluté j !
                            E_monosol = E[i]
                            indice_monosol = i
                print(
                    f"Le monosoluté le plus stable a été trouvé: en position {indice_monosol - 1}, nom : {deftype[indice_monosol]}")
                print(f"     nrj du monosoluté le plus stable: {E_monosol} eV")
                monosolute = True

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

        for i in range(len(E)):  # Seul le bulk ne doit pas être inclus dans la liste des défauts

            Eref = energie_reference[atom]
            if not ((nbY[i] == 0) and (nbV[i] == 0)):
                NRJb = E[i] + (float(nbV[i]) + float(nbY[i]) - 1.0) * Ebulk - nbV[i] * Elac - nbY[i] * E_monosol

                Etemp.append(NRJb)
                Egc_t = E[i] - Ebulk - Eref
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
                  f"Edft = {Edft[i]}, Eb = {Eb[i]}, Egc = {Egc[i]} ", f"Défaut: {dtype[i]}")
    else:
        print('Ici on lit les énergies de liaison')
        nY = nbY.copy()
        nV = nbV.copy()
        g = deg.copy()

        EFv = 2.058  # Il faut définir l'énergie de la lacune

        Eb = E.copy()
        dtype = deftype.copy()
    print("nY, nV à la fin de convNRJtoBIND:", nY, nV)
    return nY, nV, g, Eb, dtype, EFv


# ============================================================

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


# ============================================================

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


def fonction_cluster(nat, atome, DFTnrj, debug, EoREF=0.0):
    # recupérer les valeurs de convNRJ2bind après lecture du fichier
    if selection_d.get() == "DFT":
        nbY, nbV, deg, E, deftype, charge = lecture(fichier, debug, oxyde=False)
        nY, nV, g, Eb, dtype, EFv = convNRJ2bind(nbY, nbV, deg, E, deftype, nat, atome, DFTnrj, nrjREF=EoREF)
    else:
        messagebox.showwarning("Attention", "Il faut sélectionner DFT pour pouvoir continuer")
        return

    mu_O = Symbol('mu_O')
    T = Symbol('T')
    Cv = exp(- EFv / (kB * T))

    clusters = [g[i] * exp(- (Eb[i] + nY[i] * mu_O) / (kB * T)) * (Cv ** nV[i]) for i in range(len(nY))]

    return sum(clusters)


# Il y a maintenant 2 fonctions similaires "calcul_grand_canonique", une qui permet de calculer en isotherme
# (calcul_grand_canonique_isoth), l'autre en isoconcentration (calcul_grand_canonique)
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

                clusters = [g[i] * exp(- Eb[i] / (kB * T)) * (Yc ** nY[i]) * (Cv ** nV[i]) for i in range(len(nY))]
                # clusters.append(Cv)

                # équations des concentrations nominales en fonction des concentrations des différents amas: eq1 pour oxy et eq2 pour lac
                # en canonique (concentrations fixées), bizarre dans le cas de la lacune!
                eqY = sum(nY[i] * clusters[i] for i in range(len(nY))) - Ytot
                totutile = sum(clusters)  # A CONFIRMER
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
    # Cette fonction fait l'inverse de la première : elle boucle d'abord sur les températures, puis sur les concentrations
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

                clusters = [g[i] * exp(- Eb[i] / (kB * T)) * (Yc ** nY[i]) * (Cv ** nV[i]) for i in range(len(nY))]
                # clusters.append(Cv)

                # équations des concentrations nominales en fonction des concentrations des
                # différents amas: eq1 pour oxy et eq2 pour lac en canonique
                # (concentrations fixées), bizarre dans le cas de la lacune!
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


# ============================================================


def traceCONC(T, Ytot, atom, clusters, clusternum, dtype):
    """
    Tracé des concentrations en fonction de la température
    :param T:
    :param Ytot:
    :param atom:
    :param clusters:
    :param clusternum:
    :param dtype:
    :return:
    """
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

        elif "V3" in dtype[i]:
            ax1.plot(T, data, '^-', lw=2, label=dtype[i])

        elif "V" + atom in dtype[i]:
            ax1.plot(T, data, 'd-', lw=2, label=dtype[i])
        else:
            ax1.plot(T, data, 'x-', lw=2, label=dtype[i])

    plt.yscale('log')
    plt.gca().legend(loc='lower left', bbox_to_anchor=(0, 0.5), ncol=3, fontsize=8)
    plt.title("Total [" + atom + "] = " + str(Ytot), loc='left')
    ax1.set_xlabel('Temperature (in K)')
    ax1.set_ylabel('Concentrations')

    plt.ylim([1e-13, 1])
    plt.xlim([T[0], T[-1]])

    figs.append([fig1, atom, Ytot])  # Une astuce pour faire passer atom et Ytot dans les fonctions trace() et record()!


# ========================================================================

def traceCONC_isoth(T, Ytot, atom, clusters, clusternum, dtype):
    """
    Tracé des concentrations en fonction de la concentration totale
    :param T:
    :param Ytot:
    :param atom:
    :param clusters:
    :param clusternum:
    :param dtype:
    :return:
    """
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

        elif "V3" in dtype[i]:
            ax1.plot(Ytot, data, '^-', lw=2, label=dtype[i])

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


""" Ti: température initiale
    Tf: température finale
    dT: le pas de calcul
    nx, ny, nz: taille de la boîte, en nombre d'atomes
    nat: nombre d'atomes dans la boîte de simulation
"""


# ========================================================================

def trace():
    """
    Créer les graphes à partir des données de sortie
    :return:
    """
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
    cas = data[10]
    oxyde_flag = data[11]
    DFTnrj = data[12]
    isoT_flag = data[13]
    nbY = data[14]
    nbV = data[15]
    E = data[17]
    deftype = data[18]
    debug = data[20]
    print("nbY:", nbY)
    print("nbV:", nbV)
    print("degré:", nbY)
    print("E:", E)
    print("deftype:", deftype)

    nat = x * y * z * nb_atome_maille  # Maintenant, c'est adaptatif au nombre d'atomes par maille que l'on entre

    # Attention, il ne faut pas oublier d'inverser la lecture des nB et nV
    if debug:
        DFTnrj = False
        nameoutput = "C_amas_diff_T_C_canonique_C_canonique_V_gc_avc_V2.dat"
        EoREF = 0.
    else:
        nameoutput = atome + "_test.dat"
        EoREF = energie_reference[atome]  # énergie de l'atome donné
    conc = concentration[cas]  # 1 = final, 2 = test, 3 = fixe

    try:
        nbY, nbV, deg, E, deftype, charge = lecture(fichier, debug, oxyde=False)
        nY, nV, g, Eb, dtype, EFv = convNRJ2bind(nbY, nbV, deg, E, deftype, nat, atome, DFTnrj, nrjREF=EoREF)
    except Exception:
        return

    if isoT_flag:
        calcul_grand_canonique_isoth(nY, nV, g, Eb, dtype, charge, mini, maxi, pas, EFv, nameoutput, canonique, atome,
                                     oxyde_flag)
        tracer_isoth()
    else:
        calcul_grand_canonique(nY, nV, g, Eb, dtype, charge, mini, maxi, pas, conc, EFv, nameoutput, canonique, atome,
                               oxyde_flag)
        tracer()


def tracer():
    """
    Pour afficher le graphe sélectionné
    :return:
    """
    global select

    # Ouvrir une fenêtre pour afficher le graphe
    fenetreGraphe = tkt.Toplevel()
    fenetreGraphe.title("Graphes")
    fenetreGraphe.geometry("900x700")
    canvas = FigureCanvasTkAgg(figs[select % len(figs)][0], master=fenetreGraphe)
    canvas.get_tk_widget().pack(fill=tkt.BOTH, expand=True)
    canvas.draw()
    toolbar = NavigationToolbar2Tk(canvas, fenetreGraphe)
    toolbar.update()
    canvas.get_tk_widget().pack(fill=tkt.BOTH, expand=True)

    # Boutons pour changer de graphe
    def shift_left():
        global select
        select -= 1

        for widget in fenetreGraphe.winfo_children():
            widget.destroy()

        canvas = FigureCanvasTkAgg(figs[select % len(figs)][0], master=fenetreGraphe)
        canvas.get_tk_widget().pack(fill=tkt.BOTH, expand=True)
        canvas.draw()
        toolbar = NavigationToolbar2Tk(canvas, fenetreGraphe)
        toolbar.update()
        canvas.get_tk_widget().pack(fill=tkt.BOTH, expand=True)

        stockBoutons = tkt.Frame(fenetreGraphe, padx=5, pady=5)
        stockBoutons.columnconfigure(0, weight=1)
        stockBoutons.columnconfigure(1, weight=1)
        boutonGauche = tkt.Button(stockBoutons, text="Précédent", command=shift_left, padx=5, pady=5, bd=3, width=10)
        boutonGauche.grid(row=0, column=0)
        boutonDroite = tkt.Button(stockBoutons, text="Suivant", command=shift_right, padx=5, pady=5, bd=3, width=10)
        boutonDroite.grid(row=0, column=1)
        stockBoutons.pack()
        cauall = tkt.Button(fenetreGraphe, text="Enregistrer tous les graphes", command=record_all, padx=5, pady=5,
                            bd=3)
        cauall.pack()

    def shift_right():
        global select
        select += 1

        for widget in fenetreGraphe.winfo_children():
            widget.destroy()

        canvas = FigureCanvasTkAgg(figs[select % len(figs)][0], master=fenetreGraphe)
        canvas.get_tk_widget().pack(fill=tkt.BOTH, expand=True)
        canvas.draw()
        toolbar = NavigationToolbar2Tk(canvas, fenetreGraphe)
        toolbar.update()
        canvas.get_tk_widget().pack(fill=tkt.BOTH, expand=True)

        stockBoutons = tkt.Frame(fenetreGraphe, padx=5, pady=5)
        stockBoutons.columnconfigure(0, weight=1)
        stockBoutons.columnconfigure(1, weight=1)
        boutonGauche = tkt.Button(stockBoutons, text="Précédent", command=shift_left, padx=5, pady=5, bd=3, width=10)
        boutonGauche.grid(row=0, column=0)
        boutonDroite = tkt.Button(stockBoutons, text="Suivant", command=shift_right, padx=5, pady=5, bd=3, width=10)
        boutonDroite.grid(row=0, column=1)
        stockBoutons.pack()
        cauall = tkt.Button(fenetreGraphe, text="Enregistrer tous les graphes", command=record_all, padx=5, pady=5,
                            bd=3)
        cauall.pack()

    stockBoutons = tkt.Frame(fenetreGraphe, padx=5, pady=5)
    stockBoutons.columnconfigure(0, weight=1)
    stockBoutons.columnconfigure(1, weight=1)
    boutonGauche = tkt.Button(stockBoutons, text="Précédent", command=shift_left, padx=5, pady=5, bd=3, width=10)
    boutonGauche.grid(row=0, column=0)
    boutonDroite = tkt.Button(stockBoutons, text="Suivant", command=shift_right, padx=5, pady=5, bd=3, width=10)
    boutonDroite.grid(row=0, column=1)
    stockBoutons.pack()
    cauall = tkt.Button(fenetreGraphe, text="Enregistrer tous les graphes", command=record_all, padx=5, pady=5, bd=3)
    cauall.pack()

    fenetreGraphe.mainloop()


def tracer_isoth():
    """
    Pour afficher le graphe isotherme sélectionné
    :return:
    """
    global select

    # Ouvrir une fenêtre pour afficher le graphe
    fenetreGraphe = tkt.Toplevel()
    fenetreGraphe.title("Graphes isothermes")
    fenetreGraphe.geometry("900x700")
    canvas = FigureCanvasTkAgg(figs[select % len(figs)][0], master=fenetreGraphe)
    canvas.get_tk_widget().pack(fill=tkt.BOTH, expand=True)
    canvas.draw()
    toolbar = NavigationToolbar2Tk(canvas, fenetreGraphe)
    toolbar.update()
    canvas.get_tk_widget().pack(fill=tkt.BOTH, expand=True)

    # Boutons pour changer de graphe
    def shift_left():
        global select
        select -= 1

        for widget in fenetreGraphe.winfo_children():
            widget.destroy()

        canvas = FigureCanvasTkAgg(figs[select % len(figs)][0], master=fenetreGraphe)
        canvas.get_tk_widget().pack(fill=tkt.BOTH, expand=True)
        canvas.draw()
        toolbar = NavigationToolbar2Tk(canvas, fenetreGraphe)
        toolbar.update()
        canvas.get_tk_widget().pack(fill=tkt.BOTH, expand=True)

        stockBoutons = tkt.Frame(fenetreGraphe, padx=5, pady=5)
        stockBoutons.columnconfigure(0, weight=1)
        stockBoutons.columnconfigure(1, weight=1)
        boutonGauche = tkt.Button(stockBoutons, text="Précédent", command=shift_left, padx=5, pady=5, bd=3, width=10)
        boutonGauche.grid(row=0, column=0)
        boutonDroite = tkt.Button(stockBoutons, text="Suivant", command=shift_right, padx=5, pady=5, bd=3, width=10)
        boutonDroite.grid(row=0, column=1)
        stockBoutons.pack()
        cauall = tkt.Button(fenetreGraphe, text="Enregistrer tous les graphes", command=record_all, padx=5, pady=5,
                            bd=3)
        cauall.pack()

    def shift_right():
        global select
        select += 1

        for widget in fenetreGraphe.winfo_children():
            widget.destroy()

        canvas = FigureCanvasTkAgg(figs[select % len(figs)][0], master=fenetreGraphe)
        canvas.get_tk_widget().pack(fill=tkt.BOTH, expand=True)
        canvas.draw()
        toolbar = NavigationToolbar2Tk(canvas, fenetreGraphe)
        toolbar.update()
        canvas.get_tk_widget().pack(fill=tkt.BOTH, expand=True)

        stockBoutons = tkt.Frame(fenetreGraphe, padx=5, pady=5)
        stockBoutons.columnconfigure(0, weight=1)
        stockBoutons.columnconfigure(1, weight=1)
        boutonGauche = tkt.Button(stockBoutons, text="Précédent", command=shift_left, padx=5, pady=5, bd=3, width=10)
        boutonGauche.grid(row=0, column=0)
        boutonDroite = tkt.Button(stockBoutons, text="Suivant", command=shift_right, padx=5, pady=5, bd=3, width=10)
        boutonDroite.grid(row=0, column=1)
        stockBoutons.pack()
        cauall = tkt.Button(fenetreGraphe, text="Enregistrer tous les graphes", command=record_all, padx=5, pady=5,
                            bd=3)
        cauall.pack()

    stockBoutons = tkt.Frame(fenetreGraphe, padx=5, pady=5)
    stockBoutons.columnconfigure(0, weight=1)
    stockBoutons.columnconfigure(1, weight=1)
    boutonGauche = tkt.Button(stockBoutons, text="Précédent", command=shift_left, padx=5, pady=5, bd=3, width=10)
    boutonGauche.grid(row=0, column=0)
    boutonDroite = tkt.Button(stockBoutons, text="Suivant", command=shift_right, padx=5, pady=5, bd=3, width=10)
    boutonDroite.grid(row=0, column=1)
    stockBoutons.pack()
    cauall = tkt.Button(fenetreGraphe, text="Enregistrer tous les graphes", command=record_all, padx=5, pady=5, bd=3)
    cauall.pack()

    fenetreGraphe.mainloop()


def record():
    """
    Enregistre le graphe sélectionné dans le dossier courant
    :return:
    """
    global select

    figs[select][0].savefig("fig_" + figs[select][1] + "_" + str(figs[select][2]) + ".pdf", dpi=600)


def record_all():
    """
    Enregistre tous les graphes dans le dossier courant
    :return:
    """
    try:
        global select
        stock = select
        for i in range(len(figs)):
            select = i
            record()
        select = stock
        messagebox.showinfo("Enregistrement",
                            "Tous les graphes ont été enregistrés dans le dossier courant avec succès")
    except Exception:
        messagebox.showerror("Enregistrement", "Une erreur est survenue lors de l'enregistrement des graphes")


def recuperation():
    """
    Récupère les données entrées par l'utilisateur
    :return:
    """

    v1 = int(temperature_initiale.get())
    v2 = int(temperature_finale.get())
    v3 = int(temperature_pas.get())
    v4 = int(dim_x.get())
    v5 = int(dim_y.get())
    v6 = int(dim_z.get())

    v7 = int(nb_atom_maille.get())

    va = selection_y.get()
    vb = selection_c.get()
    vc = selection_r.get()
    vd = selection_o.get()
    ve = selection_x.get()
    vf = selection_d.get()
    vg = selection_type_g.get()
    debug = False

    if v1 > v2:
        v2, v1 = v1, v2
        messagebox.showwarning("Températures",
                               "La température initiale doit être plus élevée que la température initiale !")
        messagebox.showinfo("Températures", "Les valeurs ont été automatiquement inversées.")
    if v3 <= 0:
        messagebox.showerror("Températures", "Le pas de calcul doit être strictement positif !")

    if v7 <= 0:
        messagebox.showerror("Nombre atomes par maille",
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
        messagebox.showerror("Atomes", "L'atome sélectionné est incorrect!")

    if vb == "GCy: lacunaire uniquement":
        vb = "GCy"
    elif vb == "GCvy: lacunaire et interstitiel":
        vb = "GCvy"
    else:
        messagebox.showwarning("Canonique",
                               "L'ensemble canonique sélectionné est incorrect; le mode lacunaire uniquement est sélectionné par défaut.")

    if vc == "Be: Béryllium":
        vc = "Be"

    elif vc == "Ti: Titane":
        vc = "Ti"

    else:
        messagebox.showwarning("Atomes",
                               "Le cristal sélectionné est incorrect; le Béryllium est sélectionné par défaut.")

    if vd == "Exponentiel borné":
        vd = 1
    elif vd == "Test unique":
        vd = 2
    elif vd == "Valeurs linéaires":
        vd = 3
    else:
        messagebox.showerror("Concentrations",
                             "Le mode de concentration sélectionné est incorrect; les valeurs linéaires sont sélectionnées par défaut.")

    if ve == "Oxydes":
        ve = True
    elif ve == "Non chargés":
        ve = False
    else:
        ve = False
        messagebox.showwarning("Charges",
                               "L'option sélectionnée est incorrecte; les charges sont par défaut désactivées.")

    if vf == "DFT":
        vf = True
    elif vf == "Energy Binding":
        vf = False
    else:
        vf = False
        messagebox.showwarning("énergies de liaison",
                               "L'option sélectionnée est incorrecte; les énergies de liaison sont par défaut désactivées.")

    if vg == "Isotherme":
        vg = True
    elif vg == "Iso-concentration":
        vg = False
    else:
        vg = True
        messagebox.showwarning("Mode de graphiques",
                               "L'option sélectionnée est incorrecte; les graphes sont par défaut isothermes.")

    nbY, nbV, deg, E, deftype, charge = [], [], [], [], [], []
    if va in ["H", "C", "N", "O"] and vb in ["GCy", "GCvy"] and vc == ["Ti"] and vd in [1, 2, 3] and v3 > 0:
        nbY, nbV, deg, E, deftype, charge = lecture("data_Ti_" + va, False, False)

    # return [v1, v2, v3, v4, v5, v6, va, vb, vc, vd, ve, vf, nbY, nbV, deg, E, deftype, charge, debug] ##Ce qu'il y avait au début
    return [v1, v2, v3, v4, v5, v6, v7, va, vb, vc, vd, ve, vf, vg, nbY, nbV, deg, E, deftype, charge, debug]


def lire_fichier():
    """
    Ouvre une fenêtre pour sélectionner un fichier de données
    :return:
    """
    global select
    global fichier

    data = recuperation()
    atom = data[7]

    try:
        file = open(fichier)
        texte = file.read()
        file.close()
    except (FileNotFoundError, TypeError):
        messagebox.showerror("Fichier", "Le fichier demandé n'existe pas!")
        return

    texte = texte.split("\n")  # On découpe le texte en lignes
    texte = [ligne for ligne in texte if ligne != '']  # On enlève les lignes vides
    texte = [ligne for ligne in texte if 'flag' not in ligne]  # On enlève les lignes de flag
    nb_ligne = len(texte)

    texte = [ligne.replace(" ", ";") for ligne in texte]
    # Si le fichier est bien formaté, chaque ligne doit contenir 4 points virgules (i.e. 5 colonnes)
    if any(";" not in ligne for ligne in texte) or any(";" in ligne for ligne in texte if ligne.count(";") != 4):
        messagebox.showerror("Fichier", "Le fichier demandé n'est pas bien formaté!")
        return

    # Ouvrir une fenêtre pour afficher les données du fichier
    fenetreFichier = tkt.Toplevel(fenetre)
    fenetreFichier.title("Fichier de données")
    style = ttk.Style()
    style.configure("Treeview.Heading", font=('Helvetica', 14))
    style.configure("Treeview", font=('Helvetica', 13))
    tree = ttk.Treeview(fenetreFichier,
                        columns=('Nb lacunes', "Nb atomes" + atom, "Dégénéréscences", "E(DFT) [eV]", "Nom Config"),
                        show='headings')
    tree.heading('Nb lacunes', text='Nb lacunes')
    tree.heading("Nb atomes" + atom, text="Nb atomes " + atom)
    tree.heading('Dégénéréscences', text='Dégénéréscences')
    tree.heading('E(DFT) [eV]', text='E(DFT) [eV]')
    tree.heading('Nom Config', text='Nom Config')
    tree.column('Nb lacunes', anchor='center')
    tree.column("Nb atomes" + atom, anchor='center')
    tree.column('Dégénéréscences', anchor='center')
    tree.column('E(DFT) [eV]', anchor='center')
    tree.column('Nom Config', anchor='center')
    tree.pack(expand=True, fill=tkt.BOTH)

    for i in range(nb_ligne):
        ligne = texte[i].split(";")
        tree.insert("", "end", values=(ligne[0], ligne[1], ligne[2], ligne[3], ligne[4]))

    fenetreFichier.geometry("1000x" + str(25 * nb_ligne))
    fenetreFichier.mainloop()


def select_fichier():
    """
    Ouvre une fenêtre pour sélectionner un fichier de données
    :return:
    """
    global select
    global fichier

    fichier = filedialog.askopenfilename(title="Ouvrir un fichier")

    if fichier:
        file = open(fichier)
        texte = file.read()
        file.close()
        texte = texte.split("\n")

        # Si le fichier contient un flag, on récupère les valeurs pour les mettre dans les champs correspondants
        if any("flag" in ligne for ligne in texte):
            ligne = [l for l in texte if "flag" in l][0].split(";")
            ligne = [(l.split("="))[1].replace(" ", "") for l in ligne]
            t_init, t_fin, t_step, dim_supercell, nb_at, atom, ens_canon, charge, cristal, traitement = ligne

            if int(t_init) < 0 or int(t_fin) < 0 or int(t_step) < 0 or int(dim_supercell[0]) < 0 or int(
                    dim_supercell[1]) < 0 or int(dim_supercell[2]) < 0 or int(nb_at) < 0:
                messagebox.showerror("Flag", "Les valeurs du flag ne peuvent pas être négatives")
                return
            if int(t_init) > int(t_fin):
                messagebox.showerror("Flag", "La température initiale doit être inférieure à la température finale")
                return
            if int(t_step) == 0:
                messagebox.showerror("Flag", "Le pas de calcul ne peut pas être nul")
                return
            if int(dim_supercell[0]) == 0 or int(dim_supercell[1]) == 0 or int(dim_supercell[2]) == 0:
                messagebox.showerror("Flag", "Les dimensions du cristal ne peuvent pas être nulles")
                return
            if int(nb_at) == 0:
                messagebox.showerror("Flag", "Le nombre d'atomes par maille ne peut pas être nul")
                return

            temperature_initiale.delete(0, END)
            temperature_initiale.insert(0, t_init)
            temperature_finale.delete(0, END)
            temperature_finale.insert(0, t_fin)
            temperature_pas.delete(0, END)
            temperature_pas.insert(0, t_step)
            dim_x.delete(0, END)
            dim_y.delete(0, END)
            dim_z.delete(0, END)
            dim_x.insert(0, dim_supercell[0])
            dim_y.insert(0, dim_supercell[1])
            dim_z.insert(0, dim_supercell[2])
            nb_atom_maille.delete(0, END)
            nb_atom_maille.insert(0, nb_at)

            selection_y.delete(0, END)
            i = 0
            while i < len(atomes_y) and atom not in atomes_y[i]:
                i += 1
            selection_y.current(i)

            selection_c.delete(0, END)
            if ens_canon == "GCy":
                selection_c.current(0)
            else:
                selection_c.current(1)

            selection_x.delete(0, END)
            if charge == "#O":
                selection_x.current(1)
            else:
                selection_x.current(0)

            selection_r.delete(0, END)
            if cristal == "Be":
                selection_r.current(0)
            else:
                selection_r.current(1)

            selection_d.delete(0, END)
            if traitement == "DFT":
                selection_d.current(0)
            else:
                selection_d.current(1)

        # Get file name only
        textSelectData.config(text="Fichier sélectionné : " + fichier.split("/")[-1])

        messagebox.showinfo("Fichier", "Le fichier a été sélectionné avec succès")


def select_fichier_solu():
    """
    Ouvre une fenêtre pour sélectionner un ou plusieurs fichiers de données
    :return:
    """
    global fichiersSolu

    fichier = filedialog.askopenfilenames(title="Sélectionner des fichiers")

    textSelectFichiers.config(text="Fichiers sélectionnés : " + str(len(fichier)))

    # On réinitialise la liste des fichiers
    fichiersSolu = []
    listFichier.delete(0, END)

    for f in fichier:
        fichiersSolu.append(f)
        listFichier.insert(END, f.split("/")[-1])

    messagebox.showinfo("Fichiers", "Les fichiers ont été sélectionnés avec succès")


def lire_fichier_solu():
    """
    Ouvre une fenêtre pour lire le fichier de données expérimentales sélectionné
    :return:
    """

    try:
        file = open(listFichier.get(listFichier.curselection()))
        texte = file.read()
        file.close()
    except Exception:
        messagebox.showerror("Fichier", "Erreur lors de la lecture du fichier")
        return

    texte = texte.split("\n")  # On découpe le texte en lignes
    texte = [ligne for ligne in texte if ligne != '']  # On enlève les lignes vides

    # Ouvrir une fenêtre pour afficher les données du fichier
    fenetreFichier = tkt.Toplevel(fenetre)
    fenetreFichier.title("Fichier de données")
    style = ttk.Style()
    style.configure("Treeview.Heading", font=('Helvetica', 14))
    style.configure("Treeview", font=('Helvetica', 13))
    tree = ttk.Treeview(fenetreFichier, columns=("Concentration", "Température"), show='headings')
    tree.heading('Concentration', text='Concentration')
    tree.heading('Température', text='Température (K)')
    tree.column('Concentration', anchor='center')
    tree.column('Température', anchor='center')
    tree.pack(expand=True, fill=tkt.BOTH)

    for i in range(len(texte)):
        ligne = texte[i].split(";")
        tree.insert("", "end", values=(ligne[0], ligne[1]))

    fenetreFichier.geometry("800x600")
    fenetreFichier.mainloop()


def calcul_solubilite():
    """
    Calcul de la solubilité
    :return:
    """
    global fichiersSolu

    if not fichiersSolu:
        messagebox.showerror("Données", "Veuillez sélectionner au moins un fichier de données pour les solutions !")
        return

    if nomEspece.get() == "":
        messagebox.showerror("Données", "Veuillez entrer le nom de l'espèce chimique !")
        return

    if titreGraphique.get() == "":
        messagebox.showerror("Données", "Veuillez entrer le titre du graphique !")
        return

    if fichier_ftotal is None:
        messagebox.showerror("Données", "Veuillez sélectionner un fichier de données pour F[X] !")
        return

    # Calcul récupération data
    data = recuperation()
    nat = data[6]
    atome = data[7]
    DFTnrj = data[12]
    debug = data[20]

    mu_O = Symbol('mu_O')
    T = []
    F = []
    kB = 8.6173303e-5

    temp_min = float('inf')
    temp_max = float('-inf')

    solu = {}
    temp = {}

    file1 = open(fichier_ftotal)
    lines = file1.readlines()
    file1.close()
    for l in lines:
        l = l.split()
        T.append(float(l[0]))
        if float(l[0]) < temp_min:
            temp_min = float(l[0])
        if float(l[0]) > temp_max:
            temp_max = float(l[0])
        F.append(float(l[1]))

    for f in fichiersSolu:
        file = open(f, 'r')
        lines = file.readlines()
        file.close()
        sol = []
        t = []

        for l in lines:
            l = l.split(";")
            sol.append(float(l[0]))
            t.append(float(l[1]))

            if float(l[1]) < temp_min:
                temp_min = float(l[1])
            if float(l[1]) > temp_max:
                temp_max = float(l[1])

        solu[f] = sol
        temp[f] = t

    concentration = []
    mu_sol = [-0.1]
    sum_cluster = fonction_cluster(nat, atome, DFTnrj, debug)

    for i in range(50, 100):
        beta = 1 / (kB * T[i])

        ''' 
        Aide pour le calcul de la solubilité 
        On doit trouver x (qui ici est le potentiel chimique) :
        1 / 2 * beta * (x + F[i]) = somme [ g[i] * exp(- Eb[i] / (kB * T)) * (Yc ** nY[i]) * (Cv ** nV[i]) for i in range(len(nY)) ] (terme entre les crochets est la somme des clusters de la ligne 246 du code)
        A droite de l'égalité, on a pas de x, on ne peut pas résoudre l'équation telle quelle
        On la transforme comme suit :
        1 / 2 * beta * (x + F[i]) = somme [ g[i] * exp(- (Eb[i] + nY[i]*x) / (kB * T)) * (Cv ** nV[i]) for i in range(len(nY)) ] (on peut laisser les lacunes pour le moment)
        (entre les crochets de somme, c'est la variable cluster[i] de la ligne 246 du code)

        Mais sous cette forme, on peut pas utiliser les calculs de la ligne 246 du code. Il faudrait donc modifier la fonction calcul_grand_canonique pour la réécrire sous cette forme.
        C'est peut-être un peu long, donc je te propose de créer une nouvelle fonction, qui fait la même chose que calcul_grand_canonique, mais qui renvoie la somme des clusters (le terme entre les crochets de la ligne 246) pour chaque température écrit de la manière ci-dessus.

        Et après, on pourra résoudre l'équation pour x comme suit (quelque chose comme ça, à vérifier):
        value = nsolve(1 / 2 * beta * (x + F[i]) - sum[ cluster[i] for i in len(cluster) ], -0.1)

        Attention ! Le coefficient 1/2 dépend de la phase ordonnée. Dans notre cas, c'est Ti2N, d'où 1/2 (c'est TixNy, et le ratio y/x)
        Il faudra donc prévoir un endroit pour le modifier (par exemple, dans l'interface graphique, laisser un champ pour donner la steochiométrie de la phase ordonnée)
        '''

        sum_cluster = sum_cluster.subs({Symbol('T'): T[i]})
        value = nsolve(1 / 2 * beta * (mu_O + F[i]) - sum_cluster, -0.1, verify=False)
        # value = nsolve(1 / 2 * beta * (mu_O + F[i]) - sum_cluster, mu_sol[-1], verify=False)

        # On stocke la valeur de value dans mu_sol pour la prochaine itération (à travailler, car on n'obtient pas le
        # même résultat que si on met -0.1 en paramètre)
        mu_sol.append(value)
        concentration.append(1 / 2 * (value + F[i]) * beta / 100)

    fenetreSolu = tkt.Toplevel()
    fenetreSolu.title("Calcul de solubilité")
    fenetreSolu.geometry("900x700")

    fig = plt.figure()
    plot1 = fig.add_subplot(111)
    plot1.plot(concentration, T[50:100], label="DFT parametrization")
    for s in solu:
        plot1.plot(solu[s], temp[s], '.', label=s.split("/")[-1].split(".")[0])
    plot1.set_xlabel(nomEspece.get() + ' concentration')
    plot1.set_ylabel('Temperature [K]')
    plot1.set_title(titreGraphique.get())
    plot1.legend()
    plot1.yaxis.set_major_locator(ticker.MultipleLocator(100))
    plot1.grid(linestyle='--')

    canvas = FigureCanvasTkAgg(fig, master=fenetreSolu)
    canvas.get_tk_widget().pack(fill=tkt.BOTH, expand=True)
    canvas.draw()
    toolbar = NavigationToolbar2Tk(canvas, fenetreSolu)
    toolbar.update()
    canvas.get_tk_widget().pack(fill=tkt.BOTH, expand=True)

    fenetreSolu.mainloop()


def select_fichier_ftotal():
    """
    Ouvre une fenêtre pour sélectionner le fichier ftotal
    :return:
    """
    global fichier_ftotal

    fichier_ftotal = filedialog.askopenfilename(title="Ouvrir un fichier", filetypes=[('Text files', '*.txt')])

    if fichier_ftotal:
        textSelectFichierFtotal.config(text="Fichier sélectionné : " + fichier_ftotal.split("/")[-1])
        messagebox.showinfo("Fichier", "Le fichier a été sélectionné avec succès")
    else:
        messagebox.showerror("Fichier", "Le fichier demandé n'existe pas!")
    return


def lire_fichier_ftotal(v):
    """
    Ouvre une fenêtre pour lire le fichier ftotal sélectionné (affiche les valeurs des énergies libres
    en Ev ou en Kj/mol en fonction de l'unité choisie)
    :return:
    """
    global fichier_ftotal

    print(v)
    try:
        file = open(fichier_ftotal)
        texte = file.read()
        file.close()
    except (FileNotFoundError, TypeError, NameError):
        messagebox.showerror("Fichier", "Le fichier demandé n'existe pas!")
        return

    texte = texte.split("\n")  # On découpe le texte en lignes

    # Ouvrir une fenêtre pour afficher les données du fichier
    fenetreFichier = tkt.Toplevel(fenetre)
    fenetreFichier.title("Fichier de données")
    style = ttk.Style()
    style.configure("Treeview.Heading", font=('Helvetica', 14))
    style.configure("Treeview", font=('Helvetica', 13))
    tree = ttk.Treeview(fenetreFichier, columns=("Température (K)", "Energie libre totale F[X] (eV)"), show='headings')
    tree.heading('Température (K)', text='Température (K)')
    if v == "1":
        tree.heading('Energie libre totale F[X] (eV)', text='Energie libre totale F[X] (eV)')
    else:
        tree.heading('Energie libre totale F[X] (eV)', text='Energie libre totale F[X] (Kj/mol)')
    tree.column('Température (K)', anchor='center')
    tree.column('Energie libre totale F[X] (eV)', anchor='center')
    tree.pack(expand=True, fill=tkt.BOTH)

    for i in range(len(texte)):
        ligne = texte[i].split()
        if v == "1":  # On affiche les valeurs en eV
            tree.insert("", "end", values=(ligne[0], ligne[1]))
        else:  # On affiche les valeurs en Kj/mol
            tree.insert("", "end", values=(ligne[0], float(ligne[1]) * 96.485))

    fenetreFichier.geometry("800x600")
    fenetreFichier.mainloop()


def check_data():
    """
    Vérifie que les données sont bien formatées
    :return:
    """
    global fichier_ftotal
    global fichiersSolu

    if fichier_ftotal is None:
        messagebox.showerror("Données", "Veuillez sélectionner un fichier de données pour F[X] !")
        return

    if len(fichiersSolu) == 0:
        messagebox.showerror("Données", "Veuillez sélectionner au moins un fichier de données pour les solutions !")
        return

    file = open(fichier_ftotal)
    ftot = file.read()
    file.close()

    ftot = [ligne for ligne in ftot.split("\n") if ligne != '']  # On enlève les lignes vides

    if any(len(ligne.split()) != 2 for ligne in ftot):
        messagebox.showerror("Données", "Le fichier F[X] n'est pas bien formaté !")
        return

    for f in fichiersSolu:
        file = open(f, 'r')
        lines = file.readlines()
        file.close()

        lines = [ligne for ligne in lines if ligne != '']  # On enlève les lignes vides

        if any(len(l.split(";")) != 2 for l in lines):
            messagebox.showerror("Données", "Le fichier " + f + " n'est pas bien formaté !")
            return

    messagebox.showinfo("Données", "Les données sont formatées correctement !")


# ===============================================================================
# ============================ FRONT ============================================
# ===============================================================================

fichier = None  # Pour stocker le fichier de données sélectionné
fichiersSolu = []  # Pour stocker les fichiers pour le calcul de solubilité
fichier_ftotal = None  # Pour stocker le fichier ftotal sélectionné
figs = []  # Pour stocker notre série de graphes
select = 0  # Cette variable nous sera utile pour faire défiler les graphes !
fenetre = tkt.Tk()
fenetre.title("Calcul des énergies")
fenetre.geometry("1200x900")
fenetre.columnconfigure(0, weight=1)
fenetre.rowconfigure(0, weight=1)
fenetre.rowconfigure(1, weight=1)
fenetre.rowconfigure(2, weight=1)
fenetre.rowconfigure(3, weight=1)

"""y_scrollbar = ttk.Scrollbar(cadre0, orient = VERTICAL, command = my_canvas.yview)
y_scrollbar.pack(side = RIGHT, fill = Y)
my_canvas.configure(yscrollcommand = y_scrollbar.set)
my_canvas.bind("<Configure>", lambda e: my_canvas.config(scrollregion= my_canvas.bbox(ALL)))"""

# ============================ PARAMETRES DE CALCUL ============================================

cadre1 = tkt.LabelFrame(fenetre, text="Paramétrage des variations de température, en Kelvin", padx=20, pady=20,
                        font=("Helvetica", 10), bd=5, relief="groove")
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

temperature_initiale = tkt.Entry(cadre1, width=5, font=("Helvetica", 12), justify="center")
temperature_initiale.insert(0, "300")
titre_initiale = tkt.Label(cadre1, text="Entrez la température initiale :", font=("Helvetica", 11))
temperature_finale = tkt.Entry(cadre1, width=5, font=("Helvetica", 12), justify="center")
temperature_finale.insert(0, "1000")
titre_finale = tkt.Label(cadre1, text="Entrez la température maximale :", font=("Helvetica", 11))
temperature_pas = tkt.Entry(cadre1, width=5, font=("Helvetica", 12), justify="center")
temperature_pas.insert(0, "100")
titre_pas = tkt.Label(cadre1, text="Entrez le pas de calcul :", font=("Helvetica", 11))
stock_dim = tkt.Canvas(cadre1)
dim_x = tkt.Entry(stock_dim, width=2, font=("Helvetica", 12), justify="center")
dim_y = tkt.Entry(stock_dim, width=2, font=("Helvetica", 12), justify="center")
dim_z = tkt.Entry(stock_dim, width=2, font=("Helvetica", 12), justify="center")
dim_x.insert(0, "4")
dim_y.insert(0, "4")
dim_z.insert(0, "3")
titre_dim = tkt.Label(cadre1, text="Entrez les dimensions x, y et z du cristal :", font=("Helvetica", 11))

nb_atom_maille = tkt.Entry(cadre1, width=2, font=("Helvetica", 12), justify="center")
titre_atom = tkt.Label(cadre1, text="Entrez le nombre d'atomes par maille :", font=("Helvetica", 11))
nb_atom_maille.insert(0, "2")

# ============================ SELECTION DES ATOMES ============================================

cadre2 = tkt.LabelFrame(fenetre, text="Sélection des atomes", padx=20, pady=20, font=("Helvetica", 10), bd=5,
                        relief="groove")
cadre2.grid(row=1, column=0, sticky="nsew")
cadre2.rowconfigure(0, weight=1)
cadre2.rowconfigure(1, weight=1)
cadre2.rowconfigure(2, weight=1)
cadre2.rowconfigure(3, weight=1)
cadre2.rowconfigure(4, weight=1)
cadre2.rowconfigure(5, weight=1)
cadre2.rowconfigure(6, weight=1)
cadre2.columnconfigure(0, weight=1)
cadre2.columnconfigure(1, weight=1)
cadre2.columnconfigure(2, weight=1)
cadre2.columnconfigure(3, weight=1)
cadre2.columnconfigure(4, weight=1)

selection_atome = tkt.Label(cadre2, text="Molécule à faire pénétrer dans le cristal :", font=("Helvetica", 11))
atomes_y = ['H: Dihydrogène', 'C: Carbone allotrope diamant', 'N: Diazote', 'O: Dioxygène',
            'Dihydrogène, version de débugage']
selection_y = ttk.Combobox(cadre2, values=atomes_y, justify="center", state="readonly", font=("Helvetica", 11))
selection_y.current(0)

selection_cristal = tkt.Label(cadre2, text="Sélectionnez le cristal :", font=("Helvetica", 11))
atomes_c = ['Be: Béryllium', 'Ti: Titane']
selection_r = ttk.Combobox(cadre2, values=atomes_c, justify="center", state="readonly", font=("Helvetica", 11))
selection_r.current(0)

selection_canonique = tkt.Label(cadre2, text="Ensemble canonique :", font=("Helvetica", 11))
canoniques = ['GCy: lacunaire uniquement', 'GCvy: lacunaire et interstitiel']
selection_c = ttk.Combobox(cadre2, values=canoniques, justify="center", state="readonly", font=("Helvetica", 11))
selection_c.current(0)

selection_oxyde = tkt.Label(cadre2, text="Charge des atomes :", font=("Helvetica", 11))
oxydes = ['Oxydes', 'Non chargés']
selection_x = ttk.Combobox(cadre2, values=oxydes, justify="center", state="readonly", font=("Helvetica", 11))
selection_x.current(1)

selection_conc = tkt.Label(cadre2, text="Intervalle des valeurs de concentration :", font=("Helvetica", 11))
concentrations = ['Exponentiel borné', 'Test unique', 'Valeurs linéaires']
selection_o = ttk.Combobox(cadre2, values=concentrations, justify="center", state="readonly", font=("Helvetica", 11))
selection_o.current(0)

selection_dft = tkt.Label(cadre2, text="Mode de traitement des énergies de liaison :", font=("Helvetica", 11))
debugs = ['DFT', 'Energy Binding']
selection_d = ttk.Combobox(cadre2, values=debugs, justify="center", state="readonly", font=("Helvetica", 11))
selection_d.current(0)

selection_type_graphe = tkt.Label(cadre2, text="Mode de représentation des calculs :", font=("Helvetica", 11))
types_g = ['Isotherme', 'Iso-concentration']
selection_type_g = ttk.Combobox(cadre2, values=types_g, justify="center", state="readonly", font=("Helvetica", 11))
selection_type_g.current(0)

verifData = tkt.Button(cadre2, text="Vérifier les données du fichier", command=lire_fichier, bd=3)
selectData = tkt.Button(cadre2, text="Sélectionner le fichier de données", command=select_fichier, bd=3)
textSelectData = tkt.Label(cadre2, text="Aucun fichier sélectionné", font=("Helvetica", 11))
drawGraph = tkt.Button(cadre2, text="Tracer le graphe!", command=trace, bd=3)

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

titre_atom.grid(row=2, column=3, sticky=S)
nb_atom_maille.grid(row=3, column=3, sticky=N)

selection_atome.grid(row=0, column=0, sticky=NSEW)
selection_y.grid(row=0, column=1, sticky=NSEW)
selection_oxyde.grid(row=1, column=0, sticky=NSEW)
selection_x.grid(row=1, column=1, sticky=NSEW)
selection_cristal.grid(row=2, column=0, sticky=NSEW)
selection_r.grid(row=2, column=1, sticky=NSEW)
selection_canonique.grid(row=0, column=3, sticky=NSEW)
selection_c.grid(row=0, column=4, sticky=NSEW)
selection_conc.grid(row=1, column=3, sticky=NSEW)
selection_o.grid(row=1, column=4, sticky=NSEW)
selection_dft.grid(row=2, column=3, sticky=NSEW)
selection_d.grid(row=2, column=4, sticky=NSEW)
selection_type_graphe.grid(row=3, column=3, sticky=NSEW)
selection_type_g.grid(row=3, column=4, sticky=NSEW)
verifData.grid(row=5, column=2, sticky=NSEW)
selectData.grid(row=5, column=0, sticky=NS)
textSelectData.grid(row=6, column=0, sticky=NSEW)
drawGraph.grid(row=6, column=2, sticky=NSEW)

# ============================ CALCUL DE SOLUBILITE ============================================

cadre3 = tkt.LabelFrame(fenetre, text="Calcul de solubilité", padx=20, pady=20, font=("Helvetica", 10), bd=5,
                        relief="groove")
cadre3.rowconfigure(0, weight=1)
cadre3.rowconfigure(1, weight=4)
cadre3.rowconfigure(2, weight=4)
cadre3.rowconfigure(3, weight=4)
cadre3.rowconfigure(4, weight=1)
cadre3.rowconfigure(5, weight=4)
cadre3.rowconfigure(6, weight=1)
cadre3.rowconfigure(7, weight=6)

cadre3.columnconfigure(0, weight=4)
cadre3.columnconfigure(1, weight=1)
cadre3.columnconfigure(2, weight=4)
cadre3.columnconfigure(3, weight=1)
cadre3.columnconfigure(4, weight=4)
cadre3.grid(row=2, column=0, sticky=NSEW)

textSelectFichiers = tkt.Label(cadre3, text="Aucun fichier sélectionné", font=("Helvetica", 11))
boutonSelectFichiers = tkt.Button(cadre3, text="Sélectionner les fichiers de données", command=select_fichier_solu,
                                  bd=3)
listFichier = tkt.Listbox(cadre3, height=5, font=("Helvetica", 11), justify="center")
boutonLireFichierSolu = tkt.Button(cadre3, width=30, text="Lire le fichier sélectionné", command=lire_fichier_solu,
                                   bd=3)
textSelectFichiers.grid(row=0, column=0, sticky=NSEW)
boutonSelectFichiers.grid(row=1, column=0, sticky=NSEW)
listFichier.grid(row=4, column=0, sticky=NSEW)
boutonLireFichierSolu.grid(row=5, column=0, sticky=S)

textTitreGraphique = tkt.Label(cadre3, text="Titre du graphique", font=("Helvetica", 11))
titreGraphique = tkt.Entry(cadre3, width=20, font=("Helvetica", 11), justify="center")
textNomEspece = tkt.Label(cadre3, text="Nom de l'espèce", font=("Helvetica", 11))
nomEspece = tkt.Entry(cadre3, width=20, font=("Helvetica", 11), justify="center")
textTitreGraphique.grid(row=0, column=2, sticky=NSEW)
titreGraphique.grid(row=1, column=2, sticky=NSEW)
textNomEspece.grid(row=2, column=2, sticky=NSEW)
nomEspece.grid(row=3, column=2, sticky=NSEW)

textSelectFichierFtotal = tkt.Label(cadre3, text="Aucun fichier sélectionné", font=("Helvetica", 11))
textSelectFichierFtotal.grid(row=0, column=4, sticky=NSEW)

boutonSelectFtotal = tkt.Button(cadre3, text="Sélectionner Ftotal", command=select_fichier_ftotal, bd=3)
boutonSelectFtotal.grid(row=1, column=4, sticky=NSEW)

boutonLireFtotal = tkt.Button(cadre3, text="Lire Ftotal", command=lambda: lire_fichier_ftotal(v.get()), bd=3)
boutonLireFtotal.grid(row=2, column=4, sticky=NSEW)

v = tkt.StringVar()
containerRadio = tkt.Frame(cadre3)
containerRadio.rowconfigure(0, weight=1)
containerRadio.columnconfigure(0, weight=1)
containerRadio.columnconfigure(1, weight=1)
containerRadio.grid(row=3, column=4, sticky=NSEW)
radioEv = tkt.Radiobutton(containerRadio, text="Ev", value=1, variable=v)
radioKj = tkt.Radiobutton(containerRadio, text="Kj/mol", value=2, variable=v)
radioEv.select()
radioEv.grid(row=0, column=0, sticky=NSEW)
radioKj.grid(row=0, column=1, sticky=NSEW)

boutonVerifData = tkt.Button(cadre3, text="Vérifier les données", command=check_data, bd=3)
boutonVerifData.grid(row=5, column=4, sticky=NSEW)

tkt.Label(cadre3, text="", font=("Helvetica", 11)).grid(row=6, column=0, columnspan=5, sticky=NSEW)

boutonCalcul = tkt.Button(cadre3, text="Calculer", command=calcul_solubilite, bd=3)
boutonCalcul.grid(row=7, column=2, sticky=NSEW)

texte0 = tkt.Label(fenetre,
                   text="Version: 0.4.5 du 18/04/2024 par Jawad Maache \n "
                        "Versions précédentes: 0.3 du 14/12/2023 par Kevin Gautier \n "
                        "0.2 du 06/05/2023 par Gabriel Faraut \n "
                        "0.1 du 25/06/2021 par Damien Connétable \n "
                        "E_mail: damien.connetable@ensiacet.fr \n CNRS CIRIMAT")
texte0.grid(row=3, column=0, sticky=NSEW)

fenetre.mainloop()
