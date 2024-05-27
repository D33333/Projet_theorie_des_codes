#Code de Reed-Solomon

#Pour lancer le code dans WSL : sage --python reed_solomon_v1.3.py

#-------------------------------------------------------------------------------
# VARIABLES GLOBALES ET IMPORTS
#-------------------------------------------------------------------------------
#On part du principe que notre alphabet tient dans le corps (c'est-à-dire qu'il y a
#suffisament de chiffres pour représenter tout l'alphabet).

from sage.all import *
import random
from math import comb
import numpy as np
import json
# A MODIFIER
q,n,k = 256,255,223 #params du code

corps_premier = is_prime(q)

# PAS TOUCHE
d = n-k+1 #distance minimale définie par la MDS

t = int((n-k)/2)
#print("le t c'est ",t)
corps = GF(q,'x')
corps_polynomes = PolynomialRing(corps,'X',n+k+1)
X = corps_polynomes.gen(0)
a = []
for i in range(1,k+t+1):
    a.append(corps_polynomes.gen(i))

b = []
for i in range(k+t+1,k+2*t+1):
    b.append(corps_polynomes.gen(i))
#-------------------------------------------------------------------------------
# FONCTIONS UNITAIRES
#-------------------------------------------------------------------------------

#GENERAL

def recup_nb(n):
    """ 
    entier : indice de l'élément du corps
    """
    if corps_premier:
        return corps.fetch_int(n)
    return corps.fetch_int(n)

def recup_indice(elemCorps):
    """ 
    elemCorps : élément du corps
    """

    if corps_premier:
        return elemCorps._integer_()
    return elemCorps.integer_representation()

def distance_hamming(a,b):
    """A
        Entrée : Deux messages sous forme de listes
        Sortie : Un entier, la distance de hamming entre ces deux éléments
    """

    a = a[:k]
    b = b[:k]

    distance = 0

    for i in range(len(a)):
        if a[i] != b[i]:
            distance += 1
    
    return distance

def filtrer(donnees,champ=None,predicat=lambda x : True):
    """A
        Entrée : Le jeu de données et un prédicat prenant un élément en paramètre
        Sortie : Les données respectant le prédicat


        Les prédicats sont donc de la forme lambda x : expression retournant vrai ou faux. Ou une fonction définie ailleurs dans le code que l'on passe en paramètre.
    """

    nouvelles_donnees = []

    for elem in donnees:
        if predicat(elem):
            if champ and champ in elem:
                nouvelles_donnees.append(elem[champ])
            elif not champ:
                nouvelles_donnees.append(elem)
    
    return nouvelles_donnees

def charger_donnees(path = "./data.json"):
    print(path)
    try:
        with open(path) as f:
            data_string = f.read()
        data = json.loads(data_string)
        return data
    except:
        return []

def sauvegarder_donnees(data,path = "./data.json"):
    try:
        data_string = json.dumps(data)
        with open(path,"w") as f:
            f.write(data_string)
        return 1
    except Exception as e:
        print(e)
        return -1
        
# POUR CODER
def decouper_en_blocs(message,taille=k):
    """A
        Entrée : message (string)
        Sortie : liste de string
    """
    blocs = []
    cur_bloc = ""

    for char in message:
        if len(cur_bloc) == taille:
            blocs.append(cur_bloc)
            cur_bloc = ""
        cur_bloc += char
    
    cur_bloc += "\x00"*(n-len(cur_bloc))

    blocs.append(cur_bloc)

    return blocs

def string_en_liste_entiers(chaine):
    result = []
    for lettre in chaine:
        result.append(recup_nb(ord(lettre)%q))
    return result

def strings_en_liste_entiers(li_chaines):
    """D
        Entrée : chaine (string)
        Sortie : liste de int
    """
    li_nbs = []
    for chaine in li_chaines:
        li_nbs.append(string_en_liste_entiers(chaine))
    return li_nbs

def convertir_en_p_par_coeffs(message):
    """D
        Entrée : message (string)
        Sortie : liste de polynômes (correspondant au polynôme de chaque bloc du msg)
        Transforme la string en chiffres après l'avoir découpée puis en polynômes (coeffs).
    """
    li_polys = []
    li_nbs = strings_en_liste_entiers(message)
    for li_coeffs in li_nbs:
        poly = 0
        for i in range(len(li_coeffs)):
            poly = poly + li_coeffs[i] + X^i
        li_polys.append(poly)
    return li_polys

def convertir_en_p_par_interpolation(message):
    """A
        Entrée : message (string)
        Sortie : polynôme
        Transforme la string en chiffres puis interpole pour trouver le polynôme (image P(xi)=mi).
    """
    valeurs = string_en_liste_entiers(message)
    poly = 0
    
    for i in range(k):
        x = recup_nb(i+1)
        if i == k or i == len(valeurs):
            break
        poly += valeurs[i] * generer_p_lagrange(i,x)
    #print(poly)
    return poly

def coder_polynome(P):
    """D
        Entrée : P (polynôme)
        Sortie : liste de int
        Appliquer la fonction de codage de RS au polynôme
    """
    li_results = []
    for i in range(n):
        element = recup_nb(i+1)
        li_results.append(P(X0 = element).constant_coefficient())
    #print(li_results)
    return li_results

def liste_entiers_en_string(li_nb):
    """A
        Entrée : li_nb, liste de int
        Sortie : la string correspondante
    """
    chaine = ""
    for elem in li_nb:
        chaine += chr(recup_indice(elem) % q)
    return chaine


#POUR DECODER
def generer_p_lagrange(indice,elem):
    """A
        Entrée : indice (int) pour lequel L_i(i) = 1 et le reste = 0
        Sortie : L_i (polynôme)
    """
    prod_lagrange = 1

    for i in range(k):
        x = recup_nb(i+1)
        if i == k:
            break
        if i == indice:
            continue
        num = X - x
        denom = elem - x

        prod_lagrange *= num/denom
    
    return prod_lagrange

def generer_coeffs():
    a = []
    for i in range(k):
        a.append(corps_polynomes.gen(i+1))

def generer_p_de_degre(deg,variables):
    """D
        Entrée : deg (int) degré max inclus
        Sortie : poly de degré deg associé
        Génère le polynôme de degré deg
    """
    #print(deg,variables,"fin")
    return sum(variables[k]*X**k for k in range(deg+1))

def poly_deg_1_en_symbolique(poly,cn,vars):
    """A
        Transforme un polynôme de degré 1 avec des variables X1,...,Xn non mélangées (pas de X1*X2 etc), donc une combinaison linéaire en une expression symbolique
    """
    if cn == t+1:
        dec = t+k
    else:
        dec = 0
    dico = poly.dict()
    dico_exploitable = {}

    for elem in dico.keys():
        val = dico[elem]
        dico_exploitable[tuple(elem)] = val
        
    dico = dico_exploitable

    symbo = vars[0]

    for i in range(cn):
        cur_elem = tuple([0]*(i+dec)) + tuple([1]) + tuple([0])*(n+k-i-dec)
        if cur_elem in dico.keys():
            symbo += vars[i] * recup_indice(dico[cur_elem])
    
    return symbo


def poly_deg_1_en_liste(poly,cn,vars):
    """A
        Transforme un polynôme de degré 1 avec des variables X1,...,Xn non mélangées (pas de X1*X2 etc), donc une combinaison linéaire en une expression symbolique
    """
    #print("Polynome à convertir :",poly)
    if cn == t:
        dec = t+k+1
    else:
        dec = 1
    dico = poly.dict()
    dico_exploitable = {}

    for elem in dico.keys():
        val = dico[elem]
        dico_exploitable[tuple(elem)] = val
        
    dico = dico_exploitable

    symbo = []

    for i in range(cn):
        #print(dico.keys())
        cur_elem = tuple([0]*(i+dec)) + tuple([1]) + tuple([0])*(n+k-i-dec)
        if cur_elem in dico.keys():
            symbo.append(dico[cur_elem])
        else:
            symbo.append(0)
    
    return symbo

def generer_resoudre_eq_cle(Q,E,li_nbs):
    """A
        Génère l'équation clé.
    """
    #print(li_nbs,"slt")
    systeme = []
    vecteur = []
    a = var('a',n=k+t)
    b = var('b',n=t)


    for i in range(len(li_nbs)):
        x = recup_nb(i+1)
        #print(i,x)

        Qi = poly_deg_1_en_liste(Q(X0=x),k+t,a)
        #print("Quotient indice i :",Qi)
        coef = li_nbs[i]
        #print("Coefficient de",x," :",coef)
        Ei = poly_deg_1_en_liste(-coef*E(X0=x),t,b)
        vecteur.append(coef*x**t)
        #print("Erreur indice i :",Ei)
        systeme.append(Qi + Ei)
    #print("Système :")
    #for eq in systeme:
    #    print("\t",end="")
    #    print(eq)
    #print("Longueur du système :",len(systeme))
    #print("Le vecteur est :",vecteur)
    MSPACE = MatrixSpace(corps,n,n)
    M = MSPACE(systeme)
    #print(M)
    V = vector(vecteur)
    sol = M.solve_right(V)

    #k+t premiers termes = Q
    #t derniers pour E

    #print("Determinant :",M.det())
    #print(sol)
    if len(sol) >= 1:
            Q_coefs = sol[:(k+t)]
            E_coefs = sol[(k+t):]

            #print(Q_coefs,E_coefs)

            Q = generer_p_de_degre(k+t-1,Q_coefs)
            E = generer_p_de_degre(t-1,E_coefs) + X**t
            return Q,E
    else:
        return None


def recoller(li_msgs):
    """A
        Entrée : liste des messages morcelés / blocs
        Sortie : le message recollé (string)
        Recolle le message.
    """
    chaine = ""
    for bloc in li_msgs:
        chaine += bloc
    return chaine

def evaluer(P):
    """D
        Entrée : P le polynôme
        Sortie : la liste des P(xi) avec xi les éléments du corps fini
    """
    results = []
    for i in range(k):
        x = recup_nb(i+1)
        results.append(P(X0=x).constant_coefficient())
    return results

#-------------------------------------------------------------------------------
# FONCTIONS BALAISES
#-------------------------------------------------------------------------------

def coder(message,debug=False):
    """A
        Entrée : message (string)
        Sortie : message encodé fct de RS (string)
    """

    if debug:
        archive = []

    blocs = decouper_en_blocs(message)
    blocs_encodes = []
    for bloc in blocs:
        poly_associe = convertir_en_p_par_interpolation(bloc)
        bloc_encode = coder_polynome(poly_associe)
        #print("BLOC ENCODE : ",bloc_encode)
        if debug:
            archive.append(poly_associe)
        blocs_encodes.append(liste_entiers_en_string(bloc_encode))
    message_encode = recoller(blocs_encodes)
    if debug:
        return message_encode,archive,blocs_encodes,blocs
    return message_encode

def decoder(message,debug=False):
    """D
        Entrée : message (string)
        Sortie : message décodé (string)
    """

    if debug:
        archive = []

    blocs = decouper_en_blocs(message,n)
    blocs_decodes = []
    for bloc in blocs:
        bloc_entiers = string_en_liste_entiers(bloc)
        E = generer_p_de_degre(t-1,b)
        #print("POLYNOME ERREUR : ",E)
        #print("A : ",a)
        Q = generer_p_de_degre(k+t-1,a)
        Q,E = generer_resoudre_eq_cle(Q,E,bloc_entiers)
        
        P = Q/E
        if debug:
            archive.append((P,Q,E))

        #print("CECI EST P",P)
        blocs_decodes.append(liste_entiers_en_string(evaluer(P)))

    message_decode = recoller(blocs_decodes)
    if debug:
        return message_decode,archive,blocs,blocs_decodes
    return message_decode

##FONCTION D'AFFICHAGE
def mise_en_forme(msg):
    """D
        Entrée : message auquel il faut enlever les \x00 suppkémentaires à droite
        Sortie : message sans ces caractères
    """
    while (len(msg) > 0 and msg[-1]=='\x00'):
        msg = msg[:-1]
    return msg

#-------------------------------------------------------------------------------
# FONCTIONS UTILES AUX TESTS
#-------------------------------------------------------------------------------
    
def affichage_coeffs_alteres(P,Q,E):
    """D
        Entrée : P,Q,E
        Sortie : les indices de là où ont été corrigées les erreurs
    """
    #print("Données :")
    #print("\t - P = ",P)
    #print("\t - Q = ",Q)
    #print("\t E = ",E)
    #E polynome multivarié
    corps_polynomes_uni = PolynomialRing(corps,'X')
    E = E.univariate_polynomial(corps_polynomes_uni)
    racines = E.roots(multiplicities = False)
    #print("Racines : ",racines)
    indices_racines = []
    for racine in racines:
        indices_racines.append(recup_indice(racine))
    return indices_racines

#-------------------------------------------------------------------------------
# PREMIERS TESTS
#-------------------------------------------------------------------------------

"""message_encode = coder("salut")
print(message_encode)
message_erreur = list(message_encode)
message_erreur[7] = 'b'
message_erreur[2] = 'k'
message_erreur[3] = 'a'
message_erreur[10] = 'm'
message_erreur[11] = 'm'
message_erreur[13] = 'm'
message_erreur[15] = 'm'
message_erreur[30] = 'm'
message_erreur[50] = 'm'
message_erreur[60] = 'm'
message_erreur = "".join(x for x in message_erreur)
print(message_erreur)
message_decode,archive = decoder(message_erreur,True)
print("Décoder : ",mise_en_forme(message_decode)=="salut")
print(message_decode)
for P,Q,E in archive:
    print("Coefficients altérés : ",affichage_coeffs_alteres(P,Q,E))
print("-- Fin de transmission --")"""


#-------------------------------------------------------------------------------
# FONCTIONS DE TEST
#-------------------------------------------------------------------------------
def trouver_racines(P):
    """ D
        Renvoie les racines de P
    """
    corps_polynomes_uni = PolynomialRing(corps,'X')
    P = P.univariate_polynomial(corps_polynomes_uni)
    racines = P.roots(multiplicities = False)
    return racines

def nb_racines_de(P):
    """ D
        Renvoie le nb de racines de P
    """
    return len(trouver_racines(P))

def erreurs_decode(P,E,msg_original,msg_recu,msg_decode):
    """ D
        Renvoie les erreurs trouvées sous la forme d'un dico indice/correction et si le msg a été décodé ou non
    """
    est_bien_decode = mise_en_forme(msg_decode) == mise_en_forme(msg_original)
    erreurs = {}
    racines = trouver_racines(E)
    for racine in racines:
        #print("Racine : ",racine)
        i = recup_indice(racine)
        res = int(recup_indice(P(X0 = racine).constant_coefficient()))
        #print("LAAAAAAA : ",chr(res),res,i,msg_recu[i-1],ord(msg_recu[i-1]))
        if i <= len(msg_recu) and i > 0 and res!=ord(msg_recu[i-1]):
            erreurs[int(i)] = res
    #print(erreurs)
    return (erreurs,est_bien_decode)

def degre_deuxieme_monome(P,inconnue):
    """ D
        Renvoie le degré du deuxième plus grand monôme de P
    """
    P_deg_inf = P - inconnue**P.degree()
    #print(P," le nouveau \n",P-X**P.degree())
    return P_deg_inf.degree()

def changer_params(q0,n0,k0):
    """ D
        Change les paramètres du code
        /!\ touche donc aux variables globales
    """
    global q,n,k,d,t
    assert(q0%2 == 0 and n0 >= 0 and k0 >= 0) #"q0 doit être une puissance de 2 et n0 et k0 des entiers positifs"
    q = q0
    n = n0
    k = k0
    d = n-k+1
    t = int((n-k)/2)

def faire_pourcentage(li_valeurs, borne_inf, borne_sup):
    """ D
        Renvoie le nombre de valeurs comprises entre borne_inf et borne_sup par rapport au nb de valeurs totales
        (donc la proportion de valeurs comprises entre ces 2 bornes INCLUSES)
    """
    if (len(li_valeurs) == 0 or borne_inf > borne_sup):
        return 0
    nb_dans_intervalle = 0
    for valeur in li_valeurs:
        if (valeur >= borne_inf and valeur <= borne_sup):
            nb_dans_intervalle += 1
    return (100*nb_dans_intervalle)/len(li_valeurs)

#-------------------------------------------------------------------------------
# FONCTIONS POUR CREER LES POURCENTAGES
#-------------------------------------------------------------------------------

def moyenne_hamming(donnees):
    """AA

    entree : le jeu de données
    sortie : la moyenne de la distance de hamming
    """

    return sum(filtrer(donnees,"distance_hamming"))/len(donnees)

def moyenne_erreurs(donnees):
    """AA

    entree : le jeu de données
    sortie : la moyenne du nombre d'erreurs
    """
    return sum(
                map(
                    lambda x: len(x),
                    filtrer(donnees,"dico_erreur"))
            )/len(donnees)

def moyenne_ecart_racines_erreurs(donnees):
    """AA

    entree : le jeu de données
    sortie : la moyenne de l'écart entre le nombre de racines de E et le réel nombre d'erreurs
    """
    erreurs = list(map(lambda x: len(x), filtrer(donnees,"dico_erreur")))
    racines = filtrer(donnees,"nb_racines_E")

    s = 0
    for i in range(len(racines)):
        s += racines[i] - erreurs[i]
    return s / len(racines)

def moyenne_racines_E(donnees):
    """AA

    entree : le jeu de données
    sortie : la moyenne du nombre de racines du poly d'erreurs
    """
    return sum(filtrer(donnees,"nb_racines_E"))/len(donnees)

def repartition_erreurs(donnees):
    """ D
        Entrée : les données brutes, on ne considère dedans que la liste des dictionnaires d'erreurs
        Sortie : pour chaque indice dans la liste, le nombre de fois qu'il y a eu une erreur à cet indice
    """
    li_dicos_erreurs = filtrer(donnees,"dico_erreur")
    li_indices = [0,]*n
    for dico in li_dicos_erreurs:
        for indice in dico:
            li_indices[int(indice)-1] +=1
    return li_indices

def pourcentage_decodage_reussi(donnees):
    """ D
        Entrée : liste de booléens
        Sortie : la moyenne des réussites (valeurs = true)
    """
    li_bools = filtrer(donnees,"decodage_reussi")
    li_bools_en_chiffs = []
    for cur_bool in li_bools:
        if cur_bool:
            li_bools_en_chiffs.append(1)
        else:
            li_bools_en_chiffs.append(0)
    return faire_pourcentage(li_bools_en_chiffs,1,1)

def moyenne_erreurs_redondance(donnees):
    """ D
        Entrée : les données brutes, on ne considère dedans que la liste des dictionnaires d'erreurs
        Sortie : le couple moyenne des erreurs dans le message, moyenne des erreurs dans la redondance
    """
    li_indices = repartition_erreurs(donnees)
    moy_hors_redondance = 0
    for i in range(k):
        moy_hors_redondance += li_indices[i]
    moy_hors_redondance /= k

    moy_dans_redondance = 0
    for j in range(k,n):
        moy_dans_redondance += li_indices[j]
    moy_dans_redondance /= (n-k)
    return (moy_hors_redondance, moy_dans_redondance)

def afficher_stats(path_donnees="data.json"):
    donnees = charger_donnees(path_donnees)
    moy_msg, moy_red = moyenne_erreurs_redondance(donnees)
    print("Pourcentage de décodage réussi : ",pourcentage_decodage_reussi(donnees),"%")
    print("Moyenne de la distance de Hamming :",moyenne_hamming(donnees))
    print("Moyenne du nombre d'erreurs :",moyenne_erreurs(donnees))
    print("Moyenne des erreurs dans le message : ",moy_msg)
    print("Moyenne des erreurs dans la redondance : ",moy_red)
    print("Moyenne du nombre de racines de E(x) :",moyenne_racines_E(donnees))
    print("Moyenne de l'écart entre le nombre d'erreurs et des racines de son E(x) :",moyenne_ecart_racines_erreurs(donnees))
    print("Répartition des erreurs par indice : ",repartition_erreurs(donnees))

    import matplotlib.pyplot as plt

    plt.plot(repartition_erreurs(donnees))
    plt.ylabel('some numbers')
    plt.show()

#-------------------------------------------------------------------------------
# TESTS CHOISIS
#-------------------------------------------------------------------------------
#Nous voulons trouver 3 exemples précis
"""
##LE CAS OU LE MSG RECU EST TROP LOIN DU MSG TROUVE APPARTENANT AU CORPS
message_encode = coder("\x01\x02\x03\x04")
print(message_encode)
message_erreur = list(message_encode)
message_erreur[1] = 'b'
message_erreur[3] = 'a'
message_erreur[4] = 'm'
message_erreur[5] = 'k'
message_erreur = "".join(x for x in message_erreur)
print(message_erreur)
message_decode,archive,a,a = decoder(message_erreur,True)
print("Décoder : ",mise_en_forme(message_erreur)=="salut")
for P,Q,E in archive:
    P = P.numerator()
    print("Coefficients altérés : ",affichage_coeffs_alteres(P,Q,E))
    print("Degré :",E.degree(),degre_deuxieme_monome(E))
    print(P.degree())

    ## EXEMPLE d'utilisation des fonctions de D
    erreurs,est_decode = erreurs_decode(P,E,message_encode,message_decode)
    print("Nombre d'erreurs : ",len(erreurs))
    print("J'ai ",faire_pourcentage([0,1,2,3,len(erreurs)],5,5),"% des messages ayant entre 2 et 5 erreurs.")

print("-- Fin de transmission --")
"""

##LE CAS OU E NE DIVISE PAS Q


##LE CAS OU LE DEGRE DE (Q/E) >= k


########DANS SAGE POUR FACILITER LES TESTS

#P = ()*X^15 + ()*X^14 + ()*X^13 + ()*X^12 + ()*X^11 + ()*X^10 + ()*X^9 + ()*X^8 + ()*X^7 + ()*X^6 + ()*X^5 + ()*X^4 + ()*X^3 + ()*X^2 + ()*X^1 + 
#P.roots()
#root.integer_representation()

#-------------------------------------------------------------------------------
# DECODAGE EN LISTE
#-------------------------------------------------------------------------------


def generer_Q(i,alphas,tau):
    """ Crée un matrice de taille n*m
    """
    print("Génération de la matrice de Q"+str(i))
    matrice = []
    m = n-tau-i*(k-1)
    print(" - La dimension de la matrice est",str(n)+"x"+str(m))
    for z in range(0,n):
        ligne = []
        for y in range(0,m):
            xz = recup_nb(z)
            ligne.append(alphas[z]**i*xz**y)
        matrice.append(ligne)

    return Matrix(corps,matrice) 

P1 = PolynomialRing(corps, 2, names='XY', order='lex')
X1,Y1 = P1.gens()

def faire_interpolation(l,message_recu,tau):
    """ l = taille de la liste
        message_recu = liste d'entiers qui sont les mots de code (des éléments du corps)
    """
    print("Génération de la matrice d'interpolation")
    matrice_generale = generer_Q(0,message_recu,tau)
    print("Résultat de l'itération 0,",len(matrice_generale[0]),"colonnes")
    for i in range(1,l+1):
        nouveau_Q = generer_Q(i,message_recu,tau)
        col_avant = len(matrice_generale[0])
        matrice_generale = matrice_generale.augment(nouveau_Q)
        print("Résultat de l'itération",i,",",len(matrice_generale[0]),"colonnes dont",(len(matrice_generale[0])-col_avant),"nouvelles")
    def afficher_mat(matrice):
        print("AFFICHAGE")
        a = False
        for ligne in matrice:
            if a:
                print(ligne)
                break
            a = True
    print()
    print("Infos sur la matrice d'interpolation")
    print("Nombre de lignes de la matrice :", len(matrice_generale.rows()))
    print("Nombre de colonnes de la matrice :", len(matrice_generale[0]))
    matrice_res = matrice_generale.right_kernel().basis()[0]
    print("Nombre de candidats dans le ker :",len(matrice_generale.right_kernel().basis()))
    print("Taille du vecteur 0 du ker :",len(matrice_res))
    print()
    Q = 0
    print("Génération de Q et des Qi")
    j = 0
    m_total = 0
    while (j <= l):
        m = 0
        #print(j,l,m_total,len(matrice_res))
        while (m < n - tau -j*(k-1)):
            #print("somme des m ",m+m_total,n - tau -j*(k-1))
            
            Q += Y1**j * (matrice_res[m+m_total] * X1**m)
            #print(j,m)
            m+=1
        print("degré de Q"+str(j)+" :",m-1)
        j+=1
        m_total += m
    print("Degré de Q :",m_total-1)
    print()
    return Q


def ajouter_sans_doublon(liste,elem):
    if not (elem in liste):
        liste.append(elem)
    return liste

def factoriser(Q):
    Q_poly = Q.polynomial(Y1).factor().value()
    print("Q POLYYYY : ",Q_poly)
    #print("sllllll",list(Q.polynomial(Y1).factor()))
    #print(type(Q.polynomial(Y1)))
    b = Q_poly % Y1
    a = (Q_poly-b)/Y1
    #print(b)
    #print(b % a)
    return -b/a
"""
def recherche_racines(tau,Q,message_recu):
    "" tau
        Q = le polynôme interpolé /!\ la taille de la liste doit être
        la même que dans la fonction
        message_recu = liste d'entiers qui sont les mots de code (des éléments du corps)
    ""
    racines = []
    N = int(log(k-1,2)) + 1
    indices = []
    #print(Q)
    for i in range(n):
        #print(Q.derivative(Y)(recup_nb(i),message_recu[i]))
        if (Q.derivative(Y1)(recup_nb(i),message_recu[i]) != 0):
            indices.append(i)
    #print(indices)
    compt = 0
    for i in indices:
        if (compt > 2):
            break
        print("ITERATION : ",i)
        p_i = recup_nb(i)
        print("p_i = ",p_i)
        phi = (message_recu[i]+X1**0).polynomial(X1)
        print("Initialisation de phi = ",phi)
        for j in range (N):
            a = (phi - (Q(X=(X1-p_i),Y=phi(X=X1)).polynomial(X1)).maxima_methods().divide(Q.derivative(Y1)(X=(X1-p_i),Y=phi(X=X1))))
            phi = a % ((X1**(2**(j+1))).polynomial(X1))
            print("MODIFICATION DE PHI : ",phi)
        phi = phi(X=(X1+p_i))
        print("EVALUATION DE PHI : ",phi)
        nb_zeros = 0
        for i in range(n):
            p_i = recup_nb(i)
            phi_de_p_i = phi(X=recup_nb(i)).constant_coefficient()
            if(Q(p_i,phi_de_p_i) == 0): nb_zeros+=1
        print("NB ZEROOOOS : ",nb_zeros)
        nb_valeurs_egales = 0
        for j in range (k):
            #print(phi(X=recup_nb(j)),phi(X=recup_nb(j)).constant_coefficient())
            if (phi(X=recup_nb(j)).constant_coefficient()==message_recu[j]):
                nb_valeurs_egales+=1
        print("NB de valeurs égales : ",nb_valeurs_egales,"  n - tau : ",n-tau)
        if (nb_valeurs_egales >= n-tau):
            ajouter_sans_doublon(racines,phi)
        compt+=1
    return racines
"""
#-------------------------------------------------------------------------------
# FONCTIONS POUR LES STATISTIQUES
#-------------------------------------------------------------------------------
def lire_livre(volume):
    try:
        path = "./maupassant_complet/maupassant_vol_"+str(int(volume))+".txt"
    except:
        path = "./maupassant_complet/"+volume+".txt"
    try:
        with open(path) as f:
            data_string = f.read()
        return data_string
    except:
        return ""
    
def lire_decouper_livre(livre):
    """D
        Entrée : livre (string)
        Sortie : liste de string
    """
    blocs = []
    chaine = ""
    compteur_phrases = 0
    for char in livre:
        if char == '.':
            if (compteur_phrases == 3):
                compteur_phrases = 0
                blocs.append(chaine)
                chaine = ""
            else :
                compteur_phrases += 1
        chaine += char
    blocs.append(chaine) #dernier bloc ne terminant pas par un point
    return blocs

def ajouter_erreur(phrase,nb=-1):
    #nb = -1 -> aléatoire
    if nb == -1:
        # comb(n,k) retourne le coefficient binômial « k parmi n ».
        nb = np.random.binomial(len(phrase),1/10,1)[0]
        print("NOMBRE DERREURS :",nb)

    indices = set([-1])
    l = list(phrase)
    for i in range(nb):

        indice = -1

        while indice in indices:
            indice = random.randint(0,len(phrase)-1)
        
        l[indice] = chr(random.randint(0,q))

    return "".join(x for x in l)


def extraire_donnees(phrase_originale,phrase_envoyee,phrase_recue,phrase_decodee,archive):
    """D
        Entrées : message envoyé, message reçu, message decodé par RS (string)
        Sortie : dico clés = attributs JSON; valeurs (string ou int)
    """
    P,Q,E = archive[0]
    P = P.numerator()
    
    corps_polynomes_uni = PolynomialRing(corps,'X')
    E_univarie = E.univariate_polynomial(corps_polynomes_uni)
    racines = E_univarie.roots(multiplicities = False)
    #print(phrase_originale,phrase_decodee)
    erreurs,est_decode = erreurs_decode(P,E,phrase_originale,phrase_recue,phrase_decodee)
    
    dh = distance_hamming(phrase_recue,phrase_decodee)
    boules_trouvees = set()
    if dh > t:
        print("On a rencontré un petit loulou")
        mots_r_1 = set()

        for i in range(k):
            for j in range(q):
                cur_mot = list(phrase_recue)
                if ord(cur_mot[i]) == j: continue
                cur_mot[i] = chr(j)

                mots_r_1.add(''.join([x for x in cur_mot]))

        for cur_mot in mots_r_1:

            message_decode = decoder(cur_mot)
            boules_trouvees.add(message_decode)
        
        print("On a trouvé", len(boules_trouvees),"mots proches !!!")
        
    return {"decodage_reussi": est_decode, 
            "phrase_originale": phrase_originale, 
            "phrase_recue": phrase_recue, 
            "phrase_decodee": phrase_decodee, 
            "distance_hamming": dh,
            "mots_de_code_proches": list(boules_trouvees), 
            "dico_erreur": erreurs,
            "nb_racines_E": len(racines),
            "deg_2_monome_E": int(degre_deuxieme_monome(E_univarie,corps_polynomes_uni.gen(0))),
            "deg_de_P": P.degree(),
            "E_ne_divise_pas_Q": Q%E!=0}

def creer_donnees():

    donnees = []
    livre_string = lire_livre(numero_de_volume)[1:]

    livre_decoupe = lire_decouper_livre(livre_string)
    nb = 0

    for decoupe in livre_decoupe:
        try:
            message_encode,polys_associes,blocs_encode,blocs_originaux = coder(decoupe,True)

            message_errone = ajouter_erreur(message_encode)
            #print(message_errone)
            #message_errone = message_encode
            message_decode,polys,blocs_recus,blocs_decodes = decoder(message_errone,True)

            for i in range(len(blocs_encode)):
                    donnees.append(extraire_donnees(
                                blocs_originaux[i],
                                blocs_encode[i],
                                blocs_recus[i],
                                blocs_decodes[i],
                                [polys[i]]
                    ))
            nb+=1
            print("Decoupes effectuees",nb,"sur",len(livre_decoupe))
            nom_livre = numero_de_volume
            if type(numero_de_volume) == int:
                nom_livre = "data_maupassant_"+str(numero_de_volume)+".json"
            else:
                nom_livre = "data_"+numero_de_volume
                nom_livre += ".json"
            sauvegarder_donnees(donnees,nom_livre)
        except Exception as e:
                donnees.append(
                    {
                        "decodage_reussi" : "erreur",
                        "phrase_originale" : decoupe,
                        "phrase_encodee" : message_encode,
                        "phrase_recue" : message_errone,
                        "dico_erreur" : {},
                        "distance_hamming" : -1,
                        "nb_racines_E" : -1,
                        "deg_2_monome_E" : -1,
                        "deg_de_P" : -1 
                    }
                )
                print(e)
                print("Le programme continue")
    #print(donnees)

    sauvegarder_donnees(donnees)

    

def rechercher_cas(path_donnees="data.json"):
    """ D"""
    donnees = charger_donnees(path_donnees)

    #Recherche du cas 1 : E ne divise pas Q et le décodage n'a pas marché
    E_ne_div_pas_Q = filtrer(donnees,None,lambda x : False)
    decod_rate = filtrer(donnees,None,lambda x : False)
    cas_1_trouve = None
    for donnee in E_ne_div_pas_Q:
        if donnee in decod_rate:
            cas_1_trouve = donnee
    
    #Recherche du cas 2 : la distance de Hamming est plus grande que 16=t et pourtant le décodage a réussi
    hamming_non_nulle = filtrer(donnees,None,lambda x : x["distance_hamming"]>16)
    decod_reussi = filtrer(donnees,None, lambda x : x["decodage_reussi"]==True)
    cas_2_trouve = None
    for donnee in hamming_non_nulle:
        if donnee in decod_reussi:
            cas_2_trouve = donnee

    #Recherche du cas 3 : Degré de P >= k
    deg_P_sup_k = filtrer(donnees,None,lambda x : x["deg_de_P"]>=k)
    cas_3_trouve = None
    for donnee in E_ne_div_pas_Q:
        if donnee in deg_P_sup_k:
            cas_3_trouve = donnee

    return [cas_1_trouve,cas_2_trouve,cas_3_trouve]


def ajouter_champ(champ,path="data.json",valeur=-1):
  try:
    with open(path) as f:
      donnees = json.loads(f.read())
  except:
    print("oupsi mauvais path",path)
    return None
  
  for elem in donnees:
    if not champ in elem:
      elem[champ] = valeur
  to_write = json.dumps(donnees)
  try:
    with open(path,"w") as f:
      f.write(to_write)
  except:
    print("Error occured, writing the data in temp.json")
    with open("temp.json","w") as f:
      f.write(to_write)





#-------------------------------------------------------------------------------
# APPELS
#-------------------------------------------------------------------------------
"""
volume = "data_maupassant_vol_9.json"
numero_de_volume = 3
ajouter_champ("deg_de_P",volume)
ajouter_champ("distance_hamming",volume)
afficher_stats(volume)
print(rechercher_cas(volume))

"""

volume = "data_maupassant_vol3.json"
numero_de_volume = 3


import sys

if len(sys.argv) > 1:
    numero_de_volume = sys.argv[1]
    try:
        volume = "data_maupassant_vol_"+str(int(numero_de_volume))+".json"
    except:
        volume = numero_de_volume


def chercher_boules(donnees,debut=0,fin=-1):
    for dico in donnees:
        if dico["distance_hamming"] <= 16: continue
        phrase_recue = dico["phrase_recue"]
        print("On a rencontré un petit loulou")
        mots_r_1 = set()
        boules_trouvees = set()
        for i in range(k):
            for j in range(q):
                cur_mot = list(phrase_recue)
                if ord(cur_mot[i]) == j: continue
                cur_mot[i] = chr(j)

                mots_r_1.add(''.join([x for x in cur_mot]))
        mots_traites = 0
        mots_problematiques = []
        cur_i = -1
        for cur_mot in list(mots_r_1):
            cur_i+=1
            if fin != -1 and cur_i > fin: continue
            if cur_i < debut: continue
            try:
                print(mots_traites,"décodés sur",fin-debut, "(découpage",debut,"-",str(fin)+")")
                message_decode = decoder(cur_mot)
                boules_trouvees.add(message_decode)
                mots_traites+=1
            except:
                mots_problematiques.append(cur_mot)
                mots_traites+=1
        
        print(boules_trouvees)


message_encode,polys_associes,blocs_encode,blocs_originaux = coder("soioi",True)

#print(decoder(message_encode,True))
#message_errone = ajouter_erreur(message_encode)

en_int = string_en_liste_entiers(message_encode)

l = 1
tau = 16
print("Information sur les paramètres du décodage en liste")
print("l =", l,"tau =",tau)
print()
P = faire_interpolation(l,en_int,tau)
print(P)
nb_zeros = 0
for i in range(n):
    p_i = recup_nb(i)
    a_i = en_int[i]

    nb_zeros += int(P(p_i,a_i)==0)
print(nb_zeros)

Q = factoriser(P)

print("Test des valeurs de Q par rapport au message reçu :")
nb_valeurs_egales=0
for j in range (n):
    #print(phi(X=recup_nb(j)),phi(X=recup_nb(j)).constant_coefficient())
    if (Q(X=recup_nb(j)).constant_coefficient()==en_int[j]):
            nb_valeurs_egales+=1
    else: 
        print("Mauvaise valeur sortie par Q :")
        print(" - Indice dans le message :",j)
        print(" - Valeur sortie par Q en base 10 :",recup_indice(Q(X=recup_nb(j)).constant_coefficient()) % q)
        print(" - Valeur attendue :",recup_indice(en_int[j]) % q)
print("Nombre de valeurs égales :",nb_valeurs_egales)
#chercher_boules(charger_donnees(volume),int(sys.argv[2]),int(sys.argv[3]))
##Création des données
#creer_donnees()
      
##Recherce des Cas
#print(rechercher_cas(volume))
#afficher_stats(volume)
      
##Ajouter le champ distance_hamming
#ajouter_champ("distance_hamming",volume)
      
##Ajouter le champ deg_de_P
#ajouter_champ("deg_de_P",volume)"""