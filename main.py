from sympy import *
import itertools

MULABEL = "MU"
NULABEL = "NU"

#=================================================================================================
# 4-VECTORS
#=================================================================================================

class lcvec:
    def __init__(self,plus,minus,perp):
        self.__plus = plus
        self.__perp = perp
        self.__minus = minus

    def __add__(self,vec2):
        return lcvec(self.__plus + vec2.__plus,
                     self.__minus + vec2.__minus,
                     self.__perp + vec2.__perp)

    def __sub__(self,vec2):
        return lcvec(self.__plus - vec2.__plus,
                     self.__minus - vec2.__minus,
                     self.__perp - vec2.__perp)

    def getplus(self):
        return self.__plus
    def getminus(self):
        return self.__minus
    def getperp(self):
        return self.__perp
    def geti(self,i):
        return [self.__plus, self.__minus, self.__perp][i]

    def pprint(self):
        #return "<"+self.__plus+", "+self.__minus+", "+self.__perp+">"
        return "<%s, %s, %s>"%(self.__plus, self.__minus, self.__perp)

#the actual variables
p1, p1p, p1perp = symbols("p_1^+ p_1^{'+} p_{1\\perp}")
p2, y, p2perp = symbols("p_2^- y p_{2\\perp}")
lm, lperp = symbols("l^- l_{\\perp}")
kp, kperp = symbols("k^+ k_{\\perp}")
k1p, k1perp = symbols("k_1^+ k_{1\\perp}")
k2p, k2perp = symbols("k_2^+ k_{2\\perp}")

#these only exist to format the perp contractions and output text
lperpstr, kperpstr, k1perpstr, k2perpstr = symbols('l_{\\perp} k_{\\perp} k_{1\\perp} k_{2\\perp}')
lperpsq, kperpsq, k1perpsq, k2perpsq = symbols('l_{\\perp}^2 k_{\\perp}^2 k_{1\\perp}^2 k_{2\\perp}^2')
lperpdotkperp, lperpdotk1perp, lperpdotk2perp, k1perpdotk2perp = symbols('(l_{\\perp}\\cdot\\~k_{\\perp}) (l_{\\perp}\\cdot\\~k_{1\\perp}) (l_{\\perp}\\cdot\\~k_{2\\perp}) (k_{1\\perp}\\cdot\\~k_{2\\perp})')


p2v = lcvec(0,p2,p2perp)
p1v = lcvec(p1,0,p1perp)
p1pv = lcvec(p1p,0,p1perp)

lv = lcvec(lperp*lperp/(2*lm),lm,lperp)

kv = lcvec(kp,0,kperp)
k1v = lcvec(k1p,0,k1perp)
k2v = lcvec(k2p,0,k2perp)


#=================================================================================================
# HELPERS
#=================================================================================================

def isgammaplus(g):
    return g == "+"
def isgammaminus(g):
    return g == "-"
def isgammapm(g):
    return isgammaplus(g) or isgammaminus(g)
def isgammaperp(g):
    return not isgammapm(g)
def ismu(g):
	return g == MULABEL
def isnu(g):
	return g == NULABEL

def combineetas(uppers,lowers):
    def is4dindex(index):
        return index[0]=="?"
    def perpfrom4d(index):
        return index[1:]
    def issameindex(a,b):
        return a==b or (is4dindex(a) and perpfrom4d(a)==b) or (is4dindex(b) and perpfrom4d(b)==a)

    #uppers is a list of g^{a b}s
    #lowers is a list of g_{a b}s
    #in actual practical use everything is a perp index
    #but it's not hard to include normal 4-vector index structure here so I do
    #these have a ? as a prefix to show they're actual proper 4-vector labels
    #the ? is removed when applicable to "cast" them to perps
    #e.g. g^{?mu perpa} = g^{perpmu perpa}
    #the only time this WOULDN'T happen is when we have a g^{mu nu} g_{mu nu} = 4

    #need the possibility of having a g^{a}_{b}
    def makegperp(uppers,lowers):
        assert len(uppers)+len(lowers)==2

        #if one index is a full 4-vector index but the other is just perp
        #then the 4-vector index is downgraded to a perp index
        #since g^{+ perp} = g^{- perp} = 0
        has4dindex = False
        all4dindex = True
        for index in uppers:
            if is4dindex(index):
                has4dindex = True
            else:
                all4dindex = False
        for index in lowers:
            if is4dindex(index):
                has4dindex = True
            else:
                all4dindex = False
        if has4dindex and (not all4dindex):
            #replace 4d with perp
            uppers = [perpfrom4d(index) if is4dindex(index) else index for index in uppers]
            lowers = [perpfrom4d(index) if is4dindex(index) else index for index in lowers]

        return [uppers,lowers]

    #make all of them the standard form first
    allgs = []
    for gperp in uppers:
        allgs.append(makegperp([i for i in gperp],[]))
    for gperp in lowers:
        allgs.append(makegperp([],[i for i in gperp]))

    #find the repeated indices to contract
    lowerindices = [index for gperp in lowers for index in gperp]
    upperindices = [index for gperp in uppers for index in gperp]
    repeatedindices = [i for i in lowerindices for j in upperindices if issameindex(i,j)]

    coeff = 1
    while len(repeatedindices) > 0: #take care of one repeated index
        # print("ALLGS",allgs)
        # print("REPEATED",repeatedindices)

        currindex = repeatedindices[0]

        #find currindex
        uppergi = -1
        lowergi = -1
        for gi in range(len(allgs)):
            currg = allgs[gi]
            if currindex in currg[0]:
                uppergi = gi
            if currindex in currg[1]:
                lowergi = gi

        if lowergi == uppergi: #they're in the same term
            if is4dindex(currindex):
                # print("\tFOUND A G^MU_MU")
                #we have a g^{mu}_{mu} = 4
                coeff *= 4
            else:
                # print("\tFOUND A G^A_A")
                #we have a g^{perpa}_{perpa} = 2
                coeff *= 2
            allgs.pop(lowergi)
        else: #contracting two gs
            g1 = allgs[lowergi] #has a _{a}
            g2 = allgs[uppergi] #has a ^{a}
            # print("\tCONTRACTING",g1,g2)

            newg = makegperp(g1[0]+[i for i in g2[0] if not issameindex(i,currindex)] , [i for i in g1[1] if not issameindex(i,currindex)]+g2[1])

            allgs.pop(max(lowergi,uppergi))
            allgs.pop(min(lowergi,uppergi))
            allgs.append(newg)

        #recompute what we still need to do
        upperindices = [index for gperp in allgs for index in gperp[0]]
        lowerindices = [index for gperp in allgs for index in gperp[1]]
        repeatedindices = [i for i in lowerindices for j in upperindices if issameindex(i,j)]

    # print("AT THE END:",coeff,allgs)
    retguppers = []
    retglowers = []
    retguplows = []

    for eta in allgs:
        if len(eta[1]) == 0:
            retguppers.append(eta[0])
        elif len(eta[0]) == 0:
            retglowers.append(eta[1])
        else:
            retguplows.append(eta)
            print("there's a metric tensor with one upstairs and one downstairs index you should actually code that part now")

    return {
        "coeff": coeff,
        "guppers": retguppers,
        "glowers": retglowers,
        "guplows": retguplows
    }


#=================================================================================================
# TRACE OF GAMMA MATRICES
#=================================================================================================

#----------------------------------- permutation math
def allowed_permutation_generator(N):
    assert N%2 == 1
    #arrange 1 through N in some order
    #such that for each pair of adjacent numbers starting at index 1
    #the second number in the pair is greater
    #e.g. 2(34)(15) is allowed, 2(53)(14) is not
    #also...
    #in order to not overcount e.g. (ab)(cd)(ef) and (ab)(ef)(cd)
    #we can require the left values in each pair to be strictly increasing
    for perm in itertools.permutations(list(range(1,N+1))):
        isvalid = True
        maxleftsofar = -1
        for i in range(1,N,2):
            if perm[i]>perm[i+1] or perm[i]<maxleftsofar:
                isvalid = False
                break
            maxleftsofar = perm[i]
        if isvalid:
            yield [0]+list(perm) #append the fixed first element

#https://math.stackexchange.com/questions/65923/how-does-one-compute-the-sign-of-a-permutation/72221#72221
def permutation_sign(perm):
    #this finds the sign of a permutation of [0 .. N-1]
    #sign = (-1)^(# of even length cycles)
    visited = [False for _ in range(len(perm))]
    even_cycles = 0
    for i in range(len(perm)):
        if not(visited[i]): #start investigating a new cycle
            nexti = i
            cyclelen = 0
            while not(visited[nexti]):
                cyclelen += 1
                visited[nexti] = True
                nexti = perm[nexti] #list values are zero indexed
            if cyclelen%2 == 0:
                even_cycles += 1
    return (-1)**even_cycles

#----------------------------------- trace of list of gamma matrices
def trace(gammas,symmetries=[]):
    #gammas is a string or list of all indices with whatever names you give them
    # "+" is gamma^+
    # "-" is gamma^-
    # anything else is considered to be a label for a perp 
    #e.g. gammas = "-+-ab-cd+" or ["-","+","-","a","b","-","c","d","+"]

    #symmetries = e.g. [('a','c'),('b','d')] or ["ac","bd"]
    #drop terms that are antisymmetric in e.g. a <-> c and b <-> d
    #useful when contracting
    for sym in symmetries:
        assert len(sym) == 2

    if len(gammas)%2 == 1:
        return 0


    #first anticommute all the +-s to the left
    #need to count how many perps we pass by
    perpcount = [0 for _ in range(len(gammas))]
    for i in range(len(perpcount)-2,-1,-1): #start from the right end
        if isgammaperp(gammas[i]):
            perpcount[i] = perpcount[i+1] + 1
        else:
            perpcount[i] = perpcount[i+1]

    #find the overall sign from how much we had to anticommute
    pm_anticom_factors = [1 for _ in range(len(gammas))]
    for i in range(len(gammas)):
        if isgammapm(gammas[i]):
            pm_anticom_factors[i] = (-1)**(perpcount[0]-perpcount[i])

    #now multiply them all together
    pm_anticom_factor = (-1)**pm_anticom_factors.count(-1)
    


    #take the sublist of these +-s together, now anticommuted to the left
    gammapms = [g for g in gammas if isgammapm(g)]

    #if we have two + in a row or two - in a row it's zero, including at the ends
    if len(gammapms)%2 == 1:
                return 0
    for i in range(len(gammapms)):
        if isgammaplus(gammapms[i]) and isgammaplus(gammapms[(i+1)%len(gammapms)]):
            return 0
        if isgammaminus(gammapms[i]) and isgammaminus(gammapms[(i+1)%len(gammapms)]):
            return 0

    #if we get here then gammapms is +-+-+-... or -+-+-+...
    #either way (\pm \mp)^n = 2^(n-1) \pm \mp
    pm_pair_factor = 1
    if len(gammapms) >= 2:
        pm_pair_factor = 2**( len(gammapms)/2 - 1 )

    gammaperps = [g for g in gammas if isgammaperp(g)]

    if len(gammaperps) == 0:
        #no perps, everything was +/-
        #all we have left is pm_pair_factor * tr[+ -]
        return pm_pair_factor * 4


    #at this point we're left with const * \pm * \mp * all the perps in order
    #the \pm \mp only contribute a eta^{\pm \mp} = 1 in the trace; can be dropped
    #so all that matters is the perps remaining
    #at this point every term will be 4\eta{\perp \perp}, pull the 4 out now
    overall_factor = pm_anticom_factor * pm_pair_factor * 4

    #first get the symmetries to refer to the index in gammaperps, and all ways to apply the symmetries
    perpsyms = []
    for sym in symmetries:
        try:
            perpsyms.append(  sorted( [gammaperps.index(sym[i]) for i in range(len(sym))] )  )
        except ValueError: #one of the symmetric indices has been chosen to be +/- and we can't find it here
            pass
    #perpsyms = [sorted([gammaperps.index(sym[i]) for i in range(len(sym))]) for sym in symmetries]
    #https://stackoverflow.com/a/1482316
    sympowerset = list(itertools.chain.from_iterable(itertools.combinations(perpsyms, r) for r in range(len(perpsyms)+1)))


    allterms = {} #dictionary of term: factor

    #and this keeps a standard formatting for the dictionary keys
    def standard_order(perm):
        #smaller index on the left in each pair, increasing left to right
        #this is used to keep to dictionary consistent when messing around with indices due to symmetries
        pairs = [(min(perm[i],perm[i+1]),max(perm[i],perm[i+1])) for i in range(0,len(perm),2)]
        pairs.sort(key=lambda x: x[0])
        return tuple([pairs[i][j] for i in range(len(pairs)) for j in range(2)])


    #the trace of 2N perp matrices labeled 012345...
    #is a product of eta^{ij}s for all permutations of the indices WHERE
    #the first index is 0
    #each eta^{ij} has the smaller index first (otherwise we double count each one)
    #e.g. for 012345, every term is eta^{0A}eta^{BC}eta^{DE}
    #and we take all possibilities where B<C, D<E
    #the sign of each term is given as the sign of the permutation (0ABCDE)
    for allowedperm in allowed_permutation_generator(len(gammaperps)-1): #(N-1)!! in total
        #allowedperm is a list e.g. [0,3,1,2] of what indices of gammaperps to use in what order

        permsign = permutation_sign(allowedperm)

        #now account for symmetries in all possible ways using the powerset
        foundsym = False #have we found it in the dictionary yet
        for symset in sympowerset:
            symperm = [i for i in allowedperm] #deep copy

            #use some symmetries to swap elements
            for sym in symset: #has two indices to swap
                tmp = symperm.index(sym[0])
                symperm[symperm.index(sym[1])] = sym[0]
                symperm[tmp] = sym[1]

            if standard_order(symperm) in allterms: #found a symmetric term in the dictionary
                allterms[standard_order(symperm)] += permsign
                foundsym = True
                break
        if not foundsym: #new term
            allterms[standard_order(allowedperm)] = permsign

    retdict = {
        "overall_factor": overall_factor,
        "terms": []
    }
    for etaterm in allterms:
        if allterms[etaterm] != 0: #nonzero coefficient
            etas = []
            for i in range(0,len(etaterm)-1,2):
                alphi = gammaperps[etaterm[i]]
                alphj = gammaperps[etaterm[i+1]]
                etas.append((alphi,alphj))
            retdict["terms"].append((etas, allterms[etaterm]))

    return retdict


# print(trace('i+-+-+gc-d',['cg']))
# print(trace(['+','-','a','b','+','c','d','-']))
# print(trace('+-ab+cd-'))
# print(trace('+-abc-d-'))
# print(trace('+-abd+-+ef'))
# print(trace(['-','+']))
# print(trace(['-','+','-','+','-','+','a','b']))
# print(trace(['+','-','+','-','+','\\alpha','\\beta','-']))
# print(trace(['a','-','b','+','-','+','-','+','-','+']))
# print(trace('ab'))
# print(trace('abcd'))
# print(trace('abcdef'))
# print(trace('abcdefg'))
# print(trace('abcdefgh',['cf']))
# print(trace(['+','-','a','\\beta','b','\\gamma','d','e'],['ab']))
# print(trace('i+-bB+gc-d',['bc','Bg']))
# print(trace('-+gAabdB',['ab']))
# raise ValueError



#================================================================================================= 
# MAIN
#================================================================================================= 

def fulltrace(terms,symmetries=[],Gterms=[],guppers=[],glowers=[]):
    def fulltrace_helper(terms,symmetries=[],Gterms=[],guppers=[],glowers=[]):
        gammas = [] #label of the gamma matrix
        vectors = [] #name of the vector
        vectorindices = [] #associated index of the vector

        dummyname = "PERP%d"
        i = 0

        for term in terms:
            if isinstance(term, lcvec): #we have a pslash
                gammas.append(dummyname%i)
                vectors.append(term)
                vectorindices.append(dummyname%i)
                i+=1
            elif term[0] == "?": #wildcard index for gluon props
                #we keep track of the index but there's no associated vector
                gammas.append(term[1:])
                vectors.append(None)
                vectorindices.append(term[1:])
            else: #just a bare gamma matrix
                gammas.append(term)

        finalanswer = 0
        finalterms = [] #(possibleindices,result)

        #loop over every 3^n possible option for the unknown indices
        for possibleindices in itertools.product(["+","-","p"], repeat=len(vectorindices)):
            currgammas = [_ for _ in gammas]
            currsymmetries = [_ for _ in symmetries]
            currguppers = [_ for _ in guppers]
            currglowers = [_ for _ in glowers]

            vectorterms = []
            perpterms = [] #redundant but easier
            Gmunus = [[-1,-1] for _ in range(len(Gterms))] #the indices of each gluon prop

            #replace in indices in the gamma array
            for i in range(len(vectorindices)):
                coordi = ["+","-","p"].index(possibleindices[i])

                #check if this index was in a gluon prop
                for Gtermi in range(len(Gterms)):
                    Gterm = Gterms[Gtermi]
                    if vectorindices[i] == Gterm[1][0]: #first index
                        Gmunus[Gtermi][0] = coordi
                    elif vectorindices[i] == Gterm[1][1]: #second index
                        Gmunus[Gtermi][1] = coordi

                if coordi != 2: #keep the name label for the perp
                    currgammas[currgammas.index(vectorindices[i])] = possibleindices[i]

                #add the terms we're contracting with
                def innerprodi(i):
                    return [1,0,2][i]

                if vectors[i] != None:
                    vectorterms.append((innerprodi(coordi), vectors[i]))
                    if coordi == 2:
                        perpterms.append((vectorindices[i], vectors[i].geti(coordi)))

            #or maybe the G term is between non-wildcard indices of a perp
            for Gtermi in range(len(Gterms)):
                Gterm = Gterms[Gtermi]
                for label in currgammas:
                    if label == Gterm[1][0]:
                        Gmunus[Gtermi][0] = 2
                    elif label == Gterm[1][1]:
                        Gmunus[Gtermi][1] = 2

            plusminusGfactors = 1

            for Gtermi in range(len(Gterms)):
                Gterm = Gterms[Gtermi]
                assert (-1 not in Gmunus[Gtermi])
                if Gmunus[Gtermi] == [1,1]: #- -
                    val = (Gterm[0].getperp())**2
                    plusminusGfactors *= (Gterm[0].getperp())**2 / (Gterm[0].getminus())**2                    
                elif Gmunus[Gtermi] == [1,2]: #- perp
                    perpterms.append((Gterm[1][1], Gterm[0].getperp()))
                    plusminusGfactors /= Gterm[0].getminus()
                elif Gmunus[Gtermi] == [2,1]: #perp -
                    perpterms.append((Gterm[1][0], Gterm[0].getperp()))
                    plusminusGfactors /= Gterm[0].getminus()
                elif Gmunus[Gtermi] == [2,2]: #perp perp
                    currglowers.append((Gterm[1][0],Gterm[1][1]))
                    plusminusGfactors *= -1
                else: #everything else is zero
                    plusminusGfactors = 0

            #find all pairwise symmetries in the perp components
            for pair in itertools.combinations(perpterms, 2):
                comp1 = pair[0][1]
                comp2 = pair[1][1]
                if comp1-comp2 == 0 or comp1+comp2==0: #could have a minus sign
                    currsymmetries.append([pair[0][0], pair[1][0]])

            currtrace = trace(currgammas, symmetries=currsymmetries)

            if currtrace != 0 and plusminusGfactors!=0:#0 not in Gfactors:
                # print(possibleindices)
                # print(perpterms)
                # print(vectorindices)
                # print(currtrace)
                # print(currsymmetries)
                # print(currgammas)
                # print(plusterms,minusterms,Gfactors)
                # print(plusminusGfactors)
                # print(currguppers)
                # print(currglowers)
                thistermtotal = 0

                for term in vectorterms:
                    if term[0] == 0: #plus component
                        plusminusGfactors *= term[1].getplus()
                    elif term[0] == 1: #minus component
                        plusminusGfactors *= term[1].getminus()

                if isinstance(currtrace,dict):
                    #loop over each term
                    for term in currtrace['terms']:
                        #contract metric tensors as much as we can on their own
                        combinedgs = combineetas(currguppers+term[0], currglowers)
                        termguppers = [sorted(pair) for pair in combinedgs['guppers']]
                        termglowers = [sorted(pair) for pair in combinedgs['glowers']]
                        termguplows = [sorted(pair) for pair in combinedgs['guplows']]

                        #should have no free indices other than to contract with the vectors
                        assert (len(termglowers)==0) 
                        assert (len(termguplows)==0)

                        finalterm = currtrace['overall_factor'] * plusminusGfactors * term[1] * combinedgs['coeff']

                        #do the contractions
                        #each adds a minus sign since eta^{perpa perpb} v_aperp w_bperp = -v_perp cdot w_perp
                        perpcontractions = []
                        for eta in termguppers: #find the two perp indices
                            vec1 = None
                            vec2 = None
                            for vec in perpterms:
                                if vec[0] == eta[0]:
                                    vec1 = vec[1]
                                elif vec[0] == eta[1]:
                                    vec2 = vec[1]

                            if not (vec1 and vec2):
                                print("Couldn't contract index?")
                                raise SyntaxError

                            dotprod = vec1*vec2
                            dotprod = simplify(dotprod.subs([
                                (lperp*kperp, lperpdotkperp),
                                (lperp*k1perp, lperpdotk1perp),
                                (lperp*k2perp, lperpdotk2perp),
                                (k1perp*k2perp, k1perpdotk2perp),
                                # (lperp*lperp, lperpsq),
                                # (kperp*kperp, kperpsq),
                                # (k1perp*k1perp, k1perpsq),
                                # (k2perp*k2perp, k2perpsq),
                            ]))

                            finalterm *= -1
                            finalterm *= dotprod

                        thistermtotal += finalterm

                    finalanswer += thistermtotal
                    finalterms.append((possibleindices,thistermtotal))
                else:
                    #just got a float out from the trace

                    #sure hope all these metric tensors contract with themselves
                    #because there's none from the trace
                    combinedgs = combineetas(currguppers, currglowers)
                    termguppers = [sorted(pair) for pair in combinedgs['guppers']]
                    termglowers = [sorted(pair) for pair in combinedgs['glowers']]
                    termguplows = [sorted(pair) for pair in combinedgs['guplows']]

                    #should have no free indices
                    assert (len(termguppers)==0)
                    assert (len(termglowers)==0) 
                    assert (len(termguplows)==0)

                    finalterm = currtrace * plusminusGfactors * combinedgs['coeff']

                    finalanswer += finalterm
                    finalterms.append((possibleindices,finalterm))

        return {
            "allterms": finalterms,
            "answer": finalanswer
        }

    #======================================================================================

    def pprint(answer,ret=False):
        #should maybe use a regex in the future
        try:
            answer = answer.subs(lm,y*p2)
        except AttributeError:
            pass
        pretty = str(answer).replace(".0","").replace("p_2^-**","(p_2^-)**").replace("**","^").replace("*","")
        if ret:
            return pretty
        else:
            print(pretty + " \\\\")

    #check well formed args, append symmetries
    for Gterm in Gterms: # (vec, (index1, index2))
        assert (len(Gterm)==2)
        assert (len(Gterm[1])==2)
        if Gterm[1] not in symmetries:
            symmetries.append(Gterm[1])
    for eta in guppers: # (index1, index2)
        assert (len(eta)==2)
        if eta not in symmetries:
            symmetries.append(eta)
    for eta in glowers: # (index1, index2)
        assert (len(eta)==2)
        if eta not in symmetries:
            symmetries.append(eta)

    if MULABEL in terms and NULABEL in terms: #do a g_{mu nu} contraction
        muindex = terms.index(MULABEL)
        nuindex = terms.index(NULABEL)
        lindex = min(muindex,nuindex)
        rindex = max(muindex,nuindex)

        #to help show how the indices are being shuffled
        dummyindices = list(range(len(terms)-2))
        if lindex == muindex:
            dummyindices[muindex:muindex] = [MULABEL]
            dummyindices[nuindex:nuindex] = [NULABEL]
        else:
            dummyindices[nuindex:nuindex] = [NULABEL]
            dummyindices[muindex:muindex] = [MULABEL]

        print("INPUT TERMS:")
        print(dummyindices)
        print()

        nonmunuterms = [i for i in terms if i not in [MULABEL, NULABEL]]
        allsubanswers = [] #anticommute stuff around until we get mu and nu by each other

        if rindex == lindex+1:
            print("MU NU immediately contracts out as 4I")
            print("TRACE 1:")
            print(list(range(len(terms)-2)))
            result = fulltrace_helper(nonmunuterms,symmetries,Gterms,guppers,glowers)

            print("\tALLOWED TERMS:")
            for term in result['allterms']:
                print('\t',term[0], '=>', pprint(term[1],ret=True))
            print()

            allsubanswers.append((4,result))
        else:
            tracei = 1
            sign = 1
            while rindex >= lindex+2:
                #anticommute the one on the right towards the one on the left
                #the g^{A nu} gets contracted with gamma^mu g_{mu nu} to change the gamma^mu to a gamma^A
                dummyindices = list(range(len(terms)-2))
                contracted = dummyindices.pop(rindex-2) #-1 is because lindex is gone
                dummyindices[lindex:lindex] = [contracted]

                result = fulltrace_helper([nonmunuterms[i] for i in dummyindices],symmetries,Gterms,guppers,glowers)

                print("TRACE %d:"%tracei)
                print(dummyindices)
                print("\tALLOWED TERMS:")
                for term in result['allterms']:
                    print('\t',term[0], '=>', pprint(term[1],ret=True))
                print()

                #stop condition has
                #gamma^mu gamma^A gamma_mu = -2gamma^A
                if rindex == lindex+2:
                    allsubanswers.append((-2*sign,result))
                else:
                    allsubanswers.append((2*sign,result))

                sign*=-1
                rindex-=1
                tracei+=1

        print("OVERALL SUM:")
        sumstr = ""
        for i in range(len(allsubanswers)):
            subtrace = allsubanswers[i]
            if subtrace[0] > 0:
                sumstr += " + %d*[TRACE %d]"%(abs(subtrace[0]),i+1)
            else:
                sumstr += " - %d*[TRACE %d]"%(abs(subtrace[0]),i+1)
        if sumstr[:3] == " + ":
            sumstr = sumstr[3:]
        else:
            sumstr = sumstr[1:]

        print(sumstr)
        overallsum = sum(subtrace[0] * subtrace[1]['answer'] for subtrace in allsubanswers)
        pprint(overallsum)

    else:
        if MULABEL in terms or NULABEL in terms:
            raise ValueError

        result = fulltrace_helper(terms,symmetries,Gterms,guppers,glowers)
        print("ALLOWED TERMS:")
        for term in result['allterms']:
            print(term[0], '=>', pprint(term[1],ret=True))
        print()
        print("OVERALL SUM:")
        pprint(result['answer'])


#=================================================================================================
#=================================================================================================

#10a
fulltrace(["+","-",kv+p2v,"\\gamma",kv+p2v+p1pv,lv,kv+p2v+p1v,"\\beta",kv+p2v,"-"], glowers=[ ("\\gamma","\\beta") ])

#12a
#fulltrace(['+','-',kv+p2v,MULABEL,lv-p1pv,'\\gamma',lv,'\\beta',lv-p1v,NULABEL,kv+p2v,'-'], glowers=[ ("\\gamma","\\beta") ])

#12b
#fulltrace(['-','?\\beta',lv-p1v,NULABEL,p2v+kv,'-','+','-',p2v+kv,MULABEL,lv-p1pv,'?\\gamma'], Gterms=[ (lv,('\\gamma','\\beta')) ])

#15e
#fulltrace(['-',lv-kv-p2v,'?\\gamma',kv+p2v,'\\alpha','-','\\beta',kv+p2v,'?\\delta',lv-kv-p2v], Gterms=[ (lv,('\\delta','\\gamma')) ], glowers=[ ("\\alpha","\\beta") ])
