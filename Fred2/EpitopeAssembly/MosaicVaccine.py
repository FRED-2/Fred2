# This code is part of the Fred2 distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""
.. module:: EpitopeAssembly.MosaicVaccine
   :synopsis: Embedding of the mosaic vaccine design problem into the asymmetric orienteering problem.
The methods offers an exact solution for small till medium sized problems as well as heuristics based on a Matheuristic
using Tabu Search and Branch-and-Bound for large problems.

The heuristic proceeds as follows:

I: initialize solution s_best via greedy construction

s_current = s_best
WHILE convergence is not reached DO:
    I: s<-Tabu Search(s_current)
   II: s<-Intensification via local MIP(s) solution (allow only alpha arcs to change)
       if s > s_best:
          s_best = s
  III: Diversification(s) to escape local maxima
END
.. moduleauthor:: schubert

"""

from __future__ import division

import itertools as itr
import multiprocessing as mp
import math

import numpy as np
import collections
import string
import copy

from pyomo.environ import *
from pyomo.opt import SolverFactory

from Fred2.Core import EpitopePredictionResult


class TabuList(collections.MutableSet):

    def __init__(self, iterable=None, size=None):
        self.end = end = []
        end += [None, end, end]         # sentinel node for doubly linked list
        self.map = {}                   # key --> [key, prev, next]
        self.maxSize = size
        if iterable is not None:
            self |= iterable

    def __len__(self):
        return len(self.map)

    def __contains__(self, key):
        return tuple(key) in self.map

    def add(self, key):
        key = tuple(key)
        if key not in self.map:
            if len(self.map) >= self.maxSize:
                self.pop()
            end = self.end
            curr = end[1]
            curr[2] = end[1] = self.map[key] = [key, curr, end]

    def discard(self, key):
        if key in self.map:
            key, prev, next = self.map.pop(key)
            prev[2] = next
            next[1] = prev

    def __iter__(self):
        end = self.end
        curr = end[2]
        while curr is not end:
            yield curr[0]
            curr = curr[2]

    def __reversed__(self):
        end = self.end
        curr = end[1]
        while curr is not end:
            yield curr[0]
            curr = curr[1]

    def pop(self, last=True):
        if not self:
            raise KeyError('set is empty')
        key = self.end[1][0] if last else self.end[2][0]
        self.discard(key)
        return key

    def __repr__(self):
        if not self:
            return '%s()' % (self.__class__.__name__,)
        return '%s(%r)' % (self.__class__.__name__, list(self))

    def __eq__(self, other):
        if isinstance(other, OrderedSet):
            return len(self) == len(other) and list(self) == list(other)
        return set(self) == set(other)


def suffixPrefixMatch(m):
    ''' Return length of longest suffix of x of length at least k that
        matches a prefix of y.  Return 0 if there no suffix/prefix
        match has length at least k.
    '''
    x,y,k = m
    if x == y:
        return 0
    if x == "start":
        return len(y)
    if y == "start":
        return 0

    if len(x) < k or len(y) < k:
        return 0
    idx = len(y) # start at the right end of y
    # Search right-to-left in y for length-k suffix of x
    while True:
        hit = string.rfind(y, x[-k:], 0, idx)
        if hit == -1: # not found
            return len(y)
        ln = hit + k
        # See if match can be extended to include entire prefix of y
        if x[-ln:] == y[:ln]:
            return len(y)-ln # return length of prefix
        idx = hit + k - 1 # keep searching to left in Y
    return -1


def _start(args):
        return _move(*args)


def _move(func, args):
        return func(*args)


def _vertexDel(vins, vdel,  best_obj, best_sol, __imm, __arcCost):
    """
        simultaniously insert a new vertex and delete a
        existing vertex
    """
    pertuped = []
    for e in best_sol:
        if vdel == e[0]:
            pertuped.append((vins,e[1]))
        elif vdel == e[1]:
            pertuped.append((e[0],vins))
        else:
            pertuped.append(e)
    return best_obj+__imm[vins]-__imm[vdel], sum(__arcCost[i][j] for i,j in pertuped), pertuped, vins


def _vertexIns(vertex, edge,  best_obj, best_sol, __imm, __arcCost):
    """
        Insert a new vertex between vertex i and vertex j
    """
    pertuped = []
    for e in best_sol:
        if edge == e:
            pertuped.append((e[0],vertex))
            pertuped.append((vertex, e[1]))
        else:
            pertuped.append(e)
    length = sum(__arcCost[i][j] for i,j in pertuped)
    return best_obj+__imm[vertex], length, pertuped, vertex


class MosaicVaccineTS:

    def __init__(self, _results, threshold=None, k=10, solver="glpk", verbosity=0):

        #check input data
        if not isinstance(_results, EpitopePredictionResult):
            raise ValueError("first input parameter is not of type EpitopePredictionResult")

        #start constructing model
        self.__pool = mp.Pool(mp.cpu_count())
        self.__solver = SolverFactory(solver)#, solver_io = "python")
        self.__verbosity = verbosity
        self.__changed = True
        a = _results.columns.tolist()
        self.__alleleProb = self.__init_alleles(a)
        self.__k = k
        self.__result = None
        self.__thresh = {} if threshold is None else threshold
        self.__imm, self.__peps =  self.__init_imm(_results)
        self.__n = len(self.__peps)
        self.__arcCost = self.__init_arc_cost()
        self.instance = self.__init_model()

    def __init_alleles(self, _alleles, verbosity=0):
            """
                initializes allele objects an tests if they have probs
                if not copys them and uniformily distributes probability within locus
            """
            prob = []
            no_prob = []
            for a in _alleles:
                if a.prob is None:
                    no_prob.append(copy.deepcopy(a))
                else:
                    prob.append(a)

            if len(no_prob) > 0:
                #group by locus
                no_prob_grouped = {}
                prob_grouped = {}
                for a in no_prob:
                    no_prob_grouped.setdefault(a.locus, []).append(a)
                for a in prob:
                    prob_grouped.setdefault(a.locus, []).append(a)

                for g, v in no_prob_grouped.iteritems():
                    total_loc_a = len(v)
                    if g in prob_grouped:
                        remaining_mass = 1.0 - sum(a.prob for a in prob_grouped[g])
                        for a in v:
                            a.prob = remaining_mass/total_loc_a
                    else:
                        for a in v:
                            a.prob = 1.0/total_loc_a

            if verbosity:
                for a in _alleles:
                    print a.name, a.prob
            return prob+no_prob

    def __init_imm(self, _results):
            __thresh = self.__thresh
            res_df = _results.xs(_results.index.values[0][1], level="Method")
            #print "before",len(res_df)
            res_df = res_df[res_df.apply(lambda x: any(x[a] > __thresh.get(a.name, -float("inf"))
                                                       for a in res_df.columns), axis=1)]
            pep = ["start"]
            imm = [0]

            for tup in res_df.itertuples():
                p = tup[0]
                pep.append(p)
                imm.append(sum(a.prob*s for a, s in itr.izip(self.__alleleProb, tup[1:])
                        if s > __thresh.get(a.name, -float("inf"))))
            return imm, pep

    def __init_arc_cost(self):

            pool = self.__pool
            return [pool.map(suffixPrefixMatch, ((str(i), str(j), 1) for j in self.__peps)) for i in self.__peps]

    def __init_model(self):
            """
                initializes MIP model for OR with MTZ sub tour elimination
            """
            model = ConcreteModel()

            #Sets
            model.Nodes = RangeSet(0,len(self.__peps)-1)
            def arc_init(model):
                return ((i,j) for j in model.Nodes for i in model.Nodes if i != j)
            model.Arcs = Set(dimen=2,initialize=arc_init)
            def NodesOut_init(model, node):
                return [ j for (i,j) in model.Arcs if i == node]
            model.NodesOut = Set(model.Nodes, initialize=NodesOut_init)
            def NodesIn_init(model, node):
                return [ i for (i,j) in model.Arcs if j == node]
            model.NodesIn = Set(model.Nodes, initialize=NodesIn_init)

            #Params
            def i_init(model, i):
                return self.__imm[i]
            model.i = Param(model.Nodes, initialize= i_init)
            def d_init(model, i, j ):
                return self.__arcCost[i][j]
            model.d = Param(model.Arcs, initialize=d_init)
            model.TMAX = Param(initialize=self.__k, within=PositiveIntegers, mutable=True)

            #Variables
            model.x = Var(model.Arcs, domain=Binary, bounds=(0,1), initialize=0)
            model.u = Var(model.Nodes - set([0]), bounds=(1.0,len(self.__peps)-1))

            model.Obj = Objective(rule=lambda model:sum( model.x[i,j]*model.i[i] for i,j in model.Arcs),
                sense=maximize)

            #conecitfity constraint
            def Conn1_rule(model, node):
                return sum( model.x[node,j] for j in model.NodesOut[node]) <= 1.0
            model.Conn1 = Constraint(model.Nodes,rule=Conn1_rule)
            def Conn2_rule(model, node):
                return sum( model.x[i,node] for  i in model.NodesIn[node]) <= 1.0
            model.Conn2 = Constraint(model.Nodes,rule=Conn1_rule)
            #Equality constraint
            def Equal_rule(model, node):
                return sum( model.x[node,j] for j in model.NodesOut[node]) == sum( model.x[i,node] for  i in model.NodesIn[node])
            model.Equal = Constraint(model.Nodes, rule=Equal_rule)
            #Knapsack Constriant
            def Knapsack_rule(model):
                return sum( model.d[i,j]*model.x[i,j] for i,j in model.Arcs) <= self.__k
            model.Knapsack = Constraint(rule=Knapsack_rule)
            #Subout Elimination MTZ
            def Subtour_rule(model, i,j):
                return model.u[i]-model.u[j]+(len(self.__peps)-1)*model.x[i,j] <= len(self.__peps)-2
            model.SubTour = Constraint(((i,j) for i in xrange(1,len(self.__peps))
                    for j in  xrange(1,len(self.__peps)) if i != j), rule=Subtour_rule)
            model.c = ConstraintList()
            model.tabu = ConstraintList()

            return model.create()

    def solve(self, options=""):
            """
                solves the model optimally
            """
            instance = self.instance
            instance.x.reset()
            instance.u.reset()
            instance.preprocess()

            res = self.__solver.solve(instance, options=options, tee=False)
            instance.load(res)
            if self.__verbosity > 0:
                res.write(num=1)
            sol = []
            unsorted = dict([(i, j) for i, j in instance.Arcs if instance.x[i, j].value > 0])
            i=0
            while unsorted:
                j = unsorted[i]
                sol.append((i, j))
                del unsorted[i]
                i = j
            return instance.Obj(), sol

    def approximate(self, phi=0.05, options="", _greedyLP=True, _tabu=True, _intensify=True, _jump=True,
                    max_iter=10000, delta_change=1e-4, max_delta=101, seed=23478234):
            """
                Matheueristic using Tabu Search
            """

            def __greedy_init():
                __imm = self.__imm
                __arcCost = self.__arcCost
                __k = self.__k
                normalized_gain = np.divide(np.array(__imm),np.array(__arcCost))
                imm = 0
                length = 0
                possible = set(xrange(1,self.__n))
                i = 0
                result = []
                while length < __k:
                    if length == 0:
                        j = _,j = max([ (normalized_gain[0,j],j)  for j in possible])
                        result.append((0,j))
                        length += __arcCost[0][j]
                        i = j
                    else:
                        possible.discard(i)
                        _,j = max([ (normalized_gain[i,j],j)  for j in possible])
                        length += __arcCost[i][j]
                        imm += __imm[i]
                        if length > __k:
                            result.append((i,0))
                            return imm, result
                        result.append((i,j))
                        i = j
                imm += __imm[j]
                result.append((j,0))
                return imm,result

            def __lp_init():
                """
                initializes first construction based on LP relaxation of the global problem with rounding
                allows for complete constraints not only the capasity an stuff
                see TABU SEARCH FOR MIXED INTEGER PROGRAMMING Joao Pedro Pedroso
                """
                #set variables to non negativ real to obtain relaxation

                #copy is supoptimal for large instances
                #but it seams as if the change of domain can only be done once??
                #at least an error occures when solving the MIP after solving its relaxation
                instance = copy.deepcopy(self.instance)
                solver = self.__solver
                __n = self.__n
                __k = self.__k
                __imm = self.__imm
                __arcCost = self.__arcCost

                instance.x.domain = NonNegativeReals
                instance.preprocess()
                result = []
                imm = 0
                length = 0
                possible = set(xrange(1,__n))
                i = 0
                while length < __k:
                    lp_result = solver.solve(instance, options=options)
                    instance.load(lp_result)
                    if length == 0:
                        _,j = max([ (instance.x[0,j].value,j)  for j in possible])
                        instance.x[0,j].setlb(1)
                        instance.x[0,j].setub(1)
                        result.append((0,j))
                        length += __arcCost[0][j]
                        i = j
                    else:
                        possible.discard(i)
                        _,j = max([ (instance.x[i,j].value,j)  for j in possible  ])
                        length += __arcCost[i][j]
                        imm += __imm[i]
                        if length > __k:
                            result.append((i,0))
                            return imm, result
                        #possible = remove_j(possible,i)
                        result.append((i,j))
                        instance.x[i,j].setlb(1)
                        instance.x[i,j].setub(1)
                        i = j
                imm += __imm[j]
                result.append((j,0))
                del instance
                return imm, result

            def __tabu_search(best_obj, best_sol):
                """
                    performe Tabu Search according to Liang et al. IEEE 2002
                    some moves are missing -> arc swap for example. dont know if
                    needed.
                """
                __k = self.__k
                __imm = self.__imm
                __arcCost = self.__arcCost
                #generates tasks
                best_vertices = set(sum(best_sol, ()))
                remaining_vertices = set(xrange(self.__n)) - best_vertices
                tasks = []
                for v in remaining_vertices:
                    for e in best_sol:
                        if tenure.get(e[0],0) < curr_iter and tenure.get(e[1],0) < curr_iter:
                            tasks.append((_vertexIns,(v, e, best_obj, best_sol, __imm, __arcCost)))
                for r in remaining_vertices:
                    for v in best_vertices:
                        if tenure.get(v,0) < curr_iter:
                            tasks.append((_vertexDel,(r, v, best_obj, best_sol, __imm, __arcCost)))
                #make neighborhood
                neighbor = self.__pool.map(_start, tasks)
                changed_vertex = None
                best_length = 0
                b_obj = -float('inf')
                b_sol = None
                #filter neighbour with tabu list and best_obj
                for obj, length, sol, v in neighbor:
                    if sol not in tabu_list and  obj > b_obj and length <= __k:
                        b_sol = sol
                        b_obj = obj
                        best_length = length
                        changed_vertex = v
                if b_sol is None:
                    return best_obj, best_sol
                tabu_list.add(b_sol)
                #iteration dependent? the earlier the iteration the longer tabu?
                tenure[changed_vertex] = curr_iter + math.ceil(math.sqrt((len(best_vertices)+len(remaining_vertices))/curr_iter)) \
                            + math.floor((len(best_vertices)/(self.__k/9.))*best_length)
                return b_obj, b_sol

            def __mip_intensification(phi, curr_sol):
                """
                    refine solution by allowing at most phi new variables
                """
                #set warmstart
                instance = self.instance
                instance.x.reset()
                instance.u.reset()
                self.instance.x.domain = Binary
                selected_vars = set(curr_sol)
                sorted_nodes = [ i for (i,j) in curr_sol]
                for i,j in instance.Arcs:
                    instance.x[i,j] = 1 if (i,j) in selected_vars else 0
                    if i != 0:
                        try:
                            k = sorted_nodes.index(i)
                            instance.u[i] = k+1
                        except ValueError:
                            instance.u[i] = 1
                    if j != 0:
                        try:
                            k = sorted_nodes.index(j)
                            instance.u[j] = k+1
                        except ValueError:
                            instance.u[j] = 1
                #clear possible contained constraint in constraint list and add new constraint
                instance.c.clear()
                instance.tabu.clear()
                instance.c.add(sum(1-instance.x[i,j] for i,j in curr_sol) <= max(math.ceil(phi*len(curr_sol)),2))
                instance.preprocess()
                result = self.__solver.solve(instance, options=options+",timelimit=10" if options else "timelimit=10",warmstart=True)# tee=True)
                instance.load(result)
                sol = []
                unsorted = dict([(i,j) for i,j in instance.Arcs if instance.x[i,j].value > 0])
                i=0
                while unsorted:
                    j = unsorted[i]
                    sol.append((i,j))
                    del unsorted[i]
                    i=j
                return instance.Obj(), sol

            def __jump(best_sol):
                """
                    change 50% of the current solution with new nodes
                    use __lp_init with a initial starting set of vertices
                    (in random order)
                """
                #set variables to non negativ real to obtain relaxation
                instance = copy.deepcopy(self.instance)
                solver = self.__solver
                rand = np.random.rand
                __n = self.__n
                __k = self.__k
                __imm = self.__imm
                __arcCost = self.__arcCost

                instance.x.domain = NonNegativeReals
                instance.c.clear()
                instance.tabu.clear()
                for sol in tabu_list:
                    instance.tabu.add(sum(1 - instance.x[i,j] for i,j in sol) >= 1)
                #best_vertices = set(sum(best_sol, ()))
                #selected = []
                for k, (i,j) in enumerate(best_sol):
                    if rand() >= 0.5:
                        instance.x[i,j].setlb(1)
                        instance.x[i,j].setub(1)
                        #selected.append((k,i,j))
                        #if i > 0:
                        #    instance.u[i].setlb(k+1)
                        #    instance.u[i].setub(k+1)
                        #if j > 0:
                        #    instance.u[j].setlb(k+2)
                        #    instance.u[j].setub(k+2)
                instance.preprocess()
                result = []
                imm = 0
                length = 0
                possible = set(xrange(1,__n))
                i = 0
                while length < __k:
                    lp_result = solver.solve(instance, options=options+",timelimit=20" if options else "timelimit=20")
                    instance.load(lp_result)
                    if length == 0:
                        _,j = max([ (instance.x[0,j].value,j)  for j in possible])
                        instance.x[0,j].setlb(1)
                        instance.x[0,j].setub(1)
                        result.append((0,j))
                        length += __arcCost[0][j]
                        i = j
                    else:
                        possible.discard(i)
                        _,j = max([ (instance.x[i,j].value,j)  for j in possible  ])
                        length += __arcCost[i][j]
                        imm += __imm[i]
                        if length > __k:
                            result.append((i,0))
                            return imm, result
                        #possible = remove_j(possible,i)
                        result.append((i,j))
                        instance.x[i,j].setlb(1)
                        instance.x[i,j].setub(1)
                        i = j
                imm += __imm[j]
                result.append((j,0))
                return imm, result
            ##########################################

            np.random.seed(seed)

            if _greedyLP:
                best_obj, best_sol = __lp_init()
            else:
                best_obj, best_sol = __greedy_init()
            curr_obj, curr_sol = best_obj, best_sol

            if not _tabu:
                return best_obj, best_sol
            print "Start solution: ", best_obj, best_sol
            curr_iter = 1
            delta = 1
            tabu_list = TabuList(size=self.__k)
            tabu_list.add(best_sol)
            tenure = {0:max_iter}
            while (curr_iter < max_iter) and (delta < max_delta):
                curr_obj, curr_sol = __tabu_search(curr_obj, curr_sol)

                if curr_obj > best_obj:
                    best_obj = curr_obj
                    best_sol = curr_sol

                if abs(best_obj - curr_obj) <= delta_change or curr_obj < best_obj:
                    delta += 1

                if _intensify:
                    #sould not be done every time
                    #TODO: find good rule when to apply intensification
                    if not curr_iter % 50:
                        curr_obj, curr_sol = __mip_intensification(phi, curr_sol)
                        if curr_obj > best_obj:
                            best_obj = curr_obj
                            best_sol = curr_sol
                            delta = 0
                            #if _jump:
                                #diversification should be performed after each
                                #increase of objective (?)
                             #   curr_obj, curr_sol = __jump(curr_sol)
                if delta % max(math.floor(max_delta/5),10) == 0 and _jump:
                    curr_obj, curr_sol = __jump(curr_sol)
                curr_iter += 1

            self.__pool.close()
            return best_obj, sum(self.__arcCost[i][j] for i,j in best_sol), best_sol



