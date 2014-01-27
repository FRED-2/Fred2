'''
    Created on Apr 11, 2013
    
    This class implements the epitope selection functionality
    of OptiTope published by Toussaint et al. [1].
    
    This module builds upon Coopr's Pyomo, an embedded algebraic modeling
    languages [2].
    
    It allows to (de)select specific constraints of the underlying
    ILP and to solve the specific problem with a MIP solver of choice
    
    
    [1] N. C. Toussaint and O. Kohlbacher. OptiTope--a web server for the selection of
    an optimal set of peptides for epitope-based vaccines. Nucleic Acids Res,
    2009, 37, W617-W622
    [2] Pyomo - Optimization Modeling in Python. William E. Hart, Carl Laird,
    Jean-Paul Watson and David L. Woodruff. Springer, 2012.
    
    @author: Benjamin Schubert
'''

from __future__ import division
from coopr.pyomo import *
from coopr.opt import SolverFactory
from Fred2.core.Peptide import PeptideSet
from Fred2.core.Allele import AlleleSet, Allele


class EpitopeSelectionException(Exception):
    '''
        classdocs
    '''
    
    def __init__(self, expr, msg):
        '''
            Constructor
        '''
        self.expr = expr
        self.msg = msg



class EpitopeSelection(object):
    '''
        classdocs
    '''
    
    def __init__(self, peptideSet, alleleSet, k=10, solver="glpsol",verbosity=0):
        '''
            Constructor
            
            @param peptideSet: PeptideSet with all informations
            @param alleleSet: dictionary containing alleles for which the binding prediction exist and their frequency
            depending on a selected population
            @param k (int): the number of epitopes to select
            @param solver (String): the solver to be used (default glpsol)
        '''
        
        self.__solver = SolverFactory(solver)
        self.__verbosity=verbosity
        self.__changed = True
        self.__peptideSet = peptideSet
        self.__alleleProb = alleleSet
        self.__k = k
        self.__result = None
        
        #Variable, Set and Parameter preparation
        alleles_I = {}
        variations = []
        epi_var={}
        imm ={}
          
        for seq, peps in self.__peptideSet.items():
            for p in peps:
                binder=False
                #run through score list and generate imm and allele_I
                for pred,hla,score in p.get_metadata('affinity'):
                    if hla not in alleles_I:
                        alleles_I[hla] = set()
                    if score >= self.__alleleProb[hla].get_metadata('bindingThresh')[0]:
                        binder = True
                        alleles_I[hla].add(seq)
                    
                    imm[seq,hla] = score
                
                #run through the variation list and generate the variation list and epi_var dict
                for v in p.get_metadata('variation'):
                    vari = v
                    variations.append(vari)
                    epi_var[vari]=set()
                    if binder:
                        epi_var.setdefault(vari, set()).add(seq)
        
        variations = set(variations)
        model = ConcreteModel()
        
        #set definition
        model.Q = Set(initialize=variations)
        model.E = Set(initialize=self.__peptideSet.keys(),filter=
                      lambda model, epitope,peptide_set=self.__peptideSet,alleleSet=self.__alleleProb:
                      any( [True for predictor, allele, score in self.__peptideSet[epitope][0].get_metadata('affinity')
                            if score >= alleleSet[allele].get_metadata('bindingThresh')[0]] ))
        model.A = Set(initialize=self.__alleleProb.keys())
        model.E_var = Set(model.Q,initialize=lambda mode, v: epi_var[v])
        model.A_I = Set(model.A, initialize=lambda model,a: alleles_I[a])
                            
                            
        #parameter definition
        model.k = Param(initialize=self.__k, within=PositiveIntegers,mutable=True)
        model.p = Param(model.A, initialize=lambda model,a: self.__alleleProb[a].get_metadata('prob')[0])
                            
        cons = self.__calcEpitopeConservation()
        model.c = Param(model.E, initialize=lambda model,e: cons[e])
                            
        #threshold parameters
        model.i = Param(model.E, model.A, initialize=lambda model, e,a: imm[e,a])
        model.t_allele = Param(initialize=0, within=NonNegativeIntegers,mutable=True)
        model.t_var = Param(initialize=0, within=NonNegativeIntegers,mutable=True)
        model.t_c = Param(initialize=0.0, within=NonNegativeReals,mutable=True)
                            
        # Variable Definition
        model.x = Var(model.E, within=Binary)
        model.y = Var(model.A, within=Binary)
                            
        # Objective definition
        model.Obj = Objective(rule=lambda mode: sum( model.x[e] * sum( model.p[a]*model.i[e,a] for a in model.A ) for e in model.E ),sense = maximize)
                                                  
                                                  
       #Obligatory Constraint (number of selected epitopes)
       model.NofSelectedEpitopesCov = Constraint(rule=lambda model: sum(model.x[e] for e in model.E) <= model.k)
                                                  
       #optional constraints (in basic model they are disabled)
       model.IsAlleleCovConst = Constraint(model.A, rule=lambda model,a: sum(model.x[e] for e in model.A_I[a]) >= model.y[a])
       model.MinAlleleCovConst = Constraint(rule=lambda model: sum(model.y[a] for a in model.A) >= model.t_allele)
       model.AntigenCovConst = Constraint(model.Q, rule=lambda model,q: sum( model.x[e] for e in model.E_var[q] ) >= model.t_var)
       model.EpitopeConsConst = Constraint(model.E, rule=lambda model,e: (1 - model.c[e])*model.x[e] <= 1 - model.t_c)
                                                  
       #generate instance
       self.__instance = model.create()
       if self.__verbosity > 0:
           print "MODEL INSTANCE"
           self.__instance.pprint ()
                                                  
       #deactivate additional constraints and variables
       #params
       self.__instance.c.deactivate()
       self.__instance.t_c.deactivate()
       self.__instance.t_allele.deactivate()
       self.__instance.t_var.deactivate()
                                                  
       #constraints
       self.__instance.IsAlleleCovConst.deactivate()
       self.__instance.MinAlleleCovConst.deactivate()
       self.__instance.AntigenCovConst.deactivate()
       self.__instance.EpitopeConsConst.deactivate()
                                                  
       #variables
       self.__instance.y.deactivate()
    
    
    
   def __calcEpitopeConservation(self):
        '''
            Calculates the conservation of an epitopes.
            Conservation is defined as the fraction of variation the epitope can descent from
            
            @return returns a dictionary with key epitope-seq and value conservation (float)
            @rtype dict[string]=float
            @exception EpitopeSelectionException: if an peptide object does not have the metadata field 'variation'9
        '''
        var = []
        cons = {}
        for eps in self.__peptideSet.values():
            for ep in eps:
                try:
                    var.extend(ep.get_metadata('variation'))
                except:
                    raise OptiTopeException("internal function", "Peptide %s has no metadata field 'variation'."%(ep.sequence))
        total = float(len(set(var)))
        
        for seq,ep in self.__peptideSet.items():
            #print seq,len(set(ep[0].get_metadata('variation')))
            try:
                cons[seq] = len(set(ep[0].get_metadata('variation')))/total
            except:
                raise EpitopeSelectionException("internal function", "Peptide %s has no metadata field 'variation'."%(ep.sequence))
        return cons

    
    def set_k(self,k):
        '''
            sets the number of epitopes to select
            @param k: the number of epitopes
            @type k: int
            @exception OptiTopeException: if the input variable is not in the same domain as the parameter
        '''
        tmp = self.__instance.k.value
        try:
            getattr(self.__instance, str(self.__instance.k)).set_value(int(k))
            self.__changed = True
        except:
            self.__changed = False
            getattr(self.__instance, str(self.__instance.k)).set_value(int(tmp))
            raise EpitopeSelectionException('set_k', 'An error has occurred during setting parameter k. Please check if k is integer.')
    
    
    def activate_allele_coverage_const(self, minCoverage):
        '''
            enables the allele Coverage Constraint
            
            @param minCoverage (float): percentage of alleles which have to be covered
            @exception EpitopeSelectionException: if the input variable is not in the same domain as the parameter
            '''
        #parameter
        mc = self.__instance.t_allele.value
        
        try:
            self.__instance.t_allele.activate()
            getattr(self.__instance, str(self.__instance.t_allele)).set_value(
                                                                              int(len(self.__alleleProb)*minCoverage))
            #variables
            self.__instance.y.activate()
                                                                                                                                               
            #constraints
            self.__instance.IsAlleleCovConst.activate()
            self.__instance.MinAlleleCovConst.activate()
            self.__changed = True
        except:
            getattr(self.__instance, str(self.__instance.t_allele)).set_value(mc)
            self.__instance.t_allele.deactivate()
            self.__changed = False
            raise EpitopeSelectionException('activate_allele_coverage_const","An error occurred during activation of of the allele coverage constraint. '+
                                    'Please check your specified minimum coverage parameter to be in the range of 0.0 and 1.0.')
    
    def deactivate_allele_coverage_const(self):
        '''
            deactivates the allele coverage constraint
        '''
        
        #parameter
        self.__changed = True
        self.__instance.t_allele.deactivate()
        
        #variables
        self.__instance.y.deactivate()
        
        #constraints
        self.__instance.IsAlleleCovConst.deactivate()
        self.__instance.MinAlleleCovConst.deactivate()
    
    
    def activate_antigen_coverage_const(self, t_var):
        '''
            activates the variation coverage constraint
            @param t_var: the number of epitopes which have to come from each variation
            @type t_var: int
            @exception EpitopeSelectionException: if the input variable is not in the same domain as the parameter
            
        '''
        tmp = self.__instance.t_var.value
        try:
            self.__instance.t_var.activate()
            getattr(self.__instance, str(self.__instance.t_var)).set_value(int(t_var))
            self.__instance.AntigenCovConst.activate()
            self.__changed = True
        except:
            self.__instance.t_var.deactivate()
            getattr(self.__instance, str(self.__instance.t_var)).set_value(int(tmp))
            self.__instance.AntigenCovConst.activate()
            self.__changed = False
            raise EpitopeSelectionException("activate_antigen_coverage_const", "An error has occurred during activation of the coverage constraint. Please make sure your input is an integer.")
    
    def deactivate_antigen_coverage_const(self):
        '''
            deactivates the variation coverage constraint
        '''
        self.__changed = True
        self.__instance.t_var.deactivate()
        self.__instance.AntigenCovConst.deactivate()
    
    
    def activate_epitope_conservation_const(self,t_c):
        '''
            activates the epitope conservation constraint
            @param t_c: the percentage of conservation an epitope has to have.
            @type t_c: float [0.0,1.0]
            @exception EpitopeSelectionException: if the input variable is not in the same domain as the parameter
        '''
        if t_c < 0 or  t_c > 1:
            raise EpitopeSelectionException("activate_epitope_conservation_const", "The conservation threshold is out of its numerical bound. It has to be between 0.0 and 1.0.")
        
        self.__changed = True
        self.__instance.c.activate()
        self.__instance.t_c.activate()
        getattr(self.__instance, str(self.__instance.t_c)).set_value(float(t_c))
        self.__instance.EpitopeConsConst.activate()
    
    def deactivate_epitope_conservation_const(self):
        '''
            deactivates epitope conservation constraint
        '''
        self.__changed = True
        self.__instance.c.deactivate()
        self.__instance.t_c.deactivate()
        self.__instance.EpitopeConsConst.deactivate()
    
    def solve(self):
        '''
            invokes the selected solver and solves the problem.
            
            @return returns the optimal epitopes
            @rtype PeptideSet
            @exception EpitopeSelectionException: if the solver raised a problem or the solver is not accessible via the PATH environmental variable.
        '''
        if self.__changed:
            try:
                self.__instance.x.reset()
                self.__instance.y.reset()
                self.__instance.preprocess()
            
                res = self.__solver.solve(self.__instance)
                self.__instance.load(res)
                if self.__verbosity > 0:
                    res.write(num=1)
            
                if  str(res.Solution.status) != 'optimal':
                    print "Could not solve problem - "+str(res.Solution.status)+". Please check your settings"
                    sys.exit()
            
                self.__result = PeptideSet([self.__peptideSet[x][0] for x in self.__instance.x if self.__instance.x[x].value == 1.0])
                self.__result.log_metadata("obj",res.Solution.Objective.Value)
                if self.__instance.y.active:
                    self.__result.log_metadata("cov_alleles", AlleleSet([Allele(y) for y in self.__instance.y if self.__instance.y[y] == 1.0]))
            
                self.__changed = False
                return self.__result
            except:
                raise EpitopeSelectionException("solve", "An Error has occurred during solving. Please check your settings and if the solver is registered in PATH environment variable.")
        else:
            return self.__result    







