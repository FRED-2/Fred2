
def spacer_design_c1(ei, ej, k, en, cn, cl_pssm, epi_pssms, cleav_pos, allele_prob, alpha,
                    thresh, solver, beta=None, options=""):
    """
        PRIVATE:
        internal spacer design for a pre-defined spacer length between two epitopes

        :param str ei: Start epitope
        :param str ej: End epitope
        :param int k: Length of spacer
        :param str options: solver options
        :return: Tuple of ei, ej, spacer (str), cleavage score, immunogenicity score
    """

    if k <= 0:
        seq = ei+ej
        i = len(ei)-cleav_pos
        g=len(ei)+k-cleav_pos
        c1=sum(cl_pssm[j][seq[i+j]] for j in xrange(cn))
        c2=sum(cl_pssm[j][seq[g+j]] for j in xrange(cn))
        return "",(c1+c2)/2, sum(prob*sum(max(sum(epi_pssms[j,seq[i+j],a] for j in xrange(en))+epi_pssms.get((-1,"con",a),0)-thresh[a],0)
                                                                    for i in xrange(len(seq)-en))
                                                                        for a,prob in allele_prob.iteritems()),c1,c2

    def normalize_pssm(p):
        max_p = -float("inf")
        min_p = float("inf")
        norm = {}
        for i,v in p.iteritems():
            max_tmp = max(v.values())
            min_tmp = min(v.values())
            if max_tmp > max_p:
                max_p = max_tmp
            if min_tmp < min_p:
                min_p = min_tmp
        for i,v in p.iteritems():
            for a,score in v.iteritems():
                norm.setdefault(i,{}).update({a:(score-min_p)/(max_p - min_p)})
        return norm

    cl_pssm_norm = normalize_pssm(cl_pssm)
    alphabet = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
    model = ConcreteModel()
    le = len(ei)+len(ej)+k
    neg_inf = -float("inf")
    #Sets
    ep = {}
    def sequence_set(model,i):
            if i < len(ei):
                return [ei[i]]
            elif i < len(ei)+k:
                return alphabet
            else:
                return [ej[i-len(ei)-k]]
    model.A = Set(initialize=allele_prob.keys())
    model.C = Set(initialize=range(cn))
    model.EN = Set(initialize=range(en))
    model.L = Set(initialize=range(le))
    model.Sigma = Set(initialize=alphabet)
    model.S = Set(model.L, initialize=sequence_set)
    model.AUX = Set(dimen=2, initialize=lambda model: [(i,a) for i in xrange(le) for a in model.S[i]])
    model.R = Set(initialize=range(le-(en-1)))

    #param
    model.f = Param(model.C, model.Sigma, initialize=lambda model,i,a: cl_pssm_norm[i][a])
    model.ci = Param(initialize=len(ei)-cleav_pos)
    model.cj = Param(initialize=len(ei)+k-cleav_pos)
    model.p = Param(model.A, initialize=lambda model, m: allele_prob[m])
    model.bi = Param(model.A, initialize=lambda model, m: epi_pssms.get((-1,"con",m),0))
    model.bc = Param(initialize=cl_pssm_norm.get(-1,{}).get("con",0))
    #epitope part
    model.i = Param(model.EN, model.Sigma,model.A, initialize=lambda model,i,a,m: epi_pssms[i,a,m])
    model.tau_epi = Param(initialize=10**6,mutable=True)
    model.tau_cleav = Param(initialize=-10**6, mutable=True)
    model.t_a = Param(model.A, initialize=lambda model, a: thresh.get(a, neg_inf))

    #Variables
    model.x = Var(model.AUX,domain=Binary)
    model.y = Var(model.R,model.A,domain=NonNegativeReals)

    #objective linear
    model.obj_cleav = Objective(rule=lambda model: (sum( model.f[i,a]*model.x[model.ci+i,a] for i in model.C for a in model.S[model.ci+i] )
                              + 0.0*sum(model.f[j,a]*model.x[model.cj+j,a] for j in model.C for a in model.S[model.cj+j])+1*model.bc),sense=maximize)

    model.obj_epi = Objective(rule=lambda model: sum(model.y[i,a]*model.p[a] for a in model.A
                                                     for i in model.R), sense=minimize)

    #model.obj_non_cleav = Objective(rule=lambda model: sum( model.f[i+j,a]*model.x[j+i,a] for i in model.L
    #                                                        for j in model.C
    #                                                            for a in model.S[i+j]
    #                                                                if i != model.ci and i != model.cj), sense=minimize)


    #constraints
    model.cons = Constraint(model.L,rule=lambda model, i: sum(model.x[i,a] for a in model.S[i]) == 1)

    model.max_imm_c = Constraint(model.R,model.A,rule=lambda model, i, m:
                                            model.y[i,m] >= sum( model.x[i+j,a]*model.i[j,a,m]
                                                                for j in model.EN for a in model.S[i+j])+model.bi[m]-model.t_a[m])
    ##neo-epitope constraint
    model.c_epi = Constraint(rule=lambda model:sum(model.y[i,a]*model.p[a] for a in model.A
                                                   for i in model.R) <= model.tau_epi)

    #cleavage constraint
    model.c_cleavage = Constraint(rule=lambda model: (sum( model.f[i,a]*model.x[model.ci+i,a] for i in model.C for a in model.S[model.ci+i] )
                              + 0*sum(model.f[j,a]*model.x[model.cj+j,a] for j in model.C for a in model.S[model.cj+j])+1*model.bc) >= model.tau_cleav)

    instance = model.create()
    solver = SolverFactory(solver, option=options)

    instance.obj_epi.deactivate()
    instance.c_epi.deactivate()
    instance.c_cleavage.deactivate()
    #instance.obj_non_cleav.deactivate()

    instance.preprocess()
    res = solver.solve(instance)  #, tee=True)

    if (res.solver.status == SolverStatus.ok) and (res.solver.termination_condition == TerminationCondition.optimal):
        instance.load(res)

        obj_cleav = instance.obj_cleav()

        instance.obj_cleav.deactivate()
        instance.obj_epi.activate()
        instance.c_cleavage.activate()

        #set bound of now inactive objective
        getattr(instance, "tau_cleav").set_value(alpha*obj_cleav)

        instance.preprocess()
        res2 = solver.solve(instance)#, tee=True)
        if (res2.solver.status == SolverStatus.ok) and (res2.solver.termination_condition == TerminationCondition.optimal):
            instance.load(res2)

            if beta is not None:
                obj_imm = instance.obj_epi()

                instance.obj_epi.deactivate()
                instance.obj_non_cleav.activate()
                instance.c_epi.activate()

                getattr(instance, "tau_epi").set_value(beta*obj_imm)
                instance.preprocess()
                res3 = solver.solve(instance)#, tee=True)
                if (res3.solver.status == SolverStatus.ok) and (res3.solver.termination_condition == TerminationCondition.optimal):
                    instance.load(res2)
                    ci = float(sum(cl_pssm[i][a]*instance.x[model.ci+i,a] for i in instance.C for a in instance.S[instance.ci+i] ))+cl_pssm.get(-1,{}).get("con",0)
                    cj = float(sum(cl_pssm[j][a]*instance.x[model.cj+j,a] for j in instance.C for a in instance.S[instance.cj+j]))+cl_pssm.get(-1,{}).get("con",0)

                    return "".join([a for i in xrange(len(ei), len(ei) + k) for a in instance.S[i] if
                            instance.x[i, a].value]), float(ci+cj)/2, instance.obj_epi(),float(ci),float(cj)
                else:
                    raise RuntimeError("Problem could not be solved. Please check your input.")
            else:
                ci = float(sum(cl_pssm[i][a]*instance.x[model.ci+i,a] for i in instance.C for a in instance.S[instance.ci+i] ))+cl_pssm.get(-1,{}).get("con",0)
                cj = float(sum(cl_pssm[j][a]*instance.x[model.cj+j,a] for j in instance.C for a in instance.S[instance.cj+j]))+cl_pssm.get(-1,{}).get("con",0)

                return "".join([a for i in xrange(len(ei), len(ei) + k) for a in instance.S[i] if
                            instance.x[i, a].value]), float(ci+cj)/2, instance.obj_epi(),float(ci),float(cj)
        else:
            raise RuntimeError("Problem could not be solved. Please check your input.")
    else:
        raise RuntimeError("Problem could not be solved. Please check your input.")


def spacer_design_c2(ei, ej, k, en, cn, cl_pssm, epi_pssms, cleav_pos, allele_prob, weight,
                    thresh, solver, options=""):
    """
        PRIVATE:
        internal spacer design for a pre-defined spacer length between two epitopes

        :param str ei: Start epitope
        :param str ej: End epitope
        :param int k: Length of spacer
        :param str options: solver options
        :return: Tuple of ei, ej, spacer (str), cleavage score, immunogenicity score
    """
    if k <= 0:
        seq = ei+ej
        i = len(ei)-cleav_pos
        g=len(ei)+k-cleav_pos
        c1=sum(cl_pssm[j][seq[i+j]] for j in xrange(cn))
        c2=sum(cl_pssm[j][seq[g+j]] for j in xrange(cn))
        return "",(c1+c2)/2, sum(prob*sum(max(sum(epi_pssms[j,seq[i+j],a] for j in xrange(en))+epi_pssms.get((-1,"con",a),0)-thresh[a],0)
                                                                    for i in xrange(len(seq)-en))
                                                                        for a,prob in allele_prob.iteritems()),c1,c2

    def normalize_pssm(p):
        max_p = -float("inf")
        min_p = float("inf")
        norm = {}
        for i,v in p.iteritems():
            max_tmp = max(v.values())
            min_tmp = min(v.values())
            if max_tmp > max_p:
                max_p = max_tmp
            if min_tmp < min_p:
                min_p = min_tmp
        for i,v in p.iteritems():
            for a,score in v.iteritems():
                norm.setdefault(i,{}).update({a:(score-min_p)/(max_p - min_p)})
        return norm

    cl_pssm_norm = normalize_pssm(cl_pssm)
    alphabet = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
    model = ConcreteModel()
    le = len(ei)+len(ej)+k
    neg_inf = -float("inf")
    #Sets
    ep = {}
    def sequence_set(model,i):
            if i < len(ei):
                return [ei[i]]
            elif i < len(ei)+k:
                return alphabet
            else:
                return [ej[i-len(ei)-k]]
    model.A = Set(initialize=allele_prob.keys())
    model.C = Set(initialize=range(cn))
    model.EN = Set(initialize=range(en))
    model.L = Set(initialize=range(le))
    model.Sigma = Set(initialize=alphabet)
    model.S = Set(model.L, initialize=sequence_set)
    model.AUX = Set(dimen=2, initialize=lambda model: [(i,a) for i in xrange(le) for a in model.S[i]])
    model.R = Set(initialize=range(le-(en-1)))

    #param
    model.f = Param(model.C, model.Sigma, initialize=lambda model,i,a: cl_pssm_norm[i][a])
    model.ci = Param(initialize=len(ei)-cleav_pos)
    model.cj = Param(initialize=len(ei)+k-cleav_pos)
    model.p = Param(model.A, initialize=lambda model, m: allele_prob[m])
    model.bi = Param(model.A, initialize=lambda model, m: epi_pssms.get((-1,"con",m),0))
    model.bc = Param(initialize=cl_pssm_norm.get(-1,{}).get("con",0))
    #epitope part
    model.i = Param(model.EN, model.Sigma,model.A, initialize=lambda model,i,a,m: epi_pssms[i,a,m])
    model.tau_epi = Param(initialize=10**6,mutable=True)
    model.tau_cleav = Param(initialize=-10**6, mutable=True)
    model.t_a = Param(model.A, initialize=lambda model, a: thresh.get(a, neg_inf))

    #Variables
    model.x = Var(model.AUX,domain=Binary)
    model.y = Var(model.R,model.A,domain=NonNegativeReals)

    #objective linear
    model.obj_cleav = Objective(rule=lambda model: (0*sum( model.f[i,a]*model.x[model.ci+i,a] for i in model.C for a in model.S[model.ci+i] )
                              + sum(model.f[j,a]*model.x[model.cj+j,a] for j in model.C for a in model.S[model.cj+j])+1*model.bc),sense=maximize)

    model.obj_epi = Objective(rule=lambda model: sum(model.y[i,a]*model.p[a] for a in model.A
                                                     for i in model.R), sense=minimize)

    model.max_imm_c = Constraint(model.R,model.A,rule=lambda model, i, m:
                                            model.y[i,m] >= sum( model.x[i+j,a]*model.i[j,a,m]
                                                                for j in model.EN for a in model.S[i+j])+model.bi[m]-model.t_a[m])

    #constraints
    model.cons = Constraint(model.L,rule=lambda model, i: sum(model.x[i,a] for a in model.S[i]) == 1)


    ##neo-epitope constraint
    model.c_epi = Constraint(rule=lambda model:sum(model.y[i,a]*model.p[a] for a in model.A
                                                   for i in model.R) <= model.tau_epi)

    #cleavage constraint
    model.c_cleavage = Constraint(rule=lambda model: (0*sum( model.f[i,a]*model.x[model.ci+i,a] for i in model.C for a in model.S[model.ci+i] )
                              + sum(model.f[j,a]*model.x[model.cj+j,a] for j in model.C for a in model.S[model.cj+j])+1*model.bc) >= model.tau_cleav)

    instance = model.create()
    solver = SolverFactory(solver, option=options)

    instance.obj_epi.deactivate()
    instance.c_epi.deactivate()
    instance.c_cleavage.deactivate()

    instance.preprocess()
    res = solver.solve(instance)  #, tee=True)

    if (res.solver.status == SolverStatus.ok) and (res.solver.termination_condition == TerminationCondition.optimal):
        instance.load(res)

        obj_cleav = instance.obj_cleav()

        instance.obj_cleav.deactivate()
        instance.obj_epi.activate()
        instance.c_cleavage.activate()

        #set bound of now inactive objective
        getattr(instance, "tau_cleav").set_value(weight*obj_cleav)

        instance.preprocess()
        res2 = solver.solve(instance)#, tee=True)
        if (res2.solver.status == SolverStatus.ok) and (
            res2.solver.termination_condition == TerminationCondition.optimal):
            instance.load(res2)
            ci = float(sum(cl_pssm[i][a]*instance.x[model.ci+i,a] for i in instance.C for a in instance.S[instance.ci+i] ))+cl_pssm.get(-1,{}).get("con",0)
            cj = float(sum(cl_pssm[j][a]*instance.x[model.cj+j,a] for j in instance.C for a in instance.S[instance.cj+j]))+cl_pssm.get(-1,{}).get("con",0)

            return "".join([a for i in xrange(len(ei), len(ei) + k) for a in instance.S[i] if
                            instance.x[i, a].value]), float(ci+cj)/2, instance.obj_epi(),float(ci),float(cj)
        else:
            raise RuntimeError("Problem could not be solved. Please check your input.")
    else:
        raise RuntimeError("Problem could not be solved. Please check your input.")


def spacer_design_cleav(ei, ej, k, en, cn, cl_pssm, epi_pssms, cleav_pos, allele_prob, weight,
                    thresh, solver, options=""):
    """
        PRIVATE:
        internal spacer design for a pre-defined spacer length between two epitopes

        :param str ei: Start epitope
        :param str ej: End epitope
        :param int k: Length of spacer
        :param str options: solver options
        :return: Tuple of ei, ej, spacer (str), cleavage score, immunogenicity score
    """
    if k <= 0:
        seq = ei+ej
        i = len(ei)-cleav_pos
        g=len(ei)+k-cleav_pos
        c1=sum(cl_pssm[j][seq[i+j]] for j in xrange(cn))
        c2=sum(cl_pssm[j][seq[g+j]] for j in xrange(cn))
        return "",(c1+c2)/2, sum(prob*sum(max(sum(epi_pssms[j,seq[i+j],a] for j in xrange(en))+epi_pssms.get((-1,"con",a),0)-thresh[a],0)
                                                                    for i in xrange(len(seq)-en))
                                                                        for a,prob in allele_prob.iteritems()),c1,c2

    def normalize_pssm(p):
        max_p = -float("inf")
        min_p = float("inf")
        norm = {}
        for i,v in p.iteritems():
            max_tmp = max(v.values())
            min_tmp = min(v.values())
            if max_tmp > max_p:
                max_p = max_tmp
            if min_tmp < min_p:
                min_p = min_tmp
        for i,v in p.iteritems():
            for a,score in v.iteritems():
                norm.setdefault(i,{}).update({a:(score-min_p)/(max_p - min_p)})
        return norm

    cl_pssm_norm = normalize_pssm(cl_pssm)
    alphabet = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
    model = ConcreteModel()
    le = len(ei)+len(ej)+k
    neg_inf = -float("inf")
    #Sets
    ep = {}
    def sequence_set(model,i):
            if i < len(ei):
                return [ei[i]]
            elif i < len(ei)+k:
                return alphabet
            else:
                return [ej[i-len(ei)-k]]
    model.A = Set(initialize=allele_prob.keys())
    model.C = Set(initialize=range(cn))
    model.EN = Set(initialize=range(en))
    model.L = Set(initialize=range(le))
    model.Sigma = Set(initialize=alphabet)
    model.S = Set(model.L, initialize=sequence_set)
    model.AUX = Set(dimen=2, initialize=lambda model: [(i,a) for i in xrange(le) for a in model.S[i]])
    model.R = Set(initialize=range(le-(en-1)))

    #param
    model.f = Param(model.C, model.Sigma, initialize=lambda model,i,a: cl_pssm_norm[i][a])
    model.ci = Param(initialize=len(ei)-cleav_pos)
    model.cj = Param(initialize=len(ei)+k-cleav_pos)
    model.p = Param(model.A, initialize=lambda model, m: allele_prob[m])
    model.bi = Param(model.A, initialize=lambda model, m: epi_pssms.get((-1,"con",m),0))
    model.bc = Param(initialize=cl_pssm_norm.get(-1,{}).get("con",0))
    #epitope part
    model.i = Param(model.EN, model.Sigma,model.A, initialize=lambda model,i,a,m: epi_pssms[i,a,m])
    model.tau_epi = Param(initialize=10**6,mutable=True)
    model.tau_cleav = Param(initialize=-10**6, mutable=True)
    model.t_a = Param(model.A, initialize=lambda model, a: thresh.get(a, neg_inf))

    #Variables
    model.x = Var(model.AUX,domain=Binary)
    model.y = Var(model.R,model.A,domain=NonNegativeReals)

    #objective linear
    model.obj_cleav = Objective(rule=lambda model: 0.5*(sum( model.f[i,a]*model.x[model.ci+i,a] for i in model.C for a in model.S[model.ci+i] )
                              + sum(model.f[j,a]*model.x[model.cj+j,a] for j in model.C for a in model.S[model.cj+j])+2*model.bc),sense=maximize)

    model.obj_epi = Objective(rule=lambda model: sum(model.y[i,a]*model.p[a] for a in model.A
                                                     for i in model.R), sense=minimize)

    model.max_imm_c = Constraint(model.R,model.A,rule=lambda model, i, m:
                                            model.y[i,m] >= sum( model.x[i+j,a]*model.i[j,a,m]
                                                                for j in model.EN for a in model.S[i+j])+model.bi[m]-model.t_a[m])

    #constraints
    model.cons = Constraint(model.L,rule=lambda model, i: sum(model.x[i,a] for a in model.S[i]) == 1)


    ##neo-epitope constraint
    model.c_epi = Constraint(rule=lambda model:sum(model.y[i,a]*model.p[a] for a in model.A
                                                   for i in model.R) <= model.tau_epi)

    #cleavage constraint
    model.c_cleavage = Constraint(rule=lambda model: 0.5*(sum( model.f[i,a]*model.x[model.ci+i,a] for i in model.C for a in model.S[model.ci+i] )
                              + sum(model.f[j,a]*model.x[model.cj+j,a] for j in model.C for a in model.S[model.cj+j])+2*model.bc) >= model.tau_cleav)

    instance = model.create()
    solver = SolverFactory(solver, option=options)

    instance.obj_epi.deactivate()
    instance.c_epi.deactivate()
    instance.c_cleavage.deactivate()

    instance.preprocess()
    res = solver.solve(instance)  #, tee=True)
    instance.load(res)


    spacer = "".join([a for i in xrange(len(ei), len(ei) + k) for a in instance.S[i] if
                            instance.x[i, a].value])
    seq = ei+spacer+ej
    ci = float(sum(cl_pssm[i][a]*instance.x[model.ci+i,a] for i in instance.C for a in instance.S[instance.ci+i] ))+cl_pssm.get(-1,{}).get("con",0)
    cj = float(sum(cl_pssm[j][a]*instance.x[model.cj+j,a] for j in instance.C for a in instance.S[instance.cj+j]))+cl_pssm.get(-1,{}).get("con",0)
    imm = sum(prob*sum(max(sum(epi_pssms[j,seq[i+j],a] for j in xrange(en))+epi_pssms.get((-1,"con",a),0)-thresh[a],0)
                                                                    for i in xrange(len(seq)-en))
                                                                        for a,prob in allele_prob.iteritems())
    return spacer, float(ci+cj)/2, imm,float(ci),float(cj)



def spacer_design_imm(ei, ej, k, en, cn, cl_pssm, epi_pssms, cleav_pos, allele_prob, weight,
                    thresh, solver, options=""):
    """
        PRIVATE:
        internal spacer design for a pre-defined spacer length between two epitopes

        :param str ei: Start epitope
        :param str ej: End epitope
        :param int k: Length of spacer
        :param str options: solver options
        :return: Tuple of ei, ej, spacer (str), cleavage score, immunogenicity score
    """
    if k <= 0:
        seq = ei+ej
        i = len(ei)-cleav_pos
        g=len(ei)+k-cleav_pos
        c1=sum(cl_pssm[j][seq[i+j]] for j in xrange(cn))
        c2=sum(cl_pssm[j][seq[g+j]] for j in xrange(cn))
        return "",(c1+c2)/2, sum(prob*sum(max(sum(epi_pssms[j,seq[i+j],a] for j in xrange(en))+epi_pssms.get((-1,"con",a),0)-thresh[a],0)
                                                                    for i in xrange(len(seq)-en))
                                                                        for a,prob in allele_prob.iteritems()),c1,c2

    def normalize_pssm(p):
        max_p = -float("inf")
        min_p = float("inf")
        norm = {}
        for i,v in p.iteritems():
            max_tmp = max(v.values())
            min_tmp = min(v.values())
            if max_tmp > max_p:
                max_p = max_tmp
            if min_tmp < min_p:
                min_p = min_tmp
        for i,v in p.iteritems():
            for a,score in v.iteritems():
                norm.setdefault(i,{}).update({a:(score-min_p)/(max_p - min_p)})
        return norm

    cl_pssm_norm = normalize_pssm(cl_pssm)
    alphabet = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
    model = ConcreteModel()
    le = len(ei)+len(ej)+k
    neg_inf = -float("inf")
    #Sets
    ep = {}
    def sequence_set(model,i):
            if i < len(ei):
                return [ei[i]]
            elif i < len(ei)+k:
                return alphabet
            else:
                return [ej[i-len(ei)-k]]
    model.A = Set(initialize=allele_prob.keys())
    model.C = Set(initialize=range(cn))
    model.EN = Set(initialize=range(en))
    model.L = Set(initialize=range(le))
    model.Sigma = Set(initialize=alphabet)
    model.S = Set(model.L, initialize=sequence_set)
    model.AUX = Set(dimen=2, initialize=lambda model: [(i,a) for i in xrange(le) for a in model.S[i]])
    model.R = Set(initialize=range(le-(en-1)))

    #param
    model.f = Param(model.C, model.Sigma, initialize=lambda model,i,a: cl_pssm_norm[i][a])
    model.ci = Param(initialize=len(ei)-cleav_pos)
    model.cj = Param(initialize=len(ei)+k-cleav_pos)
    model.p = Param(model.A, initialize=lambda model, m: allele_prob[m])
    model.bi = Param(model.A, initialize=lambda model, m: epi_pssms.get((-1,"con",m),0))
    model.bc = Param(initialize=cl_pssm_norm.get(-1,{}).get("con",0))
    #epitope part
    model.i = Param(model.EN, model.Sigma,model.A, initialize=lambda model,i,a,m: epi_pssms[i,a,m])
    model.tau_epi = Param(initialize=10**6,mutable=True)
    model.tau_cleav = Param(initialize=-10**6, mutable=True)
    model.t_a = Param(model.A, initialize=lambda model, a: thresh.get(a, neg_inf))

    #Variables
    model.x = Var(model.AUX,domain=Binary)
    model.y = Var(model.R,model.A,domain=NonNegativeReals)

    #objective linear

    model.obj_epi = Objective(rule=lambda model: sum(model.y[i,a]*model.p[a] for a in model.A
                                                     for i in model.R), sense=minimize)

    model.max_imm_c = Constraint(model.R,model.A,rule=lambda model, i, m:
                                            model.y[i,m] >= sum( model.x[i+j,a]*model.i[j,a,m]
                                                                for j in model.EN for a in model.S[i+j])+model.bi[m]-model.t_a[m])

    #constraints
    model.cons = Constraint(model.L,rule=lambda model, i: sum(model.x[i,a] for a in model.S[i]) == 1)



    instance = model.create()
    solver = SolverFactory(solver, option=options)


    instance.preprocess()
    res = solver.solve(instance)  #, tee=True)
    instance.load(res)
    
    spacer = "".join([a for i in xrange(len(ei), len(ei) + k) for a in instance.S[i] if
                            instance.x[i, a].value])
    print "IMM, ",spacer
    seq = ei+spacer+ej
    ci = float(sum(cl_pssm[i][a]*instance.x[model.ci+i,a] for i in instance.C for a in instance.S[instance.ci+i] ))+cl_pssm.get(-1,{}).get("con",0)
    cj = float(sum(cl_pssm[j][a]*instance.x[model.cj+j,a] for j in instance.C for a in instance.S[instance.cj+j]))+cl_pssm.get(-1,{}).get("con",0)
    imm = sum(prob*sum(max(sum(epi_pssms[j,seq[i+j],a] for j in xrange(en))+epi_pssms.get((-1,"con",a),0)-thresh[a],0)
                                                                    for i in xrange(len(seq)-en))
                                                                        for a,prob in allele_prob.iteritems())
    return spacer, float(ci+cj)/2, imm,float(ci),float(cj)