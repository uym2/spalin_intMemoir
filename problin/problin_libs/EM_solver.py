from problin_libs.ML_solver import *
from math import exp,log
import cvxpy as cp
from problin_libs import min_llh, conv_eps, eps
from scipy.sparse import diags
import timeit
import numpy as np

def log_sum_exp(numlist):
    # using log-trick to compute log(sum(exp(x) for x in numlist))
    # mitigate the problem of underflow
    maxx = max(numlist)
    result = maxx + log(sum([exp(x-maxx) for x in numlist]))
    return result

def pseudo_log(x):
    return log(x) if x>0 else min_llh

class EM_solver(ML_solver):
    def __init__(self,treelist,data,prior,params={'nu':0,'phi':0,'sigma':0}):
        super(EM_solver,self).__init__(treelist,data,prior,params)
        self.has_polytomy = False
        self.__mark_polytomies__(eps_len=0)
        self.num_edges = sum([len(list(tree.traverse_postorder())) for tree in self.trees])
    
    def __mark_polytomies__(self,eps_len=0):
        # polytomy_mark and resolve all polytomies in self.tree_obj
        self.has_polytomy = False
        self.num_polytomy_mark = 0
        for tree in self.trees:
            for node in tree.traverse_preorder():
                node.polytomy_mark = False
                if len(node.children) > 2:
                    self.has_polytomy = True
            tree.resolve_polytomies()
            for node in tree.traverse_preorder():
                if not hasattr(node,'polytomy_mark'):
                    node.polytomy_mark = True
                    self.num_polytomy_mark += 1
                    node.edge_length = eps_len  
                    self.has_polytomy = True              
    
    def x2brlen(self,x):
        i = 0
        for tree in self.trees:        
            for node in tree.traverse_postorder():
                if not node.polytomy_mark:
                    node.edge_length = x[i]
                    i += 1
    
    def ultrametric_constr(self,local_brlen_opt=True):
        N = self.num_edges-self.num_polytomy_mark
        if local_brlen_opt:
            for tree in self.trees:
                N -= len([node for node in tree.traverse_postorder() if node.mark_fixed])
        M = []
        b = []
        idx = 0
        for tree in self.trees: 
            for node in tree.traverse_postorder():
                if node.is_leaf():
                    node.constraint = [0.]*N
                    node.constant = node.edge_length if node.mark_fixed else 0
                else:
                    c1,c2 = node.children
                    m = [x-y for (x,y) in zip(c1.constraint,c2.constraint)]                
                    m_compl = [-x for x in m]
                    if sum([x!=0 for x in m]) > 0 and not (m in M or m_compl in M):
                        M.append(m)
                        b.append(c2.constant-c1.constant)
                    node.constraint = c1.constraint
                    if node.mark_fixed:
                        node.constant = c1.constant + node.edge_length
                    else:
                        node.constant = c1.constant
                if not node.polytomy_mark and not (node.mark_fixed and local_brlen_opt):    
                    node.constraint[idx] = 1
                    idx += 1
        for tree in self.trees[1:]:
            m = [x-y for (x,y) in zip(self.trees[0].root.constraint,tree.root.constraint)]
            M.append(m)
            b.append(tree.root.constant-self.trees[0].root.constant)                    
        return M,b

    def Estep_in_llh(self):
        # assume az_partition has been performed so each node has the attribute node.alpha
        # compute the inside llh, store in L0 and L1 of each node
        phi = self.params.phi
        nu = self.params.nu
        for i in range(len(self.trees)):
            tree = self.trees[i]
            numsites = self.numsites[i]
            Q = self.Q[i]
            M = self.charMatrices[i]
            for node in tree.traverse_postorder():
                p = exp(-node.edge_length)
                node.L0 = [0]*numsites # L0 and L1 are stored in log-scale
                node.L1 = [0]*numsites
                for site in range(numsites):   
                    q = Q[site][node.alpha[site]] if node.alpha[site] not in ['?','z'] else 1.0
                    if node.is_leaf():
                        if node.alpha[site] == "?":
                            masked_llh = pseudo_log(1-(1-phi)*p**nu) if M[node.label][site] == '?' else pseudo_log(1-p**nu)
                            #if self.charMatrices[node.label][site] == '?':
                            #    masked_llh = log(1-(1-phi)*p**nu) if (1-(1-phi)*p**nu)>0 else min_llh
                            #else:
                            #    masked_llh = log(1-p**nu) if (1-p**nu)>0 else min_llh           
                            node.L0[site] = node.L1[site] = masked_llh
                        elif node.alpha[site] == 'z':
                            node.L0[site] = (nu+1)*(-node.edge_length) + pseudo_log(1-phi) #if 1-phi>0 else min_llh
                            node.L1[site] = min_llh
                        else:
                            node.L0[site] = nu*(-node.edge_length) + pseudo_log(1-p) + pseudo_log(q) + pseudo_log(1-phi) #if (1-p)*q*(1-phi)>0 else min_llh
                            node.L1[site] = nu*(-node.edge_length) + pseudo_log(1-phi) #if (1-phi)>0 else min_llh
                    else:
                        C = node.children
                        l0 = l1 = 0
                        for c in C:
                            l0 += c.L0[site]
                            l1 += c.L1[site]
                        # note: l0_z, l0_alpha, and l0_masked are lists    
                        l0_z = [l0 + (nu+1)*(-node.edge_length)]
                        l0_alpha = [l1 + pseudo_log(1-p) + pseudo_log(q) + nu*(-node.edge_length)] if (node.alpha[site] != 'z') else []
                        l0_masked = [pseudo_log(1-p**nu)] if (node.alpha[site] == '?') else []
                        node.L0[site] = log_sum_exp(l0_z + l0_alpha + l0_masked)
                        if node.alpha[site] == 'z':
                            node.L1[site] = min_llh
                        elif node.alpha[site] != '?' or nu == 0 or p==1:
                            node.L1[site] = l1 + nu*(-node.edge_length) 
                        else:
                            node.L1[site] = log_sum_exp([l1 + nu*(-node.edge_length), pseudo_log(1-p**nu)])
    
    def lineage_llh(self):
        # override the function of the base class
        self.Estep_in_llh()
        llh = 0
        for tree in self.trees:
            llh += sum(tree.root.L0)
        return llh

    def Estep_out_llh(self):
        # assume binary tree
        # assume az_parition and Estep_in_llh have been performed 
        # so that all nodes have `alpha`, `L0` and `L1` attribues
        # output: add the attributes `out0` and `out1` to each node
        # where v.out0 = P(~D_v,v=0) and v.out1 = P(~D_v,v=-1) 
        for i in range(len(self.trees)):
            tree = self.trees[i]
            numsites = self.numsites[i]
            Q = self.Q[i]
            for v in tree.traverse_preorder():
                if v.is_root(): # base case
                    # Auxiliary components
                    v.A = [0]*numsites
                    #v.X = [-self.params.nu*v.edge_length + log(1-exp(-v.edge_length))]*self.numsites if self.params.nu*v.edge_length > 0 else [min_llh]*self.numsites
                    v.X = [-self.params.nu*v.edge_length + pseudo_log(1-exp(-v.edge_length))]*numsites
                    v.out_alpha = [{} for j in range(numsites)]
                    # Main components    
                    v.out0 = [-(1+self.params.nu)*v.edge_length]*numsites
                    #v.out1 = [log(1-exp(-v.edge_length*self.params.nu))]*self.numsites if self.params.nu*v.edge_length > 0 else [min_llh]*self.numsites
                    v.out1 = [pseudo_log(1-exp(-v.edge_length*self.params.nu))]*numsites
                else:
                    u = v.parent
                    w = None 
                    # get the sister
                    for x in u.children:
                        if x is not v:
                            w = x
                            break
                    # Auxiliary components
                    v.A = [None]*numsites
                    v.X = [None]*numsites
                    v.out_alpha = [{} for j in range(numsites)]
                    # Main components
                    v.out0 = [None]*numsites
                    v.out1 = [None]*numsites
                    for site in range(numsites): 
                        # compute out0: P(~D_v,v=0)
                        v.out0[site] = u.out0[site] + w.L0[site] - (1+self.params.nu)*v.edge_length
                        v.A[site] = u.out0[site] + w.L0[site] 
                        v.X[site] = -self.params.nu*v.edge_length + pseudo_log(1-exp(-v.edge_length)) + v.A[site]
                        # compute out1: P(~D_v,v=-1)
                        if w.alpha[site] == 'z': # z-branch
                            v.out1[site] = pseudo_log(1-exp(-v.edge_length*self.params.nu)) + v.A[site]
                        elif w.alpha[site] == '?': # masked branch
                            v.X[site] = log_sum_exp([v.X[site],u.X[site]+w.L1[site]-self.params.nu*v.edge_length])
                            p = 1-exp(-v.edge_length*self.params.nu) # if nu=0 then p=0
                            pl = pseudo_log(p)
                            v.out1[site] = u.out1[site] if self.params.nu == 0 else log_sum_exp([pl+v.A[site],pl+u.X[site]+w.L1[site],u.out1[site]])
                        else:
                            alpha0 = w.alpha[site] 
                            if alpha0 not in u.out_alpha[site]:
                                self.__out_alpha_up__(u,Q,site,alpha0)
                            B = u.out_alpha[site][alpha0] + self.params.nu*(-v.edge_length) + w.L1[site]  
                            C = v.A[site] - self.params.nu*v.edge_length + pseudo_log(1-exp(-v.edge_length)) + pseudo_log(Q[site][alpha0])
                            v.out_alpha[site][alpha0] = log_sum_exp([B,C])
                            v.X[site] = log_sum_exp([v.X[site],B])
                            v.out1[site] = pseudo_log(1-exp(-v.edge_length*self.params.nu)) + log_sum_exp([v.A[site],w.L1[site]+u.out_alpha[site][alpha0]])

    def __out_alpha_up__(self,node,Q,site,alpha0):
        # auxiliary function, shoudn't be called outside
        v = node
        path = []   
        # going up to until reaching either a z-branch or the root
        while not v.is_root():
            u = v.parent
            for x in u.children:
                if x is not v:
                    w = x
                    break
            if w.alpha[site] != '?' and w.alpha[site] != alpha0: # the branch above u is a z-branch
                v.out_alpha[site][alpha0] = v.A[site] - self.params.nu*v.edge_length + pseudo_log(1-exp(-v.edge_length)) + pseudo_log(Q[site][alpha0])
                break
            path.append(v)
            v = u
        if v.is_root(): # no z-branch found along the way
            v.out_alpha[site][alpha0] = -self.params.nu*v.edge_length + pseudo_log(1-exp(-v.edge_length)) + pseudo_log(Q[site][alpha0])             
        # going down to compute all the out_alpha along the path
        while path:
            v = path.pop()
            u = v.parent
            for x in u.children:
                if x is not v:
                    w = x
            B = u.out_alpha[site][alpha0] + self.params.nu*(-v.edge_length) + w.L1[site]  
            C = v.A[site] - self.params.nu*v.edge_length + pseudo_log(1-exp(-v.edge_length)) + pseudo_log(Q[site][alpha0])
            v.out_alpha[site][alpha0] = log_sum_exp([B,C])

    def Estep_posterior(self):
        # assume binary tree
        # assume az_parition, Estep_in_llh, and Estep_out_llh have been performed 
        # so that all nodes have `alpha`, `L0`, `L1`, `out0`, and `out1` attribues
        # output: add the attributes `S0-4` (refer to the paper for definitions; all S are NOT stored in log-scale)
        # and `post0` and `post1` to each node where v.post0 = log P(v=0|D) and v.post1 = log P(v=-1|D) 
        for i in range(len(self.trees)):
            tree = self.trees[i]
            numsites = self.numsites[i]   
            M = self.charMatrices[i]
            full_llh = tree.root.L0
            for v in tree.traverse_preorder():
                v.post0 = [None]*numsites
                v.post1 = [None]*numsites
                v.S0 = [None]*numsites
                v.S1 = [None]*numsites
                v.S2 = [None]*numsites
                v.S3 = [None]*numsites
                v.S4 = [None]*numsites
                for site in range(numsites):
                    # compute auxiliary values: v_in1 = log P(D_v|v=-1),v_in0 = log P(D_v|v=0), v_in_alpha = log P(D_v|v=alpha0)
                    v_in1 = 0 if v.alpha[site] == '?' else None
                    if v.is_leaf():
                            c = M[v.label][site]
                            if c == 0:
                                v_in0 = pseudo_log(1-self.params.phi)
                            else:
                                v_in0 = pseudo_log(self.params.phi) if (c == '?') else min_llh 
                    else:    
                        v1,v2 = v.children
                        v_in0 = v1.L0[site] + v2.L0[site]                 
                    # compute posterior
                    v.post0[site] = v_in0 + v.out0[site] - full_llh[site] if v_in0 is not None else min_llh
                    v.post1[site] = v_in1 + v.out1[site] - full_llh[site] if v_in1 is not None else min_llh               
                    # compute S (note that all S values are NOT in log-scale)
                    if v.alpha[site] == 'z': # z-branch
                        v.S0[site] = 1.0
                        v.S1[site] = v.S2[site] = v.S3[site] = v.S4[site] = 0.0
                    elif v.is_root():
                        v.S0[site] = exp(v_in0 + (1.0+self.params.nu)*(-v.edge_length) - v.L0[site])
                        v.S2[site] = 0.0 if v.alpha[site] != '?' else (1.0-exp(-self.params.nu*v.edge_length))/exp(v.L0[site])
                        v.S1[site] = 1.0-v.S0[site]-v.S2[site]
                        v.S3[site] = v.S4[site] = 0.0
                    else:
                        u = v.parent
                        v.S0[site] = exp(u.post0[site] + v_in0 + (1.0+self.params.nu)*(-v.edge_length) - v.L0[site])
                        if v.alpha[site] != '?':
                            v.S2[site] = 0.0
                            v.S4[site] = 0.0
                        else: # masked branch
                            v.S2[site] = exp(u.post0[site]-v.L0[site])*(1.0-exp(-self.params.nu*v.edge_length))
                            v.S4[site] = (1.0-exp(u.post0[site])-exp(u.post1[site]))*(1.0-exp(-self.params.nu*v.edge_length))/exp(v.L1[site])
                        v.S1[site] = exp(u.post0[site]) - v.S0[site] - v.S2[site] 
                        v.S3[site] = 1.0-v.S0[site]-v.S1[site]-exp(v.post1[site])

    def Estep(self):
        self.Estep_in_llh()
        self.Estep_out_llh()
        self.Estep_posterior()

    def Mstep(self,optimize_phi=True,optimize_nu=True,optimize_brlens=True,verbose=1,eps_nu=1e-5,eps_s=1e-6,ultra_constr=False,local_brlen_opt=True,brlen_solver='cvxpy'):
    # assume that Estep have been performed so that all nodes have S0-S4 attributes
    # output: optimize all parameters: branch lengths, phi, and nu
    # verbose level: 1 --> show all messages; 0 --> show minimal messages; -1 --> completely silent        
        #start_time = timeit.default_timer()
        if not optimize_phi:
            if verbose > 0:
                print("Fixing phi to " + str(self.params.phi))    
            phi_star = self.params.phi
        else:       
            if verbose > 0:
                print("Optimizing phi")
            R = []
            R_tilde = []
            for x in range(len(self.trees)):
                tree = self.trees[x]
                M = self.charMatrices[x]
                numsites = self.numsites[x]
                for i,v in enumerate(tree.traverse_leaves()):
                    R.append(sum([x != '?' for x in M[v.label]]))
                    R_tilde.append(sum([1-exp(p) for (j,p) in enumerate(v.post1) if M[v.label][j] == '?'])) 
            phi_star = sum(R_tilde)/(sum(R)+sum(R_tilde))
            if abs(phi_star) < 1/(numsites*len(M)):
                phi_star = 0
        # optimize nu and all branch lengths
        N = self.num_edges-self.num_polytomy_mark
        if local_brlen_opt:
            for tree in self.trees:
                N -= len([node for node in tree.traverse_postorder() if node.mark_fixed])
        S0 = np.zeros(N)
        S1 = np.zeros(N)
        S2 = np.zeros(N)
        S3 = np.zeros(N)
        S4 = np.zeros(N)
        d_ini = np.zeros(N)
        i = 0
        for i in range(len(self.trees)):
            tree = self.trees[i]
            numsites = self.numsites[i]
            for v in tree.traverse_postorder():
                if not v.polytomy_mark and not (v.mark_fixed and local_brlen_opt):    
                    s = [sum(v.S0),sum(v.S1),sum(v.S2),sum(v.S3),sum(v.S4)]
                    s = [max(eps_s,x) for x in s]
                    s = [x/sum(s)*numsites for x in s]
                    S0[i],S1[i],S2[i],S3[i],S4[i] = s
                    d_ini[i] = v.edge_length
                    i += 1

        def __optimize_brlen__(nu,verbose=False): # nu is a single number
            var_d = cp.Variable(N,nonneg=True) # the branch length variables
            C0 = -(nu+1)*S0.T @ var_d
            C1 = -nu*S1.T @ var_d + S1.T @ cp.log(1-cp.exp(-var_d)) 
            C2 = S2.T @ cp.log(1-cp.exp(-nu*var_d)) if sum(S2) > 0 and nu > eps_nu else 0 
            C3 = -nu*S3.T @ var_d
            C4 = S4.T @ cp.log(1-cp.exp(-nu*var_d)) if sum(S4) > 0 and nu > eps_nu else 0

            objective = cp.Maximize(C0+C1+C2+C3+C4)
            constraints = [np.zeros(N)+self.dmin <= var_d, var_d <= np.zeros(N)+self.dmax] 
            if ultra_constr:
                M,b = self.ultrametric_constr(local_brlen_opt=local_brlen_opt)
                constraints += [np.array(M) @ var_d == np.array(b)]
            prob = cp.Problem(objective,constraints)
            prob.solve(verbose=False,solver=cp.MOSEK)
            return var_d.value,prob.status
       
        def __optimize_brlen_scipy__(nu,verbose=False):
            def f(x): # objective function
                C0 = -(nu+1)*np.dot(S0,x)
                C1 = -nu*np.dot(S1,x) + np.dot(S1,np.log(1-np.exp(-x)))
                C2 = np.dot(S2,np.log(1-np.exp(-nu*x))) if sum(S2) > 0 and nu > eps_nu else 0
                C3 = -nu*np.dot(S3,x)
                C4 = np.dot(S4,np.log(1-np.exp(-nu*x))) if sum(S4) > 0 and nu > eps_nu else 0
                return -(C0+C1+C2+C3+C4)
            def g(x): # gradient
                dC0 = -(nu+1)*S0
                dC1 = -nu*S1 + S1*np.exp(-x)/(1-np.exp(-x))
                dC2 = nu*S2*np.exp(-nu*x)/(1-np.exp(-nu*x)) if sum(S2) > 0 and nu > eps_nu else 0
                dC3 = -nu*S3
                dC4 = nu*S4*np.exp(-nu*x)/(1-np.exp(-nu*x)) if sum(S4) > 0 and nu > eps_nu else 0
                return -(dC0 + dC1 + dC2 + dC3 + dC4)
            def h(x): # hessian
                hC0 = 0
                hC1 = S1*np.exp(-x)/(1-np.exp(-x))**2    
                hC2 = nu*nu*S2*np.exp(-nu*x)/(1-np.exp(-nu*x))**2 if sum(S2) >0 and nu>eps_nu else 0
                hC3 = 0
                hC4 = nu*nu*S4*np.exp(-nu*x)/(1-np.exp(-nu*x))**2 if sum(S4) >0 and nu>eps_nu else 0
                return diags(hC0+hC1+hC2+hC3+hC4)
            x0 = d_ini
            bounds = optimize.Bounds(np.zeros(N)+self.dmin,np.zeros(N)+self.dmax,keep_feasible=True)
            constraints = []
            if ultra_constr:
                M,b = self.ultrametric_constr(local_brlen_opt=local_brlen_opt)
                constraints.append(optimize.LinearConstraint(csr_matrix(M),[0]*len(M),[0]*len(M),keep_feasible=False))
            #out = optimize.minimize(f, x0, method="trust-constr",jac=g,hess=h,options={'disp':True,'verbose':1,'maxiter':50000}, bounds=bounds,constraints=constraints)    
            out = optimize.minimize(f, x0, method="SLSQP",jac=g,options={'disp':True,'maxiter':50000}, bounds=bounds,constraints=constraints)    
            status = "optimal" if out.success else out.message
            return out.x,status

        def __optimize_nu__(d): # d is a vector of all branch lengths
            var_nu = cp.Variable(1,nonneg=True) # the nu variable
            C0 = -(var_nu+1)*S0.T @ d
            C1 = -var_nu*S1.T @ d + S1.T @ cp.log(1-cp.exp(-d))
            C2 = S2.T @ cp.log(1-cp.exp(-var_nu*d)) if sum(S2) > 0 else 0
            C3 = -var_nu*S3.T @ d
            C4 = S4.T @ cp.log(1-cp.exp(-var_nu*d)) if sum(S4) > 0 else 0
            objective = cp.Maximize(C0+C1+C2+C3+C4)
            prob = cp.Problem(objective)
            prob.solve(verbose=False,solver=cp.MOSEK)
            #prob.solve(verbose=False,solver=cp.ECOS,max_iters=400)
            return var_nu.value[0],prob.status

        nIters = 1
        nu_star = self.params.nu
        for r in range(nIters):
            if not optimize_brlens:
                if verbose > 0:
                    print("Fixing branch lengths to the given values.")
                d_star = d_ini
                status_d = "optimal"
            else:    
                if verbose > 0:
                    print("Optimizing branch lengths. Current phi: " + str(phi_star) + ". Current nu:" + str(nu_star))
                try:
                    opt_func =  __optimize_brlen__ if brlen_solver == 'cvxpy' else  __optimize_brlen_scipy__
                    d_star,status_d = opt_func(nu_star,verbose=False)
                except:
                    d_star = d_ini
                    status_d = "failure"
                if status_d == "infeasible": # should only happen with local EM 
                    return False,"d_infeasible"
            if not optimize_nu:
                if verbose > 0:
                    print("Fixing nu to " + str(self.params.nu))
                nu_star = self.params.nu
                status_nu = "optimal"    
            else:    
                if verbose > 0:
                    print("Optimizing nu")
                try:
                    nu_star,status_nu = __optimize_nu__(d_star)                 
                except:
                    status_nu = "failure"
                #if status_nu != "optimal":
                #    return False,status_nu
        # place the optimal value back to params
        self.params.phi = phi_star
        self.params.nu = nu_star
        i = 0
        for tree in self.trees:
            for node in tree.traverse_postorder():
                if not node.polytomy_mark and not (node.mark_fixed and local_brlen_opt):    
                    node.edge_length = d_star[i]
                    i += 1
        success = (status_d == "optimal") and (status_nu == "optimal")
        if success:
            status = "optimal"
        else:    
            status = ""
            if status_d != "optimal":
                status += "failed_d"
            if status_nu != "optimal":
                status += ",failed_nu"
        return success,status
    
    def EM_optimization(self,verbose=1,optimize_phi=True,optimize_nu=True,optimize_brlens=True,ultra_constr=False,maxIter=1000):
        # assume that az_partition has been performed
        # optimize all parameters: branch lengths, phi, and nu
        # if optimize_phi is False, it is fixed to the original value in params.phi
        # the same for optimize_nu
        # caution: this function will modify params in place!
        # verbose level: 1 --> show all messages; 0 --> show minimal messages; -1 --> completely silent
        pre_llh = self.lineage_llh()
        if verbose >= 0:
            print("Initial phi: " + str(self.params.phi) + ". Initial nu: " + str(self.params.nu) + ". Initial nllh: " + str(-pre_llh))
        em_iter = 1
        converged = False
        while em_iter <= maxIter:
            if verbose > 0:
                print("Starting EM iter: " + str(em_iter))
                print("Estep")
            self.Estep()
            if verbose > 0:
                print("Mstep")
            m_success,status=self.Mstep(optimize_phi=optimize_phi,optimize_nu=optimize_nu,optimize_brlens=optimize_brlens,verbose=verbose,ultra_constr=ultra_constr,local_brlen_opt=True)
            if not m_success:
                if status == "d_infeasible": # should only happen with local EM
                    if verbose >= 0:
                        print("Warning: EM failed to optimize parameters in one Mstep due to infeasible constraints") 
                    return -pre_llh, em_iter,status
                elif verbose >= 0:    
                    print("Warning: EM failed to optimize parameters in one Mstep.")                
            curr_llh = self.lineage_llh()
            if verbose > 0:
                print("Finished EM iter: " + str(em_iter) + ". Current nllh: " + str(-curr_llh))
            if abs((curr_llh - pre_llh)/pre_llh) < conv_eps:
                converged = True
                break
            pre_llh = curr_llh
            em_iter += 1
        if not converged and verbose >= 0:
            print("Warning: exceeded maximum number of EM iterations (" + str(maxIter) + " iters)!")
        return -curr_llh, em_iter,status    
  
    def posterior_silence(self):
        self.Estep()
        for leaf in self.tree.traverse_leaves():
            print(leaf.label,[round(exp(x),2) for x in leaf.post1])

    def optimize_one(self,randseed,fixed_phi=None,fixed_nu=None,verbose=1,ultra_constr=False,optimize_brlens=True):
        # optimize using a specific initial point identified by the input randseed
        # verbose level: 1 --> show all messages; 0 --> show minimal messages; -1 --> completely silent
        seed(a=randseed)
        param_cats,ini_params = self.ini_all(fixed_phi=fixed_phi,fixed_nu=fixed_nu)
        x0 = []
        for p in param_cats:
            if p != 'brlens' or optimize_brlens:
                x0 = x0 + ini_params[p]
        self.x2params(x0,fixed_phi=fixed_phi,fixed_nu=fixed_nu,include_brlens=optimize_brlens)
        self.az_partition()
        nllh,em_iter,status = self.EM_optimization(verbose=verbose,optimize_phi=(fixed_phi is None),optimize_nu=(fixed_nu is None),ultra_constr=ultra_constr,optimize_brlens=optimize_brlens)
        if verbose >= 0:
            print("EM finished after " + str(em_iter) + " iterations.")
        return nllh,status
