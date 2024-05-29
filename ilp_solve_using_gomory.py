import numpy as np
import copy
rfile = open("input_ilp.txt","r")
inputs = rfile.readlines()
objective = inputs[1]
i = 4
A = []
b = []
c = []
constraint_types = []
while inputs[i] != "\n":
    row = inputs[i]
    A.append([int(x) for x in row.split(", ")])
    i+=1
i += 2
while inputs[i] != "\n":
    b.append(int(inputs[i]))
    i += 1
i += 2
while inputs[i] != "\n":
    constraint_types.append(inputs[i].rstrip())
    i += 1
i += 2
c = [int(x) for x in inputs[i].split(", ")]

def to_standard(A,constraint_types):
    for m in range(len(A)):
        if constraint_types[m] == "<=":
            c.append(0)
            for i in range(len(A)):
                if (i == m):
                    A[i].append(1)
                else:
                    A[i].append(0)
        elif constraint_types[m] == ">=":
            c.append(0)
            for i in range(len(A)):
                if ( i == m):
                    A[i].append(-1)
                else:
                    A[i].append(0)
        else:
            continue

def simplex_help(A,b,c,bfs,basic_x):
    ans = {}
    tableau = np.array(A, dtype= float)
    c_b = [c[x] for x in basic_x]
    c_hat = np.array(np.array(c) - np.matrix(c_b).__matmul__(np.matrix(tableau)))[0]
    tableau_ini= copy.deepcopy(tableau)
    c_hat_in= copy.deepcopy(c_hat)
    bfs_in = []
    for e in basic_x:
        bfs_in.append(bfs[e])
    tableau_ini = np.c_[np.transpose(np.array(bfs_in)),tableau_ini]
    c_hat_in= np.insert(c_hat_in, 0, -1*(np.dot(np.transpose(np.array(c)), np.array(bfs))))
    ans["initial_tableau"] = np.r_[np.asmatrix(c_hat_in),tableau_ini]
    while True:
        pivot_col = 0
        pivot_row = 0
        for s in range(len(c_hat)):
            if c_hat[s] < 0:
                pivot_col += s
                break
        if round(c_hat[pivot_col],10) >= 0:
            tableau_fin= copy.deepcopy(tableau)
            ans["AB_invA"]= tableau
            c_hat_fin= copy.deepcopy(c_hat)
            bfs_fin=[]
            for e in basic_x:
                bfs_fin.append(bfs[e])
            tableau_fin = np.c_[np.transpose(np.array(bfs_fin)),tableau_fin]
            c_hat_fin= np.insert(c_hat_fin, 0, -1*(np.dot(np.transpose(np.array(c)), np.array(bfs))))
            ans["c_hat"] = c_hat_fin
            ans["x_b"] = bfs_fin
            ans["final_tableau"] = np.r_[np.asmatrix(c_hat_fin),tableau_fin]
            ans["solution_status"] = "optimal"
            ans["optimal_solution"] = bfs
            ans["basic_x"]= basic_x
            if(objective== "maximize\n"):
                ans["optimal_cost"]= np.dot(np.transpose(np.array(c)), np.array(bfs))*(-1)
            else:
                ans["optimal_cost"]= np.dot(np.transpose(np.array(c)), np.array(bfs))
            return ans
        
        thetas = 10**18
        
        for t in range(len(basic_x)):
            if round(tableau[t][pivot_col],10) > 0 and round(bfs[basic_x[t]]/tableau[t][pivot_col],10) < thetas:
                thetas = bfs[basic_x[t]]/tableau[t][pivot_col]
                pivot_row = t
        if round(tableau[pivot_row][pivot_col],10)<=0:
            tableau_fin= copy.deepcopy(tableau)
            ans["AB_invA"]= tableau
            c_hat_fin= copy.deepcopy(c_hat)
            bfs_fin=[]
            for e in basic_x:
                bfs_fin.append(bfs[e])
            tableau_fin = np.c_[np.transpose(np.array(bfs_fin)),tableau_fin]
            if objective == "maximize\n":
                c_hat_fin= np.insert(c_hat_fin, 0, -np.inf)
                ans["final_tableau"] = np.r_[np.asmatrix(c_hat_fin),tableau_fin]
                ans["solution_status"] = "unbounded"
                ans["optimal_cost"]= np.inf
            else:
                c_hat_fin= np.insert(c_hat_fin, 0, np.inf)
                ans["final_tableau"] = np.r_[np.asmatrix(c_hat_fin),tableau_fin]
                ans["solution_status"] = "unbounded"
                ans["optimal_cost"]= -np.inf
            return ans
        pivot = float(tableau[pivot_row][pivot_col])
        tableau[pivot_row] = tableau[pivot_row]/pivot
        bfs[basic_x[pivot_row]]= bfs[basic_x[pivot_row]]/pivot
        for h in range(len(tableau)):
            if h!=pivot_row:
                bfs[basic_x[h]] -= bfs[basic_x[pivot_row]]*tableau[h][pivot_col]
                tableau[h] = tableau[h] - tableau[h][pivot_col]*tableau[pivot_row]
        c_hat = c_hat - c_hat[pivot_col]*tableau[pivot_row]
        bfs[pivot_col]=bfs[basic_x[pivot_row]]
        bfs[basic_x[pivot_row]]=0
        basic_x[pivot_row] = pivot_col


def simplex_algo(A,b,c,constraint_types,objective):
    to_standard(A,constraint_types)
    ans = {}
    if objective == "maximize\n":
        for i in range(len(c)):
            c[i]= -1*c[i]

    aux_c = [0]*(len(A[0])) + [1]*len(A)
    for l in range(len(b)):
        if(b[l] < 0):
            b[l] *= -1
            A[l].map(lambda x: -1*x)
    aux_A = copy.deepcopy(A)
    
    for k in range(len(constraint_types)):
        for i in range(len(aux_A)):
                if ( i == k):
                    aux_A[i].append(1)
                else:
                    aux_A[i].append(0)
    basic_x = []
    bfs = [0]*(len(A[0]))
    for ind in range(len(aux_A[0])-len(A), len(aux_A[0])):
        basic_x.append(ind)
        bfs.append(b[ind-len(aux_A[0])+len(A)])
    aux_ans = simplex_help(aux_A,b,aux_c,bfs,basic_x)
    aux_bfs = aux_ans["optimal_solution"]
    cost = np.array(aux_bfs).dot(np.array(aux_c))
    if cost > 0:
        ans["solution_status"] = "infeasible"
        return ans
    else:
        tableau = aux_ans["AB_invA"]
        c_hat = aux_ans["final_tableau"][0][1:]
        aux_basic_x = aux_ans["basic_x"]
        i = 0
        while i < len(aux_basic_x):
            if aux_basic_x[i] >= len(A[0]):
                y_row = tableau[i]
                flag = False
                for j in range(len(A[0])):
                    if y_row[j] != 0:
                        flag = True
                        pivot_col = j
                        pivot_row = i
                        pivot = float(tableau[pivot_row][pivot_col])
                        tableau[pivot_row] = tableau[pivot_row]/pivot
                        aux_bfs[aux_basic_x[pivot_row]]= aux_bfs[aux_basic_x[pivot_row]]/pivot
                        for h in range(len(tableau)):
                            if h!=pivot_row:
                                aux_bfs[aux_basic_x[h]] -= aux_bfs[aux_basic_x[pivot_row]]*tableau[h][pivot_col]
                                tableau[h] = tableau[h] - tableau[h][pivot_col]*tableau[pivot_row]
                        c_hat = c_hat - c_hat[pivot_col]*tableau[pivot_row]
                        aux_bfs[pivot_col]=aux_bfs[aux_basic_x[pivot_row]]
                        aux_bfs[aux_basic_x[pivot_row]]=0
                        aux_basic_x[pivot_row] = pivot_col
                        break
                if not flag:
                    tableau = np.delete(tableau,i,axis=0)
                    aux_bfs[aux_basic_x[i]] = -5
                    aux_basic_x.pop(i)
                else:
                    i += 1
            else:
                i += 1
        for k in range(len(A[0])+len(A)-1,len(A[0])-1,-1):
            tableau = np.delete(tableau,k,axis=1)
        fin_ans = simplex_help(tableau,b,c,aux_bfs[:len(aux_bfs)-len(A)], aux_basic_x)
        return fin_ans

def dual_simplex(AB_invA, xb, chat,cost, basic_x):
    feasible= True
    while(np.min(xb)<0 and feasible==True): # When there is a negative value in the solution xb
        pivot_row= np.argmin(xb)
        max_val= -1*np.inf
        pivot_col=0
        for j in range(AB_invA.shape[1]):
            if AB_invA[pivot_row, j]<0:
                if np.round(-1*chat[j]/AB_invA[pivot_row, j],10)>np.round(max_val,10):
                    max_val= -1*chat[j]/AB_invA[pivot_row, j]
                    pivot_col= j
        if round(AB_invA[pivot_row, pivot_col], 10)>=0:
            feasible= False
        else:
            pivot_elem= AB_invA[pivot_row, pivot_col]
            count = 0
            cost= round(cost - xb[pivot_row]*round(chat[pivot_col]/pivot_elem, 10), 10)
            AB_invA[pivot_row, :]= np.round(AB_invA[pivot_row, :]/pivot_elem,10)
            chat= chat- AB_invA[pivot_row, :]*chat[pivot_col]
            for i in range(AB_invA.shape[0]):
                np.round(AB_invA,decimals=10,out=None)
                if i!=pivot_row:
                    count += 1
                    xb[i]= round(xb[i] - round(xb[pivot_row]*AB_invA[i, pivot_col]/pivot_elem, 10), 10)
                    AB_invA[i, :]= np.round(AB_invA[i, :] - AB_invA[pivot_row, :]*(round(AB_invA[i, pivot_col],10)),10)   
            xb[pivot_row]= round(xb[pivot_row]/pivot_elem,10)
            basic_x[pivot_row]= pivot_col
        #print(AB_invA)
        #print(xb)
        #print(chat)
        #print(cost)

    if feasible==False:
        return None,None,None,None,None,False
    else:
        return AB_invA,xb,chat,cost,basic_x,True

def gomory(A,b,c):
    ans = simplex_algo(A,b,c,constraint_types,objective)
    if(ans["solution_status"] == "optimal"):
        ans["optimal_solution"]=np.round(ans["optimal_solution"], 10)
        print("initial_solution: ",', '.join(map(str,ans["optimal_solution"])))
        A,b,c,cost,basic_x = ans["AB_invA"], ans["x_b"],ans["c_hat"][1:],-1*ans["optimal_cost"],ans["basic_x"]
        n = len(ans["optimal_solution"])
        #print(A)
        #print(b)
        frac=b-np.floor(b)
        num_cuts = 0
        while(round(np.max(frac),10)!=0):
            num_cuts += 1
            #print(num_cuts)
            m=A.shape[0]
            index = np.argmax(frac)
            r= A[index,:]
            r = np.floor(r)-r
            A =  np.vstack((A,r))
            col = [0 for _ in range(m+1)]
            col[m]=1
            col=np.matrix(col)
            A= np.hstack((A,np.transpose(col)))
            b = np.append(b,-frac[index])
            basic_x.append(A.shape[1])
            c= np.append(c,0)
            threshold = 1e-5
            A = np.where(np.abs(A) < threshold, 0, A)
            b= np.where(np.abs(b) < threshold,0, b)
            c= np.where(np.abs(c) < threshold,0, c)
            A,b,c,cost,basic_x,flag= dual_simplex(A,b,c,cost,basic_x)
            if not flag:
                print("final_solution:")
                print("solution_status: infeasible")
                print("number_of_cuts: ",num_cuts)
                print("optimal_value:")
                return
            frac=b-np.floor(b)
        x = [0]*n
        for i in range(len(basic_x)):
            if(basic_x[i] < n):
                x[basic_x[i]] = b[i]
        print("final_solution: ",', '.join(map(str,x)))
        print("solution_status: ",ans["solution_status"])
        print("number_of_cuts: ",num_cuts)
        print("optimal_value: ", cost)
    else:
        print("initial_solution:")
        print("final_solution:")
        print("solution_status: ",ans["solution_status"])
        print("number_of_cuts: 0")
        print("optimal_value:")
        return

gomory(A,b,c)