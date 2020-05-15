from gurobipy import *
import pandas as pd
import numpy as np

try:    
    # Alines: augmented branch-to-node matrix
    A = pd.read_csv("A.txt",header=None,sep=' ')
    # Blines: augmented node-to-node susceptance matrix
    Bd = pd.read_csv("Bd.txt",header=None,sep=' ')   
    generators = pd.read_csv("generators.txt",header=None,sep=' ')
    loads = pd.read_csv("loads.txt",header=None,sep=' ')
    costs = pd.read_csv("costs.txt",header=None,sep=' ')
    batteries = pd.read_csv("battery.txt",header=None,sep=' ')   
    uncon_load = pd.read_csv("uncontrollable_loads.txt",header=None,sep=' ') 
    uncon_gen = pd.read_csv("uncontrollable_generations.txt",header=None,sep=' ')  
    Llimits = pd.read_csv("line_limits.txt",header=None,sep=' ') 

    L=len(A)
    N=len(A.columns)
    Nt=N+1
    T=len(loads)
    timeps = np.arange(T)    
    timepsxtra = np.arange(T+1)
    timepsdef = np.arange(T-1)

    nodes = np.arange(1,Nt)
    lines = np.arange(1,L+1)

    genminuptimes = generators.iloc[1:,3]
    loadminuptimes = loads.iloc[1:,3]

    genminuptime_slack=generators.iloc[0,3]
    loadminuptime_slack=loads.iloc[0,3]

    C = np.zeros(shape=(len(timeps)+1,len(nodes)+1,3))    
    for n in nodes:
        for k in range(3*T):
            C[k//3,n,k%3] = float(costs.iloc[n,k])

    Cs = np.zeros(shape=(len(timeps)+1, 3))
    for k in range(3*T):
        Cs[k//3,k%3]=float(costs.iloc[0,k])

    BatteryLimits = batteries.iloc[:N,]
    BatteryLimits.index = range(1,Nt)

    ftemp = np.empty(shape=(len(timeps),len(lines)),dtype=gurobipy.LinExpr)

    Btemp = [np.sum(Bd.iloc[l-1,:]*A.iloc[:,n-1])for l in lines for n in nodes]
    Btemp = pd.DataFrame(np.array(Btemp).reshape(5,4))

    B = [np.sum(A.iloc[:,n-1]*Btemp.iloc[:,m-1])for n in nodes for m in nodes]
    B = pd.DataFrame(np.array(B).reshape(4,4)) 

    b0 = [-np.sum(B.iloc[i,:]) for i in B.columns]

    M=1000000000000

    m=Model('powerflow')

    y=m.addVars(2,timeps,nodes,vtype=GRB.BINARY)
    z=m.addVars(2,timeps,nodes,vtype=GRB.BINARY)
    yslack=m.addVars(2,timeps,vtype=GRB.BINARY)
    zslack=m.addVars(2,timeps,vtype=GRB.BINARY)

    powergen=m.addVars(timeps, nodes,lb=0.0, name="controllable generations")
    powerload=m.addVars(timeps, nodes, lb=0.0,ub=GRB.INFINITY, name="controllable loads")
    uncpowergen=m.addVars(timeps, nodes,lb=0.0, name="uncontrollable generations")
    uncpowerload=m.addVars(timeps, nodes, lb=0.0,ub=GRB.INFINITY, name="uncontrollable loads")
    uncpowergen_slack=m.addVars(timeps, lb=0.0, name="uncontrollable generations slack bus")
    uncpowerload_slack=m.addVars(timeps, lb=0.0,ub=GRB.INFINITY, name="uncontrollable loads slack bus")
    powerbat=m.addVars(timeps, nodes, lb=-GRB.INFINITY,ub=GRB.INFINITY, name="battery charging power")
    enerbat=m.addVars(timepsxtra, nodes, lb=0.0,ub=GRB.INFINITY, name="battery energy")
    powerbat_slack=m.addVars(timeps, lb=-GRB.INFINITY,ub=GRB.INFINITY, name="battery charging power slack bus")
    enerbat_slack=m.addVars(timepsxtra, lb=0.0,ub=GRB.INFINITY, name="battery energy slack")

    powernet=m.addVars(timeps, nodes,lb=-GRB.INFINITY,ub=GRB.INFINITY, name="net")
    phaseangles=m.addVars(timeps, nodes,lb=-40.0, ub=40.0, name="phase angles")
    powergen_slack=m.addVars(timeps, name="controllable generation slack")
    powerload_slack=m.addVars(timeps, lb=0.0,ub=GRB.INFINITY, name="controllable load slack")
    powernet_slack=m.addVars(timeps, lb=-GRB.INFINITY, ub=GRB.INFINITY, name="net slack")

    powerflow=m.addVars(timeps, lines, lb=-GRB.INFINITY, ub=GRB.INFINITY, name="powerflow")

    m.addConstrs(
        uncpowerload[i, k]==float(uncon_load[k][i]) for i in timeps for k in nodes)

    m.addConstrs(
        uncpowerload_slack[i]==float(uncon_load[0][i]) for i in timeps)

    m.addConstrs(
        uncpowergen[i, k]==float(uncon_gen[k][i]) for i in timeps for k in nodes)

    m.addConstrs(
        uncpowergen_slack[i]==float(uncon_gen[0][i]) for i in timeps)  

    m.addConstrs(
        powergen[i, k]>=float(generators[0][k]) for i in timeps for k in nodes)

    m.addConstrs(
        powergen[i, k]<=float(generators[1][k]) for i in timeps for k in nodes)

    m.addConstrs(
        powergen_slack[i]>=float(generators[0][0]) for i in timeps)

    m.addConstrs(
        powergen_slack[i]<=float(generators[1][0]) for i in timeps)

    m.addConstrs(
        powerload[i, k]>=float(loads[0][k]) for i in timeps for k in nodes)

    m.addConstrs(
        powerload[i, k]<=float(loads[1][k]) for i in timeps for k in nodes)

    m.addConstrs(
        powerload_slack[i]>=float(loads[0][0]) for i in timeps)

    m.addConstrs(
        powerload_slack[i]<=float(loads[1][0]) for i in timeps)

    m.addConstrs(
        powergen[i, k]-powergen[i+1, k]<=float(generators[2][k]) for i in timepsdef for k in nodes)    
    # Ramping, 3rd column ramping limit 
    m.addConstrs(
        powergen_slack[i]-powergen_slack[i+1]<=float(generators[2][0]) for i in timepsdef)    

    m.addConstrs(
        powerload[i, k]-powerload[i+1, k]<=float(loads[2][k]) for i in timepsdef for k in nodes)    

    m.addConstrs(
        powerload_slack[i]-powerload_slack[i+1]<=float(loads[2][0]) for i in timepsdef) 

    m.addConstrs(
        powerbat[i, k]<=float(batteries[0][k]) for i in timeps for k in nodes)

    m.addConstrs(
        powerbat_slack[i]<=float(batteries[0][0]) for i in timeps)

    m.addConstrs(
        powerbat[i, k]>=-float(batteries[1][k]) for i in timeps for k in nodes)

    m.addConstrs(
        powerbat_slack[i]>=-float(batteries[1][0]) for i in timeps)

    m.addConstrs(
        enerbat[i, k]<=float(batteries[2][k]) for i in timepsxtra for k in nodes)

    m.addConstrs(
        enerbat_slack[i]<=float(batteries[2][0]) for i in timepsxtra)

    m.addConstrs(
        enerbat[0, k]==float(batteries[3][k]) for k in nodes)

    m.addConstr(
        enerbat_slack[0]==float(batteries[3][0]))    

    m.addConstrs(
        enerbat[i+1, k]==enerbat[i, k]+powerbat[i,k] for i in timeps for k in nodes)

    m.addConstrs(
        enerbat_slack[i+1]==enerbat_slack[i]+powerbat_slack[i] for i in timeps)

    m.addConstrs(
       powernet[i, k]==powergen[i, k]+uncpowergen[i, k]-powerload[i, k]-uncpowerload[i, k]-powerbat[i, k] for i in timeps for k in nodes)

    m.addConstrs(
       powernet_slack[i]==powergen_slack[i]+uncpowergen_slack[i]-powerload_slack[i]-uncpowerload_slack[i]-powerbat_slack[i] for i in timeps)

    # min uptime downtime constraints
    for i in timeps:
        for k in nodes:
            m.addConstr(
                    powergen[i, k]+(M*y[0,i,k])+(M*y[1,i,k])>=5)
            if i+1 in timeps:
                m.addConstr(
                    powergen[i+1, k]<=(M*(1-y[0,i,k]))+(M*y[1,i,k])+0.01)
            for l in range(genminuptimes[k]):
                if genminuptimes[k]!=0:
                    if i+2+l in timeps:
                        m.addConstr(
                            powergen[i+2+l, k]+M*(1-y[1,i,k])>=5)

    for i in timeps:
        m.addConstr(
            powergen_slack[i]+(M*yslack[0,i])+(M*yslack[1,i])>=5)
        if i+1 in timeps:
            m.addConstr(
                    powergen_slack[i+1]<=(M*(1-yslack[0,i]))+(M*yslack[1,i])+0.01)
        for l in range(genminuptime_slack):
            if genminuptime_slack!=0:
                if i+2+l in timeps:
                    m.addConstr(
                        powergen_slack[i+2+l]+M*(1-yslack[1,i])>=5)

    for i in timeps:
        for k in nodes:
            m.addConstr(
                    powerload[i, k]+(M*z[0,i,k])+(M*z[1,i,k])>=5)
            if i+1 in timeps:
                m.addConstr(
                    powerload[i+1, k]<=(M*(1-z[0,i,k]))+(M* z[1,i,k])+0.01)
            for l in range(loadminuptimes[k]):
                if loadminuptimes[k]!=0:
                    if i+2+l in timeps:
                        m.addConstr(
                            powerload[i+2+l, k]+M*(1-z[1,i,k])>=5)   

    for i in timeps:
        m.addConstr(
            powerload_slack[i]+(M*zslack[0,i])+(M*zslack[1,i])>=5)
        if i+1 in timeps:
            m.addConstr(
                    powerload_slack[i+1]<=(M*(1-zslack[0,i]))+(M*zslack[1,i])+0.01)
        for l in range(loadminuptime_slack):
            if loadminuptime_slack!=0:
                if i+2+l in timeps:
                    m.addConstr(
                            powerload_slack[i+2+l]+M*(1-zslack[1,i])>=5)   

    for t in timeps:
        for n1 in nodes:
            sum=0.0
            for n2 in nodes:
                sum+=B.iloc[n1-1,n2-1]*phaseangles[t,n2]
            m.addConstr(powernet[t,n1]==sum)        

    for t in timeps:
        sum=0.0 
        for n in nodes:
            sum+=b0[n-1]*phaseangles[t,n]
        m.addConstr(powernet_slack[t]==sum)         

    for t in timeps:
        for l in lines:
            sum=0
            for n in nodes:
                sum+=A.iloc[l-1,n-1]*phaseangles[(t,n)]
            ftemp[t,l-1]=sum   

    for t in timeps:
        for l1 in lines:
            sum=0
            for l2 in lines:
                sum+=Bd.iloc[l1-1,l2-1]*ftemp[t,l2-1]           
            m.addConstr(powerflow[t,l1]==sum)          

    m.addConstrs(
       powerflow[t, l]<=Llimits[0][l-1] for t in timeps for l in lines)

    m.addConstrs(
       powerflow[t, l]>=-Llimits[0][l-1] for t in timeps for l in lines)

    cost=LinExpr()
    for t in timeps:
        cost+=((Cs[t, 0]*powergen_slack[t]*powergen_slack[t])+(Cs[t, 1]*powergen_slack[t])+Cs[t,2])
        for n in nodes:
            cost+=((C[t,n,0]*powergen[t,n]*powergen[t,n])+(C[t,n, 1]*powergen[t,n])+C[t,n, 2])

    m.setObjective(cost, GRB.MINIMIZE)

    m.optimize()

    for v in m.getVars():
        print(v.varName, "{:.4f}".format(v.x))

    print('Obj:', m.objVal)

except GurobiError:
    print('Error reported')
