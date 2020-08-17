LANG="en_US.UTF-8"
import numpy as np
from ahkab import circuit, ahkab
import cmath
##################### Functions #####################


def getRectForm(r, phi):
    # returns a the rectangular form of the number
    r = float(r)
    phi = float(phi)
    return np.complex(r * np.math.cos(phi * np.math.pi / 180), r * np.math.sin(phi * np.math.pi / 180))


#####################################################
################### Intialization ###################
circ = circuit.Circuit("M's  circuit")
File = open(input("Enter Your File Name: ") + '.txt')
# File = open('YC.txt')
Frequancy = float(File.readline())

#####################################################
################### Containers ###################


R_Names = []; R_Nodes = []; R_Values = []; R_V_I_P = []
L_Names = []; L_Nodes = []; L_Values = []; L_V_I_P = []
C_Names = []; C_Nodes = []; C_Values = []; C_V_I_P = []

V_Names = []; V_Nodes = []; V_V_I_P = []
I_Names = []; I_Nodes = []; I_Values = []; I_Types = []; I_V_I_P = []

G_Names = []; G_Nodes = []; G_Values = []; G_depNodes =[]; G_V_I_P = []
E_Names = []; E_Nodes = []; E_Values = []; E_depNodes =[]; E_V_I_P = []
H_Names = []; H_Nodes = []; H_Values = []; H_ZV_ID =[]; H_V_I_P = []
F_Names = []; F_Nodes = []; F_Values = []; F_ZV_ID =[]; F_V_I_P = []

## MOSFET
M_Names = []; M_Nodes = []; #(nd,ng,ns,nb)
M_Values = [] ; # (length,width)
M_V_I_P = []
#####################################################
################ Parsing the netlist ################


for line in File:
    # L = (File.readline()).split()  *** menna *** Something goes wrong with this line, but I couldn't find it. (line.split()) works fine.
    L=line.split()
    element = L[0][0]
    name = L[0]
    if not (L[0][0] =='M' or L[0][0] =='m'):
     nodePos = L[1]
     nodeNeg = L[2]
    else:
      nodeDrain=L[1]
      nodeGate=L[2]
      nodeSource=L[3]
      nodeBody=L[4]

    # Passive elements
    # R
    if element == 'R' or element == 'r':
        R_Names.append(name); R_Nodes.append((nodePos, nodeNeg)); R_Values.append(float(L[3]))

        circ.add_resistor(name, nodePos, nodeNeg, value=float(L[3]))
    # L
    elif element == 'L' or element == 'l':
        L_Names.append(name); L_Nodes.append((nodePos, nodeNeg)); L_Values.append(float(L[3]))
        circ.add_inductor(name, nodePos, nodeNeg, value=float(L[3]))
    # C
    elif element == 'C' or element == 'c':
        C_Names.append(name); C_Nodes.append((nodePos, nodeNeg)); C_Values.append(float(L[3]))
        circ.add_capacitor(name, nodePos, nodeNeg, value=float(L[3]))
    # independent source
    # V
    elif element == 'V' or element == 'v':
        V_Names.append(name); V_Nodes.append((nodePos, nodeNeg));
        if L[3] == 'DC' or L[3] == 'dc':
            if L[5] == 'AC' or L[5] == 'ac':
                circ.add_vsource(name, nodePos, nodeNeg, dc_value=float(L[4]), ac_value=float(L[6]))
            else:
                circ.add_vsource(name, nodePos, nodeNeg, dc_value=float(L[4]))
        elif L[3] == 'AC' or L[3] == 'ac':
            if L[5] == 'DC' or L[5] == 'dc':
                circ.add_vsource(name, nodePos, nodeNeg, dc_value=float(L[6]), ac_value=float(L[4]))
            else:
                circ.add_vsource(name, nodePos, nodeNeg, ac_value=float(L[4]))
        else:
            # the default type is AC
            circ.add_vsource(name, nodePos, nodeNeg, dc_value=0, ac_value=getRectForm(L[3], L[4]))
    # I
    elif element == 'I' or element == 'I':
        I_Names.append(name);
        I_Nodes.append((nodePos, nodeNeg));
        if L[3] == 'DC' or L[3] == 'dc':
            if L[5] == 'AC' or L[5] == 'ac':
                I_Types.append(("DC", "AC"))
                I_Values.append((L[4], L[6]))
                circ.add_isource(name, nodePos, nodeNeg, dc_value=float(L[4]), ac_value=float(L[6]))
            else:
                I_Types.append(("DC", ""))
                I_Values.append((L[4], 0))
                circ.add_isource(name, nodePos, nodeNeg, dc_value=float(L[4]))
        elif L[3] == 'AC' or L[3] == 'ac':
            if L[5] == 'DC' or L[5] == 'dc':
                I_Types.append(("DC", "AC"))
                I_Values.append((L[6], L[4]))
                circ.add_isource(name, nodePos, nodeNeg, dc_value=float(L[6]), ac_value=float(L[4]))
            else:
                I_Types.append(("", "AC"))
                I_Values.append((0, float(L[4])))
                circ.add_isource(name, nodePos, nodeNeg, dc_value=0, ac_value=float(L[4]))
        else:
            # the default type is AC
            I_Types.append(("", "AC"))
            # I_Values.append((0, getRectForm(L[3], L[4])))
            I_Values.append((L[3], L[4]))
            circ.add_isource(name, nodePos, nodeNeg, dc_value=0, ac_value=getRectForm(L[3], L[4]))

    # Dependent sources
    # ccvs --> current controlled voltage source
    elif element == 'H' or element == 'h':
        H_Names.append(name); H_Nodes.append((nodePos, nodeNeg));
        H_ZV_ID.append(L[3]); H_Values.append(float(L[4]));
        circ.add_ccvs(name, nodePos, nodeNeg,L[3],float(L[4]))

    # cccs --> current controlled current source
    elif element == 'F' or element == 'f':
        F_Names.append(name); F_Nodes.append((nodePos, nodeNeg));
        F_ZV_ID.append(L[3]); F_Values.append(float(L[4]));
        circ.add_cccs(name, nodePos, nodeNeg,L[3],float(L[4]))

    # vccs --> volt controlled current source
    elif element == 'G' or element == 'g':
        G_Names.append(name);
        G_Nodes.append((nodePos, nodeNeg));
        G_depNodes.append((L[3],L[4]));
        G_Values.append(float(L[5]));
        circ.add_vccs(name, nodePos, nodeNeg, L[3], L[4],float(L[5]))

    # vcvs --> volt controlled voltage source
    elif element == 'E' or element == 'e':
        E_Names.append(name);
        E_Nodes.append((nodePos, nodeNeg));
        E_depNodes.append((L[3],L[4]));
        E_Values.append(float(L[5]));
        circ.add_vccs(name, nodePos, nodeNeg, L[3], L[4],float(L[5]))

    # # Transistors
    # # (M --> MOSFET)
    elif element == 'M' or element == 'm':
        M_Names.append(name);
        M_Nodes.append((nodeDrain, nodeGate,nodeSource,nodeBody));
        M_Values.append((L[5],L[6]))
        circ.add_mos(name,nodeDrain,nodeGate,nodeSource,nodeBody,L[5],L[6])

#####################################################
################ Solving the netlist ################


print("The Netlist")
print(circ)
ac = ahkab.new_ac(start=Frequancy, stop=Frequancy, points=2, x0=None)
r = ahkab.run(circ, ac)

## for DC circuit      ** menna ** If our code will deal with DC circuit as well, then we need this condition and accordingly the code below will be modified to fit this.
# else:
#  op = ahkab.new_op( x0=None)
#  r = ahkab.run(circ, op)

# for i in r['ac']:
#     print(i[0] + "= ", i[1][1])
File.close()

#####################################################
################ saving the data ################

print("Solution: ID-->V-->I-->Pw")
# Passive elements
roundTO = 5
# The active resistor power values
RPS = 0
for i, node in enumerate(R_Nodes):
    v1 = 0
    v2 = 0
    if node[0] != '0':
        v1 = r['ac']['V' + node[0]][0]
    if node[1] != '0':
        v2 = r['ac']['V' + node[1]][0]
    V = v1 - v2
    # P = round((np.power(cmath.phase(V), 2)) / (2 * float(R_Values[i])), roundTO)
    P=(np.power(V, 2)) / (2 * abs(R_Values[i]))
    I = V / float(R_Values[i])
    R_V_I_P.append((R_Names[i], V, I, P))
    RPS += P
    print(R_V_I_P[i])

# The reactive inductor power values
LPS = 0
for i, node in enumerate(L_Nodes):
    v1 = 0
    v2 = 0
    if node[0] != '0':
        v1 = r['ac']['V' + node[0]][0]
    if node[1] != '0':
        v2 = r['ac']['V' + node[1]][0]
    V = v1 - v2
    print(V)
    # P = round((np.power(cmath.phase(V), 2)) / (Frequancy * 1j * 2 * np.math.pi * 2 * float(L_Values[i])), roundTO)
    P = (np.power(V, 2)) / (Frequancy * 1j * 2 * np.math.pi * 2 * float(L_Values[i]))
    I = V / ((float(L_Values[i])) * 2 * np.math.pi * Frequancy * 1j)

    L_V_I_P.append((L_Names[i], V, I, P))
    LPS += P
    print(L_V_I_P[i])

# The reactive capacitor power values
CPS = 0
for i, node in enumerate(C_Nodes):
    v1 = 0
    v2 = 0
    if node[0] != '0':
        v1 = r['ac']['V' + node[0]][0]
    if node[1] != '0':
        v2 = r['ac']['V' + node[1]][0]
    V = v1 - v2
    # P = round((np.power(cmath.phase(V), 2)) / (2 * (-1j / (2 * np.math.pi * float(C_Values[i]) * Frequancy))), roundTO)
    P = (np.power(V, 2)) / (2 * (-1j / (2 * np.math.pi * float(C_Values[i]) * Frequancy)))
    I = V / (-1j / (float(C_Values[i])) * 2 * np.math.pi * Frequancy)
    C_V_I_P.append((C_Names[i], V, I, P))
    CPS += P
    print(C_V_I_P[i])
####
####
# Delivered power by the independent voltage sources
VPS = 0
for i, node in enumerate(V_Nodes):
    v1 = 0
    v2 = 0
    if node[0] != '0':
        v1 = r['ac']['V' + node[0]][0]
    if node[1] != '0':
        v2 = r['ac']['V' + node[1]][0]
    V = v1 - v2
    P = round(V * np.conjugate(r['ac']['I' + '(' + V_Names[i] + ')'][0]) / 2, roundTO)
    I = (r['ac']['I' + '(' + V_Names[i] + ')'][0])
    V_V_I_P.append((V_Names[i], V, I, P))
    VPS += P
    print(V_V_I_P[i])
    
# Delivered power by the independent current sources sources
IPS = 0
for i, node in enumerate(I_Nodes):
    v1 = 0
    v2 = 0
    if node[0] != '0':
        v1 = r['ac']['V' + node[0]][0]
    if node[1] != '0':
        v2 = r['ac']['V' + node[1]][0]
    V = v1 - v2
    # P = round(-0.5 * V * np.conjugate(I_Values[i][1]), roundTO) ### **menna** I_Values[i] causes an error because it is an array.

    P = 0.5 * (V) *np.conjugate(getRectForm(I_Values[i][0],I_Values[i][1]))
    I = getRectForm(I_Values[i][0],I_Values[i][1])
    I_V_I_P.append((I_Names, V, I, P))
    IPS += P
    print(I_V_I_P)

# Delivered power by the vcvs
EPS = 0
for i, node in enumerate(E_Nodes):
    v1 = 0
    v2 = 0
    if node[0] != '0':
        v1 = r['ac']['V' + node[0]][0]
    if node[1] != '0':
        v2 = r['ac']['V' + node[1]][0]
    V = v1 - v2
    I=r['ac']['I(' + L[0] + ')'][0]
    P = 0.5 * (E_Values[i] * ((r['ac']['V' + E_depNodes[i][0]][0]) - r['ac']['V' +E_depNodes[i][1]][0])) * (np.conjugate(I))
    E_V_I_P.append((G_Names, V, I, P))
    EPS += P
    print(E_V_I_P)

# Delivered power by the vccs
GPS = 0
for i, node in enumerate(G_Nodes):
    v1 = 0
    v2 = 0
    dv1=0
    dv2=0
    if node[0] != '0':
        v1 = r['ac']['V' + node[0]][0]
    if node[1] != '0':
        v2 = r['ac']['V' + node[1]][0]
    V = v1 - v2
    if G_depNodes[i][0] != '0':
        dv1 = r['ac']['V' + G_depNodes[i][0]][0]
    if G_depNodes[i][1] != '0':
        dv2 = r['ac']['V' + G_depNodes[i][1]][0]
    dV=dv1-dv2
    I=G_Values[i] * dV
    P = 0.5 * (V) * (np.conjugate(I))
    G_V_I_P.append((G_Names, V, I, P))
    GPS += P
    print(G_V_I_P)


# Delivered power by the ccvs
HPS = 0
for i, node in enumerate(H_Nodes):
    v1 = 0
    v2 = 0
    if node[0] != '0':
        v1 = r['ac']['V' + node[0]][0]
    if node[1] != '0':
        v2 = r['ac']['V' + node[1]][0]
    V = v1 - v2
    I=r['ac']['I('+H_Names[i] + ')'][0]
    P = 0.5 * (H_Values[i] * (r['ac']['I('+H_ZV_ID[i] + ')'][0])) * np.conjugate(I)
    H_V_I_P.append((H_Names, V, I, P))
    HPS += P
    print(H_V_I_P)


# Delivered power by the cccs
FPS = 0
for i, node in enumerate(F_Nodes):
    v1 = 0
    v2 = 0
    if node[0] != '0':
        v1 = r['ac']['V' + node[0]][0]
    if node[1] != '0':
        v2 = r['ac']['V' + node[1]][0]
    V = v1 - v2
    I=r['ac']['I('+F_ZV_ID[i] + ')'][0]*F_Values[i]
    P = 0.5 * (V) * np.conjugate(I)
    F_V_I_P.append((F_Names, V, I, P))
    FPS += P
    print(F_V_I_P)

# MOS nodes Voltages:
MPS = 0
for i, node in enumerate(M_Nodes):
    vd = 0
    vg = 0
    vs = 0
    vb = 0
    if node[0] != '0':
        vd = r['ac']['V' + node[0]][0]
    if node[1] != '0':
        vg = r['ac']['V' + node[1]][0]
    if node[2] != '0':
        vs = r['ac']['V' + node[2]][0]
    if node[3] != '0':
        vb = r['ac']['V' + node[3]][0]
    V = []
    V.append((vd,vg,vs,vb))

    I = I_Values[i]
    M_V_I_P.append((M_Names, V))
    print(I_V_I_P)

#####################################################
################ Power Balance ################

pwdelivered = VPS+IPS
print("pwDelivered", pwdelivered)
pwconsumed = RPS
print("pwConsumed", pwconsumed)
pwreserved = LPS + CPS
print("pwReserved", pwreserved)
pwbalance = round(pwdelivered-(pwreserved-pwconsumed), roundTO)
print("pwBalance", pwbalance)

