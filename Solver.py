LANG="en_US.UTF-8"
import numpy as np
from ahkab import circuit, ahkab
import cmath
##################### Functions #####################


def isNumeric(x):
    try:
        float(x)
    except:
        return False
    return True


def getRectForm(r, phi):
    # returns a the rectangular form of the number
    r = float(r)
    phi = float(phi)
    return np.complex(r * np.math.cos(phi * np.math.pi / 180), r * np.math.sin(phi * np.math.pi / 180))

#####################################################
################### Containers ###################


R_Names = []; R_Nodes = []; R_Values = []; R_V_I_P = []
L_Names = []; L_Nodes = []; L_Values = []; L_V_I_P = []
C_Names = []; C_Nodes = []; C_Values = []; C_V_I_P = []

V_Names = []; V_Nodes = []; V_V_I_P = []
I_Names = []; I_Nodes = []; I_Values = []; I_Types = []; I_V_I_P = []

#####################################################
################### Intialization ###################
circ = circuit.Circuit(title="M's  circuit")
# File = open(input("Enter Your File Name: ") + '.txt')
File = open('YC.txt')
Frequancy = float(File.readline())
#####################################################
################ Parsing the netlist ################
for line in File:
    L = (File.readline()).split()
    element = L[0][0]
    name = L[0]
    nodePos = L[1]
    nodeNeg = L[2]
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
                I_Values.append((0, L[4]))
                circ.add_isource(name, nodePos, nodeNeg, dc_value=0, ac_value=float(L[4]))
        else:
            # the default type is AC
            I_Types.append(("", "AC"))
            I_Values.append((0, getRectForm(L[3], L[4])))
            circ.add_isource(name, nodePos, nodeNeg, dc_value=0, ac_value=getRectForm(L[3], L[4]))
#####################################################
################ Solving the netlist ################
print("The Netlist")
print(circ)
ac = ahkab.new_ac(start=Frequancy, stop=Frequancy, points=2, x0=None)
r = ahkab.run(circ, ac)
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
    P = round((np.power(cmath.phase(V), 2)) / (2 * float(R_Values[i])), roundTO)
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
    P = round((np.power(cmath.phase(V), 2)) / (Frequancy * 1j * 2 * np.math.pi * 2 * float(L_Values[i])), roundTO)
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
    P = round((np.power(cmath.phase(V), 2)) / (2 * (-1j / (2 * np.math.pi * float(C_Values[i]) * Frequancy))), roundTO)
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
    P = round(- 0.5 * V * np.conjugate(I_Values[i]), roundTO)
    I = I_Values[i]
    I_V_I_P.append((I_Names, V, I, P))
    IPS += P
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
