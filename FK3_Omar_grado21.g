# %%
LoadPackage("GBNP");
GBNP.ConfigPrint(
    "X1",  # "X12",
    "X2",  # "X13",
    "X3",  # "X23",
    "X4",  # "ext1",
    "X5"   # "ext2",
);

# %%
w:=E(3);
q1:= 1;
q2:= 1;
########################## Generators ########################################
F := FreeAssociativeAlgebraWithOne(Rationals, 5);;
gens := GeneratorsOfAlgebra(F);;
e:=gens[1];;
x1:=gens[2];;
x2:=gens[3];;
x3:=gens[4];;
x4:=gens[5];;
x5:=gens[6];;

# %%
x14 := x1*x4-x5*x1;;
x24 := -w*x5*x2+x2*x4;;
x34 := -w^2*x5*x3+x3*x4;;
x15 := -q2*x4*x1+x1*x5;;
x25  := -q2*w^2*x4*x2+x2*x5;;
x35 := -q2*w*x4*x3+x3*x5;;
Print("x114: ", "\n");;
x114 := -q1*q2*x14*x1+x1*x14;;
x124 := -q1*q2*w*x34*x1+x1*x24;;
x134 := -q1*q2*w^2*x24*x1+x1*x34;;
x214 := -q1*q2*w^2*x34*x2+x2*x14;;
Print("x224: ", "\n");
x224 := -q1*q2*x24*x2+x2*x24;;
x234 := -q1*q2*w*x14*x2+x2*x34;;
x314 := -q1*q2*w*x24*x3+x3*x14;;
x324 := -q1*q2*w^2*x14*x3+x3*x24;;
Print("x334: ", "\n");
x334 := -q1*q2*x34*x3+x3*x34;;

# %%
rels := [
    # ------------------------  grado (2,0) -------------------------- #
    x1^2,               
    x2^2,
    x3^2,
    x1*x2+x3*x1+x2*x3,
    x2*x1+x1*x3+x3*x2,
    # --------------------------  grado (0,2) -------------------------- #
    x4^2,
    x5^2,
    x4*x5+x5*x4,
    # -----------------------------  grado (2,2) ---------------------- #
    #x14^2,
    #x24^2,
    #x34^2,
    #x14*x24+x34*x14+x24*x34,
    #x24*x14+x14*x34+x34*x24,
    # -----------------------------  grado (1,1) ---------------------- #
    q1*x15+x14,
    q1*w*x25+x24,
    q1*w^2*x35+x34,
    # ------------------------------  grado (4,2) --------------------- #
    # x124^2,
    # x134^2,
    # x124*x134+x134*x124,
    # -------------------------- grado (2,1) -------------------------- #
    # x124+x234+x314,
    # x214+x134+x324,
    # w^2*x124-x234,
    # w^2*x134-x214
];;

# %%
# PrintNPList(List(rels, el-> GP2NP(el)));

# %%
# Imprimir las relaciones
I:=GP2NPList(rels);;
PrintNPList(I);

# %%
# =====================================================
# Paso 5: Calcular la base de Gröbner
# =====================================================
J := SGrobner(I);;

# %%
PrintNPList(J);

# %%
#DimQA(J, 5);

# %%
Xom1:=x124+x234+x314;
Xom2:=x214+x134+x324;
Xom3:=w^2*x124-x234;
Xom4:=w^2*x134-x214;

# %%
PrintNP(StrongNormalFormNP(GP2NP(Xom1), J)); # esta se hace cero
PrintNP(StrongNormalFormNP(GP2NP(Xom2), J)); # tambien esta 
PrintNP(StrongNormalFormNP(GP2NP(Xom3), J));
PrintNP(StrongNormalFormNP(GP2NP(Xom4), J));


