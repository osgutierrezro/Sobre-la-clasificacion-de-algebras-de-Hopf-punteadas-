LoadPackage("GBNP");
	GBNP.ConfigPrint("X1","X2","X3","X4","X5");
	
	################ Cuerpo de coeficientes: CF(3)(t1) ####################
	K  := CF(3);;
	R  := FunctionField(K, ["t1"]);;
	t1 := IndeterminatesOfFunctionField(R)[1];;
	
	################ Algebra tensorial T(HV_3) ############################
	F    := FreeAssociativeAlgebraWithOne(R, 5);;
	gens := GeneratorsOfAlgebra(F);;
	one  := One(F);;
	x1 := gens[2];; x2 := gens[3];; x3 := gens[4];;
	x4 := gens[5];; x5 := gens[6];;
	
	################ Escalares ############################################
	w  := (E(3)*One(R)) * one;;      # raiz cubica primitiva de la unidad
	q1 := t1 * one;;                 # parametro q1
	q2 := t1^-2 * one;;              # q2 = q1^(-2), asi q1^2*q2 = 1
	
	Print("q1^2*q2 = 1 : ", q1^2*q2 = one, "\n");
	Print("1+w+w^2 = 0 : ", one + w + w^2 = Zero(F), "\n");
	
	################ x_{ij} := ad_c(x_i)(x_j) #############################
	x14 := x1*x4 - x5*x1;;
	x24 := x2*x4 - w*x5*x2;;
	x34 := x3*x4 - w^2*x5*x3;;
	x15 := x1*x5 - q2*x4*x1;;
	x25 := x2*x5 - q2*w^2*x4*x2;;
	x35 := x3*x5 - q2*w*x4*x3;;
	
	################ x_{1j4} := ad_c(x_1)(x_{j4}) #########################
	x124 := x1*x24 - q1*q2*w*x34*x1;;
	x134 := x1*x34 - q1*q2*w^2*x24*x1;;
	
	################ Conjunto G: relaciones (rel1)-(rel8) #################
	G := [ x1^2, x2^2, x3^2, x4^2, x5^2,
	x1*x2 + x3*x1 + x2*x3,
	x2*x1 + x1*x3 + x3*x2,
	x4*x5 + x5*x4,
	x14 + q1*x15,
	x24 + w*q1*x25,
	x34 + w^2*q1*x35,
	x14^2, x24^2, x34^2,
	x124*x134 + x134*x124 ];;
	
	nombres := [ "x1^2", "x2^2", "x3^2", "x4^2", "x5^2",
	"x1x2+x3x1+x2x3", "x2x1+x1x3+x3x2", "x4x5+x5x4",
	"x14+q1x15", "x24+wq1x25", "x34+w^2q1x35",
	"x14^2", "x24^2", "x34^2", "x124x134+x134x124" ];;
	
	grados := [ 2,2,2,2,2, 2,2,2, 2,2,2, 4,4,4, 6 ];;
	
	################ 1) Dimension de T(HV_3)/<G> ##########################
	GB := SGrobner(GP2NPList(G));;
	Print("dim T(HV_3)/<G> = ", DimQA(GB, 5), "\n");     # 2304
	
	################ 2) Minimalidad: r no esta en <G \ {r}> ###############
	for k in [1..Length(G)] do
	resto := G{Difference([1..Length(G)], [k])};;
	GBr := SGrobnerTrunc(GP2NPList(resto), grados[k], [1,1,1,1,1]);;
	if IsRecord(GBr) then GBr := GBr.G; fi;
	if StrongNormalFormNP(GP2NP(G[k]), GBr) = [[],[]] then
	estado := "REDUNDANTE";
	else
	estado := "esencial";
	fi;
	Print(nombres[k], " (grado ", grados[k], "): ", estado, "\n");
	od;