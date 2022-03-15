ParamHermite:=module()
description "Parametric Hermite matrix computation";
option package;
export DRLMatrix, NonSpecializationPolynomial, DRLMatrixNoDenom,
    RealRootClassification, LEXMatrix:

local cleanGB, computeEntry, naiveHermite, specHM, liftMatrix, HermiteParams,
QuotientBasis, NumberOfRealSolutions, CleanFactors,
DRLMatrixQuotientRing,DRLMatrix_inter, DRLMatrix_old, computeMatrix,
computeP, RRC_core,
indexDV,SamplePoints,GenericDimension,
multMatrices, matrixBasis, naiveHermite_modp, directHermite_modp,
directHermite, multMatrices_modp, matrixBasis_modp, specHM_modp, liftMatrix_modp, HermiteParams_modp,
sign_variant,np_factor,naiveHermite_nodenom,naiveHermite_quotient,euler_candy,signPatterns,
ppc1,ppc2,isol,classify,mipoint,
AssumptionTest,ZeroDimReduce,
mulCoeff, GenerateData, mInterp, mycoeffs, rmvDup,computeMatrix_quo, matrixBasis_quo, MonomialDiv, allMinors:

$include "supfun.m":
$include "main_Q.m":
$include "main_modp.m":
$include "main_quotient.m":
$include "sampling.m":
$include "LEXMatrix.m":
$include "RRC.m":

DRLMatrix:=proc(F,vars,params,char,Q:=1)
    description "Direct implementation without removing denominators";
    local H:
    if char = 0 then:
        H:=naiveHermite(F,vars,params,Q):
    else:
        H:=naiveHermite_modp(F,vars,params,char):
    end if:
    return H:
end proc:

DRLMatrixNoDenom:=proc(F,vars,params,char,Q:=1)
    description "direct implementation";
    local H:
    if char = 0 then:
        H:=naiveHermite_nodenom(F,vars,params,Q):
    else:
        H:=naiveHermite_modp(F,vars,params,char):
    end if:
    return H:
end proc:

DRLMatrix_inter:=proc(F,vars,params,char)
    description "evaluation / interpolation implementation";
    local H:
    if char = 0 then:
        H:=HermiteParams(F,vars,params):
    else:
        H:=HermiteParams_modp(F,vars,params,char):
    end if:
    return H:
end proc:

DRLMatrix_old:=proc(F,vars,params,char)
    description "old implementation : bad performance";
    local H:
    if char = 0 then:
        H:=directHermite(F,vars,params):
    else:
        H:=directHermite_modp(F,vars,params,char):
    end if:
    return H:
end proc:
end module: