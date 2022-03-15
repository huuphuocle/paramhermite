# new implementation of parametric Hermite matrices for quotient ring Q[t]/w(t)

# params := [t] : # only one parameter case for now

computeMatrix_quo := proc(m,B,gb,vars,params,char,ord:=lexdeg(vars,params))
    local i,j,k,d,M,s,tmp,B2:
    B2:=map(primpart,B):
    d:=nops(B):
    M:=Matrix(d):
    for i from 1 to d do:
        tmp := expand(m*B[i]):
        if member(primpart(tmp),B2,'j') then:
            M[i,j]:=1:
        else:
            s:=Groebner[NormalForm](tmp,gb,ord,characteristic=char): 
            for k from 1 to d do:
                M[i,k]:=mulCoeff(s,B2[k],vars)/content(B[k]):
            end do:
        end if:
    end do:
    return M:
end proc:

matrixBasis_quo := proc(multmat,B,w,vars,varw,n)
    local d,matB,i,var,t1,t2:
    d:=nops(B):
    matB := Array(1..d):
    matB[1] := IdentityMatrix(d):
    for i from 2 to d do:
        var := convert(indets(B[i]),list):
        member(var[1],vars,'t1'):
        member(normal(B[i]/var[1]),B,'t2'):
        matB[i] := map(p->normal(rem(p,w,varw)),Multiply(multmat[t1],matB[t2])):
    end do:
    return matB:
end proc:

DRLMatrixQuotientRing:=proc(F,w,vars,params)
    description "Compute the DRL Hermite matrix for an ideal over a quotient ring Q[t]/w(t)";
    local n,B,d,gb,multmat,matB,H,i,j,ind,Vtr1,distinct_index,indexMat,gblm:
    n:=nops(vars):
    
    B:=QuotientBasis([op(F),w],vars,params,0):
    d:=nops(B):

    gblm:=FGb[fgb_gbasis_lm]([op(F),w],0,vars,params):

    gb:=gblm[1]:
    
    multmat:=Array(1..n,i->computeMatrix_quo(vars[i],B,gb,vars,params,0)):
    matB := matrixBasis_quo(multmat,B,w,vars,varw,n):
    H:=Matrix(d):
    
    Vtr1:=Vector(d):
    for i from 1 to d do:
        Vtr1[i]:=normal(expand(Trace(matB[i]))): 
        H[1,i]:=Vtr1[i]:
    end do:
    
    distinct_index, indexMat:=rmvDup(B):

    distinct_index:=distinct_index[(d+1)..-1]:

    for ind in distinct_index do:
        for j from 1 to d do:
            H[op(ind)] := H[op(ind)] + matB[ind[1]][ind[2],j]*Vtr1[j]:
            H[op(ind)]:=normal(rem(expand(H[op(ind)]),w,varw)):
        end do:
    end do:
    
    # fill the repeated entries

    for i from 1 to d do:
        for j from i to d do:
            if indexMat[i,j] <> 0 then:
                H[i,j]:=H[op(indexMat[i,j])]:
            end if:
        end do:
    end do:
    for i from 1 to d do:
        for j from 1 to i-1 do: 
            H[i,j] := H[j,i]:
        end do:
    end do:

    return H:
end proc:

# naiveHermite_quotient:=proc(F,vars,params,w,t)
#     local n,B,d,gb,multmat,matB,H,i,j,ind,Vtr1,
#         distinct_index,indexMat,
#         gblm:
#     n:=nops(vars):
    
#     B:=QuotientBasis([op(F),w],vars,[op(params),t],0):
#     d:=nops(B):

#     gblm:=FGb[fgb_gbasis_lm]([op(F),w],0,vars,[op(params),t]):

#     gb:=gblm[1][1..-2]:
    
#     multmat:=Array(1..n,i->computeMatrix(vars[i],B,gb,vars,[op(params),t],0)):
#     matB := matrixBasis(multmat,B,d,vars,n):
#     H:=Matrix(d):
    
#     Vtr1:=Vector(d):
#     for i from 1 to d do:
#         Vtr1[i]:=normal(expand(Trace(matB[i]))): 
#         H[1,i]:=Vtr1[i]:
#     end do:
    
#     distinct_index, indexMat:=rmvDup(B):

#     distinct_index:=distinct_index[(d+1)..-1]:

#     for ind in distinct_index do:
#         for j from 1 to d do:
#             H[op(ind)] := H[op(ind)] + matB[ind[1]][ind[2],j]*Vtr1[j]:
#             H[op(ind)]:=normal(expand(H[op(ind)])):
#         end do:
#     end do:

#     for i from 1 to d do:
#         for j from i to d do:
#             if indexMat[i,j] <> 0 then:
#                 H[i,j]:=H[op(indexMat[i,j])]:
#             end if:
#         end do:
#     end do:
#     for i from 1 to d do:
#         for j from 1 to i-1 do: 
#             H[i,j] := H[j,i]:
#         end do:
#     end do:

#     return H:
# end proc: