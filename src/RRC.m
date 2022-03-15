# Function for Real Root Classification

RealRootClassification:=proc(F,vars,params)
    return RRC_core(F,vars,params):
end proc:

RRC_core:=proc(F,vars,params)
    local n,lm,lmr,q,B,d,gb,multmat,matB,H,i,j,ind,Vtr1,ldenom,
        distinct_index,indexMat,gblm, NFQ,coeffsQ,termsQ,matQ,
        r,k,checked,b,nonproper,ld,wH,w,sp,nrealsols,nr,formulas,Dmax:

    n:=nops(vars):
    gblm:=FGb[fgb_gbasis_lm](F,0,vars,params):
    
    # Construct the basis of the quotient ring
    lmr:=subs([seq(params[i]=1,i=1..nops(params))],map(primpart,gblm[2])):
    lm:=FGb[fgb_gbasis](lmr,2,[],vars):

    nonproper:=[]:
    lmr:=ListTools[Reverse](lmr):
    for i from 1 to nops(lm) do:
        if member(lm[i],lmr,'k') then:
            nonproper:=[op(nonproper),mulCoeff(gblm[1][nops(lmr)+1-k],lmr[k],vars)]:
        end if:
    end do:
    nonproper:=remove(p->degree(p)=0,nonproper):
    nonproper:=np_factor(nonproper):

    q:=queue[new](1):
    B:=[]:
    checked:=[]:
    while not queue[empty](q) do:
        b:=queue[dequeue](q):
        if not member(b,checked) then:
            r:=Groebner[NormalForm](b,lm,tdeg(op(vars))):
            checked:=[op(checked),b]:
            if r <> 0 then:
                B:=[op(B),b]:
                for i from 1 to n do:
                    queue[enqueue](q,expand(vars[i]*b)):
                end do:
            end if:
        end if:
    end do:
    B := sort(B,(a,b)->Groebner[TestOrder](a,b,tdeg(op(vars)))):
    d:=nops(B):

    if nops(convert(indets(gblm[2]),list)) > nops(vars) then:
        gb:=cleanGB(gblm,vars,params):
    else:
        gb:=gblm[1]:
    end if:
    
    multmat:=Array(1..n,i->computeMatrix(vars[i],B,gb,vars,params,0)):
    matB := matrixBasis(multmat,B,d,vars,n):
    H:=Matrix(d):
    
    Vtr1:=Vector(d):
    for i from 1 to d do:
        Vtr1[i]:=normal(expand(Trace(matB[i]))): 
        H[1,i]:=Vtr1[i]:
    end do:

    ldenom:=map(denom,Vtr1):
    
    distinct_index, indexMat:=rmvDup(B):

    distinct_index:=distinct_index[(d+1)..-1]:

    for ind in distinct_index do:
        for j from 1 to d do:
            H[op(ind)] := H[op(ind)] + matB[ind[1]][ind[2],j]*H[1,j]:
            H[op(ind)]:=normal(expand(H[op(ind)])):
        end do:
        H[op(ind)]:=normal(expand(H[op(ind)]*ldenom[ind[1]]*ldenom[ind[2]])):
    end do:

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

    ld:=[seq(Determinant(SubMatrix(H,1..i,1..i)),i=2..d)]:
    if ld=[] then:
        return [H,ld,[[1,true]]]:
    end if:
    wH:=CleanFactors(numer(ld[-1]),nonproper):
    w:=expand(mul(wH)*mul(nonproper)):
    Dmax:=max(map(deg,F)):
    if 4*8^t*Dmax^(4*nops(vars)*nops(t)) > 2^d then:
        sp:=SamplePoints(w,params):
        nrealsols:=[seq(NumberOfRealSolutions(sp[i],H),i=1..nops(sp))]:
        nr:=convert(convert(nrealsols,set),list):
	    if nops(nr) = 1 then:
            if nr[1] = 0 then: return [H,ld,[[0,false]]]:
            else: return [H,ld,[[nr[1],true]]]: end if:
        else:
            formulas:=[seq([nr[i],signPatterns(d,nr[i])],i=1..nops(nr))]:
            return [H,ld,formulas]:
	    end if:
    else:
        sp:=RAG[PointsPerComponents]([seq(ld[i]<>0,i=1..nops(ld))]):
        nrealsols:=[seq(NumberOfRealSolutions(sp[i],H),i=1..nops(sp))]:
        nr:=convert(convert(nrealsols,set),list):
	    if nops(nr) = 1 then:
            if nr[1] = 0 then: return [H,ld,[[0,false]]]:
            else: return [H,ld,[[nr[1],true]]]: end if:
	    else:
            formulas:=[seq([nr[i],signPatterns(d,nr[i])],i=1..nops(nr))]:
            return [H,ld,formulas]:
	    end if:
    end if:
end proc:
