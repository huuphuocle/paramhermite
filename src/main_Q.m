#  NonSpecializationPolynomial:=proc(F,vars,params,char)
#      local gblm,lm,lm_red,lc,i,k:
#      gblm:=FGb[fgb_gbasis_lm](F,char,vars,params,{"verb"=3}):
#      lm:=map(primpart,subs([seq(params[i]=1,i=1..nops(params))],gblm[2])):
#      lm_red:=FGb[fgb_gbasis](lm,2,[],vars):
#      lc:=[]:
#      lm:=ListTools[Reverse](lm):
#      for i from 1 to nops(lm_red) do:
#          if member(lm_red[i],lm,'k') then:
#              lc:=[op(lc),mulCoeff(gblm[1][nops(lm)+1-k],lm[k],vars)]:
#          end if:
#      end do:
#      lc:=remove(p->degree(p)=0,lc):
#      lc:=np_factor(lc):
#      return lc:
#  end proc:

NonSpecializationPolynomial:=proc(F,vars,params,char)
    local m,gblm,lm,lm_red,lc,i,k:
    gblm:=FGb[fgb_gbasis_lm](F,char,vars,params,{"verb"=3}):
    lm:=map(primpart,subs([seq(params[i]=1,i=1..nops(params))],gblm[2])):
    lm_red:=FGb[fgb_gbasis](lm,2,[],vars):
    lc:=[seq([m,[]],m in lm_red)]:
    for i from 1 to nops(lm) do:
        if member(lm[i],lm_red,'k') then:
            lc[k][2]:=[op(lc[k][2]),mulCoeff(gblm[1][i],lm[i],vars)]:
        end if:
    end do:
    return lc:
end proc:

### General implementation without removing denominators

naiveHermite:=proc(F,vars,params,Q:=1)
    local n,lm,q,B,d,gb,multmat,matB,H,i,j,ind,Vtr1,
        distinct_index,indexMat,gblm, NFQ,coeffsQ,termsQ,matQ,
        r,k,checked,b:

    n:=nops(vars):
    gblm:=FGb[fgb_gbasis_lm](F,0,vars,params):
    
    # Construct the basis of the quotient ring
    lm:=subs([seq(params[i]=1,i=1..nops(params))],map(primpart,gblm[2])):
    lm:=FGb[fgb_gbasis](lm,2,[],vars):

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

    if Q <> 1 then:
        
        # entries when Q <> 1

        NFQ := Groebner[NormalForm](Q,gb,tdeg(op(vars))):
        # we write a normal form of Q, then we add the terms w.r.t. this normal form to obtain the matrix for Q
        coeffsQ := [coeffs(NFQ,vars,'termsQ')]: 
        termsQ:=[termsQ]:
        matQ:=Matrix(d):
        for i from 1 to nops(coeffsQ) do:
            member(termsQ[i],B,'j'):
            matQ:=matQ + coeffsQ[i]*matB[j]:
        end do:
        matQ:=map(expand,matQ): # multiplication matrix by Q
        j:='j':

        for i from 1 to d do:
            H[1,i]:=0:
            for j from 1 to d do:
                for k from 1 to d do:
                    H[1,i]:=H[1,i]+matB[i][j,k]*matQ[k,j]:
                end do:
            end do:
        end do:
    else:

        # Q = 1
        
        Vtr1:=Vector(d):
        for i from 1 to d do:
            Vtr1[i]:=normal(expand(LinearAlgebra[Trace](matB[i]))): 
            H[1,i]:=Vtr1[i]:
        end do:
    end if:
    
    distinct_index, indexMat:=rmvDup(B):

    distinct_index:=distinct_index[(d+1)..-1]:

    for ind in distinct_index do:
        for j from 1 to d do:
            H[op(ind)] := H[op(ind)] + matB[ind[1]][ind[2],j]*H[1,j]:
            H[op(ind)]:=normal(expand(H[op(ind)])):
        end do:
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

    return H:
end proc:

# General implementation with denominators removed

naiveHermite_nodenom:=proc(F,vars,params,Q:=1)
    local n,B,d,gb,multmat,matB,H,i,j,ind,Vtr1,
        distinct_index,indexMat,ldenom,
        gblm,proper, NFQ,coeffsQ,termsQ,matQ,k,
        lm,q,checked,b,r:

    n:=nops(vars):
    gblm:=FGb[fgb_gbasis_lm](F,0,vars,params):
    
    # Construct the basis of the quotient ring
    lm:=subs([seq(params[i]=1,i=1..nops(params))],map(primpart,gblm[2])):
    lm:=FGb[fgb_gbasis](lm,2,[],vars):

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
        proper:=false:
        gb:=cleanGB(gblm,vars,params):
    else:
        proper:=true:
        gb:=gblm[1]:
    end if:
    
    multmat:=Array(1..n,i->computeMatrix(vars[i],B,gb,vars,params,0)):
    matB := matrixBasis(multmat,B,d,vars,n):
    H:=Matrix(d):

    if Q <> 1 then:
        
        # entries when Q <> 1 (there is one inequality)

        NFQ := Groebner[NormalForm](Q,gb,tdeg(op(vars))):
        coeffsQ := [coeffs(NFQ,vars,'termsQ')]:
        termsQ:=[termsQ]:
        matQ:=Matrix(d):
        for i from 1 to nops(coeffsQ) do:
            member(termsQ[i],B,'j'):
            matQ:=matQ + coeffsQ[i]*matB[j]:
        end do:
        matQ:=map(expand,matQ):
        j:='j':

        for i from 1 to d do:
            H[1,i]:=0:
            for j from 1 to d do:
                for k from 1 to d do:
                    H[1,i]:=H[1,i]+matB[i][j,k]*matQ[k,j]:
                end do:
            end do:
        end do:
    else:

        # Q = 1 (there is no inequality)
        
        Vtr1:=Vector(d):
        for i from 1 to d do:
            Vtr1[i]:=normal(expand(LinearAlgebra[Trace](matB[i]))):
            H[1,i]:=Vtr1[i]:
            # H[1,i]:=numer(Vtr1[i]):
        end do:
    end if:

    ldenom:=map(denom,Vtr1):

    distinct_index, indexMat:=rmvDup(B):

    distinct_index:=distinct_index[(d+1)..-1]:

    for ind in distinct_index do:
        for j from 1 to d do:
            H[op(ind)]:=H[op(ind)] + matB[ind[1]][ind[2],j]*Vtr1[j]:
            H[op(ind)]:=normal(expand(H[op(ind)])):
        end do:
        H[op(ind)]:=normal(expand(H[op(ind)])):
        # H[op(ind)]:=normal(expand(H[op(ind)]*ldenom[ind[1]]*ldenom[ind[2]])):
    end do:

    for i from 1 to d do:
        for j from i to d do:
            if indexMat[i,j] <> 0 then:
                H[i,j]:=H[op(indexMat[i,j])]:
            end if:
        end do:
    end do:
    for i from 1 to d do:
        for j from 1 to i do:
            H[j,i]:=normal(expand(H[j,i]*ldenom[i]*ldenom[j])):
            H[i,j] := H[j,i]:
        end do:
    end do:

    return H:
end proc:

### Direct implemetation of Hermite matrices. 
### It computes each entry separately and is not efficient.
### Should be used for testing other implementations.

directHermite:=proc(F,vars,params)
    local st,gb,B,M,i,j,l,computed,computed_index,t,bb: 
    st := time():
    B := QuotientBasis(F,vars,params,0,tdeg(op(vars))):
    l := nops(B):
    printf("Basis : %d ", nops(B));
    print(B);
    gb := FGb:-fgb_gbasis(F,0,vars,params):

    M := Matrix(l):
    
    computed:=B:
    computed_index:=[seq([1,i],i=1..l)]:

    st := time():

    M[1,1]:=l:
    for j from 2 to l do:
        M[1,j]:=computeEntry(B[j],B,gb,vars,params,0,tdeg(op(vars))):
    end do:

    for i from 2 to l do:
        for j from i to l do:
            bb:=expand(B[i]*B[j]):
            if member(bb,computed,'t') then:
                M[i,j] := M[computed_index[t][1],computed_index[t][2]]:
            else:
                M[i,j] := computeEntry(bb,B,gb,vars,params,0,tdeg(op(vars))):
                computed := [op(computed),bb]:
                computed_index:=[op(computed_index),[i,j]]:
            end if:
        end do:
    end do:

    for i from 1 to l do:
        for j from 1 to i-1 do:
            M[i,j] := M[j,i]:
        end do:
    end do:
    printf("Entries: %d/%d \n",nops(computed),l*(l+1)/2);
    printf("Elapsed time: %f \n",time()-st);
    return M:
end proc:

### Compute multiplication matrices of x_i using row reduction
multMatrices := proc(gb,lt,B,D,vars,n,C,S,lC)
    local ll, NF, dmin, dmax, T, i, j, k, ind, ind2, ind3, M, d:
    dmin := degree(C[1]):
    dmax := degree(C[-1]):
    ll := [op(C),op(B)]:
    NF := Array(1..(D+S)):
    for i from 1 to D do:
        NF[i] := B[i]:
    end do:
    
    M := Matrix(S,D+S):
    ind := 1:
    
    for d from 1 to dmax-dmin+1 do:
        ind2 := ind:
        for i from 1 to nops(lC[d]) do:
            if member(lC[d][i],lt,'t') then:
                for j from 1 to D+S do:
                    M[ind,j] := mulCoeff(gb[t],ll[j],vars):
                end do:
            else:
                for k from 1 to n do:
                    if member(normal(lC[d][i]/vars[k]),lC[d-1],'t') then:
                        for j from 1 to D+S do:
                            M[ind,j] := mulCoeff(expand(lC[d][i]-vars[k]*NF[ind3+t-1+D]),ll[j],vars):
                        end do:
                        break;
                    end if:
                end do:
            end if:
            ind := ind + 1:
        end do:
        #M := MTM:-rref(M):
        M:=LinearAlgebra:-ReducedRowEchelonForm(M):

        for i from ind2 to ind-1 do:
            for j from 1 to D do:
                NF[i+D] := NF[i+D] - M[i,S+j]*B[j]: 
            end do:
        end do:
        ind3 := ind2:
    end do:
    T := Array(1..n,i->Matrix(D)):
    ll := [op(B),op(C)]:
    
    # build the matrix T(xi)
    for j from 1 to D+S do:
        for i from 1 to n do:
            if member(vars[i],indets(ll[j])) and member(normal(ll[j]/vars[i]),B,'t') then:
                for k from 1 to D do:
                    T[i][k,t] := mulCoeff(NF[j],B[k],vars):
                end do:
            end if:
        end do:
    end do:
    
    return T:
end proc:

# compute multiplication tables of B: L_1,...,L_d

matrixBasis := proc(multmat,B,d,vars,n)
    local matB,i,var,t1,t2:
    matB := Array(1..d):
    matB[1] := IdentityMatrix(d):
    for i from 2 to d do:
        var := convert(indets(B[i]),list):
        member(var[1],vars,'t1'):
        member(normal(B[i]/var[1]),B,'t2'):
        matB[i] := map(normal,LinearAlgebra[Multiply](multmat[t1],matB[t2])):
    end do:
    return matB:
end proc:

# compute each specialized Hermite matrix

specHM := proc(F,B,D,vars,n,C,S,listC,distinct_index)
    local gb,lt,H,i,j,k,multmat,matB,Vtr1,st:
    
    gb:=FGb:-fgb_gbasis_lm(F,0,[],vars):
    lt:=map(p->primpart(p),gb[2]):
    gb:=gb[1]:
    #st:=time():
    multmat:=multMatrices(gb,lt,B,D,vars,n,C,S,listC):
    #st:=time():
    matB := matrixBasis(multmat,B,D,vars,n):
    H:=Matrix(D):
    
    #st:=time():
    for j from 1 to D do:
        H[1,j]:=LinearAlgebra[Trace](matB[j]):
    end do:

    Vtr1 := [seq(H[1,k],k=1..D)]:
    
    #st:=time():
    for i in distinct_index do:
        for k from 1 to D do:
            H[op(i)] := H[op(i)] + matB[i[1]][k,i[2]]*Vtr1[k]:
        end do:
    end do:
    #printf("entries : ");
    #print(time()-st);
    return H:
end proc:

# interpolation step

liftMatrix:=proc(F,npoints,params,B,D,vars,n,C,S,listC,distinct_index)
	local xdata, poldata, H, x, i, j, k, st, index_short:
    
    # initialize H (Hermite matrix)
	H:=Matrix(D):

    # random points
    #xdata:=GenerateData(npoints,0):
    xdata:=[seq(i,i=1..npoints)]:
	poldata:=[]:
    index_short:=distinct_index[(D+1)..-1]:

    if nops(params) = 0 then:
        H:=specHM(F,B,D,vars,n,C,S,listC,index_short):
        return H:
    end if:
    
    # evaluation/ interpolation
    if nops(params)=1 then:
        # one variable left 
        
        # evaluation : call to specHM
        #st:=time():
        poldata:=[seq(specHM(subs(params[1]=x,F),B,D,vars,n,C,S,listC,index_short),x in xdata)]:
        #printf("specHM : ");
		#print(time()-st);
        # interpolate 
        for i in distinct_index do:
                H[op(i)]:=interp(xdata,[seq(poldata[k][op(i)],k=1..npoints)],params[1]):
        end do:
    else:
        # more than one variable : recursive call to liftMatrix
        for x in xdata do:
            poldata:=[op(poldata),liftMatrix(subs(params[1]=x,F),npoints,params[2..-1],B,D,vars,n,C,S,listC,distinct_index)]:
        end do:
        for i in distinct_index do:
            H[op(i)]:=mInterp(xdata,[seq(poldata[k][op(i)],k=1..npoints)],params[2..-1],params[1],0):
        end do:
    end if:
    return H:
end proc:

HermiteParams:=proc(F,vars,params)
    local C,listC,dmin,dmax,S,D,n,H,i,j,npoints,herdata,Herdata,p1data,p2data,a1,a2,st,tmp,m,B,distinct_index,indexMat:
    n:=nops(vars):
    m:=nops(params):
    npoints := 2 * n * (degree(F[1])-1) + 1:
    printf("Interpolate with %d points \n",npoints);
    
    # associated basis B
    B := QuotientBasis(F,vars,params,0,tdeg(op(vars))):
    D := nops(B):
    printf("Basis : ");
    print(B);

    # entries to be computed
    distinct_index, indexMat:=rmvDup(B):
    printf("Compute %d entries \n", nops(distinct_index));

    C := [seq(seq(vars[i]*B[j],i=1..n),j=1..D)]:
    C := convert(convert(C,set) minus convert(B,set),list): ### construct C
    C := sortBasis(C,tdeg(op(vars))): ### sort C
    dmin := degree(C[1]):
    dmax := degree(C[-1]):
    S := nops(C):

    # separate C by degree
    listC := [seq([],i=1..(dmax-dmin+1))]:
    for i from 1 to S do:
        listC[degree(C[i])-dmin+1] := [op(listC[degree(C[i])-dmin+1]),C[i]]:
    end do:

    # compute parametric Hermite matrix (only distinct entries)
    st:=time():
    H:=liftMatrix(F,npoints,params,B,D,vars,n,C,S,listC,distinct_index):
    print(time()-st);

    # fill in other entries
    for i from 1 to D do:
        for j from i to D do:
            if indexMat[i,j] <> 0 then:
                H[i,j]:=H[op(indexMat[i,j])]:
            end if:
        end do:
    end do:
    for i from 1 to D do:
        for j from 1 to i-1 do: 
            H[i,j] := H[j,i]:
        end do:
    end do:

    return H:
end proc: