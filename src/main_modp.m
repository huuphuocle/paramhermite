naiveHermite_modp:=proc(F,vars,params,char)
    local n,lm,q,B,d,gb,multmat,matB,H,i,j,ind,Vtr1,
        distinct_index,indexMat,gblm,r,k,checked,b:
    n:=nops(vars):
    gblm:=FGb[fgb_gbasis_lm](F,char,vars,params):
    
    # Construct the basis of the quotient ring
    lm:=subs([seq(params[i]=1,i=1..nops(params))],map(primpart,gblm[2])):
    lm:=FGb[fgb_gbasis](lm,2,[],vars):
    q:=queue[new](1):
    B:=[]:
    checked:=[]:
    while not queue[empty](q) do:
        b:=queue[dequeue](q):
        if not member(b,checked) then:
            r:=MonomialDiv(b,lm):
            checked:=[op(checked),b]:
            if r <> 0 then:
                B:=[op(B),b]:
                for i from 1 to n do:
                    queue[enqueue](q,expand(vars[i]*b)):
                end do:
            end if:
        end if:
    end do:
    B:=sort(B,(a,b)->Groebner[TestOrder](a,b,tdeg(op(vars)))):
    d:=nops(B):

    gb:=gblm[1]:
    multmat:=[seq(computeMatrix(vars[i],B,gb,vars,params,char),i=1..n)]:
    matB := matrixBasis_modp(multmat,B,d,vars,n,char):
    # matB:=[seq(map(p->Expand(p) mod char,matB[i]),i=1..d)]:
    H:=Matrix(d):
    
    for j from 1 to d do:
        H[1,j]:=LinearAlgebra[Trace](matB[j]):
    end do:

    Vtr1 := [seq(H[1,i],i=1..d)]:
    
    distinct_index, indexMat:=rmvDup(B):

    distinct_index:=distinct_index[(d+1)..-1]:

    for ind in distinct_index do:
        for j from 1 to d do:
            H[op(ind)] := H[op(ind)] + matB[ind[1]][ind[2],j]*Vtr1[j]:
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

directHermite_modp:=proc(F,vars,params,char)
    local st,gb,lm,B,M,i,j,l,computed,computed_index,t,bb: 
    
    st := time():
    B := QuotientBasis(F,vars,params,char,tdeg(op(vars))):
    st:=time()-st:
    l:=nops(B):
    printf("Basis ( %d elements ): ",l);
    printf("Elapsed time : %f\n",time()-st);
    
    st := time():
    gb := FGb:-fgb_gbasis(F,char,vars,params):
    printf("GBasis: %f \n",time()-st);
    
    M := Matrix(l):

    computed:=B:
    computed_index:=[seq([1,i],i=1..l)]:


    M[1,1]:=l:
    st := time():
    for j from 2 to l do:
        M[1,j]:=computeEntry(B[j],B,gb,vars,params,char,tdeg(op(vars))) mod char:
    end do:

    for i from 2 to l do:
        for j from i to l do:
            bb:=expand(B[i]*B[j]):
            if member(bb,computed,'t') then:
                M[i,j] := M[computed_index[t][1],computed_index[t][2]]:
            else:
                M[i,j] := computeEntry(bb,B,gb,vars,params,char,tdeg(op(vars))) mod char:
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
    st:=time()-st:
    printf("Distinct entries: %d/%d \n",nops(computed),l*(l+1)/2);
    printf("Total time : %f \n",time()-st);
    return M:
end proc:


### Compute multiplication matrices of x_i using row reduction
multMatrices_modp:=proc(gb,lt,B,D,vars,n,F,S,lF,char)
    local ll, NF, dmin, dmax, T, i, j, k, ind, ind2, ind3, M, d:
    dmin := degree(F[1]):
    dmax := degree(F[-1]):
    ll := [op(F),op(B)]:
    NF := [seq(0,i=1..(D+S))]:
    for i from 1 to D do:
        NF[i] := B[i]:
    end do:
    M := Matrix(S,D+S):
    ind := 1:
    for d from 1 to dmax-dmin+1 do:
        ind2 := ind:
        for i from 1 to nops(lF[d]) do:
            if member(lF[d][i],lt,'t') then:
                for j from 1 to D+S do:
                    M[ind,j] := mulCoeff(gb[t],ll[j],vars):
                end do:
            else:
                for k from 1 to n do:
                    if member(normal(lF[d][i]/vars[k]),lF[d-1],'t') then:
                        for j from 1 to D+S do:
                            M[ind,j] := mulCoeff(expand(lF[d][i]-vars[k]*NF[ind3+t-1+D]),ll[j],vars):
                        end do:
                        break;
                    end if:
                end do:
            end if:
            ind := ind + 1:
        end do:
        M := LinearAlgebra[Modular]:-Mod(char, M, integer[8]):
        LinearAlgebra[Modular]:-RowReduce(char,M,S,D+S,S,0,0,0,0,0,true):
        for i from ind2 to ind-1 do:
            for j from 1 to D do:
                NF[i+D] := NF[i+D] - M[i,S+j]*B[j]:  ### need to change this
            end do:
        end do:
        ind3 := ind2:
    end do:
    T := [seq(Matrix(D),k=1..n)]:
    ll := [op(B),op(F)]:

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

# Compute multiplication tables of B: L_1,...,L_d

matrixBasis_modp:=proc(multmat,B,d,vars,n,char)
    local matB,i,var,t,t1,t2:
    matB := Array(1..d):
    matB[1] := LinearAlgebra[Modular]:-Identity(char,d,integer):
    
    for i from 2 to d do:
        var := convert(indets(B[i]),list):
        member(var[1],vars,'t1'):
        member(normal(B[i]/var[1]),B,'t2'):
        matB[i] := LinearAlgebra[Modular]:-Multiply(char,multmat[t1],matB[t2]):
    end do:
    
    return matB:
end proc:

# Compute each specialized Hermite matrix

specHM_modp:=proc(F,B,D,vars,n,C,S,listC,char,distinct_index)
    local gb,lt,H,i,j,k,multmat,matB,Vtr1,st:
    
    gb:=FGb:-fgb_gbasis_lm(F,char,[],vars):
    lt:=gb[2]:
    gb:=gb[1]:
    #st:=time():
    multmat:=multMatrices_modp(gb,lt,B,D,vars,n,C,S,listC,char):
    #st:=time():
    matB := matrixBasis_modp(multmat,B,D,vars,n,char):
    H:=Matrix(D):
    
    for j from 1 to D do:
        H[1,j]:=LinearAlgebra[Trace](matB[j]):
    end do:

    Vtr1 := [seq(H[1,k],k=1..D)]:
    for i in distinct_index do:
        for k from 1 to D do:
            H[op(i)] := H[op(i)] + matB[i[1]][k,i[2]]*Vtr1[k]:
        end do:
        H[op(i)]:=H[op(i)] mod char:
    end do:
    
    return H:
end proc:

# Interpolation

liftMatrix_modp:=proc(F,npoints,params,char,B,D,vars,n,C,S,listC,distinct_index)
	local xdata, poldata, H, x, i, j, k, st, index_short:
    
    # initialize H (Hermite matrix)
	H:=Matrix(D):

    # random points
    #xdata:=GenerateData(npoints,char):
    xdata:=[seq(i,i=1..npoints)]:
	poldata:=[]:
    index_short:=distinct_index[(D+1)..-1]:

    if nops(params) = 0 then:
        H:=specHM_modp(F,B,D,vars,n,C,S,listC,char,index_short):
        return H:
    end if:
    
    # evaluation/ interpolation
    if nops(params)=1 then:
        # one variable left 
        
        # evaluation : call to specHM
        #st:=time():
        poldata:=[seq(specHM_modp(subs(params[1]=x,F),B,D,vars,n,C,S,listC,char,index_short),x in xdata)]:
        #printf("specHM : ");
		#print(time()-st);

        # interpolate 

        for i in distinct_index do:
            H[op(i)]:=Interp(xdata,[seq(poldata[k][op(i)],k=1..npoints)],params[1]) mod char:
        end do:
    else:
        # more than one variable : recursive call to liftMatrix
        for x in xdata do:
            poldata:=[op(poldata),liftMatrix_modp(subs(params[1]=x,F),npoints,params[2..-1],char,B,D,vars,n,C,S,listC,distinct_index)]:
        end do:
        
        for i in distinct_index do:
            H[op(i)]:=mInterp(xdata,[seq(poldata[k][op(i)],k=1..npoints)],params[2..-1],params[1],char):
        end do:
    end if:

    return H:
end proc:

HermiteParams_modp:=proc(F,vars,params,char)
    local C,listC,dmin,dmax,S,D,n,H,i,j,npoints,herdata,Herdata,p1data,p2data,a1,a2,st,tmp,m,B,distinct_index,indexMat:
    n:=nops(vars):
    m:=nops(params):
    npoints := 2 * n * (degree(F[1])-1) + 1:
    printf("Interpolate with %d points \n",npoints);
    
    # associated basis B
    B := QuotientBasis(F,vars,params,char,tdeg(op(vars))):
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
    H:=liftMatrix_modp(F,npoints,params,char,B,D,vars,n,C,S,listC,distinct_index):
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