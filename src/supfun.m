# function to reduce in a monomial basis (change it using HASH?)
MonomialDiv:=proc(m,lm)
    local b:
    # return Groebner[NormalForm](b,lm,tdeg(op(vars)),characteristic=2):
    for b in gb do:
        if divide(m,b) then:
            return true:
        end if:
    end do:
    return false:
end proc:

# list of coefficients of a multivariate polynomial
mycoeffs := proc(pol::{polynom},var::{list})
	local deg,j:
	deg := degree(pol,var):
	return [seq(coeff(pol,var,deg+1-j),j=1..(deg+1))]:
end proc:

# return the coefficient of m in f

mulCoeff := proc(f,m,vars)
    local c,i,k;
    c:=[coeffs(f,vars,k)];
    if member(m,[k],i) then:
		return c[i]:
	else:
		return 0:
	end if:
end:

# generate a list of l different elements of GF(char) for evaluation
GenerateData:=proc(npoints::{list},char::{integer})
	local xdata,roll,t,i:
	xdata := []:
	if char = 0 then:
		i := 1:
		while i <= npoints do:
			t := rand():
			if not (t in xdata) then:
				xdata := [op(xdata),t]:
				i := i+1:
			end if:
		end do:
	else:
		roll := rand(char):
		i := 1:
		while i <= npoints do:
			t := roll():
			if not (t in xdata) then:
				xdata := [op(xdata),t]:
				i := i+1:
			end if:
		end do:
	end if:
	return xdata:
end proc:

# compute a basis of the vector space K[X]/I
QuotientBasis:=proc(F::{list},vars::{list},params::{list},char::{integer},ord:=tdeg(op(vars)))
        local n,q,B,b,v,checked,flag,m,rvars:
    rvars:=ListTools[Reverse](vars):
    n:=nops(vars):
    #build all the monomials here
    q:=queue[new](1):
    B:=[]:
    checked:=[]:
    while not queue[empty](q) do:
        b:=queue[dequeue](q):
        if not member(b,checked) then:
            checked:=[op(checked),b]:
            flag:=true:
            for m in lm do:
                if divide(b,m) then:
                    flag:=false:
                    break:
                end if:
            end do:
            if flag then:
                B:=[op(B),b]:
                for v in rvars do:
                    queue[enqueue](q,v*b):
                end do:
            end if:
        end if:
    end do:
    # B := sort(B,(a,b)->Groebner[TestOrder](a,b,ord)):
    return B:
end proc:

# multivariate interpolation
mInterp:=proc(xdata,poldata,vars,var,char)
    local npoints,k,terms,ydata,f,i,j:
    npoints:=nops(xdata):
    (*if npoints<>nops(poldata) then:
        error:
        return:
    end if:*)
    k:=coeffs(expand(poldata[1]),vars,'terms'):
	terms:=[terms]:
    ydata:=[[k]]:
    for i from 2 to npoints do:
        k:=coeffs(expand(poldata[i]),vars):
        ydata:=[op(ydata),[k]]:
    end do:
    f:=0:
	if char=0 then:
		for i from 1 to nops(terms) do:
			f:=f+terms[i]*interp(xdata,[seq(ydata[j][i],j=1..npoints)],var):
		end do:
	else:
		for i from 1 to nops(terms) do:
			f:=f+terms[i]*(Interp(xdata,[seq(ydata[j][i],j=1..npoints)],var) mod char):
		end do:
		f:=f mod char:
	end if:
    return f:
end proc:

# compute the entry Trace(L_m)
computeEntry := proc(m,B,gb,vars,params,char,ord:=tdeg(op(vars)))
    local s,b,tmp,i:
    s := subs([seq(vars[i]=0,i=1..nops(vars))],Groebner[NormalForm](m,gb,ord,characteristic=char)):
    for b in B[2..-1] do:
        tmp := expand(m*b):
        if not tmp in B then:
            s := s + mulCoeff(Groebner[NormalForm](tmp,gb,ord,characteristic=char),b,vars):
        end if:
    end do:
    return s:
end proc:

# compute the matrix of the map L_m

computeMatrix := proc(m,B,gb,vars,params,char,ord:=tdeg(op(vars)))
    local i,j,k,d,M,s,tmp,B2:
    B2:=map(primpart,B):
    d:=nops(B):
    M:=Matrix(d):
    for i from 1 to d do:
        tmp:=m*B[i]:
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

# remove duplicated entries in the matrix
rmvDup:=proc(B::{list})
    local computed,distinct_index,d,i,j,tmp,t,M:
    d:=nops(B):
    computed := B:
    distinct_index:=[seq([1,j],j=1..d)]: # store indices
    M:=Matrix(d):
    for i from 2 to d do:
        for j from i to d do:
            tmp:=expand(B[i]*B[j]):
            if member(tmp,computed,'t') then:
                M[i,j]:=distinct_index[t]:
            else:
                computed:=[op(computed),tmp]:
                distinct_index := [op(distinct_index),[i,j]]:
            end if:
        end do:
    end do:
    return distinct_index, M:
end proc:

# obtain the index of the principal minor to be computed
indexDV:=proc(H,params:=convert(indets(H),list))
    local Hx, d, dx, tmp, p:
    Hx:=subs([seq(p=rand()*t+rand(), p in params)],H):
    d:=LinearAlgebra[RowDimension](Hx)-1:
    dx:=Determinant(Hx):
    dx:=quo(dx,gcd(dx,diff(dx,t)),t):
    while d > 1 do:
        tmp:=Determinant(SubMatrix(Hx,1..d,1..d),method=multivar):
        if not divide(tmp,dx) then:
            return d+1:
        end if:
    end do:
end proc:

# remove redundant elements in GB of DRL(x) > DRL(y)
cleanGB:=proc(gblm,vars,params)
    local lm, lmt, gb, i, v:
    lm:=map(primpart,subs([seq(v = 1, v in params)],gblm[2])):
    lmt:=fgb_gbasis(lm,2,[],vars):
    gb:=[]:
    for i from 1 to nops(lm) do:
        if member(lm[i],lmt) then:
            gb:=[op(gb),gblm[1][i]]:
        end if:
    end do:
    return gb:
end proc:

np_factor:=proc(np)
    local res, fnp, p, flag, i, j:
    res:=[]:
    fnp:=map(factors,np):
    for p in fnp do:
        for i from 1 to nops(p[2]) do:
            flag:=true:
            for j from 1 to nops(res) do:
                if divide(res[j],p[2][i][1]) then:
                    flag:=false:
                    break;
                end if:
            end do: 
            if flag then:
                res:=[op(res),p[2][i][1]]:
            end if:
        end do:
    end do:
    return res:
end proc:

# compute the number of sign changes in a list of integers
sign_variant:=proc(ll::{list})
    local cmpt,i:
    cmpt:=0:
    for i from 1 to nops(ll)-1 do:
        if(ll[i+1]*ll[i] < 0) then:
            cmpt:=cmpt+1:            
        end if:
    end do:
    return cmpt:
end proc:

# count the number of real solutions at the point sp
NumberOfRealSolutions:=proc(sp,H)
    local n,k,ld,i,M:
    M:=subs(sp,H):
    n:=LinearAlgebra[RowDimension](H):
    ld:=map(signum,[seq(LinearAlgebra[Determinant](LinearAlgebra[SubMatrix](M,1..i,1..i)),i=1..n)]):
    if member(0,ld) then:
        error "Sign sequence contains zero(s)!";
    else:
        k:=sign_variant(ld):
        return n-2*k:
    end if:
end proc:

CleanFactors:=proc(f,np)
    description "Take a polynomial f and a list of polynomials np. Return factors of f with np removed";
    local fac_f, p, i, j:
    fac_f:=[seq(p[1], p in factors(f)[2])]:
    for i from 1 to nops(fac_f) do:
        for j from 1 to nops(np) do:
            if divide(np[j],fac_f[i]) then:
                fac_f[i]:=1: 
            end if:
        end do:
    end do:
    return fac_f:
end proc:

euler_candy:=proc(d,k)
    local i,tmp,res,v:
    if k = 1 then:
        return [[d]]:
    fi:
    res:=[]:
    for i from 1 to d-k+1 do:
        tmp:=euler_candy(d-i,k-1):
        res:=[op(res),seq([i,op(v)], v in tmp)]:
    od:
    return res:
end proc:

# sign patterns correspond to r real roots
signPatterns:=proc(d,r)
    local k,i,j,p,part,res:
    k:=(d-r)/2:
    part:=euler_candy(d,k+1):
    res:=[seq([seq(seq((-1)^(i-1), j=1..p[i]), i=1..nops(p))][2..-1], p in part)]:
end proc:

GenericDimension:=proc(sys,vars,params)
    local sysf, v, gb, hil:
    sysf:=subs([seq(v=rand(),v in params)],sys):
    gb:=FGb[fgb_gbasis](sysf,65519,[],vars):
    hil:=FGb[fgb_hilbert](gb,65519,[],vars,'u'):
    return hil[2]:
end proc:

allMinors := proc(M,c)
    local lminors,lc1,lc2,i,j:
    lminors := []:
    lc1 := combinat:-choose(LinearAlgebra[RowDimension](M),c):
    lc2 := combinat:-choose(LinearAlgebra[ColumnDimension](M),c):
    for i in lc1 do:
        for j in lc2 do:
            lminors := [op(lminors),LinearAlgebra[Determinant](LinearAlgebra[SubMatrix](M,i,j))]:
        end do:
    end do:
	return lminors:
end proc:

# every input sorted in the order ord
allMatrices := proc(B,gb,lm,vars,params,char,ord:=tdeg(op(vars)))
    local n,d,m,F,M,s,tmp,i,j,k,normalForms,list_NormalForms:
    n:=nops(vars):
    d:=nops(B):
    F:=[seq(seq(vars[i]*B[j],i=1..n),j=1..d)]:
    F:=sort(F):
    normalForms:=table();
    for m in F do:
        if member(m,B) then:
            normalForms[m]:=1:
        else:
            if member(m,lm,'j') then:
                list_NormalForms[m]:=[seq(k=1..d)]:
                M[i]:=Vector([seq(mulCoeff(s,B2[k],vars)/content(B[k]),k=1..d)]):
            else:

            end if:
        end if:
    end do:
    for i from 1 to n do:
        M:=Matrix(d):
        for j from 1 to d do:
            tmp := expand(vars[i]*B[j]):
            normalForms[tmp]:
        end do:
    end do:
end proc: