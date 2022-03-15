LEXMatrix:=proc(F,vars,params,char)
    local g,B,M,i,j,l: 
    l := mul(degree(F[i]),i=1..nops(F)):
    B := [seq(vars[-1]^i,i=0..(l-1))]:
    g := FGb[fgb_gbasis_elim](F,char,vars[1..-2],[vars[-1],op(params)])[1]:
    M := Matrix(l):
    if char=0 then:
        for i from 1 to l do:
            M[1,i]:=expand(add(coeff(rem(vars[-1]^(j+i-2),g,vars[-1]),vars[-1],j-1),j=1..l)):
            M[i,l]:=expand(add(coeff(rem(vars[-1]^(j+l+i-3),g,vars[-1]),vars[-1],j-1),j=1..l)):
        end do:
    else:
        for i from 1 to l do:
            M[1,i]:=Expand(add(coeff(Rem(vars[-1]^(j+i-2),g,vars[-1]) mod char,vars[-1],j-1),j=1..l)) mod char:
            M[i,l]:=Expand(add(coeff(Rem(vars[-1]^(j+l+i-3),g,vars[-1]) mod char,vars[-1],j-1),j=1..l)) mod char:
        end do:
    end if:
    for i from 2 to l do:
		for j from 1 to l-i+1 do:
			M[i,j]:=M[1,i+j-1]:
		end do:
		for j from l+2-i to l-1 do:
			M[i,j]:=M[i+j-l,l]:
		end do:
	end do:
    return M:
end proc:

#P * H * P^t = Han
computeP := proc(F,vars,params,char)
    local gb,lm,B2,B1,l,M,i,j,tmp,p:
    gb := fgb_gbasis_lm(subs([seq(p=1,p in params)],F),char,[],vars):
    lm := map(p->normal(p/lcoeff(p)),gb[2]):
    B2 := QuotientBasis(F,vars,params,char):
    l := nops(B2):
    B1 := [seq(vars[-1]^i,i=0..(l-1))]:
    gb := fgb_gbasis(F,char,vars,params):
    M := Matrix(l):
    for i from 1 to l do:
        tmp := Groebner:-NormalForm(B1[i],gb,tdeg(op(vars)),characteristic=char):
        for j from 1 to l do:
            M[i,j] := mulCoeff(tmp,B2[j],vars):
            # printf(".");
        end do:
    end do:
    return M:
end proc: