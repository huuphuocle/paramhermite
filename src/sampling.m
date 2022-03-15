mipoint:=proc(r, s)
    if ceil(r) < s and r < ceil(r) then:
        return ceil(r):
    else:
        if floor(s) > r and floor(s) < s then:
            return floor(s):
        else:
            return (r+s)/2:
        end if:
    end if:
end proc:

isol:=proc(pol,var,mth:='ABND')
    local points, rpoints, i:
    points := map(rhs,RootFinding[Isolate](primpart(pol),var,output='interval',method=mth)):
    if points = [] then:
        rpoints := [0]:
    else:
        # rpoints := [floor(points[1][1])-1,seq((points[i][2]+points[i+1][1])/2, i=1..(nops(points)-1)),floor(points[-1][2])+1]:
        rpoints := [floor(points[1][1])-1,seq(mipoint(points[i][2],points[i+1][1]), i=1..(nops(points)-1)),floor(points[-1][2])+1]:
    end if:
    return rpoints:
end proc:

# sample points for f <> 0 with one parameter

ppc1:=proc(f,var)
    local rpoints:
    rpoints:=isol(f,var,'RS'):
    return map(p->[var = p],rpoints):
end proc:

# sample points for f <> 0 with two parameters

ppc2:=proc(w,params)
    local fw,res,w0,rpoints,lpoints,v,t,x,tmp,p,i,j,st:

    fw:=map(p->p[1],factors(w)[2]):
    res:=1:
    for i from 1 to nops(fw) do:
        tmp:=resultant(fw[i],diff(fw[i],params[1]),params[1]):
        res:=res*mul(v[1], v in sqrfree(tmp)[2]):
    end do:
    for i from 1 to nops(fw)-1 do:
        for j from i+1 to nops(fw) do:
            tmp:=resultant(fw[i],fw[j],params[1]):
            res:=res*mul(v[1], v in sqrfree(tmp)[2]):
        end do:
    end do:
    w0:=mul(v[1], v in sqrfree(res)[2]):
    rpoints:=ppc1(w0,params[2]):
    lpoints:=[]:
    for x in rpoints do:
        tmp := ppc1(subs(x,w),params[1]):
        lpoints:=[op(lpoints),seq([op(p),op(x)], p in tmp)]:
    end do:
    return lpoints:
end proc:

classify:=proc(w,params,g:=params[-1])
    local gb,v,l,u,rpoints:
    gb:=FGb[fgb_matrixn_radical2]([seq(diff(l*w-g,v), v in params)],[w,u-g],0,[l],[op(params),u],0,{"verb"=3}):
    rpoints:=isol(gb[-1],u):
    return rpoints:
end proc:

# sample points for w <> 0

SamplePoints:=proc(w,params)
    local n,rpoints,lpoints,tmp,p,i,lf,v,u,roll,wu,u0:
    n:=nops(params):
    if n = 1 then:
        return ppc1(w,params[1]):
    elif n = 2 then:
        return ppc2(w,params):
    else:
        roll:=rand(1..4): # choice of linear form
        lf:=add(roll()*v, v in params[1..-2]):
        rpoints:=classify(w,params,lf+params[-1]):
        printf("\n %a : %d \n",params[-1],nops(rpoints));
        
        wu:=subs(params[-1]=u-lf,w):
        lpoints:=[]:
        for u0 in rpoints do:
            tmp := SamplePoints(subs(u=u0,wu),params[1..-2]):
            lpoints:=[op(lpoints),seq([op(p), params[-1] = u0-subs(p,lf)], p in tmp)]:
        end do:

    end if:
    return lpoints:
end proc: