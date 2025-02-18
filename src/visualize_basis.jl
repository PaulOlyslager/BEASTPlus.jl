using Plotly

function combine_doubles(pos,val)
    pos_new = []
    val_new = []
    for (i,posi,vali) in zip(1:length(pos),pos,val)
        if approx_in(posi,pos_new)
            o = approx_equal.(pos_new ,Ref(posi))
      
            val_new[o] .+= Ref(vali)
        else
            push!(pos_new,posi)
            push!(val_new,vali)
        end
    end
    return pos_new,val_new
end
function approx_equal(v1,v2;tol=100*sqrt(eps()))
    return norm(v1-v2) < tol
end

function approx_in(p,vect;tol =100* sqrt(eps()))
    
    for v in vect
        
        if norm(v-p) < tol
            return true
        end
    end
    return false
end


function showfn(space,i;size=1,N = 1, points=[i for i in 0:N]/N)
    geo = geometry(space)
    T = coordtype(geo)
    X = T[]
    Y = T[]
    Z = T[]
    U = T[]
    V = T[]
    W = T[]
    c = []
    vals = []
    barpoints = [(@SVector [i,j]) for i in points for j in points if i+j<(1+sqrt(eps()))]
    npoints=length(barpoints)
    for tr in 1:numcells(geo)
        chrt = chart(geo,tr)
        mps = neighborhood.(Ref(chrt),barpoints)
        valsi = [(@SVector [0.0,0.0,0.0]) for i in 1:npoints]
        
        for sh in space.fns[i]
            sh.cellid != tr && continue
            for (j,mp) in enumerate(mps)
                valsi[j] += sh.coeff*refspace(space)(mp)[sh.refid].value
            end

        end
        vals = [vals;valsi]
        ci = cartesian.(mps)
        c = [c;ci]
        #vals *= size

    end

    c,vals = combine_doubles(c,vals)
    for v in vals
        push!(U,v[1])
        push!(V,v[2])
        push!(W,v[3])
    end
    
    for x in c
        push!(X,x[1])
        push!(Y,x[2])
        push!(Z,x[3])
    end
    Plotly.cone(x=X,y=Y,z=Z,u=U,v=V,w=W,sizemode="absolute",sizeref=size)
end


function show_from_fns(space,fns;size=1,N = 1, points=[i for i in 0:N]/N)
    geo = geometry(space)
    T = coordtype(geo)
    X = T[]
    Y = T[]
    Z = T[]
    U = T[]
    V = T[]
    W = T[]
    barpoints = [(@SVector [i,j]) for i in points for j in points if i+j<(1+sqrt(eps()))]
    npoints=length(barpoints)
    for tr in 1:numcells(geo)
        chrt = chart(geo,tr)
        mps = neighborhood.(Ref(chrt),barpoints)
        vals = [(@SVector [0.0,0.0,0.0]) for i in 1:npoints]
        
        for sh in fns
            sh.cellid != tr && continue
            for (j,mp) in enumerate(mps)
                vals[j] += sh.coeff*refspace(space)(mp)[sh.refid].value
            end

        end
        #vals *= size
        for v in vals
            push!(U,v[1])
            push!(V,v[2])
            push!(W,v[3])
        end
        c = cartesian.(mps)
        for x in c
            push!(X,x[1])
            push!(Y,x[2])
            push!(Z,x[3])
        end
    end
    Plotly.cone(x=X,y=Y,z=Z,u=U,v=V,w=W,sizemode="absolute",sizeref=size)
end

function show_div_from_fns(space,fns;size=1,N = 1, points=[i for i in 0:N]/N)
    geo = geometry(space)
    T = coordtype(geo)
    X = T[]
    Y = T[]
    Z = T[]
    U = T[]
    V = T[]
    W = T[]
    barpoints = [(@SVector [i,j]) for i in points for j in points if i+j<(1+sqrt(eps()))]
    npoints=length(barpoints)
    for tr in 1:numcells(geo)
        chrt = chart(geo,tr)
        mps = neighborhood.(Ref(chrt),barpoints)
        vals = [(@SVector [0.0,0.0,0.0]) for i in 1:npoints]
        
        for sh in fns
            sh.cellid != tr && continue
            for (j,mp) in enumerate(mps)
                vals[j] += sh.coeff*refspace(space)(mp)[sh.refid].divergence*normal(mp)
            end

        end
        vals *= size
        for v in vals
            norm(v) < 1e-10 && continue
            push!(U,v[1])
            push!(V,v[2])
            push!(W,v[3])
        end
        c = cartesian.(mps)
        for (i,x) in enumerate(c)
            norm(vals[i]) < 1e-10 && continue
            push!(X,x[1])
            push!(Y,x[2])
            push!(Z,x[3])
        end
    end
    #Plotly.cone(x=X,y=Y,z=Z,u=U,v=V,w=W,sizemode="absolute",sizeref=size)
    Plotly.scatter(x = X.+U,y=Y.+V,z=Z.+W,mode="markers",marker=attr(size=2),type="scatter3d")
end
function show_scallar_from_fns(space,fns;size=1,N = 1, points=[i for i in 0:N]/N)
    geo = geometry(space)
    T = coordtype(geo)
    X = T[]
    Y = T[]
    Z = T[]
    U = T[]
    V = T[]
    W = T[]
    barpoints = [(@SVector [i,j]) for i in points for j in points if i+j<(1+sqrt(eps()))]
    npoints=length(barpoints)
    for tr in 1:numcells(geo)
        chrt = chart(geo,tr)
        mps = neighborhood.(Ref(chrt),barpoints)
        vals = [(@SVector [0.0,0.0,0.0]) for i in 1:npoints]
        
        for sh in fns
            sh.cellid != tr && continue
            for (j,mp) in enumerate(mps)
                vals[j] += sh.coeff*refspace(space)(mp)[sh.refid].value*normal(mp)
            end

        end
        #vals *= size
        for v in vals
            push!(U,v[1])
            push!(V,v[2])
            push!(W,v[3])
        end
        c = cartesian.(mps)
        for x in c
            push!(X,x[1])
            push!(Y,x[2])
            push!(Z,x[3])
        end
    end
    #Plotly.cone(x=X,y=Y,z=Z,u=U,v=V,w=W,sizemode="absolute",sizeref=size)
    Plotly.scatter(x = X.+U,y=Y.+V,z=Z.+W,mode="markers",marker=attr(size=2),type="scatter3d")
end

function show_face_indices(space)
    geo = geometry(space)
    T = coordtype(geo)
    X = T[]
    Y = T[]
    Z = T[]
    names = []
    for i in 1:numcells(geo)
        c = cartesian(center(chart(geo,i)))
        push!(X,c[1])
        push!(Y,c[2])
        push!(Z,c[3])
        push!(names,string(i))
    end
    
    Plotly.scatter(x = X,y=Y,z=Z,mode="markers+text",marker=attr(size=2),type="scatter3d",text=names)
end
function show_vertex_indices(space)
    geo = geometry(space)
    T = coordtype(geo)
    X = T[]
    Y = T[]
    Z = T[]
    names = []
    for (i,c) in enumerate(vertices(geo))
        #c = cartesian(center(chart(geo,i)))
        push!(X,c[1])
        push!(Y,c[2])
        push!(Z,c[3])
        push!(names,string(i))
    end
    
    Plotly.scatter(x = X,y=Y,z=Z,mode="markers+text",marker=attr(size=2),type="scatter3d",text=names)
end

function show_basis_indices(space;sep=0.5)
    geo = geometry(space)
    T = coordtype(geo)
    X = T[0]
    Y = T[0]
    Z = T[5]
    names = ["none"]
    for i in 1:numfunctions(space)
        c = space.pos[i]
        push!(X,c[1])
        push!(Y,c[2])
        push!(Z,c[3]+sep*(-1)^i)
        push!(names,string(i))
    end
    
    Plotly.scatter(x = X,y=Y,z=Z,mode="markers+text",marker=attr(size=2),type="scatter3d",text=names)
end

function show_tangents(space;size=1)
    geo = geometry(space)
    T = coordtype(geo)
    X = T[]
    Y = T[]
    Z = T[]
    U = T[]
    V = T[]
    W = T[]
    #barpoints = [(@SVector [i,j]) for i in points for j in points if i+j<(1+sqrt(eps()))]
    barpoint = (@SVector [0.2,0.2])

    for tr in 1:numcells(geo)
        chrt = chart(geo,tr)
        t1,t2 = tangents.(Ref(neighborhood(chrt,barpoint)),1:2)
        x = cartesian(neighborhood(chrt,barpoint))
        
        for v in [2*normalize(t1),normalize(t2)]
            push!(U,v[1])
            push!(V,v[2])
            push!(W,v[3])
        end

        for i in 1:2
            push!(X,x[1])
            push!(Y,x[2])
            push!(Z,x[3])
        end
    end
    Plotly.cone(x=X,y=Y,z=Z,u=U,v=V,w=W,sizemode="absolute",sizeref=size)

end

function make_div_continuous(space,leftchart,centerchart,fns,edge,functions)
    geo = space.geo
    t1 = vertices(geo)[edge[2]].-vertices(geo)[edge[1]]
    points = [i/(length(functions)+1)*t1.+vertices(geo)[edge[1]] for i in 1:length(functions)]
    mpcenter = carttobary.(Ref(chart(space.geo,centerchart)),points)
    mpleft = carttobary.(Ref(chart(space.geo,leftchart)),points)
    rfs = refspace(space)
    centervals = zeros(length(functions))
    leftvals = zeros(length(functions))
    crfs = rfs.(mpcenter)
    lrfs = rfs.(mpleft)
    for sh in fns
        sh.cellid != leftchart && continue
        ldivs = sh.coeff*getindex.(getindex.(lrfs,sh.refid),2)
        leftvals .+= ldivs
    end
    for sh in fns
        sh.cellid != centerchart && continue
        cdivs = sh.coeff*getindex.(getindex.(crfs,sh.refid),2)
        centervals .+= cdivs
    end
    m = zeros((length(functions),length(functions)))
    for (i,mp) in enumerate(mpcenter)
        for (j,f) in enumerate(functions)
            for sh in f
                display(sh.cellid)
                sh.cellid != centerchart && continue
                println("added")
                m[i,j] += sh.coeff*rfs(mp)[sh.refid].divergence
            end
        end
    end
    y = leftvals-centervals
    display(m)
    c = m\y
    return c
end