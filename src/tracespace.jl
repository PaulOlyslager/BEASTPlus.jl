### definition of the normal trace of the nedelecc3d space
### definition of the tangential trace of the nedelecd3d space 

using BEAST
using SparseArrays
function ntrace(X::BEAST.NDLCCBasis{T}) where {T}
    facemesh = skeleton(X.geo,2)
    tetmesh = X.geo
    fns_new = typeof(X.fns[1])[]
    U = sparse(transpose(connectivity(facemesh,tetmesh,x->x)))
    rows = rowvals(U)
    vals = nonzeros(U)
    for fn in X.fns
        fn_new = typeof(fn[1])[]
        for shape in fn
            for i in nzrange(U,shape.cellid)
                faceid = rows[i]
                face_orient = sign(vals[i])
                local_index = abs(vals[i])
                ϕ = refspace(X)
                Q = normal_interpolate(BEAST.LagrangeRefSpace{T,1,3,0}(),chart(facemesh,faceid),ϕ,chart(tetmesh,shape.cellid))
              
                for (j,c) in enumerate(Q[shape.refid,1:end])
                    if abs(c) > sqrt(eps())
                        push!(fn_new,BEAST.Shape(faceid,j,face_orient*c*shape.coeff))
                    end

                end

            end
        end
        push!(fns_new,fn_new)
    end
    return BEAST.LagrangeBasis{1,0,3}(facemesh,unique_fn.(fns_new),X.pos)
end

function ttrace(X::BEAST.NDLCDBasis{T}) where {T}
    facemesh = skeleton(X.geo,2)
    tetmesh = X.geo
    fns_new = typeof(X.fns[1])[]
    U = sparse(transpose(connectivity(facemesh,tetmesh,x->x)))
    rows = rowvals(U)
    vals = nonzeros(U)
    for fn in X.fns
        fn_new = typeof(fn[1])[]
        for shape in fn
            for i in nzrange(U,shape.cellid)
                faceid = rows[i]
                face_orient = sign(vals[i])
                local_index = abs(vals[i])
                ϕ = refspace(X)
                Q = tangential_interpolate(BEAST.RTRefSpace{T}(),chart(facemesh,faceid),ϕ,chart(tetmesh,shape.cellid))
                for (j,c) in enumerate(Q[shape.refid,1:end])
                    if abs(c) > sqrt(eps())
                        push!(fn_new,BEAST.Shape(faceid,j,c*shape.coeff))
                    end

                end

            end
        end
        push!(fns_new,fn_new)
    end
    return BEAST.RTBasis(facemesh,unique_fn.(fns_new),X.pos)
end
function ttrace(X::BEAST.NDLCCBasis{T}) where {T}
    facemesh = skeleton(X.geo,2)
    tetmesh = X.geo
    fns_new = typeof(X.fns[1])[]
    U = sparse(transpose(connectivity(facemesh,tetmesh,x->x)))
    rows = rowvals(U)
    vals = nonzeros(U)
    for fn in X.fns
        fn_new = typeof(fn[1])[]
        for shape in fn
            for i in nzrange(U,shape.cellid)
                faceid = rows[i]
                face_orient = sign(vals[i])
                local_index = abs(vals[i])
                ϕ = refspace(X)
                Q = tangential_interpolate(BEAST.RTRefSpace{T}(),chart(facemesh,faceid),ϕ,chart(tetmesh,shape.cellid))
                for (j,c) in enumerate(Q[shape.refid,1:end])
                    if abs(c) > sqrt(eps())
                        push!(fn_new,BEAST.Shape(faceid,j,c*shape.coeff))
                    end

                end

            end
        end
        push!(fns_new,fn_new)
    end
    return BEAST.RTBasis(facemesh,unique_fn.(fns_new),X.pos)
end
function strace(X::BEAST.NDLCDBasis{T}) where {T}
    facemesh = skeleton(X.geo,2)
    tetmesh = X.geo
    fns_new = typeof(X.fns[1])[]
    U = sparse(transpose(connectivity(facemesh,tetmesh,x->x)))
    rows = rowvals(U)
    vals = nonzeros(U)
    for fn in X.fns
        fn_new = typeof(fn[1])[]
        for shape in fn
            for i in nzrange(U,shape.cellid)
                @assert length(nzrange(U,shape.cellid)) == 4
                faceid = rows[i]
                face_orient = sign(vals[i])
                local_index = abs(vals[i])
                ϕ = refspace(X)
                Q = ntangential_interpolate(BEAST.BDMRefSpace{T}(),chart(facemesh,faceid),ϕ,chart(tetmesh,shape.cellid))
                for (j,c) in enumerate(Q[shape.refid,1:end])
                    if abs(c) > sqrt(eps())
                        push!(fn_new,BEAST.Shape(faceid,j,face_orient*c*shape.coeff))
                    end
                end

            end
        end
        push!(fns_new,fn_new)
    end
    return BEAST.BDMBasis(facemesh,unique_fn.(fns_new),X.pos)
end
function strace(X::BEAST.NDLCCBasis{T}) where {T}
    facemesh = skeleton(X.geo,2)
    tetmesh = X.geo
    fns_new = typeof(X.fns[1])[]
    U = sparse(transpose(connectivity(facemesh,tetmesh,x->x)))
    rows = rowvals(U)
    vals = nonzeros(U)
    for fn in X.fns
        fn_new = typeof(fn[1])[]
        for shape in fn
            for i in nzrange(U,shape.cellid)
                @assert length(nzrange(U,shape.cellid)) == 4
                faceid = rows[i]
                face_orient = sign(vals[i])
                local_index = abs(vals[i])
                ϕ = refspace(X)
                Q = ntangential_interpolate(BEAST.BDMRefSpace{T}(),chart(facemesh,faceid),ϕ,chart(tetmesh,shape.cellid))
                for (j,c) in enumerate(Q[shape.refid,1:end])
                    if abs(c) > sqrt(eps())
                        push!(fn_new,BEAST.Shape(faceid,j,face_orient*c*shape.coeff))
                    end
                end

            end
        end
        push!(fns_new,fn_new)
    end
    return BEAST.BDMBasis(facemesh,unique_fn.(fns_new),X.pos)
end
function _trace(X::BEAST.LagrangeBasis{D,C,M,T,NF}) where {D,C,M,T,NF}
    facemesh = skeleton(X.geo,2)
    tetmesh = X.geo
    fns_new = typeof(X.fns[1])[]
    U = sparse(transpose(connectivity(facemesh,tetmesh,x->x)))
    rows = rowvals(U)
    vals = nonzeros(U)
    for fn in X.fns
        fn_new = typeof(fn[1])[]
        for shape in fn
            for i in nzrange(U,shape.cellid)
                faceid = rows[i]
                face_orient = sign(vals[i])
                local_index = abs(vals[i])
                ϕ = refspace(X)
                Q = BEAST.interpolate(BEAST.LagrangeRefSpace{T,1,3,0}(),chart(facemesh,faceid),ϕ,chart(tetmesh,shape.cellid))
                for (j,c) in enumerate(Q[shape.refid,1:end])
                    if abs(c) > sqrt(eps())
                        push!(fn_new,BEAST.Shape(faceid,j,c*shape.coeff))
                    end

                end

            end
        end
        push!(fns_new,fn_new)
    end
    return BEAST.LagrangeBasis{1,-1,length(X.fns)}(facemesh,unique_fn.(fns_new),X.pos)
end
function trace(X::BEAST.LagrangeBasis{D,C,M,T,NF}) where {D,C,M,T,NF}
    facemesh = skeleton(X.geo,2)
    tetmesh = X.geo
    fns_new = typeof(X.fns[1])[]
    U = sparse(transpose(connectivity(facemesh,tetmesh,x->x)))
    rows = rowvals(U)
    vals = nonzeros(U)
    for fn in X.fns
        fn_new = typeof(fn[1])[]
        for shape in fn
            for i in nzrange(U,shape.cellid)
                faceid = rows[i]
                face_orient = sign(vals[i])
                local_index = abs(vals[i])
                ϕ = refspace(X)
                Q = BEAST.interpolate(BEAST.LagrangeRefSpace{T,1,3,0}(),chart(facemesh,faceid),ϕ,chart(tetmesh,shape.cellid))
                for (j,c) in enumerate(Q[shape.refid,1:end])
                    if abs(c) > sqrt(eps())
                        push!(fn_new,BEAST.Shape(faceid,j,face_orient*c*shape.coeff))
                    end

                end

            end
        end
        push!(fns_new,fn_new)
    end
    return BEAST.LagrangeBasis{1,-1,length(X.fns)}(facemesh,unique_fn.(fns_new),X.pos)
end
ntimestrace(X) = _BasisTimes(n,trace(X))
# ### coordinate function accepts meshpoint and evaluates refspace in chart
# function evalbar(ϕ::BEAST.RefSpace,chart)
#     return function (mp)
#         x = cartesian(mp)
#         b = CompScienceMeshes.carttobary(chart,x)
#         return ϕ(b)
#     end
# end

function normal_interpolate(interpolant::BEAST.RefSpace, chart1, interpolee::BEAST.RefSpace, chart2)
    function fields(p)
        x = cartesian(p)
        v = carttobary(chart2, x)
        r = neighborhood(chart2, v)
        no = normal(chart1)
        fieldvals = [dot(no,f.value) for f in interpolee(r)]
    end

    BEAST.interpolate(fields, interpolant, chart1)
end
function tangential_interpolate(interpolant::BEAST.RefSpace, chart1, interpolee::BEAST.RefSpace, chart2)
    function fields(p)
        x = cartesian(p)
        v = carttobary(chart2, x)
        r = neighborhood(chart2, v)
        no = normal(chart1)
        fieldvals = [-no × (no× f.value) for f in interpolee(r)]
    end

    BEAST.interpolate(fields, interpolant, chart1)
end
function ntangential_interpolate(interpolant::BEAST.RefSpace, chart1, interpolee::BEAST.RefSpace, chart2)
    function fields(p)
        x = cartesian(p)
        v = carttobary(chart2, x)
        r = neighborhood(chart2, v)
        no = normal(chart1)
        fieldvals = [no × f.value for f in interpolee(r)]
    end

    BEAST.interpolate(fields, interpolant, chart1)
end

function BEAST.interpolate(fields, interpolant::BEAST.BDMRefSpace, chart)
    Q = map(faces(chart)) do face
        p1 = vertices(face)[1]
        p2 = vertices(face)[2]
        
        u1 = carttobary(chart, p1)
        u2 = carttobary(chart, p2)
        q1 = neighborhood(chart, u1)
        q2 = neighborhood(chart, u2)
        n = normal(q1)

        # minus because in CSM the tangent points towards vertex[1]
        t = -tangents(center(face),1)
        m = cross(t,n)

        fieldvals1 = fields(q1)
        q11 = [dot(fv,m) for fv in fieldvals1]
        fieldvals2 = fields(q2)
        q22 = [dot(fv,m) for fv in fieldvals2]
        q = hcat(q11,q22)
    end

    return hcat(Q...)
end


function BEAST.lagrangec0d1(mesh, vertexlist::Vector, ::Type{Val{4}})

    T = coordtype(mesh)
    U = universedimension(mesh)

    cellids, ncells = vertextocellmap(mesh)

    Cells = cells(mesh)
    Verts = vertices(mesh)

    # create the local shapes
    fns = Vector{BEAST.Shape{T}}[]
    pos = Vector{vertextype(mesh)}()

    sizehint!(fns, length(vertexlist))
    sizehint!(pos, length(vertexlist))
    for v in vertexlist

        numshapes = ncells[v]
        numshapes == 0 && continue

        shapes = Vector{BEAST.Shape{T}}(undef,numshapes)
        for s in 1: numshapes
            c = cellids[v,s]
            # cell = mesh.faces[c]
            cell = Cells[c]

            localid = something(findfirst(isequal(v), cell),0)
            @assert localid != 0

            shapes[s] = BEAST.Shape(c, localid, T(1.0))
        end

        push!(fns, shapes)
        push!(pos, Verts[v])
    end

    NF = 4
    BEAST.LagrangeBasis{1,0,NF}(mesh, fns, pos)
end

function unique_fn(fn)
    fn_new = typeof(fn[1])[]
    #display(fn)
    tups_new = []
    for shape in fn
        if (shape.cellid,shape.refid) in tups_new
            i = findfirst(x->x==(shape.cellid,shape.refid),tups_new)
            s = BEAST.Shape(shape.cellid,shape.refid,fn_new[i].coeff+shape.coeff)
            fn_new[i] = s
           # println("shorter")
        else
            push!(tups_new,(shape.cellid,shape.refid))
            push!(fn_new,shape)
        end
    end
    #@assert length(fn_new) != length(fn)
    return fn_new
end