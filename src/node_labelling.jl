function rcmpermute(A::SparseMatrixCSC)
    adj = coladjacency(A)
    sdeg = sortperm(coldegrees(A))
    p = symrcm(adj, sdeg)
    return A[p, p]
end


function symrcm(A::SparseMatrixCSC, args...)
    adj = coladjacency(A)
    sdeg = sortperm(coldegrees(A))
    return symrcm(adj, sdeg, args...)
end


function symrcm(
    adjac::Array{Array{Int64,1}},
    sdegr::Array{Int64,1}, rev=true, invert=false, warnunconnected=false
)
    # Based on the description of the RCM algorithm by Ciprian Zavoianu.
    # See: http://ciprian-zavoianu.blogspot.com/2009/01/project-bandwidth-reduction.html
    nverts = length(adjac)
    checkbounds(sdegr, nverts)
    sd = copy(sdegr)
    R = Array{Int}(undef, 0)
    Q = Array{Int}(undef, 0)
    P = 0
    C = 0
    unconnected = 0
    inserted = falses(nverts)
    while !isempty(sd)
        found = false
        while !isempty(sd)
            P = popfirst!(sd)
            if !inserted[P]
                push!(R, P)
                inserted[P] = true
                found = true
                break
            end
        end
        if found
            empty!(Q)
            append!(Q, adjac[P])
            while !isempty(Q)
                while !isempty(Q)
                    C = popfirst!(Q)
                    if !inserted[C]
                        push!(R, C)
                        inserted[C] = true
                        append!(Q, adjac[C])
                        break
                    end
                end
            end
        end
        if length(R) == nverts
            break
        else
            unconnected += 1
        end
    end
    if unconnected > 0 && warnunconnected
        warn("There are $(unconnected) unconnected regions.")
    end
    if rev
        reverse!(R)
    end
    if invert
        L = zeros(Int64, nverts)
        for i = 1:nverts
            L[R[i]] = i
        end
        return L
    else
        return R
    end
end

function coldegrees(A::SparseMatrixCSC)
    # assumes A is structurally symmetric
    cptr = A.colptr
    ncols = length(cptr) - 1
    degr = zeros(Int, ncols)
    @inbounds for j = 1:ncols
        degr[j] = coldeg(A, j)
    end
    return degr
end

function coldeg(A::SparseMatrixCSC, j::Int)
    # assumes A is structurally symmetric
    cptr = A.colptr
    return cptr[j+1] - cptr[j]
end

function coladjacency(A::SparseMatrixCSC)
    # assumes A is structurally symmetric
    cptr = A.colptr
    rval = A.rowval
    ncols = length(cptr) - 1
    adjac = Vector{Vector{Int}}(undef, ncols)
    sbyf = let A = A
        j -> coldeg(A, j)
    end
    for j = 1:ncols
        strt = cptr[j]
        jdeg = coldeg(A, j)
        jadj = Vector{Int}(undef, jdeg)
        for i = 1:jdeg
            jadj[i] = rval[strt+i-1]
        end
        sort!(jadj, by=sbyf)
        adjac[j] = jadj
    end
    return adjac
end


"""
Returns
-------
a Dictionary of node neighbours with the following key-value
pairs

**key**: Integer, node tag

**value**: Vector{Integer}, tags of the neighbouring nodes 

Args:
----
**`element_connectivity`**: Dictionary with key-value pairs

   *key*: Element type ID (as of now gmsh id is used) \n
   *value*: Matrix{Integer} whose each column contains element information
with tag at first position followed by nodal connectivity

"""
function get_node_neighbours(
    element_connectivity::Dict{Int64,Matrix{Int64}}
)::Dict{Int64,Vector{Int64}}
    neihbours::Dict{Int64,Vector{Int64}} = Dict()
    for (gmsh_ele_type, a_eltype_conn) in element_connectivity
        neighbour_local_ids = adjacent_local_node_indices(gmsh_ele_type)
        for a_ele_conn in eachcol(a_eltype_conn)
            current_ele_node_tags = a_ele_conn[2:end]
            for (k, node_tag) in enumerate(current_ele_node_tags)  # as the first term is ele_tag
                node_neighbours = current_ele_node_tags[neighbour_local_ids[:, k]]
                # here these neighbors of node are on the current
                # element
                #
                if haskey(neihbours, node_tag)
                    for a_node in node_neighbours
                        if !(a_node in neihbours[node_tag])
                            push!(neihbours[node_tag], a_node)
                        end
                    end
                else
                    neihbours[node_tag] = node_neighbours
                end
            end
        end
    end
    return neihbours
end


function merge_ele_conn(ele_conn_tables::Vector{Dict{Int, Matrix{Int}}})
    all_ele_conn = Dict{Int,Matrix{Int}}()
    for a_ec_table in ele_conn_tables
        for (aele_typ, aele_conn) in a_ec_table
            if aele_typ in keys(all_ele_conn)
                all_ele_conn[aele_typ] = hcat(all_ele_conn[aele_typ], aele_conn)
            else
                all_ele_conn[aele_typ] = aele_conn
            end
        end
    end
    return all_ele_conn
end


function get_rcm_node_labels(
    all_ele_conn::Dict{Int,Matrix{Int}},
)
    # get the node neighbours
    node_neighbours = get_node_neighbours(all_ele_conn)
    #
    num_nodes = length(node_neighbours)
    II = Int[]
    JJ = Int[]
    for (ant, an_nhbs) in node_neighbours
        push!(II, ant)
        push!(JJ, ant)
        for ann in an_nhbs
            push!(II, ant)
            push!(JJ, ann)
        end
    end
    old_perm_mat = sparse(II, JJ, 1.0, num_nodes, num_nodes)
    return symrcm(old_perm_mat)
end


"""

Returns bandwidth from matrix

"""
function bandwidth(
    A::AbstractMatrix;
    tol::Float64=0.0
)::Int64
    non_zero_indices = findall(abs.(A) .> tol)
    return maximum([abs(k[1] - k[2]) for k in non_zero_indices]) + 1
end



"""
Returns bandwidth of the global connectivity matrix from element connectivity

Args:
-----
*element_connectivity*
*dof*


```
element_connectivity = Dict(
    16 => [
        1 2 3 4 5 6
        21 23 25 27 29 31
        22 24 26 28 30 32
        23 25 27 29 31 33
        15 16 17 18 19 20
        3 5 7 9 11 13
        2 4 6 8 10 12
        1 3 5 7 9 11
        14 15 16 17 18 19
    ]
)

>julia bandwidth(element_connectivity)
91
```

"""
function bandwidth(
    element_connectivity::Dict{T,Matrix{T}};
    dof::Int64=1,
    half_bw=true
) where {T}
    lbw::T = 0
    ubw::T = 0
    for (_, a_eltype_conn) in element_connectivity
        for (_, a_ele_node_tags...) in eachcol(a_eltype_conn)
            for nt_i in a_ele_node_tags
                for nt_j in a_ele_node_tags
                    ubw = max(ubw, nt_j - nt_i)
                    lbw = max(lbw, nt_i - nt_j)
                end
            end
        end
    end
    lbw = ((lbw + 1) * dof) - 1  # dof: degrees of freedom on a node
    ubw = ((ubw + 1) * dof) - 1  # dof: degrees of freedom on a node
    if half_bw
        return (lbw, ubw)
    else
        return lbw + ubw + 1
    end
end


function update_element_connectivity(
    old_ele_connectivity::Dict{Int,Matrix{Int}},
    new_order::Vector{Int}
)::Dict{Int,Matrix{Int}}
    ntag_map = Dict((i => k) for (k, i) in enumerate(new_order))

    new_ele_connectivity::Dict{Int,Matrix{Int}} = Dict{Int,Matrix{Int}}()

    for (ael_type, aelt_conn_table) in old_ele_connectivity
        # updated element_ntags table
        m, n = size(aelt_conn_table)
        for i in 1:n
            for j in 2:m
                aelt_conn_table[j, i] = ntag_map[aelt_conn_table[j, i]]
            end
        end
        if ael_type in keys(new_ele_connectivity)
            new_ele_connectivity[ael_type] = hcat(new_ele_connectivity[ael_type], aelt_conn_table)
        else
            new_ele_connectivity[ael_type] = aelt_conn_table
        end
    end
    return new_ele_connectivity
end

