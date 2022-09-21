export coxeter_matrix_from_group_type

"""
    coxeter_matrix_from_group_type(groupType::String)

Create a coxeter matrix for a given group type. To denote the direct product of group types,
enter the group types seperated by the 'x' character. Affine groups are denoted by including a '~' character
after a letter. 

#Examples
```julia-repl
julia> coxeter_matrix_from_group_type("A4")
4×4 Matrix{Int64}:
 1  3  2  2
 3  1  3  2
 2  3  1  3
 2  2  3  1
```
```julia-repl
julia> coxeter_matrix_from_group_type("D~4")
5×5 Matrix{Int64}:
 1  2  3  2  2
 2  1  3  2  2
 3  3  1  3  3
 2  2  3  1  2
 2  2  3  2  1
```
```julia-repl
julia> coxeter_matrix_from_group_type("A3 x B4")
SubString{String}["A3", "B4"]
7×7 Matrix{Int64}:
 1  3  2  2  2  2  2
 3  1  3  2  2  2  2
 2  3  1  2  2  2  2
 2  2  2  1  3  2  2
 2  2  2  3  1  3  2
 2  2  2  2  3  1  4
 2  2  2  2  2  4  1
```
```julia-repl
julia> coxeter_matrix_from_group_type("I2(4)")
2×2 Matrix{Int64}:
 1  4
 4  1
```
"""
function coxeter_matrix_from_group_type(groupType::String)
    parsedGroupTypes = validate_and_parse_coxeter_group_type_string(groupType)

    matrix_array = Array{Matrix{Int64},1}(undef,length(parsedGroupTypes))

    for i in 1:length(parsedGroupTypes)
        matrix_array[i] = build_coxeter_matrix_from_group_type(parsedGroupTypes[i]...)
    end

    M = matrix_array[1]
    for i in 2:length(matrix_array)
        M = hcat(vcat(M,fill(2,size(matrix_array[i])[1],size(M)[1])),vcat(fill(2,size(M)[1],size(matrix_array[i])[1]),matrix_array[i]))
    end
    return M
end

"""
    build_coxeter_matrix_from_group_type(C::Char, N::Int, isAffine::Bool)

Take a given group type letter (A through H, or I2), group type integer, and boolean flag if the group is affine. Return the 
corresponding Coxeter Matrix, if it is possible to construct.

Matrices are constructed based on the diagrams in Kac's Infinite Dimensional Lie Algebras, pg. 43
"""
function build_coxeter_matrix_from_group_type(C::Char, N::Int, isAffine::Bool)
    if C == 'A'
        if isAffine
            if N==0
                error("A~ does not have size 0")
            elseif N==1
                M = build_coxeter_matrix_from_group_type('I',0,false)
            else
                M = build_coxeter_matrix_from_group_type('A',N+1,false)
                M[1,N+1] = 3
                M[N+1,1] = 3
            end
        else
            M = fill(2,N,N)
            M[diagind(M)] .= 1
            M[diagind(M, 1)] .= 3
            M[diagind(M, -1)] .= 3
        end
    elseif C =='B' || C == 'C'
        if isAffine
            M = build_coxeter_matrix_from_group_type('B',N+1,false)
            if C == 'B' && N>2
                M[1,2]=2
                M[2,1]=2
                M[1,3]=3
                M[3,1]=3
            else
                M[1,2]=4
                M[2,1]=4
            end
        else
            M = build_coxeter_matrix_from_group_type('A',N,false)
            if N>1
                M[N-1,N] =4
                M[N,N-1]=4
            else
                error("B and C types must be of size >=2")
            end
        end
    elseif C == 'D'
        if isAffine
            if N == 3
                M = build_coxeter_matrix_from_group_type('A',N,true)
            elseif N>3
                M = build_coxeter_matrix_from_group_type('D',N+1,false)
                M[1,2]=2
                M[2,1]=2
                M[1,3]=3
                M[3,1]=3
            else
                error("D~ types must be of size >=3")
            end
        else
            if N == 2
                M = build_coxeter_matrix_from_group_type('I',2,false)
            elseif N>2
                M = build_coxeter_matrix_from_group_type('A',N,false)
                M[N-1,N]=2
                M[N,N-1]=2
                M[N-2,N]=3
                M[N,N-2]=3
            else
                error("D types must be of size >=2")
            end
        end
    elseif C == 'E'
        if isAffine
            M = build_coxeter_matrix_from_group_type('A',N+1,false)
            if N == 6
                M[5,6]=2
                M[6,5]=2
                M[3,6]=3
                M[6,3]=3
            elseif N == 7
                M[7,8]=2
                M[8,7]=2
                M[4,8]=3
                M[8,4]=3
            elseif N == 8
                M[8,9]=2
                M[9,8]=2
                M[6,9]=3
                M[9,6]=3
            else
                error("E~ types must be of size 6,7 or 8")
            end
        else
            M = build_coxeter_matrix_from_group_type('A',N,false)
            if N == 6
                M[5,6]=2
                M[6,5]=2
                M[3,6]=3
                M[6,3]=3
            elseif N ==7
                M[6,7]=2
                M[7,6]=2
                M[3,7]=3
                M[7,3]=3
            elseif N == 8
                M[7,8]=2
                M[8,7]=2
                M[5,8]=3
                M[8,5]=3
            else
                error("E types must be of size 6,7 or 8")
            end
        end
    elseif C == 'F'
        if N == 4
            if isAffine
                M = build_coxeter_matrix_from_group_type('A',5,false)
            else
                M = build_coxeter_matrix_from_group_type('A',4,false)
            end
            M[N-2,N-1] = 4
            M[N-1,N-2] = 4
        else
            error("F and F~ types must be of size 4")
        end
    elseif C == 'G'
        if N==2
            if isAffine
                M = build_coxeter_matrix_from_group_type('A',3,false)
            else
                M = build_coxeter_matrix_from_group_type('A',2,false)
            end
            M[N-1,N]=6
            M[N,N-1]=6
        else
            error("G and G~ types must be of size 2")
        end
    elseif C == 'H'
        if !isAffine
            if N > 1
                M = build_coxeter_matrix_from_group_type('A',N,false)
                M[1,2]=5
                M[2,1]=5
            else
                error("H types must be of size >= 2")
            end
        else
            error("H~ type not supported")
        end
    elseif C == 'I'
        if isAffine
            error("I~ type not supported")
        else
            if N >= 2 || N==0
                M = build_coxeter_matrix_from_group_type('A',2,false)
                M[1,2]=N
                M[2,1]=N
            else
                error("I types must be of size 0 or >= 2")
            end
            
        end
    end

    return M
end

""" 
    parseCoxeterGroupType(groupType::String)

Validate and parse a group type string (For example, "C5" or "A~7"). Return the parsed elements of the string as a list.
"""
function validate_and_parse_coxeter_group_type_string(groupType::String)
    groupTypeList = split(filter(x -> !isspace(x),groupType), "x")
    parsedGroupTypes = []
    for element in groupTypeList
        groupTypeRegex = r"([A-H]~?[0-9]+)"
        m = match(groupTypeRegex,element)
        if m !== nothing && length(m.match)==length(element)
            C = element[1]
            if element[2] == '~'
                N =  parse(Int, element[3:end])
                isAffine = true
            else
                N =  parse(Int, element[2:end])
                isAffine = false
            end
            push!(parsedGroupTypes,[C,N,isAffine])
            continue   
        end
        groupTypeRegex = r"I2\([0-9]+\)"
        m = match(groupTypeRegex,element)
        if m !== nothing && length(m.match)==length(element)
            C='I'
            N = parse(Int,element[4:end-1])
            isAffine=false
            push!(parsedGroupTypes,[C,N,isAffine])
            continue
        end
        error("Unrecognized group type $(element)")

    end
    return parsedGroupTypes
end




