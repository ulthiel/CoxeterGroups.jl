export CoxeterMatrix

function CoxeterMatrix(groupType::String)
    parsedGroupType = parseCoxeterGroupType(groupType)
    return buildCoxeterMatrixFromGroupType(parsedGroupType[1], parsedGroupType[2], parsedGroupType[3])
end

# Matrices constructed by diagrams in Kac's Infinite Dimensional Lie Algebras, pg. 43
function buildCoxeterMatrixFromGroupType(C::Char, N::Int, isAffine::Bool)
    if C == 'A'
        if !isAffine
            M = fill(2,N,N)
            for i = 1:N
                M[i,i] = 1
                if i != N
                    M[i,i+1] = 3
                    M[i+1,i] = 3
                end
            end
        else
            M = buildCoxeterMatrixFromGroupType('A',N+1,false)
            M[1,N+1] = 3
            M[N+1,1] = 3
        end
    elseif C =='B' || C == 'C'
        if !isAffine
            M = buildCoxeterMatrixFromGroupType('A',N,false)
            if N>1
                M[N-1,N] =4
                M[N,N-1]=4
            else
                error("Cannot create matrix of this size")
            end
        else
            M = buildCoxeterMatrixFromGroupType('B',N+1,false)
            if C == 'B'
                M[1,2]=2
                M[2,1]=2
                M[1,3]=3
                M[3,1]=3
            else
                M[1,2]=4
                M[2,1]=4
            end

        end
    elseif C == 'D'
        if !isAffine 
            if N<4 
                error("Cannot create matrix of this size")
            end
            M = buildCoxeterMatrixFromGroupType('A',N,false)
            M[N-1,N]=2
            M[N,N-1]=2
            M[N-2,N]=3
            M[N,N-2]=3
        else
            M = buildCoxeterMatrixFromGroupType('D',N+1,false)
            M[1,2]=2
            M[2,1]=2
            M[1,3]=3
            M[3,1]=3
        end
    elseif C == 'E'
        if !isAffine
            M = buildCoxeterMatrixFromGroupType('A',N,false)
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
                error("N must be 6,7 or 8")
            end
        else
            M = buildCoxeterMatrixFromGroupType('A',N+1,false)
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
                error("N must be 6,7 or 8")
            end
        end
    elseif C == 'F'
        if !isAffine
            if N < 4
                M = buildCoxeterMatrixFromGroupType('A',N,false)
                M[N-2,N-1] = 4
                M[N-1,N-2] = 4
            else
                error("Cannot create matrix of this size")
            end
        else
            M = buildCoxeterMatrixFromGroupType('F',N+1,false)
        end 
    elseif C == 'G'
        if !isAffine
            if N<2
                M = buildCoxeterMatrixFromGroupType('A',N,false)
                M[N-1,N]=6
                M[N,N-1]=6
            else
                error("Cannot create matrix of this size")
            end
        else
            M = buildCoxeterMatrixFromGroupType('G',N+1,false)
        end
    elseif C == 'H'
        if !isAffine
            if N > 1
                M = buildCoxeterMatrixFromGroupType('A',N,false)
                M[1,2]=5
                M[2,1]=5
            else
                error("N must be >= 2")
            end
        else
            error("H affine type not supported")
        end
    elseif C == 'I'
        if !isAffine
            if N >= 2
                M = buildCoxeterMatrixFromGroupType('A',2,false)
                M[1,2]=N
                M[2,1]=N
            else
                error("N must be >= 2")
            end
        else
            error("I affine type not supported")
        end
    end

    return M
end

""" 
    parseCoxeterGroupType


"""
function parseCoxeterGroupType(groupType::String)
    groupTypeRegex = r"[A-I]~?[0-9]+"
    m = match(groupTypeRegex,groupType)
    if m === nothing || length(m.match)!=length(groupType)
        error("Unrecognized group type $(groupType)")
    end
    C = groupType[1]
    if groupType[2] == '~'
        N =  parse(Int, groupType[3:end])
        isAffine = true
    else
        N =  parse(Int, groupType[2:end])
        isAffine = false
    end
    return [C,N,isAffine]
end




