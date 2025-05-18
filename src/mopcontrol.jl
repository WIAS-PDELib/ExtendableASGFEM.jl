"""
$(TYPEDSIGNATURES)

generates multi-indices for `M` dimensions up to degree `deg`
"""
function generate_multiindices(M, deg)
    multi_indices = Array{Array{Int, 1}, 1}([[j] for j in 0:deg])
    L::Int = length(multi_indices)
    for m in 1:(M - 1)
        for i in 1:L
            for j in 0:deg
                new_mi = deepcopy(multi_indices[i])
                push!(new_mi, j)
                push!(multi_indices, new_mi)
            end
        end
        multi_indices = multi_indices[(L + 1):end]
        L *= (deg + 1)
    end

    return multi_indices
end

"""
$(TYPEDSIGNATURES)

ensures that all multi-indices have the same length
"""
function prepare_multi_indices!(multi_indices; minimal_length = 0)
    new_length = max(minimal_length, maximum(length.(multi_indices)))
    for j in 1:length(multi_indices)
        while length(multi_indices[j]) < new_length
            push!(multi_indices[j], 0)
        end
    end
    return
end

"""
$(TYPEDSIGNATURES)

adds new stochastic modes:
- for each existing mode all neighbouring modes are added (= all possible copies of that mode where one dimension is increased by 1)
- `p_extension` additionally ensures that the polynomial degree of the first dimension is increased by this amount
- `tail_extension` activates `tail_expansion[1]` many new stochastic modes with order 1,
   for each existing mode also the next `tail_expansion[2]` many higher stochastic modes are activated or increased
"""
function add_boundary_modes(multi_indices; p_extension = 1, tail_extension = [10, 2])
    last_nonzero = 0 # = support
    maxlength = 1
    maxdegree1 = 0
    nmodes = length(multi_indices)
    for j in 2:length(multi_indices)
        maxlength = max(maxlength, length(multi_indices[j]))
        maxdegree1 = max(maxdegree1, multi_indices[j][1])
        for k in length(multi_indices[j]):-1:(last_nonzero + 1)
            if multi_indices[j][k] != 0
                last_nonzero = k
            end
        end
    end
    sleep(1)
    prepare_multi_indices!(multi_indices; minimal_length = last_nonzero + tail_extension[1])
    extended_multiindices = deepcopy(multi_indices)
    maxlength2 = maximum(length.(extended_multiindices))

    ## always push the next tail_extension many [0 ... 0 1] mode first
    for k in 1:maxlength2
        new_multi_index = deepcopy(multi_indices[1])
        new_multi_index[k] = 1
        if new_multi_index in extended_multiindices
        else
            push!(extended_multiindices, new_multi_index)
        end
    end

    ## always push the next p_extension degrees in first mode
    for k in (maxdegree1 + 1):(maxdegree1 + p_extension)
        new_multi_index = deepcopy(multi_indices[1])
        new_multi_index[1] = k
        if new_multi_index in extended_multiindices
        else
            push!(extended_multiindices, new_multi_index)
        end
    end

    last_nonzero_pos::Int = 1
    for j in 1:nmodes
        ## find last nonzero mode
        last_nonzero_pos = 1
        for k in length(multi_indices[j]):-1:1
            if multi_indices[j][k] != 0
                last_nonzero_pos = k
                break
            end
        end

        for k in 1:(last_nonzero_pos + tail_extension[2])
            if k > last_nonzero + tail_extension[2]
                break
            end
            new_multi_index = deepcopy(multi_indices[j])
            new_multi_index[k] += 1
            if new_multi_index in extended_multiindices
            else
                push!(extended_multiindices, deepcopy(new_multi_index))
            end
            # if add_second_layer
            #     for k2 = 1 : last_nonzero_pos + 1
            #         new_multi_index[k2] += 1
            #         if new_multi_index in extended_multiindices
            #         else
            #             push!(extended_multiindices, deepcopy(new_multi_index))
            #         end
            #     end
            # end
        end
    end
    return extended_multiindices
end

function get_neighbourhood_relation_matrix(multi_indices, multi_indices_extended)
    NRM = ExtendableSparseMatrix{Int, Int}(length(multi_indices), length(multi_indices_extended))
    distance::Int = 0
    for j in 1:length(multi_indices)
        for k in 1:length(multi_indices_extended)
            distance = sum(abs.(multi_indices[j] .- multi_indices_extended[k]))
            if distance < 3
                NRM[j, k] = distance
            end
        end
    end
    flush!(NRM)
    return NRM
end

"""
$(TYPEDSIGNATURES)

classifies all multi indices into these categories:
- `active_int` = mode and all its direct neighbours (+1 in all active and the next dimensions) are in active_modes
- `active_bnd` = mode but not all its direct neighbours are in active_modes
- `inactive_bnd` = inactive mode that has an active neighbour
- `inactive_bnd2` = inactive mode that has an neighbour in inactive_bnd 
- `inactive_else` = all other
"""
function classify_modes(multi_indices, active_modes = multi_indices)
    active_bnd = Int[] # have an inactive neighbour
    active_int = Int[] # all neighbour modes are active
    inactive_bnd = Int[] # first layer of inactive modes
    inactive_bnd2 = Int[] # second layer of inactive modes
    inactive_else = Int[] # other inactive modes
    last_nonzero_pos::Int = 1
    bnd_level::Int = 0
    new_multi_index = deepcopy(multi_indices[1])
    active::Bool = true
    for j in 1:length(multi_indices)
        last_nonzero_pos = 1
        active = true
        for k in length(multi_indices[j]):-1:1
            if multi_indices[j][k] != 0
                last_nonzero_pos = k
                break
            end
        end
        if multi_indices[j] in active_modes # active mode
            if last_nonzero_pos == length(multi_indices[j])
                push!(active_bnd, j)
            else
                for k in 1:(last_nonzero_pos + 1)
                    new_multi_index .= multi_indices[j]
                    new_multi_index[k] += 1
                    if new_multi_index in active_modes
                    else
                        active = false
                        break
                    end
                end
                if active == true
                    push!(active_int, j)
                else
                    push!(active_bnd, j)
                end
            end
        else # inactive mode
            bnd_level = 0
            for k in 1:min(last_nonzero_pos + 1, length(new_multi_index))
                new_multi_index .= multi_indices[j]
                if new_multi_index[k] > 0
                    new_multi_index[k] -= 1
                    if new_multi_index in active_modes
                        bnd_level = 1
                        push!(inactive_bnd, j)
                        break
                    else
                    end
                end
            end
            if bnd_level == 0
                for k in 1:min(last_nonzero_pos + 1, length(new_multi_index)), k2 in 1:min(last_nonzero_pos + 1, length(new_multi_index))
                    new_multi_index .= multi_indices[j]
                    if new_multi_index[k] > 0 && new_multi_index[k2] > 0
                        new_multi_index[k] -= 1
                        new_multi_index[k2] -= 1
                        if new_multi_index in active_modes
                            bnd_level = 2
                            push!(inactive_bnd2, j)
                            break
                        else
                        end
                    end
                end
                if bnd_level == 0
                    push!(inactive_else, j)
                end
            end
        end
    end
    return inactive_else, inactive_bnd, inactive_bnd2, active_bnd, active_int
end
