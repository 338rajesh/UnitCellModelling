using UnitCellModelling

test_case_1 = false
test_case_2 = false
test_case_3 = false
test_case_4 = true


#
if test_case_1
    """

    """
    element_connectivity_1 = Dict(
    3 => [1 2 3;9 1 8; 1 3 2; 8 2 7; 4 8 5],
    2 => reshape([1; 2; 6; 7], 4, 1)
    )
    println("\n============= SET-1 ========")
    println(
        "Bandwidth before RCM renumbering: ",
        UnitCellModelling.bandwidth(element_connectivity_1, dof=1)
    )
    new_nt_order = UnitCellModelling.update_node_labels([element_connectivity_1, ])
    # RVE.FEPrep.renumber_node_tags_rcm!(element_connectivity_1)

    println("Before: element_connectivity_1: ", element_connectivity_1)
    println("new_nt_order: ", new_nt_order)
    UnitCellModelling.update_element_connectivity!(element_connectivity_1, new_nt_order)
    println("After: element_connectivity_1: ", element_connectivity_1)
    println(
        "Bandwidth after RCM renumbering: ",
        UnitCellModelling.bandwidth(element_connectivity_1, dof=1)
    )
end

if test_case_2
    """

    """
    element_connectivity_2 = Dict(
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
    println("\n============= SET-2 ========")
    println(
        "Bandwidth before RCM renumbering: ",
        UnitCellModelling.bandwidth(element_connectivity_2, dof=1)
    )
    new_nt_order = UnitCellModelling.update_node_labels([element_connectivity_2, ])
    UnitCellModelling.update_element_connectivity!(element_connectivity_2, new_nt_order)
    println(
        "Bandwidth after RCM renumbering: ",
        UnitCellModelling.bandwidth(element_connectivity_2, dof=1)
    )
end

if test_case_3
    """

    """
    element_connectivity_3 = Dict(
        3=>[
            1 2 3 4 5 6 8
            15 1 3 4 8 2 7
            1 3 11 8 2 13 6
            8 2 13 5 7 6 14
            4 8 2 10 5 7 12
        ],
        2 => [
            7 9
            5 12
            7 14
            12 9
        ],
    )
    println("\n============= SET-3 ========")
    println(
        "Bandwidth before RCM renumbering: ",
        UnitCellModelling.bandwidth(element_connectivity_3, dof=1)
    )
    new_nt_order = UnitCellModelling.update_node_labels([element_connectivity_3, ])
    UnitCellModelling.update_element_connectivity!(element_connectivity_3, new_nt_order)
    println(
        "Bandwidth after RCM renumbering: ",
        UnitCellModelling.bandwidth(element_connectivity_3, dof=1)
    )
end

if test_case_4

    raw"""
    source: https://people.sc.fsu.edu/~jburkardt/cpp_src/rcm/rcm.html

    expected bandwidth reduction from 17-> 7

    """
    element_connectivity_4 = Dict(
        3=>[
            1 2 3
            9 1 8
            1 3 2
            8 2 7
            4 8 5
        ],
        2 => reshape([4; 2; 6; 7], :, 1)
    )
    #
    println("\n============= SET-4 ========")
    println(
        "Bandwidth before RCM renumbering: ",
        UnitCellModelling.bandwidth(element_connectivity_4, dof=1)
    )
    new_nt_order = UnitCellModelling.get_rcm_node_labels(element_connectivity_4, )
    element_connectivity_4 = UnitCellModelling.update_element_connectivity(element_connectivity_4, new_nt_order)
    println(
        "Bandwidth after RCM renumbering: ",
        UnitCellModelling.bandwidth(element_connectivity_4, dof=1)
    )
    for (k, v) in element_connectivity_4
        println(k)
        display(v)
        println()
    end

end

