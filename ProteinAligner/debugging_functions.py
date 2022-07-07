def print_dp_table(alignment_obj):
    column = [" ", " "]
    for s in [alignment_obj.graph.nodes[x].seq for x in alignment_obj.graph.sorted]:
        for c in s:
            column.append(c)

    print("\t".join(column))
    counter = 0
    for l in alignment_obj.dp_table:
        if counter > 0:
            to_print = [alignment_obj.seq[counter - 1]] + [str(int(x)) for x in l]
            # print(alignment_obj.seq[counter - 1], l)
            print("\t".join(to_print))
        else:
            to_print = [" "] + [str(int(x)) for x in l]
            print("\t".join(to_print))
        counter += 1
