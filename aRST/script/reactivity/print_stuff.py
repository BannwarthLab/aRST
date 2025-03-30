from os.path import join, exists


def printafo(status=0, printtype=None, printinfo=None):
    totw = 80
    border_width = 70
    col1 = 5
    col2 = 5
    col3 = 5
    col4 = 20
    side_padding = (totw - border_width) // 2
    side = " " * side_padding

    if status == 0:
        # head
        print(side + "=" * border_width)
        title = f"{printtype} Reactivity Summary"
        print(title.center(totw))
        print(side + "-" * border_width)
        print(side + f"{'struc':^{col1}} {'at':^{col2}} {'ele':^{col3}} {f'{printtype}_a (eV)':^{col4}}|"
                     f"{'struc':^{col1}} {'at':^{col2}} {'ele':^{col3}} {f'{printtype}_b (eV)':^{col4}}")
        print(side + "-" * border_width)
        return
    elif status == 1:
        # print value

        struc_a = printinfo["struc_a"]
        at_a = printinfo["at_a"]
        ele_a = printinfo["ele_a"]
        afo_a = printinfo["afo_a"]

        struc_b = printinfo["struc_b"]
        at_b = printinfo["at_b"]
        ele_b = printinfo["ele_b"]
        afo_b = printinfo["afo_b"]
        print(side + f"{str(struc_a).center(col1)}"
                     f"{str(at_a).center(col2)}"
                     f"{str(ele_a).center(col3)}"
                     f"{str(afo_a).center(col4)}|"
              + f"{str(struc_b).center(col1)}"
                f"{str(at_b).center(col2)}"
                f"{str(ele_b).center(col3)}"
                f"{str(afo_b).center(col4)}")
        return
    elif status == 2:
        print(side + "=" * border_width)
        print()
        return
    return

