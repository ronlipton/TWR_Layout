import LayoutScript
from LayoutScript import *
import datetime

def editGDS(gdsv):

    from pprint import pprint
    # import os
    from pathlib import Path
    home_directory = Path.home()
    l = project.newLayout()  # open new instance of layout class
    global dr
    dr = l.drawing
    # pointer to the main drawing
    #
    #
    now = datetime.datetime.now()
    filetime = now.strftime("_%Y_%m_%d_%H_%M")
    #  input main gds file
    ingdsfile = str(home_directory) + "/Dropbox/Programming/TWR_layout/TWR_111V26_r2.gds"
    # file containing cells to update
    modgdsfile1 = str(home_directory) + "/Dropbox/Programming/TWR_layout/V1_Via_4x4.gds"
    modgdsfile2 = str(home_directory) + "/Dropbox/Programming/TWR_layout/V2_Via_4x4.gds"
    modgdsfile3 = str(home_directory) + "/Dropbox/Programming/TWR_layout/Contact_4x4.gds"
    modlist = [modgdsfile1, modgdsfile2, modgdsfile3]

    # time stamped updated file
    # Input Reticule
    dr.importFile(ingdsfile)
    outgdsfile = str(home_directory) + "/Dropbox/Programming/TWR_layout/Edit_" + gdsv + filetime + ".gds"

    for i in range(len(modlist)):
        dr.updateFile(modlist[i])

    cnam = "Reticule"
    dr.setCell(cnam)
    dr.stripUnneeded()  # remove unneeded calls
    print("Writing " + outgdsfile)
    l.drawing.saveFile(outgdsfile)

    return

def main():
    global l
    l = project.newLayout()  # open new instance of layout class
    global dr
    dr = l.drawing
    # pointer to the main drawing

    SetUp = setup()  # work around as static string variables are not handled correctly

    # import os
    editGDS("V26_r3")


if __name__ == "__main__":
    main()