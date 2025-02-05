import LayoutScript
from LayoutScript import *
import datetime

def editGDS(gdslist):

    from pprint import pprint

    l = project.newLayout()  # open new instance of layout class
    global dr
    dr = l.drawing
    # pointer to the main drawing
    #
    #
    # Input Reticule
    dr.importFile(gdslist[0])

    update = [False, False, True]
    # Update borders
    #  - note that ZA fill is not executed
    if update[0]: dr.updateFile(gdslist[1])

    cnam = "Reticule"
    dr.setCell(cnam)
    dr.stripUnneeded()  # remove unneeded calls

    print("Writing " + gdslist[2])
    l.drawing.saveFile(gdslist[2])
    return

def main():
    global l
    l = project.newLayout()  # open new instance of layout class
    global dr
    dr = l.drawing
    # pointer to the main drawing
    setup.gdsTextToPolygon = True
    setup.gdsTextToPolygonDefaultWidth = 200000
    setup.defaultTextWidth = 200000

    SetUp = setup()  # work around as static string variables are not handled correctly

    # import os
    from pathlib import Path
    home_directory = Path.home()
    gdsversion = "V26_r0"
    now = datetime.datetime.now()
    filetime = now.strftime("_%Y_%m_%d_%H_%M")
    #  input main gds file
    ingdsfile = str(home_directory) + "/Dropbox/Programming/TWR_layout/TWR_" + gdsversion + ".gds"
    # file containing cells to update
    modgdsfile = str(home_directory) + "/Dropbox/Programming/TWR_layout/TWR_" +  "compile_border_v17.gds"
    # time stamped updated file
    outgdsfile = str(home_directory) + "/Dropbox/Programming/TWR_layout/Edit_" + gdsversion + filetime + ".gds"
    gds_list = [ingdsfile, modgdsfile, outgdsfile]
    editGDS(gds_list)
    return

if __name__ == "__main__":
    main()