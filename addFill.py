import LayoutScript
from LayoutScript import *
from pprint import pprint
import time
import sys
# from TWR_Layout import ZA_FillCell


def addMZFill(Assy, exclude, tlayer, fillCell, csize, osize):
    #   add Mx Fill
    #   Tlayer[1] = outline, Tlayer[2] - fill, mask TLayer[3] - areas to exclude
    #   csize - dimension to expand metal
    #   osuize - dimension to compress frame
    start_time = time.time()
    dr.setCell(Assy)
    # merge excluded layers into tlayer[2]
    for ilyr in range(len(exclude)):
        l.booleanTool.boolOnLayer(tlayer[2], exclude[ilyr], tlayer[2], "A+B", 0, 0, 2)
    dr.deselectAll()
    Assy.selectLayer(tlayer[2])
    dr.currentCell.sizeAdjustSelect(csize, 0)
    l.booleanTool.boolOnLayer(tlayer[0], tlayer[2], tlayer[1], "A-B", 0, 0, 2)
    dr.deselectAll()
    Assy.selectLayer(tlayer[1])
    dr.currentCell.sizeAdjustSelect(osize, 0)
    dr.fillSelectedShapes(fillCell, 0)
    #   debug
    debug = True
    if debug:
        test = dr.currentCell.cellName
        l.drawing.saveFile("/Users/lipton/" + test + ".gds")
    else:
        Assy.deleteLayer(tlayer[1])
        Assy.deleteLayer(tlayer[2])
    print("MZ FIll time elapsed: {:.2f}s".format(time.time() - start_time))
    return True


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
    gdsversion = "V26_r1"
    ingdsfile = str(home_directory) + "/Dropbox/Programming/TWR_layout/TWR_000" + gdsversion + ".gds"
    outgdsfile = str(home_directory) + "/Dropbox/Programming/TWR_layout/Fill_" + gdsversion + ".gds"
    # Input Reticule
    dr.importFile(ingdsfile)
    #	l.drawing.saveFile("/Users/lipton/Test" +  ".gds")
    # print(dr.existCellname("Za_Fill"))
    CCell = dr.findCell("Fill_25pct")
    ZAfill = dr.findCell("Za_Fill")
    fCell = [CCell, ZAfill]
    # print(fCell)
    # layers for both Mx and ZA
    mlayer = [[43, 47, 49],[57, 57, 57]]
    tlayer = [[201, 202, 203],[201, 202, 203]]
    # dimension to grow target layer (close up holes)
    csize = [2000, 2000]
    # dimension to modify the resulting inverted layer
    osize = [0, 0]

    # Target single cell, ZA only
    Str_100_Cell = dr.findCell("Assy_str_3mm_100")
    ZA_exclude = [228, 204, 57]
    tlyr = [201, 202, 203]
    addMZFill(Str_100_Cell, ZA_exclude, tlyr, ZAfill, csize[0], osize[0])

    # cmod = [["DJ_50", "DJ_100"],["Str100_AC", "Str100_DJ"]]
    # for j in range(len(mlayer)):
    #     for i in range(len(cmod[j])):
    #         cnam = cmod[j][i]
    #         dr.setCell(cnam)
    #         Assy = dr.currentCell
    #         addMZFill(Assy, mlayer[i], tlayer[i], fCell[i], csize[i], osize[i])

	print(outgdsfile)
    l.drawing.saveFile(outgdsfile)
    print("Writing Version " + gdsversion +  " Python script completed")

if __name__ == "__main__":
    main()