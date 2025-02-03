import LayoutScript
from LayoutScript import *
from pprint import pprint


def addMxFill(Assy, exclude, tlayer, fillCell, csize, osize):
    #   add Mx Fill
    #   Tlayer[1] = outline, Tlayer[2] - fill, mask TLayer[3] - areas to exclude
    #   csize - dimension to expand metal
    #   osuize - dimension to compress frame
    #	start_time = time.time()
    dr.setCell(Assy)
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
    #	print("Mx FIll time elapsed: {:.2f}s".format(time.time() - start_time))
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
    gdsversion = "V26_r0"
    ingdsfile = str(home_directory) + "/Dropbox/Programming/TWR_layout/TWR_" + gdsversion + ".gds"

    # Input Reticule
    dr.importFile(ingdsfile)
    #	l.drawing.saveFile("/Users/lipton/Test" +  ".gds")
    print(dr.existCellname("Fill_25pct"))
    fCell = dr.findCell("Fill_25pct")
    print(fCell)
    mlayer = [43, 47, 49]
    tlayer = [201, 202, 203]

    cmod = ["DJ_50Base", "DJ_100Base"]
    for i in range(len(cmod)):
        cnam = cmod[i]
        dr.setCell(cnam)
        Assy = dr.currentCell
        csize = 2000
        osize = 0
        addMxFill(Assy, mlayer, tlayer, fCell, csize, osize)


#	print(outgdsfile)

#	l.drawing.saveFile(outgdsfile)

#	print("Writing Version " + gdsversion +  " Python script completed")

if __name__ == "__main__":
    main()