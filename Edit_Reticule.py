# -*- codin
import LayoutScript
from LayoutScript import *
from datetime import datetime

from pprint import pprint


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
now = datetime.now()
filetime = now.strftime("_%Y_%m_%d_%H_%M")
ingdsfile = str(home_directory) +  "/Dropbox/Programming/TWR_layout/TWR_" + gdsversion + ".gds"
outgdsfile = str(home_directory) +  "/Dropbox/Programming/TWR_layout/Edit_" + gdsversion + filetime + ".gds"

#
#
# Input Reticule
dr.importFile(ingdsfile)

update = [False, False, True]
# Update borders
#  - note that ZA fill is not executed
if update[0]: dr.updateFile("/Users/lipton/Dropbox/Programming/TWR_layout/SLAC_layouts/compile_border_v17.gds")
# Update Test structures
if update[1]: dr.updateFile("/Users/lipton/Dropbox/Programming/TWR_layout/SLAC_layouts/TWR_Test3x3_v4nool100.GDS")
# other update
if update[2]: dr.updateFile("/Users/lipton/Dropbox/Programming/TWR_layout/Bond_Pad_80_mod.gds")

cnam = "Reticule"
dr.setCell(cnam)
dr.stripUnneeded()  # remove unneeded calls

print("Writing " + outgdsfile)

l.drawing.saveFile(outgdsfile)
