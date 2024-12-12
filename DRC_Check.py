import LayoutScript
from LayoutScript import *

# from TWR_functions.TWR_fns import *

l = project.newLayout();  # open new instance of layout class
global dr
dr = l.drawing
# pointer to the main drawing
setup.gdsTextToPolygon = True
setup.gdsTextToPolygonDefaultWidth = 200000
setup.defaultTextWidth = 200000

SetUp = setup()  # work around as static string variables are not handled correctly

LogStat = True
# Layer definitions

#
#   Import Layout
#
dr.importFile("/Users/lipton/Dropbox/Programming/TWR_layout/TWR_23_r2.gds")

# Layer definitions
OTL = 201  # outline for drawing
OD = 1  # Defines active window
JTE = 116  # Junction termination extension IMPLANT (NP-JTE)
PGN = 117  # boron gain layer IMPLANT (NC)
NC = 117
ACN = 55  # AC phos layer IMPLANT (NQ)
NPL = 66  # n+ IMPLANT (NP)
NP = 66
MET = 44  # original design Metal
M1 = 43  # METAL 1
M2 = 47  # Metal 2
M3 = 49  # Metal 3
V1 = 46  # Via from M1 to M2
V2 = 48  # Hole to connect M3 to M2
ZG = 56  # Hole to connect ZA to M3
ZA = 57  # pad Metal (Al,Top Metal)
ZP = 58  # top passivation openings
CON = 25  # Contact (CA)
CA = 25  # contact
# ACC = 9
DJP = 18  # Deep junction p IMPLANT (PX)
DJN = 16  # Deep junction n IMPLANT (DW?)
PSB = 21  # p substrate contact IMPLANT (PD??)
PD = 21  # p contact implant
# PST = 117  # p stop (NC)
ND = 19  # n contact IMPLANT
PW = 78  # P-well IMPLANT
OF = 228  # p well
PWD = 228 # proper name
NWD = 227  # NOT NWD is PW
#
WLayer1 = 201
WLayer2 = 202
WLayer3 = 203

# List of layers to check
DRCLYRList = [JTE, PGN, ACN, ND, NPL, M1, M2, M3, DJP, DJN, PWD, NWD]
# DRCLYRList = [M1]

DRC_Cell_Names = ["6mm_100um_pitch_bumps"]

for cell in DRC_Cell_Names:
    dr.setCell(cell)
    if(LogStat): print("Cell " + cell)
    for layer in DRCLYRList:
        l.drcTool.ruleName= "Ang_45_deg " + str(layer)
        l.drcTool.angle45OnLayer(layer,True)
        if(LogStat): print("Layer " + str(layer) + " done")
l.drcTool.saveViolationList("/Users/lipton/DRC_list.layout")
nerr = l.drcTool.error
# print("errors = " + str(nerr))
for i in range(nerr):
    Type = l.drcTool.getViolationType(i)
    P1 = l.drcTool.getViolationPoint1(i)
    P2 = l.drcTool.getViolationPoint2(i)
    PX = P1.x()
    PY = P1.y()
    VValue = l.drcTool.getViolationValue(i)
    VName = l.drcTool.getViolationName(i)
    print("error " + str(i) + " Name " + str(VName)  + " Type " + str(Type) + " X = " + str(PX) + " Y " + str(PY) + " Value " +str(VValue))
