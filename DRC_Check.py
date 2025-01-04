import LayoutScript
import sys
import time
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

LogStat = False
# Layer definitions
#
#   Import Layout
#
start_time = time.time()
dr.importFile("/Users/lipton/Dropbox/Programming/TWR_layout/TWR_24_r4.gds")

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

WLayers = [WLayer1, WLayer2, WLayer3]
# List of layers to check
DRCLYRList = [JTE, PGN, ACN, ND, NPL, M1, M2, M3, DJP, DJN, PWD, NWD, ZA, ZP]
DRCMList = [M1, M2, M3]
nlayer = 2
Filen = "/Users/lipton/Dropbox/Programming/TWR_layout/DRC/DenDRC_L" + str(DRCMList[nlayer]) + ".txt"
DRCfile = open(Filen, "w")
DRC_Cell_Names = ["DJ_100", "DJ_PN", "JTE", "DJ_50", "ACC-3mm", "ACP_3mm", "St100AC_40", "St100AC_60", "St100AC_80", "St100AC_20"]
# DRC_Cell_Names = ['Assy_AC_50', 'Assy_AC_100', 'AssyDJ_50', 'AssyDJ_100', 'AssyPX_50_NG', 'AssyPX_100_NG',
#                  'Assy_Str125_Arr', 'AssyDJ_NB_100', 'Assy_str_3mm_50', 'Assy_str_3mm_100', 'Assy_ACStr_Elec', 'Assy_RTPixel', 'Assy_RT_100']
# DRC_Cell_Names = ['Assy_AC_50']
#Border_Cell_Names = ["6mm_100um_pitch_bumps", "6mm_100um_pitch_bumps_DJ_ASIL", "6mm_50um_pitch_bumps", "6mm_with_pads",  "3mm_with_pads"]
Border_Cell_Names = []
# Checks = ["density", "inset" , "grid", "angle", "dimension"]
Checks = ["inset"]
#
#  Metal density parameters [M1,M2, M3] window sizes MDensity
#
# MDensity  =  [[2500000, 1000000, 70000, 14000], [2500000, 1000000, 70000, 14000], [2500000, 1000000, 20000]]
MDensity  =  [[2500000, 1000000, 70000], [2500000, 1000000, 70000], [2500000, 1000000, 20000]]
# MDensity  =  [[2500000], [2500000], [2500000, 1000000]]
LowDensity  =  [[30., 25., 20., 5.], [25., 25., 15., 5.], [25., 15., 5.]]
HighDensity  =  [[50., 55., 70., 85.], [50., 55., 70., 90.], [50., 70., 90.]]
for cell in DRC_Cell_Names:
    exists = dr.setCell(cell)
    if (not exists):
        print("Cell " + cell + " notfound")
        sys.exit(2)
    if "inset" in Checks:
        l.drcTool.ruleName = "Inside CA and M1"
        l.drcTool.minimumInside(50, CA, M1, -1, -1)
        l.drcTool.ruleName = "Inside V1 and M1"
        l.drcTool.minimumInside(50, V1, M1, -1, -1)
        l.drcTool.ruleName = "Inside V1 and M1"
        l.drcTool.minimumInside(50, V1, M1, -1, -1)
        l.drcTool.ruleName = "Inside V2 and M3"
        l.drcTool.minimumInside(50, V2, M3, -1, -1)
        l.drcTool.ruleName = "Inside ZG and ZA"
        l.drcTool.minimumInside(2000, ZG, ZA, -1, -1)


    if "density" in Checks:
        ncell = dr.currentCell
        # Define checking region
        l.drcTool.setRegionMode()
        minval = point(-2500000, 0)
        maxval = point(0, 2500000)
        l.drcTool.setCheckRegion(minval, maxval)
        #for nlayer in range(1):
        lyr = DRCMList[nlayer]
        print(" Layer " + str(lyr))
        # WL = WLayers[nlayer]
        l.booleanTool.boolOnLayer(lyr, 0, WLayer3, "A merge", 0, 0, 2)
        for edge in range(len(MDensity[nlayer])):
            MDen  = MDensity[nlayer][edge]
            LDen = LowDensity[nlayer][edge]
            HDen = HighDensity[nlayer][edge]
            l.drcTool.ruleName= "Cell " + cell + "Density layer " + str(lyr) + " area " + str(MDen)
            l.drcTool.densityOnLayer(WLayer3 , MDen , LDen, HDen)
            print(" Density check region " + str(MDen) + " layer " + str(lyr) + " time elapsed: {:.2f}s".format(time.time() - start_time) )
        # ncell.deleteLayer(WLayer3)

    for layer in DRCLYRList:
        if "angle" in Checks:
             l.drcTool.ruleName= "Ang_45_deg " + str(layer)
             l.drcTool.angle45OnLayer(layer,True)
        if "grid" in Checks:
             l.drcTool.ruleName = "On Grid " + str(layer)
             l.drcTool.onGrid(5, layer)
        if "dimension" in Checks:
            for layer in DRCMList:
                l.drcTool.ruleName = "Max Dimension " + str(layer)
                l.drcTool.maximumDimensionOnLayer(100000000, 1000, layer, False)

    print("Cell " + cell + " done")
# for cell in Border_Cell_Names:
#     exists = dr.setCell(cell)
#     if (not exists):
#         print("Cell " + cell + " notfound")
#         sys.exit(2)
#     if(LogStat): print("Cell " + cell)
#     for layer in DRCLYRList:
#         l.drcTool.ruleName= "Ang_45_deg " + str(layer)
#         l.drcTool.angle45OnLayer(layer,True)
#         l.drcTool.ruleName = "On Grid " + str(layer)
#         l.drcTool.onGrid(5, layer);
#         if(LogStat): print("Layer " + str(layer) + " done")
#     print("Cell " + cell + " done time elapsed: {:.2f}s".format(time.time() - start_time))

l.drcTool.saveViolationList("/Users/lipton/DRC_list_V25_r0.layout")
nerr = l.drcTool.error
print("errors = " + str(nerr))
# ys.exit(2)
for i in range(min(100,nerr)):
    Type = l.drcTool.getViolationType(i)
    P1 = l.drcTool.getViolationPoint1(i)
    P2 = l.drcTool.getViolationPoint2(i)
    PX = P1.x()
    PY = P1.y()
    VValue = l.drcTool.getViolationValue(i)
    VName = l.drcTool.getViolationName(i)
    DRCfile.write("error " + str(i) + " Name " + str(VName)  + " Type " + str(Type) + " X = " + str(PX) + " Y " + str(PY) + " Value " +str(VValue) +  "\n")
