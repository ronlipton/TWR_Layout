import LayoutScript
from LayoutScript import *
import datetime
import sys

# from TWR_Layout import TLayer

def addLD(chip):
    from pprint import pprint
    # import os
    from pathlib import Path
    import time
    home_directory = Path.home()
    l = project.newLayout()  # open new instance of layout class
    global dr
    dr = l.drawing
    # pointer to the main drawing
    #
    #
    start_time = time.time()
    #  input main gds file
    ingdsfile = str(home_directory) + "/Dropbox/Programming/TWR_layout/Chips_V31_r2/" + chip + ".gds"
    outgdsfile = str(home_directory) + "/Dropbox/Programming/TWR_layout/Chips_V31_r2/" + chip + "_LD.gds"

    dr.importFile(ingdsfile)
    exists = dr.setCell(chip)

    if (not exists):
        print("Cell " + chip + " notfound")
        sys.exit(2)
    #   exclude regions ND(19), PD(21), RS_PAC(55)
    #`  temp layer is 203
    TLayer3 = 203
    TLayer4 = 204
    PD = 21
    ND = 19
    RS_PAC = 55
    LD =119
    GASO = 226
    l.booleanTool.boolOnLayer(PD,ND, TLayer3, "A+B")
    l.booleanTool.boolOnLayer(ND, RS_PAC, TLayer4, "A+B")
    l.booleanTool.boolOnLayer(GASO, TLayer4, LD, "A-B", 0, 0, 2)
    l.drawing.activeLayer = LD
    l.drawing.selectActiveLayer()
    l.drawing.currentCell.sizeAdjustSelect(50, 0)

    dr.deleteLayer(TLayer3)
    dr.deleteLayer(TLayer4)
    l.filename = "/Users/lipton/Library/CloudStorage/Dropbox/Programming/TWR_layout/Chips_V31_r2/Assy_Str125_Arr_LD.gds"
    print("Writing " + outgdsfile +" time elapsed: {:.2f}s".format(time.time() - start_time))
    l.drawing.saveFile(outgdsfile)
    dr.deleteAllCell()
    return


def main():
    global l
    l = project.newLayout()  # open new instance of layout class
    global dr
    dr = l.drawing
    # pointer to the main drawing

    SetUp = setup()  # work around as static string variables are not handled correctly
    # Cells for NWD generation (InvertOF)
    NWDFill_name = ["Assy_AC_50", "Assy_AC_100", "AssyDJ_50", "AssyDJ_100", "AssyPX_50_NG", "AssyPX_100_NG",
                    "Assy_Str125_Arr",
                    "AssyDJ_NB_100", "Assy_RTPixel", "Assy_RT_100", "AssyDJ_100_NoPS"]

    subNWD_list = ["Str100_AC_Arr3mm", "Str100_DJNPS_Arr3mm", "Str100_DJ_Arr3mm", "Str100_NOGN_Arr3mm",
                   "Str50_AC_Arr3mm", "Str50_DJNPS_Arr3mm", "Str50_DJ_Arr3mm", "Str50_NOGN_Arr3mm",
                   "Str100AC_20_Arr3mm", "Str100AC_40_Arr3mm", "Str100AC_60_Arr3mm", "Str100AC_80_Arr3mm"]

    SLAC_chips = ["AC_T3x3Arry_Mid_Gain", "AC_T3x3Arry_Wide_Gain", "AC_TPad_Gain_MidGap", "AC_TPad_NoGain_MidGap",
                  "CEP_DJ_T3x3_Arry_gain", "CEP_DJ_TPAD_gain", "CEP_RT_T3x3_Arry_Gain_JTE_pStop", "CEP_RT_TPAD_gain",
                  "DJ_T3x3_Arry_gain", "DJ_T3x3_Arry_nogain", "DJ_TPAD_gain", "DJ_TPAD_nogain",
                  "RT_T3x3_Arry_Gain_JTE_pStop", "RT_T3x3_Arry_Gain_noJTE_pstop", "RT_TPAD_gain", "RT_TPAD_nogain",
                  "JTE_0ring", "JTE_1ring", "JTE_2ring", "JTE_3ring",
                  "DJ_0ring", "DJ_1ring", "JTE_5ring", "JTE_5ring_smaller_space"]

    All_chips = SLAC_chips + NWDFill_name + subNWD_list
    All_chips= ["Assy_AC_50"]
    for chip in All_chips:
        addLD(chip)


if __name__ == "__main__":
    main()