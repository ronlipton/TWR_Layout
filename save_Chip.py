import LayoutScript
from LayoutScript import *
import datetime
import sys

def saveGDS(gdsfile):
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
    ingdsfile = str(home_directory) + "/Dropbox/Programming/TWR_layout/" + gdsfile

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
    # All_chips = ["AC_T3x3Arry_Mid_Gain", "AC_T3x3Arry_Wide_Gain"]

    # time stamped updated file
    # Input Reticule
    dr.importFile(ingdsfile)

    for chip in All_chips:
        dr.importFile(ingdsfile)
        outgdsfile = str(home_directory) + "/Dropbox/Programming/TWR_layout/Chips_V31_r2/" + chip + ".gds"
        exists = dr.setCell(chip)

        if (not exists):
            print("Cell " + chip + " notfound")
            sys.exit(2)
        dr.stripUnneeded()  # remove unneeded calls
        print("Writing " + outgdsfile)
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

    # import os
    saveGDS("TWR_11165V31_r2.gds")


if __name__ == "__main__":
    main()