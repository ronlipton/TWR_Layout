# -*- coding: utf-8 -*-
from LayoutScript import *
from pprint import pprint
def BoxDraw(c,xgr,ygr,radius,whigh,wlow, layer):
#
#       Draw Top, bottom, left, right boxes to enclosure a pixel or strip area
#   ie y of top =- ygr(0) + radius + whigh ...
#   c - target cell
#   xgr - x values (4) of guide box
#   ygr - y values of guide box
#   radius - radius of offset reference
#   whigh - high offset of box from radius
#   wlow - low offset from radius
#   layer - drawing layer
#
    ind = [0,1,2,3,0]
    for  j in range(4):
        pa = pointArray()
        il = ind[j]
        ih = ind[j+1]
        pa = pointArray()
        x = xgr[0]
        y = ygr[0]
        pa.attach(x, y + radius - wlow)
        pa.attach(x, y + radius + whigh)
        xn = xgr[1]
        yn = ygr[1]
        pa.attach(xn, yn + radius + whigh)
        pa.attach(xn, yn + radius - wlow)
        c.addBox(pointArray(pa), layer)

        pa = pointArray()
        x = xgr[1]
        y = ygr[1]
        pa.attach(x - radius + wlow, y)
        pa.attach(x - radius - whigh, y)
        xn = xgr[2]
        yn = ygr[2]
        pa.attach(xn - radius - whigh, yn)
        pa.attach(xn - radius + wlow, yn)
        c.addBox(pointArray(pa), layer)

        pa = pointArray()
        x = xgr[2]
        y = ygr[2]
        pa.attach(x, y - radius + wlow)
        pa.attach(x, y - radius - whigh)
        xn = xgr[3]
        yn = ygr[3]
        pa.attach(xn, yn - radius - whigh)
        pa.attach(xn, yn - radius + wlow)
        c.addBox(pointArray(pa), layer)

        pa = pointArray()
        x = xgr[3]
        y = ygr[3]
        pa.attach(x + radius - wlow, y)
        pa.attach(x + radius + whigh, y)
        xn = xgr[0]
        yn = ygr[0]
        pa.attach(xn + radius + whigh, yn)
        pa.attach(xn + radius - wlow, yn)
        c.addBox(pa, layer)
#        pa.clear()
    return 1


def sign(num):
    return -1 if num < 0 else 1


def erdraw(c, xgr, ygr, wlow, whigh, layer, radius):
    #
    #  Draw guard rings at radius
    #
    #       c - cell name
    #       xgr - array o x coord of GR arc centers
    #       ygr - array o y coord of GR arc centers
    #		radius - radius of reference arc
    #		wlow - low offset from reference
    #		whigh - high offset from reference
    #
    #   Uses boxdraw to draw enclosing horrizontal and vertical sides
    angle = 0
    #   edge polygons
    for ind in range(len(xgr)):
        #    horiz = bool(True)
        x = xgr[ind]
        y = ygr[ind]
        sx = sign(x)
        sy = sign(y)
        #    print(str(angle), str(layer), str(radius))
        #    c.addPolygonArc(point(x, y), int(radius-wlow), int(radius+whigh), angle, angle+90, layer)
        pa = pointArray()
        pa.attach(x + sx * (radius + whigh), y)
        pa.attach(x + sx * (radius - wlow), y)
        pa.attach(x, y + sy * (radius - wlow))
        pa.attach(x, y + sy * (radius + whigh))
        c.addPolygon(pointArray(pa), layer)
    #
    #      connect arcs
    #
    rval = BoxDraw(c, xgr, ygr, radius, whigh, wlow, layer)


def adddrBox(c, xb, yb, xl, yl, rad, layer):
    #
    #  Add box with chamfered edges
    pa = pointArray()
    pa.attach(xb + rad, yb)
    pa.attach(xb, yb + rad)
    xc = xb
    yc = yb + yl
    pa.attach(xc, yc - rad)
    pa.attach(xc + rad, yc)
    xc = xc + xl
    pa.attach(xc - rad, yc)
    pa.attach(xc, yc - rad)
    yc = yc - yl
    pa.attach(xc, yc + rad)
    pa.attach(xc - rad, yc)
    xc = xc - xl
    pa.attach(xb + rad, yb)
    c.addPolygon(pointArray(pa), layer)
    return 1


def grarray(c, nr, xg, yg, lyr, low, high, rad):
    #
    #   Draw guard ring array of nr guard rings
    #
    for ind in range(nr):
        grdraw(c, xg, yg, low[ind], high[ind], lyr, rad[ind])


def drcorner(xarr, yarr, npts, layer):
    #
    #   draw 4 corners
    #   points assumed as upper right
    #
    dr.clearPoints()
    dr.activeLayer = layer
    for ind in range(npts):
        dr.point(xarr[ind], yarr[ind])
    dr.polygon()
    # Upper left
    dr.clearPoints()
    dr.activeLayer = layer
    for ind in range(npts):
        dr.point(-xarr[ind], yarr[ind])
    dr.polygon()
    # lower left
    dr.clearPoints()
    dr.activeLayer = layer
    for ind in range(npts):
        dr.point(-xarr[ind], -yarr[ind])
    dr.polygon()
    # lower right
    dr.clearPoints()
    dr.activeLayer = layer
    for ind in range(npts):
        dr.point(xarr[ind], -yarr[ind])
    dr.polygon()

def makeAssy(cel, celllist):
    #
    #  make an assembly of multiple cells, all at 0,0
    #  names contained in celllist
    #
    for ind in range(len(celllist)):
        cnew = l.drawing.findCell(clist[ind])
        p = point(0, 0)
        #    print(cnew, " - ", clist[ind])
        cel.addCellref(cnew, p)


def bpArray(celtgt, bpcell, xur, yur):
    #
    #  make an 2x2 array of corner bond pad of bpcell
    #
    p = point(xur, yur)
    e = celtgt.addCellref(bpcell, p)
    p = point(xur, -yur)
    e = celtgt.addCellref(bpcell, p)
    p = point(-xur, -yur)
    e = celtgt.addCellref(bpcell, p)
    p = point(-xur, yur)
    e = celtgt.addCellref(bpcell, p)


def padArray(homecell, padcell, npad, xstart, xpitch, ypos, smwidth):
    xpad = xstart
    for ind in range(npad):
        p = point(xpad, ypos)
        e = homecell.addCellref(padcell, p)
        #    print(p,"-",)
        #    pa = pointArray()
        #    pa.attach(xpad-smwidth, 0)
        #    pa.attach(xpad+smwidth, 0)
        #    pa.attach(xpad+smwidth, ypos)
        #    pa.attach(xpad-smwidth, ypos)
        #    pa.attach(xpad-smwidth, 0)
        #    e = homecell.addBox(pa, 6)
        xpad = xpad + xpitch


def CSpace(W, S, L):
    # Calculate spaces for uniform mesh (microns) with line width W and total lenght L
    nex = (W - S) / (L + S)
    nfl = int(nex)
    Spac = (W - nfl * L) / (nfl - 1)
    return Spac


def NewCell(CName):
    ec = l.drawing.addCell()
    ec.thisCell.cellName = CName
    ep = ec.thisCell
    return ep


import math
def DrawBump(BP, radius, layer):
    angle = 0.785398163
    sangle = angle / 2
    vertices = [(round(radius * math.cos(i * angle - sangle)), round(radius * math.sin(i * angle - sangle))) for i in
                range(8)]
    vertices.append(vertices[0])
    pts = pointArray()
    x, y = zip(*vertices)
    for i in range(len(vertices)):
        vt = vertices[i]
        pts.attach(x[i], y[i])
    BP.addPolygon(pts, layer)
    return

#	print(vertices)

#  define spacing list 1d mesh with defined limit and  spacing
def Space_1d(width, spacing):
    NCY = max(width // spacing, 1)
    points = [j * spacing for j in range(NCY)]
    mean = sum(points) / len(points)
    out_points = [point - mean for point in points]
    return out_points

def Make_M1M2_Mesh(Prefix, Pad_Layers, Pad_Widths, M12_Pitch, Len_XY, CA_Pitch, CA, V1):
    #
    #   Make mesh contact with Metal 1 and 2
    cd = NewCell(Prefix + "_M1M2_Mesh")
    M1Str = NewCell(Prefix + "_M1_strip")
    e = M1Str.addBox(-Pad_Widths[0], -Len_XY[1] // 2, 2 * Pad_Widths[0], Len_XY[1], Pad_Layers[0])
    ypoints = Space_1d(Len_XY[1], CA_Pitch)
    for k in range(len(ypoints)):
        e = M1Str.addCellref(CA, point(0, ypoints[k]))
    M1_strips = Space_1d(Len_XY[0], M12_Pitch[0])
    for k in range(len(M1_strips)):
        e = cd.addCellref(M1Str, point(M1_strips[k], 0))
    #  Metal two connection
    M2Str = NewCell(Prefix + "_M2_Strip")
    e = M2Str.addBox(-Len_XY[0] // 2, -Pad_Widths[1], Len_XY[0], 2 * Pad_Widths[1], M2)
    for k in range(len(M1_strips)):
        e = M2Str.addCellref(V1, point(M1_strips[k], 0))
    #
    #	M2_Strips = Space_1d(Len_XY[1], M12_Pitch[1])
    #	print(M2_Strips)
    for k in range(len(ypoints) - 1):
        e = cd.addCellref(M2Str, point(0, ypoints[k] + CA_Pitch // 2))
    return cd

def Make_M1M2M3_Mesh(Prefix, Pad_Layers, Pad_Widths, M12_Pitch, Len_XY, CA_Pitch, CA, Via_list):
    #
    #   Make mesh contact with Metal 1,2 and 3
    V1 = Via_list[0]
    V2 = Via_list[1]
    cd = NewCell(Prefix + "_M1M2M3_Mesh")
    M1Str = NewCell(Prefix + "_M1_strip")
    e = M1Str.addBox(-Pad_Widths[0], -Len_XY[1] // 2, 2 * Pad_Widths[0], Len_XY[1], Pad_Layers[0])
    ypoints = Space_1d(Len_XY[1], CA_Pitch)
    for k in range(len(ypoints)):
        e = M1Str.addCellref(CA, point(0, ypoints[k]))
    M1_strips = Space_1d(Len_XY[0], M12_Pitch[0])
    for k in range(len(M1_strips)):
        e = cd.addCellref(M1Str, point(M1_strips[k], 0))

    #  Metal two connection
    M2Str = NewCell(Prefix + "_M2_Strip")
    e = M2Str.addBox(-Len_XY[0] // 2, -Pad_Widths[1], Len_XY[0], 2 * Pad_Widths[1], M2)
    for k in range(len(M1_strips)):
        e = M2Str.addCellref(V1, point(M1_strips[k], 0))
    for k in range(len(ypoints) - 1):
        e = cd.addCellref(M2Str, point(0, ypoints[k] + CA_Pitch // 2))

        # Metal 3 connection
    M3Str = NewCell(Prefix + "_M3_Strip")
    e = M3Str.addBox(-Len_XY[0] // 2, -Pad_Widths[2], Len_XY[0], 2 * Pad_Widths[2], M3)
    for k in range(len(M1_strips)):
        e = M3Str.addCellref(V2, point(M1_strips[k], 0))
    for k in range(len(ypoints) - 1):
        e = cd.addCellref(M3Str, point(0, ypoints[k] + CA_Pitch // 2))
    return cd


def Place_Pad(cd, Name, SXY_Active, SPitch, Slength ):
    astr = NewCell(Name + "_Arr")
    #	Strip_Arrays.append(Strip_name[i] + "_Arr")
    NstX = SXY_Active // SPitch
    NstY = SXY_Active // Slength
    xoff = -(((NstX) - 1) * SPitch) // 2  # bottom left
    yoff = -(((NstY) - 1) * Slength) // 2  # bottom left
    pref = point(xoff, yoff)
    poff = point(xoff + SPitch, yoff + Slength)
    e = astr.addCellrefArray(cd, pref, poff, NstX, NstY)
    return astr

def erdrawOD(c, xgr, ygr, wlow, whigh, layer, radius, OD, OD_inset):
    # draw implant ring shape weith inset active layer (OD)
    erdraw(c, xgr, ygr, wlow, whigh, layer, radius)
    erdraw(c, xgr, ygr, wlow-OD_inset, whigh-OD_inset, OD, radius)

def adddrBoxOD(c, xb, yb, xl, yl, rad, layer, OD, OD_inset):
    # draw implant chamfered box shape weith inset active layer (OD)
    adddrBox(c, xb, yb, xl, yl, rad, layer)
    adddrBox(c, xb+OD_inset, yb+OD_inset, xl-2*OD_inset, yl-2*OD_inset, rad, OD)

def make_EdgeArray(Cname, BCell, ECell, NX, NY, DX, DY):
    #  make an NX by NY array of BCells with the ECell at the 4 edges.
    # at some point need to add mirroring
    # NX, NY
    dcell = NewCell(Cname)
    xecell = (DX*NX/2) - DX/2
    yecell = (DY*NY/2) - DY/2
    dcell.addCellref(ECell, point(xecell, yecell))
    dcell.addCellref(ECell, point(-xecell, yecell))
    dcell.addCellref(ECell, point(-xecell, -yecell))
    dcell.addCellref(ECell, point(xecell, -yecell))
    p1 = point(-xecell + DX, -yecell)
    p2 = point(-xecell + 2*DX, -yecell+DY)
    dcell.addCellrefArray(BCell, p1, p2, NX-2, NY)
    p1 = point(-xecell, -yecell+DY)
    p2 = point(-xecell, -yecell+2*DY)
    dcell.addCellrefArray(BCell, p1, p2, 1, NY-2)
    p1 = point(xecell, -yecell+DY)
    p2 = point(xecell, -yecell+2*DY)
    dcell.addCellrefArray(BCell, p1, p2, 1, NY-2)
    return dcell

import LayoutScript
from LayoutScript import *
#from TWR_functions.TWR_fns import *

l = project.newLayout();  # open new instance of layout class
global dr
dr = l.drawing
# pointer to the main drawing
setup.gdsTextToPolygon = True
setup.gdsTextToPolygonDefaultWidth = 200000
setup.defaultTextWidth = 200000

SetUp = setup()  # work around as static string variables are not handled correctly
# Layer definitions

#
#   Import SLAC portions
#
dr.importFile("/Users/ronlipton/Dropbox/Programming/TWR_layout/SLAC_layouts/compile_border_v6.gds")

OTL = 0  # outline for drawing
OD = 1 # Defines active window
JTE = 116  # Junction termination extension IMPLANT (NP-JTE)
PGN = 117  # boron gain layer IMPLANT (NC)
NC = 117
ACN = 55  # AC phos layer IMPLANT (NQ)
ZP = 58  # top passivation openings
NPL = 66  # n+ IMPLANT (NP)
NP = 66
MET = 44  # original design Metal
M1 = 43  # METAL 1
M2 = 47  # Metal 2
M3 = 49  # Metal 3
V1 = 46  # Via from M1 to M2
V2 = 48  # Hole to connect M3 to M2
ZG = 56  # Hole to connect ZA to M3
ZA = 57  # 4th Metal (Al,Top Metal)
CON = 25  # Contact (CA)
CA = 25  # contact
# ACC = 9
DJP = 18  # Deep junction p IMPLANT (PX)
DJN = 16  # Deep junction n IMPLANT (DW?)
PSB = 21  # p substrate contact IMPLANT (PD??)
PD = 21
# PST = 117  # p stop (NC)
ND = 19 # n contact IMPLANT
PW = 78 # P-well IMPLANT

# inset of active
OD_inset = 100
# Number of guard rings
nrings = 5
# Device pixel pitches
Pitch = [50000, 100000]
# Pixel Rows
NPXRow = [100, 50]
# Pixel Columns
NPXCol = [100, 50]
# Cell size
XYCell = 6000000
XYCell_2 = XYCell // 2
# Reticule subcells
NCellx = 4
NCelly = 5

# Active dimension of pixel arays
Length_2 = 2500000  # half length of chip
Lng = 2 * Length_2
length = [Lng]
# Reach through LGAD
# RT LGAD rows, columns
RTRow = 8
RTCol = 8
# RT LGAD Pitch
RTPitch = [600000, 600000]  # need to replace x, y (below) by an index
RTPitchx = 600000
RTPitchy = 600000
# RT inter pixel gap
RTGap = 100000
# RT Active length
RTLenx = RTPitchx - RTGap  # active lengths
RTLeny = RTPitchy - RTGap
RTRound = 5000  # corner rounding radius
RTOXRad = 6000
RTMRad = RTLenx // 4

# JTE parameters
JTERound = 25000
JTEInset = 2000
JTEWidth = 20000

# Pixel isolation half width
PIWid_2 = 1500

#	contact half-width
cwidth = 6000  # width of contact
CA_Width = 40  # half width of Tower contact
CA_Space = 200  # Spacing for CA contact array

V1_Width = 50  # V1 half width
V1_Space = 200

V2_Width = 50  # V2 half width
V2_Space = 200

V3_Width = 50  # V3 half width
V3_Space = 300

ZG_Width = 1500  # M3 to top metal half-width
ZG_Space = 9000

# Make Via list
Via_Widths = [V1_Width, V2_Width, ZG_Width]
Via_Spaces = [V1_Space, V2_Space, ZG_Space]
Via_Layers = [V1, V2, ZG]

# metal ine half-widths
M1_Width_2 = 500
M2_Width_2 = 2000
M3_Width_2 = 3000
M1_Width = 2 * M1_Width_2
M2_Width = 2 * M2_Width_2
M3_Width = 2 * M3_Width_2

Pad_Widths = [M1_Width_2, M2_Width_2, M3_Width_2]
# Max area square pads
Pad_Cells = ["M1_Pad", "M2_Pad", "M3_Pad"]
Via_Cells = ["V1_Via", "V2_Via", "ZG_Via"]
Pad_Layers = [M1, M2, M3]
Pad_Ref = []
# Make Metal Pad cells
Pad_List = []
for i in range(len(Pad_Cells)):
    m1p = NewCell(Pad_Cells[i])
    m1p.addBox(-Pad_Widths[i], -Pad_Widths[i], 2 * Pad_Widths[i], 2 * Pad_Widths[i], Pad_Layers[i])
    Pad_List.append(m1p)

M1Pad = l.drawing.findCell(Pad_Cells[0])
#  CA contact Cells
#
ml = NewCell("Contact")
ml.addBox(-CA_Width, -CA_Width, 2 * CA_Width, 2 * CA_Width, CA)
#  4x4 contact Cell
CA4x4 = NewCell("Contact_4x4")
e = CA4x4.addCellrefArray(ml, point(-CA_Space // 2, -CA_Space // 2), point(CA_Space // 2, CA_Space // 2), 2, 2)
e = CA4x4.addCellref(M1Pad, point(0, 0))
#  Add ohmic n contact region for CA cell
CA4x4.addBox(-M1_Width_2+1, -M1_Width_2+1, M1_Width-2, M1_Width-2, ND)
CA4x4.addBox(-M1_Width_2+1+OD_inset, -M1_Width_2+1+OD_inset, M1_Width-2-2*OD_inset, M1_Width-2-2*OD_inset, OD)
#   Distinguish p and n
CAN4x4 = CA4x4


CAP4x4 = NewCell("PContact_4x4")
e = CAP4x4.addCellrefArray(ml, point(-CA_Space // 2, -CA_Space // 2), point(CA_Space // 2, CA_Space // 2), 2, 2)
e = CAP4x4.addCellref(M1Pad, point(0, 0))
#  Add ohmic n contact region for CA cell
CAP4x4.addBox(-M1_Width_2+1, -M1_Width_2+1, M1_Width-2, M1_Width-2, NP)
CAP4x4.addBox(-M1_Width_2+1+OD_inset, -M1_Width_2+1+OD_inset, M1_Width-2-2*OD_inset, M1_Width-2-2*OD_inset, OD)
#   Distinguish p and n

#  Via Cell list
Via_List = []
for i in range(len(Via_Cells)):
    ml = NewCell(Via_Cells[i])
    ml.addBox(-Via_Widths[i], -Via_Widths[i], 2 * Via_Widths[i], 2 * Via_Widths[i], Via_Layers[i])
    #  4x4 contact Cell
    ml4x4 = NewCell(Via_Cells[i] + "_4x4")
    VS = Via_Spaces[i] // 2
    e = ml4x4.addCellrefArray(ml, point(-VS, -VS), point(VS, VS), 2, 2)
    if i == 2:
        Via_List.append(ml)
    else:
        Via_List.append(ml4x4)
#    print(Via_List[i])

BP_M3_Via_80 = NewCell("BP_M3_via80")
BPM3Width = 73000
BPM3Length = 3000
e = adddrBox(BP_M3_Via_80, -BPM3Width//2, -BPM3Length//2, BPM3Width, BPM3Length, 0, ZG)
BP_M3_Via_60 = NewCell("BP_M3_via60")
BPM3Width = 42000
BPM3Length = 3000
e = adddrBox(BP_M3_Via_60, -BPM3Width//2, -BPM3Length//2, BPM3Width, BPM3Length, 0, ZG)

# 80 micron bond pad
BP_80 = NewCell("Bond_Pad_80")
padWidth_80 = 77000
padLength_80 = 205000
OXLength_80 = 190000
OXWidth_80 = 71000
e = adddrBox(BP_80, -OXWidth_80// 2, -OXLength_80//2, OXWidth_80, OXLength_80, 0, ZA)
e = adddrBox(BP_80, -padWidth_80 // 2, -padLength_80 // 2, padWidth_80, padLength_80, 0, ZP)
e = BP_80.addCellref(BP_M3_Via_80, point(0, padLength_80//2-2500))
e = BP_80.addCellref(BP_M3_Via_80, point(0, -padLength_80//2+2500))

#  60 micron bond pad
BP_60 = NewCell("Bond_Pad_60")
padWidth_60 = 57000
padLength_60 = 205000
OXLength_60 = 190000
OXWidth_60 = 51000
e = adddrBox(BP_60, -OXWidth_60// 2, -OXLength_60//2, OXWidth_60, OXLength_60, 0, ZA)
e = adddrBox(BP_60, -padWidth_60 // 2, -padLength_60 // 2, padWidth_60, padLength_60, 0, ZP)
e = BP_60.addCellref(BP_M3_Via_60, point(0, padLength_60//2-2500))
e = BP_60.addCellref(BP_M3_Via_60, point(0, -padLength_60//2+2500))

#   Bump Pad_Cell
BCell = NewCell("BumpPad")
DrawBump(BCell, 8300, ZA)  # Standard bump pad - check dimensions
DrawBump(BCell, 6500, ZP)


# Metal-Via stack for pads not including CA
M23V23Z = NewCell("M1M2V1V2ZG_Pad")
for Layer in range(3):
    M23V23Z.addCellref(Pad_List[Layer], point(0, 0))
    M23V23Z.addCellref(Via_List[Layer], point(0, 0))
    M23V23Z.addCellref(BCell, point(0, 0))

Letters = []
CNames = "abcdefghijklmn"
for j in range(13):
    num = l.drawing.addCell()
    cname = "c_" + CNames[j]
    num.thisCell.cellName = cname
    ##	e = l.drawing.setCell(cname)
    cCell = num.thisCell
    e = cCell.addText(MET, point(0, 0), CNames[j])
    e.setWidth(200000)
    l.drawing.textSelect()
    l.drawing.toPolygon()
    Letters.append(cCell)

ntype = 3
CA_Contact = l.drawing.findCell("Contact_4x4")

##############################################
#  make DC strxcel
#  def makeStrip(cl, activeLength, sgap, siwidth, scwidth, moffs, LP, LM, LC):
#
#   Active dimensions for 3x3 strip arrays
#
Str_Length = 2000000
Str_Length_2 = Str_Length//2
empty_cell = NewCell("empty")
Strip_Pitch = [12500, 50000, 100000, 50000, 100000]
Strip_Length = [50000, Str_Length, Str_Length, Str_Length, Str_Length]
Strip_name = ["Strp125", "Strp50", "Strp100", "AC50","AC100"]
Strip_Contacty = ["Strp125CY", "Strp50CY", "Strp100CY", "Strp50CY", "Strp100CY"]
ConCell = [CA_Contact, CA_Contact, CA_Contact, empty_cell, empty_cell]
Imp_Lyr = [NPL, NPL, NPL, 0, 0]
Strip_M1_Lines = []
Strip_M2 = []

STInset = 5000  # Strip implant inset from edge
ST_CA_Pitch = 10000  # Strip contact pitch
ST_M1_Pitch = 2000  # strip M1 contact strip Pitch
SXY_Active = 5000000
STXY_Active = 2000000
offset_3mm = 1500000
Strip_imp_width = []
Strip_Arrays = []

gainsurr = 40000  # gain surround of AC
gainround = 25000  # radius of reference arc for GR edges
gsurr = 10000

# DJ Implants
DJNinset = 15000  # inset of deep N (phos) from active length
DJPinset = 5000  # inset of deep P (Boron) from active length
DJRound = 5000  # edge rounding radius
DJRound = 2000  # edge rounding radius
DJInset = 10000


name = "DJ_P_3mm"
activeLength = STXY_Active
cpad = NewCell(name)
e = adddrBoxOD(cpad, -(activeLength + DJPinset) // 2, -(activeLength + DJPinset) // 2, \
             activeLength + DJPinset, activeLength + DJPinset, DJRound, DJP, OD, OD_inset)

name = "DJ_N_3mm"
dpad  = NewCell(name)
e = adddrBoxOD(dpad, -(activeLength + DJNinset) // 2, -(activeLength + DJNinset) // 2, \
             activeLength + DJNinset, activeLength + DJNinset, DJRound+4000, DJN, OD, OD_inset)
name = "DJ_3mm_2x2"
bstr = NewCell(name)
pref = point(-offset_3mm, -offset_3mm)
poff = point(offset_3mm, offset_3mm)
bstr.addCellrefArray(cpad, pref, poff, 2, 2)
bstr.addCellrefArray(dpad, pref, poff, 2, 2)

gainsurr = 40000  # gain surround of AC
gainround = 5000  # radius of reference arc for GR edges
gsurr = 10000
#  AC gain and casthode implants
name = "ACP_3mm"
activeLength = STXY_Active
cpad = NewCell(name)
e = adddrBoxOD(cpad, -(activeLength + gsurr) // 2, -(activeLength + gsurr) // 2, \
             activeLength + gsurr, activeLength + gsurr, gainround, PGN, OD, OD_inset)
# AC layer
name = "ACC-3mm"
activeLength = STXY_Active
wbox = (activeLength + gainsurr) // 2
bpad = NewCell(name)
e = adddrBoxOD(bpad, -(activeLength + gainsurr) // 2, -(activeLength + gainsurr) // 2, \
                 activeLength + gainsurr, activeLength + gainsurr, gainround, ACN, OD, OD_inset)
# Draw JTE for 3mm
cjte_AC = NewCell("JTE_3mm")
wbox = (activeLength + gainsurr) // 2 -JTERound
wbox = (activeLength + gainsurr) // 2 -gainround -2500
xpm = wbox
xm = [xpm, -xpm, -xpm, xpm]
ym = [xpm, xpm, -xpm, -xpm]
# erdrawOD(cjte_AC, xm, ym, JTEInset - 2000, JTEWidth, JTE, JTERound, OD, OD_inset)
erdrawOD(cjte_AC, xm, ym, JTEInset - 2000, JTEWidth, JTE, JTERound, OD, OD_inset)
#  Add ND
# erdrawOD(cjte_AC, xm, ym, JTEInset - 2000, JTEWidth, ND, JTERound, OD, OD_inset)
erdrawOD(cjte_AC, xm, ym, JTEInset - 2000, JTEWidth, ND, JTERound, OD, OD_inset)

name = "ACP_3mm_2x2"
bstr = NewCell(name)
pref = point(-offset_3mm, -offset_3mm)
poff = point(offset_3mm, offset_3mm)
bstr.addCellrefArray(cpad, pref, poff, 2, 2)
bstr.addCellrefArray(bpad, pref, poff, 2, 2)
bstr.addCellrefArray(cjte_AC, pref, poff, 2, 2)

for i in range(len(Strip_Pitch)):
    Strip_imp_width.append(Strip_Pitch[i] - STInset)
# print(Strip_imp_width)
for i in range(len(Strip_Pitch)):
    sgap = 4000
    SPitch = Strip_Pitch[i]
    Slength = Strip_Length[i]
    DCStripn = Strip_name[i]

    #	siwidth = Strip_Pitch[i] - STInset
    siwidth = Strip_imp_width[i]
    #	scwidth = 2000
    smwidth = M1_Width
    cd = NewCell(DCStripn)
    lng = Slength - sgap
    wid = siwidth
    #	lox = wid
    STRound = 2000  # edge rounding parameter for strips
    STMSurr = 1000
    e = adddrBoxOD(cd, -wid // 2, -lng // 2, wid, lng, STRound, Imp_Lyr[i], OD, OD_inset)
    #
    #   Add edge field plate
    #
    ypm = lng // 2 - STRound
    xpm = wid // 2 - STRound
    xm = [xpm, -xpm, -xpm, xpm]
    ym = [ypm, ypm, -ypm, -ypm]
    erdraw(cd, xm, ym, 1000, 1000, M2, STRound)
#
#   mesh subroutine
    Met_Pitch = [ST_CA_Pitch//5, ST_M1_Pitch]
    lxy = [siwidth, lng]
#   ce = Make_M1M2_Mesh(Strip_name[i], Pad_Layers, Pad_Widths, Met_Pitch, lxy, ST_CA_Pitch, ConCell[i], Via_List[1] )
    ce = Make_M1M2M3_Mesh(Strip_name[i], Pad_Layers, Pad_Widths, Met_Pitch, lxy, ST_CA_Pitch, ConCell[i], Via_List)
    cd.addCellref(ce, point(0, 0))


# P-stop for DC coupled
    if Strip_name[i][0:2] == 'AC':
        pass
    # Add p-stop
    else:
        rad = 0
        pswidth = 500
        psx = SPitch // 2
        pslen = Slength // 2
        xps = [psx, -psx, -psx, psx]
        yps = [pslen, pslen, -pslen, -pslen]
        erdrawOD(cd, xps, yps, pswidth, 0, PW, rad, OD, OD_inset)

    #      Add bond pads
    BPinset = 255000
    Row2inset = BPinset + 260000
    if (SPitch == 100000):
        cd.addCellref(BP_80, point(0, Slength // 2 - BPinset))
        cd.addCellref(BP_80, point(0, -Slength // 2 + BPinset))

        bstr = NewCell(Strip_name[i] + "_Arr")
        Strip_Arrays.append(Strip_name[i] + "_Arr")

        astr = NewCell(Strip_name[i] + "_Arr3mm")
#       e = ST_add_Array(astr, cd, SPitch, Slength, STXY_Active)
        NstX = STXY_Active // SPitch
        NstY = STXY_Active // Slength
        xoff = -(((NstX) - 1) * SPitch) // 2  # bottom left
        yoff = -(((NstY) - 1) * Slength) // 2  # bottom left
        pref = point(xoff, yoff)
        poff = point(xoff + SPitch, yoff + Slength)
        e = astr.addCellrefArray(cd, pref, poff, NstX, NstY)

        pref = point(-offset_3mm, -offset_3mm)
        poff = point(offset_3mm, offset_3mm)
        bstr.addCellrefArray(astr, pref, poff, 2, 2)

    if (SPitch == 50000):
    #   make two copies of strip for staggered pads
        dr.setCell(DCStripn)
        dr.selectAll()
        dr.point(-25000,0)
        dr.move()
        dr.point(50000,0)
        dr.copy()
        cd.addCellref(BP_60, point(-25000, Slength // 2 - BPinset))
        cd.addCellref(BP_60, point(-25000, -Slength // 2 + BPinset))
        cd.addCellref(BP_60, point(25000, Slength // 2 - Row2inset))
        cd.addCellref(BP_60, point(25000, -Slength // 2 + Row2inset))

        bstr = NewCell(Strip_name[i] + "_Arr")
        Strip_Arrays.append(Strip_name[i] + "_Arr")

        astr = NewCell(Strip_name[i] + "_Arr3mm")
        NstX = STXY_Active // (SPitch*2)
        NstY = STXY_Active // Slength
        xoff = -(((NstX) - 1) * (SPitch*2)) // 2    # bottom left
        yoff = -(((NstY) - 1) * Slength) // 2    # bottom left
        pref = point(xoff, yoff)
        poff = point(xoff + (SPitch*2), yoff + Slength)
        e = astr.addCellrefArray(cd, pref, poff, NstX, NstY)

        pref = point(-offset_3mm, -offset_3mm)
        poff = point(offset_3mm, offset_3mm)
        bstr.addCellrefArray(astr, pref, poff, 2, 2)


    if (SPitch == 12500):
        astr = NewCell(Strip_name[i] + "_Arr")
        Strip_Arrays.append(Strip_name[i] + "_Arr")
        NstX = SXY_Active // SPitch
        NstY = SXY_Active // Slength
        xoff = -(((NstX) - 1) * SPitch) // 2  # bottom left
        yoff = -(((NstY) - 1) * Slength) // 2  # bottom left
        pref = point(xoff, yoff)
        poff = point(xoff + SPitch, yoff + Slength)
        e = astr.addCellrefArray(cd, pref, poff, NstX, NstY)

#
#   Stqggered 125 um bump array cell
#       Current design - left - straight array of bumps, left - staggered pads
#
indx_125 = 0
Bump_125 = NewCell("Bump_125_Arr")
# Bump_125_stgr = NewCell("Bump_125_Stg")  # staggered bumps
SPitch = Strip_Pitch[indx_125]
Slength = Strip_Length[indx_125]
NstX = SXY_Active // SPitch
NstY = SXY_Active // Slength
xoff = -(((NstX) - 1) * SPitch) // 2
yoff = -(((NstY) - 1) * Slength) // 2
# pref = point(xoff, yoff)
# poff = point(xoff + SPitch, y = point(xoff, yoff)
# # poff = point(xoff + SPitch, yoff + Slength)
# # e = Bump_125.addCellrefArray(M23V23Z, pref, poff, NstX//2, NstY)off + Slength)
# e = Bump_125.addCellrefArray(M23V23Z, pref, poff, NstX//2, NstY)
#
#   Staggered pad section - MAKE IT ALL STAGGERED 9/27/24
#
# Lower pad
# pref = point(SPitch//2, yoff-ST_CA_Pitch)
pref = point(xoff, yoff-ST_CA_Pitch)
# poff = point(SPitch//2 + SPitch*2 , yoff + Slength - ST_CA_Pitch)
poff = point(xoff + SPitch*2, yoff + Slength - ST_CA_Pitch)
e = Bump_125.addCellrefArray(M23V23Z, pref, poff, NstX//2, NstY)
# upper pad
pref = point(xoff + SPitch, yoff+ST_CA_Pitch)
poff = point(xoff + SPitch*3 , yoff + Slength + ST_CA_Pitch)
e = Bump_125.addCellrefArray(M23V23Z, pref, poff, NstX//2, NstY)

st125cell = l.drawing.findCell(Strip_name[indx_125] + "_Arr")
e = st125cell.addCellref(Bump_125, point(0, 0))

##############################################
#
#  DJ Pixel
#
DJName = ["DJ_50", "DJ_100"]
DJAName = ["DJA_50", "DJA_100"]
DJMRad = 22000  # Radius of m1 layer
DJMSur = 5000  # 2 x metal surround of implant

## moff = 80000  # for strips
DJIWid_2 = 500
DJM2_Wid_2 = 1000
Len_XY = [Pitch[0] - DJInset, Pitch[0] - DJInset]
M12_Pitch = [2 * M1_Width, 2 * M2_Width]

for i in range(len(Pitch)):

    cpad = NewCell(DJName[i] + "Base")
    lng = Pitch[i] - DJInset
    e = adddrBoxOD(cpad, -lng // 2, -lng // 2, lng, lng, DJRound, NPL, OD, OD_inset)
    crad = 12000
    #  Mesh Contact
    Len_XY = [Pitch[i] - DJInset, Pitch[i] - DJInset]
    M12_Pitch = [2 * M1_Width, 2 * M2_Width]
#    Ref = Make_M1M2_Mesh(DJName[i], Pad_Layers, Pad_Widths, M12_Pitch, Len_XY, ST_CA_Pitch, CA4x4, Via_List[0])
    Ref = Make_M1M2M3_Mesh(DJName[i], Pad_Layers, Pad_Widths, M12_Pitch, Len_XY, ST_CA_Pitch, CA4x4, Via_List)
    e = cpad.addCellref(Ref, point(0, 0))

#    xpm = lng//2 - DJInset//2
    xpm = lng//2 - DJRound
    xm = [xpm, -xpm, -xpm, xpm]
    ym = [xpm, xpm, -xpm, -xpm]
    erdraw(cpad, xm, ym, DJM2_Wid_2, DJM2_Wid_2, M2, DJRound)
    # Add isolation - assume square
    wbox = Pitch[i] // 2
    xpm = wbox
    xm = [xpm, -xpm, -xpm, xpm]
    ym = [xpm, xpm, -xpm, -xpm]
    erdrawOD(cpad, xm, ym, DJIWid_2, 0, PW, 0, OD, OD_inset)
# cell without central bump
    dpad = NewCell(DJName[i]+"_NoBump")
    e = dpad.addCellref(cpad, point(0, 0))
# cell with bump
    epad = NewCell(DJName[i])
    e = epad.addCellref(cpad, point(0, 0))

    e = epad.addCellref(M23V23Z, point(0, 0))  # bump in cell  center
    #
    # pixel array cell
    apad = NewCell(DJAName[i])
    xoff = -(((NPXRow[i]) - 1) * Pitch[i]) // 2
    yoff = -(((NPXCol[i]) - 1) * Pitch[i]) // 2
    pref = point(xoff, yoff)
    poff = point(xoff + Pitch[i], yoff + Pitch[i])
    e = apad.addCellrefArray(epad, pref, poff, NPXRow[i], NPXCol[i])

    fpad = NewCell(DJAName[i] + "_noBump")
    e = fpad.addCellrefArray(dpad, pref, poff, NPXRow[i], NPXCol[i])
#    epad = NewCell("Edge_Array")
#    Edge_Cell = DJName[i] + "_Edge"
#    print(Edge_Cell)
#    dcell = make_EdgeArray(Edge_Cell, dpad, cpad, NPXRow[i], NPXCol[i], Pitch[i], Pitch[i])

# DJ Implants
name = "DJ_PN"
activeLength = Lng
cpad = NewCell(name)
dinset = 1500
rinset = dinset
djcorner = activeLength + DJPinset
e = adddrBoxOD(cpad, -(djcorner) // 2, -(djcorner) // 2, \
             djcorner, djcorner, DJRound+dinset, DJP, OD, OD_inset)

djcorner = activeLength + DJNinset
e = adddrBoxOD(cpad, -(djcorner) // 2, -(djcorner) // 2, \
             djcorner, djcorner, DJRound+dinset+3500, DJN, OD, OD_inset)


##############################################
#
# Reach through pixel and array
#
RTname = "RTPixel"
## RTMSurr = -10000 # overlap with cathode
RTPad_inset = 25000
## RTOin = 30000
RTGin = 8000  # gain layer inset
RTCname = "RTArray"
RT_M1_Lines = []
RT_M1_Pitch = 2 * M1_Width
RT_M2_Pitch = 2 * M2_Width

rtpad = NewCell(RTname)
RTPx_2 = RTPitchx // 2
RTPy_2 = RTPitchy // 2
rtpad.addBox(-RTPx_2, -RTPy_2, RTPitchx, RTPitchy, OTL)
lng = RTPitchx - RTGap
e = adddrBoxOD(rtpad, -lng // 2, -lng // 2, lng, lng, RTRound, NPL, OD, OD_inset)
e = adddrBoxOD(rtpad, -lng // 2 + RTGin, -lng // 2 + RTGin, lng - 2 * RTGin, lng + -2 * RTGin, RTRound, PGN, OD, OD_inset)
# add JTE
wbox = lng // 2 - RTGin
xpm = wbox
xm = [xpm, -xpm, -xpm, xpm]
ym = [xpm, xpm, -xpm, -xpm]
erdrawOD(rtpad, xm, ym, JTEInset, JTEWidth, JTE, RTRound, OD, OD_inset)
erdrawOD(rtpad, xm, ym, JTEInset, JTEWidth, JTE, RTRound, ND, OD_inset)
# erdraw(rtpad, xm, ym, JTEInset-OD_inset, JTEWidth-OD_inset, OD, RTRound)

RT_Pad_len = lng - 20000
#
#   Add contact mesh
#
Len_XY = [RT_Pad_len, RT_Pad_len]
Ref = Make_M1M2M3_Mesh("RT", Pad_Layers, Pad_Widths, M12_Pitch, Len_XY, ST_CA_Pitch, CA4x4, Via_List)
rtpad.addCellref(Ref, point(0, 0))
wbox = RTPitchx // 2 - RTRound - PIWid_2

# Pad at center
rtpad.addCellref(BP_80, point(0,0))

xpm = wbox
xm = [xpm, -xpm, -xpm, xpm]
ym = [xpm, xpm, -xpm, -xpm]
erdrawOD(rtpad, xm, ym, PIWid_2, PIWid_2, PW, RTRound, OD, OD_inset)
# erdraw(rtpad, xm, ym, PIWid_2-OD_inset, PIWid_2-OD_inset, OD, RTRound)
#
# Array cell
#
# RTpd = l.drawing.addCell()
# RTpd.thisCell.cellName =  RTCname
RTpad = NewCell(RTCname)
xoff = -(((RTRow) - 1) * RTPitchx) // 2
yoff = -(((RTCol) - 1) * RTPitchy) // 2
pref = point(xoff, yoff)
poff = point(xoff + RTPitchx, yoff + RTPitchy)
e = RTpad.addCellrefArray(rtpad, pref, poff, RTRow, RTCol)

##############################################
# AC Pads
# metal dimensions
ACMetx = [30000, 60000]  # 50 and 100 micron pad size
ACMety = [30000, 60000]
ACPitch_M1X = [2000, 2000]  # metal pitch (for 50% coverage)
ACPitch_M2Y = [8000, 8000]
# M1_Strips = []

name = ["ACPad_50", "ACPad_100"]
aname = ["AC_Array_50", "AC_Array_100"]

for i in range(len(ACMetx)):
    cd = NewCell(name[i])
    Len_XY = [ACMetx[i], ACMety[i]]
#    print("Lenxy=",Len_XY)
    Ref = Make_M1M2M3_Mesh(name[i], Pad_Layers, Pad_Widths, M12_Pitch, Len_XY, ST_CA_Pitch, CA4x4, Via_List)
    cd.addCellref(Ref, point(0, 0))

    # Outline
    x = -Pitch[i] // 2
    y = -Pitch[i] // 2
    lx = Pitch[i]
    ly = Pitch[i]
    e = cd.addBox(x, y, lx, ly, OTL)

    apad = NewCell(aname[i])
    xoff = -(((NPXRow[i]) - 1) * Pitch[i]) // 2
    yoff = -(((NPXCol[i]) - 1) * Pitch[i]) // 2
    pref = point(xoff, yoff)
    poff = point(xoff + Pitch[i], yoff + Pitch[i])
    e = apad.addCellrefArray(cd, pref, poff, NPXRow[i], NPXCol[i])
    e = cd.addCellref(M23V23Z, point(0, 0))

# Gain Layer
#name = "Gain_Layer"
for i in range(len(length)):
    activeLength = length[i]
    #	wbox = (activeLength+gainsurr)//2
    #epd = l.drawing.addCell()
    #epd.thisCell.cellName = name
    cpad = NewCell("Gain_Layer")
    #cpad = epd.thisCell
    e = adddrBoxOD(cpad, -(activeLength + gsurr) // 2, -(activeLength + gsurr) // 2, \
                 activeLength + gsurr, activeLength + gsurr, gainround, PGN, OD, OD_inset)
# AC gain Layer 3mm


# AC layer
name = "AC_Layer"

for i in range(len(length)):
    activeLength = length[i]
    wbox = (activeLength + gainsurr) // 2
    #	epd = l.drawing.addCell()
    #	epd.thisCell.cellName = name
    #	cpad = epd.thisCell
    cpad = NewCell(name)
    e = adddrBoxOD(cpad, -(activeLength + gainsurr) // 2, -(activeLength + gainsurr) // 2, \
                 activeLength + gainsurr, activeLength + gainsurr, gainround, ACN, OD, OD_inset)

# JTE
cjte = NewCell("JTE")
wbox = (activeLength + gainsurr) // 2 -JTERound
xpm = wbox
xm = [xpm, -xpm, -xpm, xpm]
ym = [xpm, xpm, -xpm, -xpm]
erdrawOD(cjte, xm, ym, JTEInset, JTEWidth, JTE, JTERound, OD, OD_inset)
erdrawOD(cjte, xm, ym, JTEInset, JTEWidth, ND, JTERound, OD, OD_inset)

#
#	Psubstrate Contact
#
PSWid = 50000  # width of contact in cell
PSubS = NewCell("PSubC")
PSubS.addBox(-XYCell_2, -XYCell_2, PSWid, XYCell, PSB)
PSubS.addBox(-XYCell_2, XYCell_2 - PSWid, XYCell, PSWid, PSB)
PSubS.addBox(XYCell_2 - PSWid, -XYCell_2, PSWid, XYCell, PSB)
PSubS.addBox(-XYCell_2, -XYCell_2, XYCell, PSWid, PSB)
#
#	Cell Outline
#
cotl = "Cell_Outline"
Outline = NewCell(cotl)
Outline.addBox(-XYCell_2, -XYCell_2, XYCell, XYCell, OTL)

CellList = []
#
#	 AC LGAD Cells
#
# ACL_list = ["Assy_AC_", "Assy_AC_"]
Border_List = ["6mm_50um_pitch_bumps","6mm_100um_pitch_bumps"]
AC_Type = ["50", "100"]
for i in range(len(AC_Type)):
    cname = "Assy_AC_" + AC_Type[i]
    Assy_AC = NewCell(cname)
    ACPName = "AC_Array_" + AC_Type[i]
    #	print(ACPName)
    clist = ["AC_Layer", ACPName, cotl, "JTE", "PSubC", "Gain_Layer",Border_List[i]]
    makeAssy(Assy_AC, clist)
    CellList.append(cname)
#
#	 DJ LGAD Cells
#
DJ_Type = ["50", "100"]
for i in range(len(DJ_Type)):
    cname = "AssyDJ_" + DJ_Type[i]
    Assy_DJ = NewCell(cname)
    DJPName = "DJA_" + DJ_Type[i]
    clist = [DJPName, cotl, "JTE", "PSubC", "DJ_PN",Border_List[i]]
    makeAssy(Assy_DJ, clist)
    CellList.append(cname)

AC_Type = ["50", "100"]
for i in range(len(AC_Type)):
    cname = "AssyAC_Str" + AC_Type[i]
    Assy_AC = NewCell(cname)
    ACCName = "AC" +AC_Type[i]+"_Arr"
    clist = [ACCName, cotl, "PSubC", "ACP_3mm_2x2","3mm_quad"]
    makeAssy(Assy_AC, clist)
    CellList.append(cname)

for i in range(3):
    cname = "Assy_" + Strip_Arrays[i]
    if Strip_Arrays[i].startswith("Strp125"):
        border = "6mm_with_pads"
        djimp = "DJ_PN"
    else:
        border = "3mm_quad"
        djimp = "DJ_3mm_2x2"
    STX_Ass = NewCell(cname)
    clist = ("PSubC", Strip_Arrays[i], border)
    makeAssy(STX_Ass, clist)
    CellList.append(cname)
    #
    #	DJ Version
    cname = "AssyDJ" + Strip_Arrays[i]
    STDJ_Ass = NewCell(cname)
    clist = ("PSubC", Strip_Arrays[i], djimp ,border)
    makeAssy(STDJ_Ass, clist)
    CellList.append(cname)

cnam = "AssyRT"
RTAss = NewCell(cnam)
clist = (RTCname, "PSubC","6mm_with_pads")
makeAssy(RTAss, clist)
CellList.append(cnam)

for i in range(len(DJ_Type)):
    cname = "AssyDJ_NB_" + DJ_Type[i]
    Assy_DJ_NB = NewCell(cname)
    DJPName = "DJA_" + DJ_Type[i] + "_noBump"
    clist = [DJPName, cotl, "JTE", "PSubC", "DJ_PN",Border_List[i]]
    makeAssy(Assy_DJ_NB, clist)
    CellList.append(cname)
#
#	Make Reticule
#
Ret_W = 25500000  # reticule width
Ret_H = 32500000  # reticule height
Ret_W_2 = Ret_W // 2
Ret_H_2 = Ret_H // 2

X0 = -XYCell * (NCellx - 1) // 2
Y0 = -XYCell * (NCelly - 1) // 2

c = l.drawing.addCell()
cnam = "Reticule"
c.thisCell.cellName = cnam
Ret = c.thisCell
Ret.addBox(-Ret_W_2, -Ret_H_2, Ret_W, Ret_H, OTL)

for i in range(len(CellList)):
    RRow = i // NCellx
    RCol = i % NCellx
    Rx = X0 + RCol * XYCell
    Ry = Y0 + RRow * XYCell
    cnew = l.drawing.findCell(CellList[i])
    p = point(Rx, Ry)
#   print(cnew, CellList[i])
#    print(RRow)
#    print(RCol)
    Ret.addCellref(cnew, p)

print(CellList)
# import os
from pathlib import Path
home_directory = Path.home()
#print(home_directory)
gdsversion = "19.3"
gdsfile = str(home_directory) + "/Dropbox/Programming/TWR_layout/TWR_" + gdsversion + ".gds"
#print(gdsfile)
l.drawing.saveFile(gdsfile)

print("Python script completed")
