import LayoutScript
from LayoutScript import *

def BoxDraw(c, xgr, ygr, radius, whigh, wlow, layer):
    #
    #       Draw Top, bottom, left, right boxes to enclose a pixel or strip area
    #       usded with edger draw to complete a chamfered box
    #   ie y of top =- ygr(0) + radius + whigh ...
    #   c - target cell
    #   xgr - x values (4) of guide box
    #   ygr - y values of guide box
    #   radius - radius of offset reference
    #   whigh - high offset of box from radius
    #   wlow - low offset from radius
    #   layer - drawing layer
    #
    ind = [0, 1, 2, 3, 0]
    for j in range(4):
        pa = pointArray()
        il = ind[j]
        ih = ind[j + 1]
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
    #  Draw rings at radius
    #
    #       c - cell pointer
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



def findCell_CK(Cell_name):
#
#   find cell Cell_name - stop program if not found
#
    cnew = l.drawing.findCell(Cell_name)
    if cnew is None:
        print("Cell " + Cell_name + " Not Found")
        exit(1)
    return cnew


def makeAssy(cel, celllist):
    #
    #  make an assembly of multiple cells, all at 0,0
    #  names contained in celllist
    #
    for ind in range(len(celllist)):
        cnew = findCell_CK(clist[ind])
        p = point(0, 0)
        #    print(cnew, " - ", clist[ind])
        cel.addCellref(cnew, p)



def NewCell(CName):
    ec = l.drawing.addCell()
    ec.thisCell.cellName = CName
    ep = ec.thisCell
    return ep


import math


def DrawBump(BP, len, layer):
    angle = 0.785398
    sangle = angle / 2
    grid = 5
    radius = float(len) / math.cos(sangle)
    #    radius = len//math.cos(sangle)
    vertices = [(roundGrid(radius * math.cos(i * angle - sangle), grid), \
                 roundGrid(radius * math.sin(i * angle - sangle), grid)) for i in range(8)]
    vertices.append(vertices[0])
    # print(vertices)
    pts = pointArray()
    x, y = zip(*vertices)
    for i in range(9):
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


def make_2dmesh(cellname, Metal_list, mlayer, via_list):
#   make a metal mesh with layer numbers defined on mlayer, vias on via_list and spacings and numbers in Metal_list
#   Meash is assumed to be a square
    ep = NewCell(cellname + "_mesh")
    for i in range(len(mlayer)):
        line_width = Metal_list[i][1]
        width = Metal_list[i][0]
        space = Metal_list[i][2]
        nline = Metal_list[i][3]
        pitch = line_width + space
        yoff = -((nline - 1) // 2) * pitch
        xoff = width // 2
        origin_via = point(-xoff + line_width // 2, yoff)
        ref_via = point(-xoff + line_width // 2 + pitch, yoff + pitch)
        for j in range(len(via_list)):
            e = ep.addCellrefArray(via_list[j], origin_via, ref_via, nline, nline)
        #   Horizontal lines
        for j in range(nline):
            ep.addBox(-xoff, yoff - line_width // 2, width, line_width, mlayer[i])
            yoff = yoff + pitch
        xoff = -((nline - 1) // 2) * pitch
        yoff = width // 2
    #   Vertical lines
        for j in range(nline):
            ep.addBox(-xoff - line_width // 2, -yoff, line_width, width, mlayer[i])
            xoff = xoff + pitch
    return ep


def erdrawOD(c, xgr, ygr, wlow, whigh, layer, radius, OD, OD_inset):
    # draw implant ring shape with inset active layer (OD)
    erdraw(c, xgr, ygr, wlow, whigh, layer, radius)  # draw the ring
    erdraw(c, xgr, ygr, wlow - OD_inset, whigh - OD_inset, OD, radius)  # draw the inset


def adddrBoxOD(c, xb, yb, xl, yl, rad, layer, OD, OD_inset):
    # draw implant chamfered box shape weith inset active layer (OD)
    adddrBox(c, xb, yb, xl, yl, rad, layer)
    adddrBox(c, xb + OD_inset, yb + OD_inset, xl - 2 * OD_inset, yl - 2 * OD_inset, rad, OD)
    return 1


def roundGrid(x, grid):
    ax = (int(math.ceil((abs(x) - (grid // 2)) / float(grid))) * grid)
    # print(abs(x))
    # print(ax)
    ax = ax*abs(x)/x
    return ax


# import numpy as np
def makeFillCell(FCname, NFill, Pad_Metal, Pad_Width, BSide, Bot_point, NXY, FCpitch):
    #
    # Make a fill cell of 2x2 metal squares
    # FCname - final cell name
    # Nx, Ny number of columns and rows
    # FCpitch - pitch of parent array (fill is x 2 of this)
    # Pad_Metal list of m1-3 pads
    # Bot_point - anchor point of parent array
    rcell = NewCell(FCname)
    fcell = NewCell(FCname + "_square")
    fcell.addBox(-BSide // 2, -BSide // 2, BSide, BSide, OTL)
    for i in range(len(Pad_Metal)):
        spacing = roundGrid((BSide - NFill[i] * 2 * Pad_Width[i]) / (NFill[i] + 1), 50)
        fpitch = 2 * Pad_Width[i] + spacing
        # p1 = point(-fpitch*NFill[i]/2, -fpitch*NFill[i]/2)
        p1 = point(-BSide // 2 + Pad_Width[i] + spacing, -BSide // 2 + Pad_Width[i] + spacing)
        p2 = point(p1.x() + fpitch, p1.y() + fpitch)
        fcell.addCellrefArray(Pad_Metal[i], p1, p2, NFill[i], NFill[i])
    V0 = point(Bot_point.x() - FCpitch // 2, Bot_point.y() - FCpitch // 2)
    V1 = point(V0.x() + FCpitch, V0.y() + FCpitch // 2)

    H0 = point(Bot_point.x() - FCpitch // 2, Bot_point.y() - FCpitch // 2)
    H1 = point(H0.x() + FCpitch // 2, H0.y() + FCpitch)
    rcell.addCellrefArray(fcell, V0, V1, NXY + 1, 2 * NXY + 1)
    rcell.addCellrefArray(fcell, H0, H1, 2 * NXY + 1, NXY + 1)
    return rcell


# Fill = makeACFill(FCname, NFill[i], Pad_List, Pad_Widths, BSide, Pitch[i])
def makeACFill(FCname, NFill, Pad_List, Pad_Widths, BSide, Pitch):
    fcell = NewCell(FCname)
    for i in range(len(NFill)):
        width = Pad_Widths[i] * 2
        space = roundGrid(Pitch / NFill[i], 50)
        for jx in range(NFill[i]):
            for jy in range(NFill[i]):
                xbox = -Pitch // 2 + space // 2 + jx * space
                ybox = -Pitch // 2 + space // 2 + jy * space
                ref = point(xbox, ybox)
                if is_outside(xbox, ybox, width, BSide + 2000):
                    fcell.addCellref(Pad_List[i], ref)
    return fcell


def is_outside(x, y, z, a):
    # Half-sides of each box
    half_z = z / 2
    half_a = a / 2

    # Calculate the smaller box's boundary positions
    min_x = x - half_z
    max_x = x + half_z
    min_y = y - half_z
    max_y = y + half_z

    # Check if any boundary exceeds the larger box's boundaries
    if (min_x > half_a) or (max_x < -half_a) or (min_y > half_a) or (max_y < -half_a):
        return True  # The box is outside
    else:
        return False  # The box is inside


def M1M2M3Fill(Cellname, Layers, Density, SSpace, MWidth, MSpace, FSize):
    dr.setCell(Cellname)
    print(" Filling " + Cellname)
    for i in range(len(Layers)):
        dr.densityFill(Layers[i], Density[i], SSpace[i], MWidth[i], MSpace[i], FSize[i])
    dr.deselectAll()
    return True


def makeMeshContact(Cell_name, DXY, LList, SList, LayerList, CList):
    #
    #	Makse an XY mesh
    #	Cell - mesh cell pointer
    #	DXY - X and Y widths of the contact (modified if necessary top aling with 0,0)
    #	L - X and Y width of metal
    #	S - x and Y space
    #	L - Layer
    #	CList - list of contact cells (Sub-M1, M1-M2, M2-M3)
    #
    Cell = NewCell(Cell_name)
    for m in range(len(LayerList)):
        L = LList[m]
        S = SList[m]
        Layer = LayerList[m]
        P = []
        NXY = []
        OffXY = []
        DXY_mod = []
        for i in range(2):
            P.append(L[i] + S[i])
            NMod = (DXY[i] - L[i]) // P[i]
            NMod = NMod if NMod % 2 == 0 else NMod - 1
            NXY.append(NMod)
            OffXY.append((P[i] * NXY[i] + L[i]) // 2)

        Off = -OffXY[0]
        for i in range(NXY[0] + 1):
            D = P[1] * (NXY[1] // 2) + L[1] // 2
            Cell.addBox(Off, -D, L[0], 2 * D, Layer)
            Off = Off + P[0]

        Off = -OffXY[1]
        for i in range(NXY[1] + 1):
            D = P[0] * (NXY[0] // 2) + L[0] // 2
            Cell.addBox(-D, Off, 2 * D, L[1], Layer)
            Off = Off + P[1]

        V0 = point(-OffXY[0] + L[0] // 2, -OffXY[1] + L[1] // 2)
        V1 = point(-OffXY[0] + L[0] // 2 + P[0], -OffXY[1] + L[0] // 2 + P[1])
        Cell.addCellrefArray(CList[m], V0, V1, NXY[0] + 1, NXY[1] + 1)

    return Cell


def rotate_list(list, k):
    #
    # Circular rotate list by k steps
    # Ensure k is within bounds of the list length
    return list[-k % len(list):] + list[:-k % len(list)]


def Edge_Polygon(Outer, inner):
    #
    #	Provide a list of points for an enclosed polygon
    #	defined by a set of inner and outer points on top
    #	left corner of a square
    #
    Factor = [[1, 1], [-1, 1], [-1, -1], [1, -1]]
    TList = []
    ilist = Outer
    # print(Outer)
    for j in range(len(Factor)):
        LFact = Factor[j]
        for i in range(len(Outer)):
            TLst = [ilist[i][0] * LFact[0], ilist[i][1] * LFact[1]]
            TList.append(TLst)
        #			print(TLst)
        ilist = rotate_list(Outer, j + 1)
    TList.append(Outer[0])

    jlist = inner
    for j in range(len(Factor)):
        LFact = Factor[j]
        for i in range(len(inner)):
            TLst = [jlist[i][0] * LFact[0], jlist[i][1] * LFact[1]]
            TList.append(TLst)
        jlist = rotate_list(inner, j + 1)
    TList.append(inner[0])
    return TList


def make_filled_cell(fill_cell, fcell_name, inner, outer, layer):
    #
    #	Make a fill cell defined by outer and inner top left corner points
    #	fill_cell - cell reference to the fill
    #	fcell_name - - name of result cell
    fcell = NewCell(fcell_name)
    OFrame = Edge_Polygon(outer, inner)
    # print(OFrame)
    dr.activeLayer = layer
    pa = pointArray()
    for i in range(len(OFrame)):
        XF = OFrame[i][0]
        YF = OFrame[i][1]
        pa.append(point(XF, YF))

    pa.append(point(OFrame[0][0], OFrame[0][1]))
    fcell.addPolygon(pa, layer)
    dr.setCell(fcell)
    dr.selectAll()
    dr.fillSelectedShapes(fill_cell, 0)
    return fcell

def makeFrame(cell, corner, width, layer):
#
# make a square frame with edge "corner" and width "width"
    fct = [[-1,1],[1,1],[1,-1],[-1,-1],[-1,1]]
    inout = [corner, corner-width]
    pa = pointArray()
    for i in range(2):
        for j in range(len(fct)):
            pa.attach(inout[i]*fct[j][0] , inout[i]*fct[j][1])
    cell.addPolygon(pointArray(pa), layer)
    return 0

def makeChamferedFrame(cell, corner, width, inset, layer):
    fct = [[-1,1],[1,1],[1,-1],[-1,-1],[-1,1]]
    inmul = [[[0,-1],[-1,0]], [[-1, 0],[0, -1]], [[0, -1],[-1, 0]],
             [[-1, 0], [0, -1]], [[0,-1],[0,-1]]]
    inout = [corner, corner-width]
    pa = pointArray()
    for i in range(2):
        for j in range(len(fct)):
            for facet in range(2):
                inx = inout[i] + inset*inmul[j][facet][0]
                iny =  inout[i] + inset*inmul[j][facet][1]
                pa.attach(inx*fct[j][0] , iny*fct[j][1])
    cell.addPolygon(pointArray(pa), layer)
    return 0

def cArray(pitch, length, active):
#
#	calculate array size and reference points given
#	pitch, length and enclosure size
#
    NstX = active // pitch
    NstY = active // length
    xoff = -((NstX - 1) * pitch) // 2  # bottom left
    yoff = -((NstY - 1) * length) // 2  # bottom left
    pref = point(xoff, yoff)
    poff = point(xoff + pitch, yoff + length)
    return NstX, NstY, pref, poff

import time
import sys

def addZAFill(Assy, zalyr, templyr, fillCell):
#   add ZA Fill

    start_time = time.time()
    dr.setCell(Assy)
    l.booleanTool.boolOnLayer(ZA, 0, templyr, 'A invert', 0, 0, 2)
    Assy.selectLayer(templyr)
    dr.currentCell.sizeAdjustSelect(-7000,0)
    dr.fillSelectedShapes(fillCell, 0)
#   debug
#    l.drawing.saveFile("/Users/lipton/Test_ZA4.gds")
#    sys.exit()
#
#   Fix lower left edge issues - removed after size adjust
#
    Assy.deleteLayer(templyr)
#   delete edge
    mult = [[-1,1], [1,1], [1, -1], [-1,-1]]
    for imul in range(len(mult)):
        dr.point(mult[imul][0]*2894598, mult[imul][1]*2895459)
        dr.point(mult[imul][0]*2955159, mult[imul][1]*2956250)
        dr.cSelect()
        dr.deleteSelect()
    print("ZA FIll time elapsed: {:.2f}s".format(time.time() - start_time))
    return True

def addMxFill(Assy, exclude,  tlayer, fillCell, csize, osize):
#   add Mx Fill
#   Tlayer[1] = outline, Tlayer[2] - fill, mask TLayer[3] - areas to exclude (inverted for mask)
#   csize - dimension to expand metal to eliminate voids
#   osize - dimension to compress frame
    start_time = time.time()
    dr.setCell(Assy)
    for ilyr in range(len(exclude)):
    	l.booleanTool.boolOnLayer(tlayer[2],exclude[ilyr],tlayer[2],"A+B",0,0,2)
    dr.deselectAll()
    Assy.selectLayer(tlayer[2])
    dr.currentCell.sizeAdjustSelect(csize,0)
    l.booleanTool.boolOnLayer(tlayer[0],tlayer[2],tlayer[1],"A-B",0,0,2)
    dr.deselectAll()
    Assy.selectLayer(tlayer[1])
    dr.currentCell.sizeAdjustSelect(osize,0)
    dr.fillSelectedShapes(fillCell, 0)
#   debug
    debug = False
    if debug:
        test = dr.currentCell.cellName
        l.drawing.saveFile("/Users/lipton/" + test + ".gds")
    else:
        Assy.deleteLayer(tlayer[1])
        Assy.deleteLayer(tlayer[2])
    print("Mx FIll time elapsed: {:.2f}s".format(time.time() - start_time))
    return True

def addMZFill(Assy, exclude, tlayer, fillCell, csize, osize):
    #   add Mx, ZA  Fill
    #   Tlayer[1] = outline, Tlayer[2] - fill, mask TLayer[3] - areas to exclude
    #   csize - dimension to expand metal
    #   osuize - dimension to compress frame
    start_time = time.time()
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
    debug = False
    if debug:
        test = dr.currentCell.cellName
        l.drawing.saveFile("/Users/lipton/" + test + ".gds")
    else:
        Assy.deleteLayer(tlayer[1])
        Assy.deleteLayer(tlayer[2])
    print("Mx FIll time elapsed: {:.2f}s".format(time.time() - start_time))
    return True
