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


# def grarray(c, nr, xg, yg, lyr, low, high, rad):
#     #
#     #   Draw guard ring array of nr guard rings
#     #
#     for ind in range(nr):
#         grdraw(c, xg, yg, low[ind], high[ind], lyr, rad[ind])
#
#
# def drcorner(xarr, yarr, npts, layer):
#     #
#     #   draw 4 corners
#     #   points assumed as upper right
#     #
#     dr.clearPoints()
#     dr.activeLayer = layer
#     for ind in range(npts):
#         dr.point(xarr[ind], yarr[ind])
#     dr.polygon()
#     # Upper left
#     dr.clearPoints()
#     dr.activeLayer = layer
#     for ind in range(npts):
#         dr.point(-xarr[ind], yarr[ind])
#     dr.polygon()
#     # lower left
#     dr.clearPoints()
#     dr.activeLayer = layer
#     for ind in range(npts):
#         dr.point(-xarr[ind], -yarr[ind])
#     dr.polygon()
#     # lower right
#     dr.clearPoints()
#     dr.activeLayer = layer
#     for ind in range(npts):
#         dr.point(xarr[ind], -yarr[ind])
#     dr.polygon()

def findCell_CK(Cell_name):
#
#   find cell Cell_name - stop program if not found
#
    cnew = l.drawing.findCell(Cell_name)
    if cnew is None:
        print("Cell " + clist[ind] + " Not Found")
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


# def bpArray(celtgt, bpcell, xur, yur):
#     #
#     #  make an 2x2 array of corner bond pad of bpcell
#     #
#     p = point(xur, yur)
#     e = celtgt.addCellref(bpcell, p)
#     p = point(xur, -yur)
#     e = celtgt.addCellref(bpcell, p)
#     p = point(-xur, -yur)
#     e = celtgt.addCellref(bpcell, p)
#     p = point(-xur, yur)
#     e = celtgt.addCellref(bpcell, p)
#
#
# def padArray(homecell, padcell, npad, xstart, xpitch, ypos, smwidth):
#     #
#     # x array of npad pad cells
#     xpad = xstart
#     for ind in range(npad):
#         p = point(xpad, ypos)
#         e = homecell.addCellref(padcell, p)
#         #    print(p,"-",)
#         #    pa = pointArray()
#         #    pa.attach(xpad-smwidth, 0)
#         #    pa.attach(xpad+smwidth, 0)
#         #    pa.attach(xpad+smwidth, ypos)
#         #    pa.attach(xpad-smwidth, ypos)
#         #    pa.attach(xpad-smwidth, 0)
#         #    e = homecell.addBox(pa, 6)
#         xpad = xpad + xpitch
#
#
# def CSpace(W, S, L):
#     # Calculate spaces for uniform mesh (microns) with line width W and total lenght L
#     nex = (W - S) / (L + S)
#     nfl = int(nex)
#     Spac = (W - nfl * L) / (nfl - 1)
#     return Spac


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


# def Make_M1M2_Mesh(Prefix, Pad_Layers, Pad_Widths, M12_Pitch, Len_XY, CA_Pitch, CA, V1):
#     #
#     #   Make mesh contact with Metal 1 and 2
#     cd = NewCell(Prefix + "_M1M2_Mesh")
#     M1Str = NewCell(Prefix + "_M1_strip")
#     e = M1Str.addBox(-Pad_Widths[0], -Len_XY[1] // 2, 2 * Pad_Widths[0], Len_XY[1], Pad_Layers[0])
#     ypoints = Space_1d(Len_XY[1], CA_Pitch)
#     for k in range(len(ypoints)):
#         e = M1Str.addCellref(CA, point(0, ypoints[k]))
#     M1_strips = Space_1d(Len_XY[0], M12_Pitch[0])
#     for k in range(len(M1_strips)):
#         e = cd.addCellref(M1Str, point(M1_strips[k], 0))
#     #  Metal two connection
#     M2Str = NewCell(Prefix + "_M2_Strip")
#     e = M2Str.addBox(-Len_XY[0] // 2, -Pad_Widths[1], Len_XY[0], 2 * Pad_Widths[1], M2)
#     for k in range(len(M1_strips)):
#         e = M2Str.addCellref(V1, point(M1_strips[k], 0))
#     #
#     #	M2_Strips = Space_1d(Len_XY[1], M12_Pitch[1])
#     #	print(M2_Strips)
#     for k in range(len(ypoints) - 1):
#         e = cd.addCellref(M2Str, point(0, ypoints[k] + CA_Pitch // 2))
#     return cd
#
#
# def Make_M1M2M3_Mesh(Prefix, Pad_Layers, Pad_Widths, M12_Pitch, Len_XY, CA_Pitch, CA, Via_list):
#     #
#     #   Make mesh contact with Metal 1,2 and 3
#     V1 = Via_list[0]
#     V2 = Via_list[1]
#     cd = NewCell(Prefix + "_M1M2M3_Mesh")
#     M1Str = NewCell(Prefix + "_M1_strip")
#     e = M1Str.addBox(-Pad_Widths[0], -Len_XY[1] // 2, 2 * Pad_Widths[0], Len_XY[1], Pad_Layers[0])
#     ypoints = Space_1d(Len_XY[1], CA_Pitch)
#     for k in range(len(ypoints)):
#         e = M1Str.addCellref(CA, point(0, ypoints[k]))
#     M1_strips = Space_1d(Len_XY[0], M12_Pitch[0])
#     for k in range(len(M1_strips)):
#         e = cd.addCellref(M1Str, point(M1_strips[k], 0))
#
#     #  Metal two connection
#     M2Str = NewCell(Prefix + "_M2_Strip")
#     e = M2Str.addBox(-Len_XY[0] // 2, -Pad_Widths[1], Len_XY[0], 2 * Pad_Widths[1], M2)
#     for k in range(len(M1_strips)):
#         e = M2Str.addCellref(V1, point(M1_strips[k], 0))
#     for k in range(len(ypoints) - 1):
#         e = cd.addCellref(M2Str, point(0, ypoints[k] + CA_Pitch // 2))
#
#         # Metal 3 connection
#     M3Str = NewCell(Prefix + "_M3_Strip")
#     e = M3Str.addBox(-Len_XY[0] // 2, -Pad_Widths[2], Len_XY[0], 2 * Pad_Widths[2], M3)
#     for k in range(len(M1_strips)):
#         e = M3Str.addCellref(V2, point(M1_strips[k], 0))
#     for k in range(len(ypoints) - 1):
#         e = cd.addCellref(M3Str, point(0, ypoints[k] + CA_Pitch // 2))
#     return cd


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


# def Place_Pad(cd, Name, SXY_Active, SPitch, Slength):
#     astr = NewCell(Name + "_Arr")
#     #	Strip_Arrays.append(Strip_name[i] + "_Arr")
#     NstX = SXY_Active // SPitch
#     NstY = SXY_Active // Slength
#     xoff = -((NstX - 1) * SPitch) // 2  # bottom left
#     yoff = -((NstY - 1) * Slength) // 2  # bottom left
#     pref = point(xoff, yoff)
#     poff = point(xoff + SPitch, yoff + Slength)
#     e = astr.addCellrefArray(cd, pref, poff, NstX, NstY)
#     return astr


def erdrawOD(c, xgr, ygr, wlow, whigh, layer, radius, OD, OD_inset):
    # draw implant ring shape with inset active layer (OD)
    erdraw(c, xgr, ygr, wlow, whigh, layer, radius)  # draw the ring
    erdraw(c, xgr, ygr, wlow - OD_inset, whigh - OD_inset, OD, radius)  # draw the inset


def adddrBoxOD(c, xb, yb, xl, yl, rad, layer, OD, OD_inset):
    # draw implant chamfered box shape weith inset active layer (OD)
    adddrBox(c, xb, yb, xl, yl, rad, layer)
    adddrBox(c, xb + OD_inset, yb + OD_inset, xl - 2 * OD_inset, yl - 2 * OD_inset, rad, OD)
    return 1


# def make_EdgeArray(Cname, BCell, ECell, NX, NY, DX, DY):
#     #  make an NX by NY array of BCells with the ECell at the 4 edges.
#     # at some point need to add mirroring
#     # NX, NY
#     dcell = NewCell(Cname)
#     xecell = (DX * NX / 2) - DX / 2
#     yecell = (DY * NY / 2) - DY / 2
#     dcell.addCellref(ECell, point(xecell, yecell))
#     dcell.addCellref(ECell, point(-xecell, yecell))
#     dcell.addCellref(ECell, point(-xecell, -yecell))
#     dcell.addCellref(ECell, point(xecell, -yecell))
#     p1 = point(-xecell + DX, -yecell)
#     p2 = point(-xecell + 2 * DX, -yecell + DY)
#     dcell.addCellrefArray(BCell, p1, p2, NX - 2, NY)
#     p1 = point(-xecell, -yecell + DY)
#     p2 = point(-xecell, -yecell + 2 * DY)
#     dcell.addCellrefArray(BCell, p1, p2, 1, NY - 2)
#     p1 = point(xecell, -yecell + DY)
#     p2 = point(xecell, -yecell + 2 * DY)
#     dcell.addCellrefArray(BCell, p1, p2, 1, NY - 2)
#     return dcell


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
def addZAFill(Assy, zalyr, templyr, fillCell):
#   add ZA Fill

    start_time = time.time()
    dr.setCell(Assy)
    l.booleanTool.boolOnLayer(ZA, 0, templyr, "A invert", 0, 0, 2)
    Assy.selectLayer(templyr)
    dr.fillSelectedShapes(fillCell, 0)
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


# -*- codin
import LayoutScript
from LayoutScript import *
from pprint import pprint


l = project.newLayout()  # open new instance of layout class
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
#dr.importFile("/Users/lipton/Dropbox/Programming/TWR_layout/SLAC_layouts/compile_border_v12_r1.gds")
dr.importFile("/Users/lipton/Dropbox/Programming/TWR_layout/SLAC_layouts/compile_border_v13.gds")

CellFill = False # turn on/off the layout editor fill algorithm
InvertOF = False # turn on inversion of OF to NWD
ZA_Fill = False

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

WLayer1 = 201
WLayer2 = 202
WLayer3 = 203

# inset of active
OD_inset = 120
# Number of guard rings
nrings = 5
# Device pixel pitches
Pitch = [50000, 100000]
# Pixel Rows
NPXRow = [100, 50, 50]
# Pixel Columns
NPXCol = [100, 50, 50]
# Cell size
XYCell = 6000000
XYCell_2 = XYCell // 2
# Reticule subcells
NCellx = 4
NCelly = 5

# define outlines
OTL_names = ["OTL_50", "OTL_100"]
PCell_Outline = []
for i in range(len(OTL_names)):
    e = NewCell(OTL_names[i])
    e.addBox(-Pitch[i] // 2, -Pitch[i] // 2, Pitch[i], Pitch[i], OTL)
    PCell_Outline.append(e)

# Active dimension of pixel arays
Length_2 = 2500000  # half length of chip
Lng = 2 * Length_2
length = [Lng]
# Reach through LGAD
# RT LGAD rows, columns
RTRow = 8
RTCol = 8

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

ZG_Width = 750  # M3 to top metal half-width
ZG_Space = 9000

# Make Via list
Via_Widths = [V1_Width, V2_Width, ZG_Width]
Via_Spaces = [V1_Space, V2_Space, ZG_Space]
Via_Layers = [V1, V2, ZG]

# metal ine half-widths
M1_Width_2 = 500
M2_Width_2 = 2000
M3_Width_2 = 2000
M1_Width = 2 * M1_Width_2
M2_Width = 2 * M2_Width_2
M3_Width = 2 * M3_Width_2

#   Default fill parameters
FLayers = [43, 47, 49]
FDensity = [25, 25, 25]
FOffset = [4000, 4000, 4000]
FWidth = [500, 2000, 2000]
FSpace = [2000, 4000, 4000]
FFrame = [10000, 10000, 10000]

# Standard 25% fill cell
FCell_25 = NewCell("Fill_25pct")
FCell_25.addBox(-500, -500, 1000, 1000, M1)
FCell_25.addBox(-500, -500, 1000, 1000, M2)
FCell_25.addBox(-500, -500, 1000, 1000, M3)
FCell_25.addBox(-1000, -1000, 2000, 2000, OTL)

# Cell for ZA fill
ZA_FillCell = NewCell("Za_Fill")
ZA_FillCell.addBox(-2000, -2000, 4000, 4000, ZA)
ZA_FillCell.addBox(-4000, -4000, 8000, 8000, OTL)


Pad_Widths = [M1_Width_2, M2_Width_2, M3_Width_2]
# Max area square pads
Pad_Cells = ["M1_Pad", "M2_Pad", "M3_Pad"]
Via_Cells = ["V1_Via", "V2_Via", "ZG_Via"]

Pad_Layers = [M1, M2, M3]
# Pad_Ref = []
# Make Metal Pad cells
Pad_List = []
for i in range(len(Pad_Cells)):
    m1p = NewCell(Pad_Cells[i])
    m1p.addBox(-Pad_Widths[i], -Pad_Widths[i], 2 * Pad_Widths[i], 2 * Pad_Widths[i], Pad_Layers[i])
    Pad_List.append(m1p)

M1Pad = findCell_CK(Pad_Cells[0])
#  CA contact Cells
#
ml = NewCell("Contact")
ml.addBox(-CA_Width, -CA_Width, 2 * CA_Width, 2 * CA_Width, CA)
#  4x4 contact Cell
CA4x4 = NewCell("Contact_4x4")
e = CA4x4.addCellrefArray(ml, point(-CA_Space // 2, -CA_Space // 2), point(CA_Space // 2, CA_Space // 2), 2, 2)
e = CA4x4.addCellref(M1Pad, point(0, 0))
#  Add ohmic n contact region for CA cell
CA4x4.addBox(-M1_Width_2 + 1, -M1_Width_2 + 1, M1_Width - 2, M1_Width - 2, ND)
CA4x4.addBox(-M1_Width_2 + 1 + OD_inset, -M1_Width_2 + 1 + OD_inset, M1_Width - 2 - 2 * OD_inset,
             M1_Width - 2 - 2 * OD_inset, OD)
#   Distinguish p and n
CAN4x4 = CA4x4

CAP4x4 = NewCell("PContact_4x4")
e = CAP4x4.addCellrefArray(ml, point(-CA_Space // 2, -CA_Space // 2), point(CA_Space // 2, CA_Space // 2), 2, 2)
e = CAP4x4.addCellref(M1Pad, point(0, 0))
#  Add ohmic n contact region for CA cell
CAP4x4.addBox(-M1_Width_2 + 1, -M1_Width_2 + 1, M1_Width - 2, M1_Width - 2, NP)
CAP4x4.addBox(-M1_Width_2 + 1 + OD_inset, -M1_Width_2 + 1 + OD_inset, M1_Width - 2 - 2 * OD_inset,
              M1_Width - 2 - 2 * OD_inset, OD)
#   Distinguish p and n

#  Via Cell list
Via_List = []
for i in range(len(Via_Cells)):
    ml = NewCell(Via_Cells[i])
    ml.addBox(-Via_Widths[i], -Via_Widths[i], 2 * Via_Widths[i], 2 * Via_Widths[i], Via_Layers[i])
    #  4x4 contact Cell
    ml4x4 = NewCell(Via_Cells[i] + "_4x4")
    VS = (Via_Spaces[i] + 2 * Via_Widths[i])
    e = ml4x4.addCellrefArray(ml, point(-VS, -VS), point(0, 0), 3, 3)
    if i == 2:
        Via_List.append(ml)
    else:
        Via_List.append(ml4x4)

BP_M3_Via_80 = NewCell("BP_M3_via80")
BPM3Width = 73000
# 10/24/24 modified for 2 um zg m3 overlap rule
BPM3Length = 2000
ZGM3_ovr = 2000
e = adddrBox(BP_M3_Via_80, -BPM3Width // 2, -BPM3Length // 2, BPM3Width, BPM3Length, 0, ZG)
e = adddrBox(BP_M3_Via_80, -BPM3Width // 2 - ZGM3_ovr, -BPM3Length // 2 - ZGM3_ovr, BPM3Width + 2 * ZGM3_ovr,
             BPM3Length + 2 * ZGM3_ovr, 0, M3)
BP_M3_Via_60 = NewCell("BP_M3_via60")
BPM3Width = 42000
# BPM3Length = 3000
e = adddrBox(BP_M3_Via_60, -BPM3Width // 2, -BPM3Length // 2, BPM3Width, BPM3Length, 0, ZG)
e = adddrBox(BP_M3_Via_60, -BPM3Width // 2 - ZGM3_ovr, -BPM3Length // 2 - ZGM3_ovr, BPM3Width + 2 * ZGM3_ovr,
             BPM3Length + 2 * ZGM3_ovr, 0, M3)

# 80 micron bond pad
BP_80 = NewCell("Bond_Pad_80")
padWidth_80 = 77000
padLength_80 = 205000
OXLength_80 = 190000
OXWidth_80 = 71000
e = adddrBox(BP_80, -OXWidth_80 // 2, -OXLength_80 // 2, OXWidth_80, OXLength_80, 0, ZP)
e = adddrBox(BP_80, -padWidth_80 // 2, -padLength_80 // 2, padWidth_80, padLength_80, 0, ZA)
e = BP_80.addCellref(BP_M3_Via_80, point(0, padLength_80 // 2 - 3000))
e = BP_80.addCellref(BP_M3_Via_80, point(0, -padLength_80 // 2 + 3000))

#  60 micron bond pad
BP_60 = NewCell("Bond_Pad_60")
padWidth_60 = 57000
padLength_60 = 205000
OXLength_60 = 190000
OXWidth_60 = 51000
e = adddrBox(BP_60, -OXWidth_60 // 2, -OXLength_60 // 2, OXWidth_60, OXLength_60, 0, ZP)
e = adddrBox(BP_60, -padWidth_60 // 2, -padLength_60 // 2, padWidth_60, padLength_60, 0, ZA)
e = BP_60.addCellref(BP_M3_Via_60, point(0, padLength_60 // 2 - 3000))
e = BP_60.addCellref(BP_M3_Via_60, point(0, -padLength_60 // 2 + 3000))

#   Bump Pad_Cell
BCell = NewCell("BumpPad")
DrawBump(BCell, 9500, ZA)  # Standard bump pad - check dimensions
DrawBump(BCell, 6500, ZP)

# Metal-Via stack for pads not including CA
M23V23Z = NewCell("M1M2V1V2ZG_Pad")
for Layer in range(3):
    M23V23Z.addCellref(Pad_List[Layer], point(0, 0))
    M23V23Z.addCellref(Via_List[Layer], point(0, 0))
    M23V23Z.addCellref(BCell, point(0, 0))

M3ZG = NewCell("ZG_Pad")
M3ZG.addCellref(Pad_List[2], point(0, 0))
M3ZG.addCellref(Via_List[2], point(0, 0))
M3ZG.addCellref(BCell, point(0, 0))

#   Fill Cells
FillCell = NewCell("Fill_Cell")
FillCell.addBox(-500, -500, 1000, 1000, M1)
FillCell.addBox(-500, -500, 1000, 1000, M2)
FillCell.addBox(-500, -500, 1000, 1000, M3)
FillCell.addBox(-1250, -1250, 2500, 2500, OTL)

#   RT Fill
# Fill cells for edge of array
EFill_Cells = ["RT_Fill", "Str_Fill", "ACPxl_Fill", "DJPxl_Fill"]
Inner = [[-2402000, 2402000]], \
    [[-996000, 1000000], [-1000000, 996000]], \
    [[-2502000, 2498000], [-2498000, 2502000]], \
    [[-2498000, 2492000], [-2492000, 2498000]]
Outer = [[-2506000, 2506000]], \
    [[-1010000, 1002000], [-1002000, 1010000]], \
    [[-2518000, 2508000], [-2508000, 2518000]], \
    [[-2518000, 2508000], [-2508000, 2518000]]
for i in range(len(EFill_Cells)):
    fcell_RT = make_filled_cell("Fill_Cell", EFill_Cells[i], Inner[i], Outer[i], OTL)

#   exclude edge cell
corner = 3000000
width = 117000
exedge_PWD = NewCell("Exclude_edge_PWD")
result = makeFrame(exedge_PWD, corner, width, PWD)

ntype = 3
CA_Contact = findCell_CK("Contact_4x4")

##############################################
#  make DC strxcel
#  def makeStrip(cl, activeLength, sgap, siwidth, scwidth, moffs, LP, LM, LC):
#
#   Active dimensions for 3x3 strip arrays
#
Str_Length = 2000000
Str_Length_2 = Str_Length // 2
empty_cell = NewCell("empty")
offset_3mm = 1500000

#	strip parameters for 6 and 3mm cells
Strip_Pitch = [12500, 50000, 50000, 50000, 50000, 100000, 100000, 100000, 100000,
               100000, 100000, 100000, 100000]

Strip_Length = [50000, Str_Length, Str_Length, Str_Length, Str_Length,
                Str_Length, Str_Length, Str_Length, Str_Length,
                Str_Length, Str_Length, Str_Length, Str_Length]

Strip_name = ["Str125", "Str50_DJ", "Str50_DJNPS", "Str50_AC", "Str50_NOGN",
              "Str100_DJ", "Str100_DJNPS", "Str100_AC", "Str100_NOGN",
              "St100AC_20", "St100AC_40", "St100AC_60", "St100AC_80"]

Strip_PS = [True, True, False, False, True, True, False, False, True,
            False, False, False, False]
Strip_Contacty = ["Strp125CY", "Strp50CY", "Strp50CY", "Strp50CY", "Strp50CY",
                  "Strp100CY", "Strp100CY", "Strp100CY", "Strp100CY",
                  "Strp100CY", "Strp100CY", "Strp100CY", "Strp100CY"]

VM1M2 = findCell_CK("V1_Via_4x4")
VM2M3 = findCell_CK("V2_Via_4x4")
ConCell = [CA_Contact, CA_Contact, CA_Contact, CA_Contact, CA_Contact,
           CA_Contact, CA_Contact, CA_Contact, CA_Contact,
           CA_Contact, CA_Contact, CA_Contact, CA_Contact]
V1Cell = [VM1M2, VM1M2, VM1M2, empty_cell, VM1M2,
           VM1M2, VM1M2, empty_cell, VM1M2,
           empty_cell, empty_cell, empty_cell, empty_cell]
V2Cell = [VM2M3, VM2M3, VM2M3, VM2M3, VM2M3,
           VM2M3, VM2M3, VM2M3, VM2M3,
           VM2M3, VM2M3, VM2M3, VM2M3]

Imp_Lyr = [NPL, NPL, NPL, 0, NPL, NPL, NPL, 0, NPL,
           0, 0, 0, 0]

Border3mm = ["empty", "3mm_with_pads", "3mm_with_pads", "3mm_with_pads", "3mm_with_pads",
             "3mm_with_pads", "3mm_with_pads", "3mm_with_pads", "3mm_with_pads",
             "3mm_with_pads", "3mm_with_pads", "3mm_with_pads", "3mm_with_pads"]

AC_Strip_Index = [3, 7, 9, 10, 11, 12]
AC_Str_Width = [25000, 75000, 20000, 40000, 60000, 80000]

array_3mm = ["Str50_DJ", "Str50_DJNPS", "Str50_AC", "Str50_NOGN"]
offx_3mm = [-offset_3mm, -offset_3mm, offset_3mm, offset_3mm]
offy_3mm = [-offset_3mm, offset_3mm, -offset_3mm, offset_3mm]

Strip_M1_Lines = []
Strip_M2 = []

STInset = 5000  # Strip implant inset from edge
ST_CA_Pitch = 10000  # Strip contact pitch
ST_M1_Pitch = 2000  # strip M1 contact strip Pitch
SXY_Active = 5000000
STXY_Active = 2000000
#   X,Y line width and space for metal layers 1, 2, 3
# LList = [[750, 750], [4000, 4000], [4000, 4000]]
# SList = [[1250, 1250], [6000, 6000], [6000, 6000]]
LList = [[750, 750], [750, 750], [4000, 4000]]
SList = [[1250, 1250], [1250, 1250], [6000, 6000]]


STRound = 2000  # edge rounding parameter for strips
STMSurr = 1000
Strip_imp_width = []
Strip_Arrays = []

# DJ Implants
DJNinset = 15000  # inset of deep N (phos) from active length
DJPinset = 5000  # inset of deep P (Boron) from active length
DJRound = 2000  # edge rounding radius
DJInset = 10000

name = "DJ_P_3mm"
activeLength = STXY_Active
cpad = NewCell(name)
e = adddrBoxOD(cpad, -(activeLength + DJPinset) // 2, -(activeLength + DJPinset) // 2, \
               activeLength + DJPinset, activeLength + DJPinset, DJRound, DJP, OD, OD_inset)

name = "DJ_N_3mm"
dpad = NewCell(name)
e = adddrBoxOD(dpad, -(activeLength + DJNinset) // 2, -(activeLength + DJNinset) // 2, \
               activeLength + DJNinset, activeLength + DJNinset, DJRound + 4000, DJN, OD, OD_inset)

name = "DJ_PN_3mm"
cstr = NewCell(name)
cstr.addCellref(cpad, point(0, 0))
cstr.addCellref(dpad, point(0, 0))

gainsurr = 40000  # gain surround of AC
gainround = 5000  # radius of reference arc for GR edges
gsurr = 10000
#  AC gain and casthode implants

name = "ACP_3mm"
activeLength = STXY_Active
cacp = NewCell(name)
e = adddrBoxOD(cacp, -(activeLength + gsurr) // 2, -(activeLength + gsurr) // 2, \
               activeLength + gsurr, activeLength + gsurr, gainround, PGN, OD, OD_inset)
# AC layer
name = "ACC-3mm"
activeLength = STXY_Active
wbox = (activeLength + gainsurr) // 2
cacc = NewCell(name)
e = adddrBoxOD(cacc, -(activeLength + gainsurr) // 2, -(activeLength + gainsurr) // 2, \
               activeLength + gainsurr, activeLength + gainsurr, gainround, ACN, OD, OD_inset)

# Draw JTE for 3mm
cjte_AC = NewCell("JTE_3mm")
wbox = (activeLength + gainsurr) // 2 - JTERound
wbox = (activeLength + gainsurr) // 2 - gainround - 2500
xpm = wbox
xm = [xpm, -xpm, -xpm, xpm]
ym = [xpm, xpm, -xpm, -xpm]
erdrawOD(cjte_AC, xm, ym, JTEInset - 2000, JTEWidth, JTE, JTERound, OD, OD_inset)
#  Add ND
erdrawOD(cjte_AC, xm, ym, JTEInset - 2000, JTEWidth, ND, JTERound, OD, OD_inset)

sgap = 4000

for i in range(len(Strip_Pitch)):
    Strip_imp_width.append(Strip_Pitch[i] - STInset)
# modify for AC
for i in range(len(AC_Strip_Index)):
    j = AC_Strip_Index[i]
    Strip_imp_width[j] = AC_Str_Width[i]

cellnames_3mm = []
smwidth = M1_Width
for i in range(len(Strip_Pitch)):

    SPitch = Strip_Pitch[i]
    Slength = Strip_Length[i]
    DCStripn = Strip_name[i]
    siwidth = Strip_imp_width[i]
    cd = NewCell(DCStripn)
    lng = Slength - sgap
    wid = siwidth
    # check for OD in AC
    e = adddrBoxOD(cd, -wid // 2, -lng // 2, wid, lng, STRound, Imp_Lyr[i], OD, OD_inset)
    if Imp_Lyr[i] == NPL:
        e = adddrBox(cd, -wid // 2, -lng // 2, wid, lng, STRound, ND)
        #
        #   Add edge field plate
        #
        ypm = lng // 2 - STRound
        xpm = wid // 2 - STRound
        xm = [xpm, -xpm, -xpm, xpm]
        ym = [ypm, ypm, -ypm, -ypm]
        erdraw(cd, xm, ym, 1000, 1000, M2, STRound)
        #   Connect FP to mesh 4 micron wide at y=0
        cd.addBox(-wid // 2, -2000, wid, 4000, M2)
    #
    #   Cell boundry for fill
    e = cd.addBox(-SPitch // 2, -Slength // 2, SPitch, Slength, OTL)
    #   mesh subroutine
    lxy = [siwidth, lng]
    # cd.addCellref(ce, point(0, 0))
    MCell_name = Strip_name[i] + "_Electrode"
    CList = [ConCell[i], V1Cell[i], V2Cell[i]]
    ce = makeMeshContact(MCell_name, lxy, LList, SList, Pad_Layers, CList)
    cd.addCellref(ce, point(0, 0))

    # P-stop for DC coupled
    if not Strip_PS[i]:
        pass
    # Add p-stop
    else:
        rad = 0
        pswidth = 500
        psx = SPitch // 2
        pslen = Slength // 2
        xps = [psx, -psx, -psx, psx]
        yps = [pslen, pslen, -pslen, -pslen]
        erdraw(cd, xps, yps, pswidth, 0, PWD, rad)
    #        cd.addCellref(bpad, point(0, 0))
    #      Add bond pads, make array
    BPinset = 249500
    Row2inset = BPinset + 260000
    # Add fill
    if CellFill:
        e = M1M2M3Fill(DCStripn, FLayers, FDensity, FOffset, FWidth, FSpace, FFrame)

    if SPitch == 100000:
        cd.addCellref(BP_80, point(0, Slength // 2 - BPinset))
        cd.addCellref(BP_80, point(0, -Slength // 2 + BPinset))

        #        bstr = NewCell(Strip_name[i] + "_Arr")
        Strip_Arrays.append(Strip_name[i] + "_Arr3mm")
        name = Strip_name[i] + "_Arr3mm"
        astr = NewCell(name)
        cellnames_3mm.append(name)

        #       e = ST_add_Array(astr, cd, SPitch, Slength, STXY_Active)
        NstX, NstY, pref, poff = cArray(SPitch, Slength, STXY_Active)

        e = astr.addCellrefArray(cd, pref, poff, NstX, NstY)
        # Add fill - not needed - border metal is close
        # fref = l.drawing.findCell("Str_Fill")
        # e = astr.addCellref(fref, point(0, 0))
        # Add border
        bref = findCell_CK(Border3mm[i])
        e = astr.addCellref(bref, point(0, 0))
        # Add DJ implants
        if "DJ" in Strip_name[i]:
            astr.addCellref(cstr, point(0, 0))
        # Add AC Layer
        if "AC" in Strip_name[i]:
            astr.addCellref(cacp, point(0, 0))
            astr.addCellref(cacc, point(0, 0))

    if SPitch == 50000:
        #   make two copies of strip for staggered pads
        dr.setCell(DCStripn)
        dr.selectAll()
        dr.point(-25000, 0)
        dr.move()
        dr.point(50000, 0)
        dr.copy()
        cd.addCellref(BP_60, point(-25000, Slength // 2 - BPinset))
        cd.addCellref(BP_60, point(-25000, -Slength // 2 + BPinset))
        cd.addCellref(BP_60, point(25000, Slength // 2 - Row2inset))
        cd.addCellref(BP_60, point(25000, -Slength // 2 + Row2inset))
        #   build array
        bstr = NewCell(Strip_name[i] + "_Arr")
        Strip_Arrays.append(Strip_name[i] + "_Arr")

        name = Strip_name[i] + "_Arr3mm"
        astr = NewCell(name)
        cellnames_3mm.append(name)

        NstX, NstY, pref, poff = cArray(SPitch*2, Slength, STXY_Active)
        e = astr.addCellrefArray(cd, pref, poff, NstX, NstY)
        # Add fill
        # Add border
        bref = findCell_CK(Border3mm[i])
        e = astr.addCellref(bref, point(0, 0))
        # Add DJ implants
        if "DJ" in Strip_name[i]:
            astr.addCellref(cstr, point(0, 0))
        # Add AC Layer
        if "AC" in Strip_name[i]:
            astr.addCellref(cacp, point(0, 0))
            astr.addCellref(cacc, point(0, 0))

    if SPitch == 12500:
        astr = NewCell(Strip_name[i] + "_Arr")
        Strip_Arrays.append(Strip_name[i] + "_Arr")
        # NstX = SXY_Active // SPitch
        # NstY = SXY_Active // Slength
        # xoff = -((NstX - 1) * SPitch) // 2  # bottom left
        # yoff = -((NstY - 1) * Slength) // 2  # bottom left
        # pref = point(xoff, yoff)
        # poff = point(xoff + SPitch, yoff + Slength)
        NstX, NstY, pref, poff = cArray(SPitch, Slength, SXY_Active)
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
xoff = -((NstX - 1) * SPitch) // 2
yoff = -((NstY - 1) * Slength) // 2
#
#   Staggered pad section - MAKE IT ALL STAGGERED 9/27/24
#
# Lower pad
pref = point(xoff, yoff - ST_CA_Pitch)
poff = point(xoff + SPitch * 2, yoff + Slength - ST_CA_Pitch)
e = Bump_125.addCellrefArray(M23V23Z, pref, poff, NstX // 2, NstY)
# upper pad
pref = point(xoff + SPitch, yoff + ST_CA_Pitch)
poff = point(xoff + SPitch * 3, yoff + Slength + ST_CA_Pitch)
e = Bump_125.addCellrefArray(M23V23Z, pref, poff, NstX // 2, NstY)

st125cell = findCell_CK(Strip_name[indx_125] + "_Arr")
e = st125cell.addCellref(Bump_125, point(0, 0))

##############################################
#
#  DJ Pixel
#
DJName = ["DJ_50", "DJ_100"]
DJAName = ["DJA_50", "DJA_100"]
DJMRad = 22000  # Radius of m1 layer
DJMSur = 5000  # 2 x metal surround of implant

Metal_list = [[[37000, 1000, 2000, 13], [37000, 1000, 2000, 13], [37000, 1000, 2000, 13]],
              [[85000, 1000, 2000, 29], [85000, 1000, 2000, 29], [85000, 1000, 2000, 29]]]
Via_list = [Via_List[0], Via_List[1], CA4x4]
## moff = 80000  # for strips
DJIWid_2 = 500
DJM2_Wid_2 = 1000
Len_XY = [Pitch[0] - DJInset, Pitch[0] - DJInset]
M12_Pitch = [2 * M1_Width, 2 * M2_Width]

for i in range(len(Pitch)):
    cpad_name = DJName[i] + "Base"
    cpad = NewCell(cpad_name)
    lng = Pitch[i] - DJInset
    e = adddrBoxOD(cpad, -lng // 2, -lng // 2, lng, lng, DJRound, NPL, OD, OD_inset)
    e = adddrBox(cpad, -lng // 2, -lng // 2, lng, lng, DJRound, ND)
    crad = 12000
    #  Mesh Contact
    Len_XY = [Pitch[i] - DJInset, Pitch[i] - DJInset]
    #    M12_Pitch = [2 * M1_Width, 2 * M2_Width]
    #    Ref = Make_M1M2_Mesh(DJName[i], Pad_Layers, Pad_Widths, M12_Pitch, Len_XY, ST_CA_Pitch, CA4x4, Via_List[0])
    #    Ref = Make_M1M2M3_Mesh(DJName[i], Pad_Layers, Pad_Widths, M12_Pitch, Len_XY, ST_CA_Pitch, CA4x4, Via_List)
    Ref = make_2dmesh(DJName[i], Metal_list[i], Pad_Layers, Via_list)
    e = cpad.addCellref(Ref, point(0, 0))

    #   add field plates
    xpm = lng // 2 - DJRound  # -*- coding: utf-8 -*-
    xm = [xpm, -xpm, -xpm, xpm]
    ym = [xpm, xpm, -xpm, -xpm]
    erdraw(cpad, xm, ym, DJM2_Wid_2 + 2000, DJM2_Wid_2, M2, DJRound)
    # Add isolation - assume square
    wbox = Pitch[i] // 2
    xpm = wbox
    xm = [xpm, -xpm, -xpm, xpm]
    ym = [xpm, xpm, -xpm, -xpm]
    erdraw(cpad, xm, ym, DJIWid_2, 0, PWD, 0)
    #
    e = cpad.addCellref(PCell_Outline[i], point(0, 0))
    # add fill to cell
    if CellFill:
        e = M1M2M3Fill(cpad_name, FLayers, FDensity, FOffset, FWidth, FSpace, FFrame)
    # cell without central bump
    dpad = NewCell(DJName[i] + "_NoBump")
    e = dpad.addCellref(cpad, point(0, 0))
    # cell with bump
    epad = NewCell(DJName[i])
    e = epad.addCellref(cpad, point(0, 0))
    e = epad.addCellref(M3ZG, point(0, 0))  # bump in cell  center
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
e = adddrBoxOD(cpad, -djcorner // 2, -djcorner // 2, \
               djcorner, djcorner, DJRound + dinset, DJP, OD, OD_inset)

djcorner = activeLength + DJNinset
e = adddrBoxOD(cpad, -djcorner // 2, -djcorner // 2, \
               djcorner, djcorner, DJRound + dinset + 3500, DJN, OD, OD_inset)

##############################################
#
# Reach through pixel and array
#
##############################################
#
# Reach through pixel and array
#
# RT inter pixel gap
RTGap = [80000, 16000]
RTRound = [5000, 4000]  # corner rounding radius
# RTOXRad = 6000
# RTMRad = RTLenx // 4

# Reach through LGAD
# RT LGAD rows, columns
RTRow = [8, 50]
RTCol = [8, 50]
# RT LGAD Pitch
# RTPitch = [600000, 600000]  # need to replace x, y (below) by an index

RTPitchx = [625000, 100000]
RTPitchy = [625000, 100000]
RTname = ["RTPixel", "RT_100"]
## RTMSurr = -10000 # overlap with cathode
RTPad_inset = 25000
RTGin = 8000  # gain layer inset
RTCname = ["RTArray", "RTP100_Array"]
# RT_M1_Lines = []
# RT_M1_Pitch = 2 * M1_Width
# RT_M2_Pitch = 2 * M2_Width

Via_list = [Via_List[0], Via_List[1], CA4x4]
Metal_list = [[[499000, 1000, 2000, 167], [499000, 1000, 2000, 167], [499000, 1000, 2000, 167]], \
              [[61000, 1000, 2000, 21], [61000, 1000, 2000, 21], [61000, 1000, 2000, 21]]]

# RT JTE parameters
RT_JTERound = [7000, 7000]
RT_JTEInset = [2000, 4000]
#RT_JTEWidth = [28000, 6000]
RT_JTEWidth = [38000, 6000]

for i in range(len(RTname)):
    # RT Active length
    RTLenx = RTPitchx[i] - RTGap[i]  # active lengths
    RTLeny = RTPitchy[i] - RTGap[i]
    rtpad = NewCell(RTname[i])
    RTPx_2 = RTPitchx[i] // 2
    RTPy_2 = RTPitchy[i] // 2
    rtpad.addBox(-RTPx_2, -RTPy_2, RTPitchx[i], RTPitchy[i], OTL)
    lng = RTPitchx[i] - RTGap[i]
    e = adddrBoxOD(rtpad, -lng // 2, -lng // 2, lng, lng, RTRound[i], NPL, OD, OD_inset)
    e = adddrBox(rtpad, -lng // 2, -lng // 2, lng, lng, RTRound[i], ND)
    e = adddrBoxOD(rtpad, -lng // 2 + RTGin, -lng // 2 + RTGin, lng - 2 * RTGin, lng + -2 * RTGin, RTRound[i]-2000, PGN, OD,
                   OD_inset)
    # add JTE
    wbox = lng // 2 - RTGin
    xpm = wbox
    xm = [xpm, -xpm, -xpm, xpm]
    ym = [xpm, xpm, -xpm, -xpm]
    # print(str(i) + "  - " + str(xpm))
#    erdrawOD(rtpad, xm, ym, RT_JTEInset[i], RT_JTEWidth[i], JTE, RT_JTERound[i], OD, OD_inset)
    JTEEdge  = xpm+RT_JTERound[i]+RT_JTEWidth[i]
    JTEWid = RT_JTEWidth[i]+RT_JTEInset[i]
    e = makeChamferedFrame(rtpad, JTEEdge, JTEWid, 3000, JTE)
    e = makeChamferedFrame(rtpad, JTEEdge-OD_inset, JTEWid-2*OD_inset, 3000, OD)
    RT_Pad_len = lng - 20000

    #
    #   Add contact mesh
    #
    Len_XY = [RT_Pad_len, RT_Pad_len]

    Ref = make_2dmesh(RTname[i], Metal_list[i], Pad_Layers, Via_list)
    rtpad.addCellref(Ref, point(0, 0))
    wbox = RTPitchx[i] // 2 - RTRound[i] - PIWid_2

    # Pad at center
    if RTname[i] == "RTPixel":
        rtpad.addCellref(BP_80, point(0, 0))
    else:
        rtpad.addCellref(BCell, point(0, 0))
        rtpad.addCellref(M3ZG, point(0, 0))
    # Pixel isolation

    e = makeFrame(rtpad, RTPx_2, 500, PWD)

    #   add field plates
    FP_round = RT_JTERound[i] + RT_JTEWidth[i]
    FP_metalw = [500, 1000, 1000]  # metal half width
    for k in range(len(Pad_Layers)):
#        erdraw(rtpad, xm, ym, FP_metalw[k], FP_metalw[k], Pad_Layers[k], FP_round)
        e = makeChamferedFrame(rtpad, JTEEdge+FP_metalw[k]//2, FP_metalw[k], 3000, Pad_Layers[k])
        # Contact electrode
        left_right = JTEEdge
        rtpad.addBox(-left_right, -FP_metalw[k], 2 * left_right, 2 * FP_metalw[k], Pad_Layers[k])
    # add fill
    #
    dr.setCell(RTname[i])
    # dr.activeLayer = PGN  # was ND
    rtpad.selectLayer(JTE)
    l.booleanTool.setA()
    rtpad.selectLayer(M2)
    l.booleanTool.setB()
    l.booleanTool.aMinusB()
    l.drawing.currentCell.moveToLayerSelect(WLayer1)
    dr.fillSelectedShapes("Fill_25pct", 0)
    #
    # Array cell
    #
    RTpad = NewCell(RTname[i] + "_Arry")
    xoff = -(((RTRow[i]) - 1) * RTPitchx[i]) // 2
    yoff = -(((RTCol[i]) - 1) * RTPitchy[i]) // 2
    pref = point(xoff, yoff)
    poff = point(xoff + RTPitchx[i], yoff + RTPitchy[i])
    e = RTpad.addCellrefArray(rtpad, pref, poff, RTRow[i], RTCol[i])

##############################################
# AC Pads
# list of electrode metal [50um][100um] [electrode width, M1 width, M2-3 width, rows]
# note that we could use makeMeshContact for this now 11/10/24
ACPitch = Pitch
Metal_list = [[[33000, 1000, 1000, 17], [33000, 1000, 1000, 17], [33000, 1000, 1000, 17]],
              [[61000, 1000, 1000, 31], [61000, 1000, 1000, 31], [61000, 1000, 1000, 31]],
              [[33000, 1000, 1000, 17], [33000, 1000, 1000, 17], [33000, 1000, 1000, 17]]]
# Metal layers (M1, M2, M3)
mlayer = [43, 47, 49]
ACPitch.append(100000)
name = ["ACPad_50", "ACPad_100", "ACPad_100D33"]
aname = ["AC_Array_50", "AC_Array_100", "AC_Array_100D33"]
# Vias [ no CA, V1, V2]
#   Only vias between M1-M2 and M2-M3
# 1/2/25 modify to couple M1 to
Via_list = [CA4x4, empty_cell, Via_List[1]]
#  fill rows, columns for 50, 100 micron
NFill = [[25, 6, 6], [50, 12, 12], [50, 12, 12]]
# Nfill = [8, 4, 4]
for i in range(len(name)):
    BSide = Metal_list[i][0][0]
    # print(BSide)
    cd = NewCell(name[i])
    #    Len_XY = [ACMetx[i], ACMety[i]]

    Ref = make_2dmesh(name[i], Metal_list[i], Pad_Layers, Via_list)
    cd.addCellref(Ref, point(0, 0))

    # Outline
    x = -ACPitch[i] // 2
    y = -ACPitch[i] // 2
    lx = ACPitch[i]
    ly = ACPitch[i]
    e = cd.addBox(x, y, lx, ly, OTL)

    # add fill
    FCname = name[i] + "_Fill"
    Fill = makeACFill(FCname, NFill[i], Pad_List, Pad_Widths, BSide, ACPitch[i])
    e = cd.addCellref(Fill, point(0, 0))

    apad = NewCell(aname[i])
    xoff = -(((NPXRow[i]) - 1) * ACPitch[i]) // 2
    yoff = -(((NPXCol[i]) - 1) * ACPitch[i]) // 2
    pref = point(xoff, yoff)
    poff = point(xoff + Pitch[i], yoff + ACPitch[i])
    e = apad.addCellrefArray(cd, pref, poff, NPXRow[i], NPXCol[i])
    e = cd.addCellref(M3ZG, point(0, 0))

# Gain Layer
# name = "Gain_Layer"
for i in range(len(length)):
    activeLength = length[i]
    cpad = NewCell("Gain_Layer")
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
wbox = (activeLength + gainsurr) // 2 - JTERound
xpm = wbox
xm = [xpm, -xpm, -xpm, xpm]
ym = [xpm, xpm, -xpm, -xpm]
erdrawOD(cjte, xm, ym, JTEInset, JTEWidth, JTE, JTERound, OD, OD_inset)
erdrawOD(cjte, xm, ym, JTEInset, JTEWidth, ND, JTERound, OD, OD_inset)

#
#	Psubstrate Contact
#
PSWid = 50000  # width of contact in cell
#
#	Cell Outline
#
cotl = "Cell_Outline"
Outline = NewCell(cotl)
Outline.addBox(-XYCell_2, -XYCell_2, XYCell, XYCell, OTL)
#
#  Make assemblies
#
CellList = []
#
#	 AC LGAD Cells
# add third option for small electrode cell
#
# ACL_list = ["Assy_AC_", "Assy_AC_"]
Border_List = ["6mm_50um_pitch_bumps", "6mm_100um_pitch_bumps",  "6mm_100um_pitch_bumps"]
AC_Type = ["50", "100"]
for i in range(len(AC_Type)):
    cname = "Assy_AC_" + AC_Type[i]
    Assy_AC = NewCell(cname)
    ACPName = "AC_Array_" + AC_Type[i]
    #	print(ACPName)
    # clist = ["AC_Layer", ACPName, cotl, "JTE", "Gain_Layer", "ACPxl_Fill", Border_List[i]]
    clist = ["AC_Layer", ACPName, cotl, "JTE", "Gain_Layer", Border_List[i], "Exclude_edge_PWD"]
    makeAssy(Assy_AC, clist)
#   add ZA Fill
    if (ZA_Fill):
        addZAFill(Assy_AC, ZA, WLayer2, ZA_FillCell)
        print(cname + " ZA Fill")
    CellList.append(cname)
#
#	 DJ LGAD Cells
#
DJ_Type = ["50", "100"]
for i in range(len(DJ_Type)):
    cname = "AssyDJ_" + DJ_Type[i]
    Assy_DJ = NewCell(cname)
    DJPName = "DJA_" + DJ_Type[i]
    # clist = [DJPName, cotl, "JTE", "DJ_PN", "DJPxl_Fill", Border_List[i]]
    clist = [DJPName, cotl, "JTE", "DJ_PN", Border_List[i], "Exclude_edge_PWD"]
    makeAssy(Assy_DJ, clist)
    #   add ZA Fill
    if (ZA_Fill):
        addZAFill(Assy_DJ, ZA, WLayer2, ZA_FillCell)
        print(cname + " ZA Fill")
    CellList.append(cname)
    #
    #

#
#	 no gain Cells
#
NG_Type = ["50_NG", "100_NG"]
for i in range(len(NG_Type)):
    cname = "AssyPX_" + NG_Type[i]
    Assy_NGPX = NewCell(cname)
    NGName = "DJA_" + DJ_Type[i]
    # clist = [DJPName, cotl, "JTE", "DJ_PN", "DJPxl_Fill", Border_List[i]]
    clist = [NGName, cotl, "JTE", Border_List[i], "Exclude_edge_PWD"] # leave out DJ implants
    makeAssy(Assy_NGPX, clist)
    #   add ZA Fill
    if (ZA_Fill):
        addZAFill(Assy_NGPX, ZA, WLayer2, ZA_FillCell)
        print(cname + " ZA Fill")
    #
    CellList.append(cname)

for i in range(1):
    cname = "Assy_" + Strip_Arrays[i]
    # 6mm border for 12.5 micon pitch
    if Strip_Arrays[i].startswith("Str125"):
        border = "6mm_with_pads"
        #    djimp = "DJ_PN"
        #    else:
        #        border = "empty"
        djimp = "empty"
    STX_Ass = NewCell(cname)
    clist = (Strip_Arrays[i], cotl, border, "Exclude_edge_PWD")
    makeAssy(STX_Ass, clist)
    CellList.append(cname)
######## no 50?
No_Pad_border = ["6mm_100um_pitch_bumps_DJ_ASIL", "6mm_100um_pitch_bumps_DJ_ASIL"]
# Just one of these for the moment
#for i in range(len(DJ_Type)):
cname = "AssyDJ_NB_" + DJ_Type[1]
Assy_DJ_NB = NewCell(cname)
DJPName = "DJA_" + DJ_Type[1] + "_noBump"
# clist = [DJPName, cotl, "JTE", "DJ_PN", "DJPxl_Fill", No_Pad_border[i]]
clist = [DJPName, cotl, "JTE", "DJ_PN", No_Pad_border[1], "Exclude_edge_PWD"]
#
makeAssy(Assy_DJ_NB, clist)
if (ZA_Fill):
    addZAFill(Assy_DJ_NB, ZA, WLayer2, ZA_FillCell)
    print(cname + " ZA Fill")
CellList.append(cname)
#
#    Assemble 3mm strip cells
cells_3mm = ["Assy_str_3mm_50", "Assy_str_3mm_100", "Assy_ACStr_Elec"]
cell_3mm_50 = NewCell("Assy_str_3mm_50")
cell_3mm_100 = NewCell("Assy_str_3mm_100")
cell_3mm_AC = NewCell(cells_3mm[2])
cell3mm_list = [cell_3mm_50, cell_3mm_100, cell_3mm_AC]
for i in range(2):
    cname = cells_3mm[i]
    c_3mm = cell3mm_list[i]
    for j in range(4):
        ccell = findCell_CK(cellnames_3mm[4 * i + j])
        # print(cellnames_3mm[4 * i + j], i, j)
        px = offx_3mm[j]
        py = offy_3mm[j]
        c_3mm.addCellref(ccell, point(px, py))
    c_3mm.addCellref(Outline, point(0,0))
    # exclude edge from NWD
    c_3mm.addCellref(exedge_PWD, point(0, 0))

    #   add ZA Fill
    if (ZA_Fill):
        addZAFill(c_3mm, ZA, WLayer2, ZA_FillCell)
        print(cname + " ZA Fill")
    CellList.append(cname)
    CellList.append(cname)
#
#   Add AC pitch variants
indx = 2
cname = cells_3mm[2]
c_3mm = cell3mm_list[2]
for j in range(4):
    ccell = findCell_CK(cellnames_3mm[4 * indx + j])
    #    print(cellnames_3mm[4 * indx + j], i, j)
    px = offx_3mm[j]
    py = offy_3mm[j]
    c_3mm.addCellref(ccell, point(px, py))
c_3mm.addCellref(Outline, point(0,0))
# exclude edge from NWD
c_3mm.addCellref(exedge_PWD, point(0, 0))
if (ZA_Fill):
    addZAFill(c_3mm, ZA, WLayer2, ZA_FillCell)
    print(cname + " ZA Fill")
CellList.append(cname)

RTF_list = ["empty", "empty"]
for i in range(len(RTCname)):
    cnam = "Assy_" + RTname[i]
    RTAss = NewCell(cnam)
#    clist = (RTname[i] + "_Arry", "6mm_with_pads", "RT_Fill")
    clist = (RTname[i] + "_Arry", cotl, "6mm_with_pads", RTF_list[i], "Exclude_edge_PWD")
    makeAssy(RTAss, clist)
    if (ZA_Fill):
        addZAFill(RTAss, ZA, WLayer2, ZA_FillCell)
        print(cnam + " ZA Fill")
    CellList.append(cnam)

#   invert OF/PWD for NWD layer
if InvertOF:
    for i in range(len(CellList)):
        #  Try to use OF to define NWD
        print( " Form NWD for " + CellList[i])
        dr.setCell(CellList[i])
        l.booleanTool.boolOnLayer(PWD, 0, NWD, "A invert")
    dr.deleteLayer(PWD)
#   remove PWD

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
    cnew = findCell_CK(CellList[i])
    p = point(Rx, Ry)
    Ret.addCellref(cnew, p)


print(CellList)
# import os
from pathlib import Path

home_directory = Path.home()

gdsversion = "25_r0"
gdsfile = str(home_directory) +  "/Dropbox/Programming/TWR_layout/TWR_" + gdsversion + ".gds"
print(gdsfile)

l.drawing.saveFile(gdsfile)

print("Writing Version " + gdsversion +  " Python script completed")
