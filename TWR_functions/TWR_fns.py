from LayoutScript import *
def BoxDraw(c,xgr,ygr,radius,whigh,wlow, layer):
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

def NewCell( CName):
    ec = dr.addCell()
    ec.thisCell.cellName = CName
    ep = ec.thisCell
    return ep

import math
def grdraw(c, xgr, ygr, wlow, whigh, layer, radius):
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
    angle = 0
    for ind in range(len(xgr)):
        horiz = bool(True)
        x = xgr[ind]
        y = ygr[ind]
        #    print(str(angle), str(layer), str(radius))
        c.addPolygonArc(point(x, y), int(radius - wlow), int(radius + whigh), angle, angle + 90, layer)
        angle = angle + 90
        horiz = not horiz
    #
    #      connect arcs
    #
    rval = BoxDraw(c, xgr, ygr, radius, whigh, wlow, layer)


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
    angle = 0
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

#
def makeAssy(cel, celllist):
    #
    #  make an assembly of multiple cells, all at 0,0
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


def DrawBump(radius, name, layer):
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
    BP = NewCell(name)
    BP.addPolygon(pts, layer)
    return BP


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