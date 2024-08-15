# -*- coding: utf-8 -*-
#
from pprint import pprint

def ecdraw(c, xgr, ygr, wlow, whigh, layer, radius):
    #
    #  Draw edge contacts for x y location tests
    #
    angle = 0
    whigh = wlow
    for ind in range(len(xgr)):
        horiz = bool(True)
        x = xgr[ind]
        y = ygr[ind]
        #    print(str(angle), str(layer), str(radius))
        #    c.addPolygonArc(point(x, y), int(radius-wlow), int(radius+whigh), angle, angle+90, layer)
        angle = angle + 90
        horiz = not horiz
        #
        #      connect arcs
        #
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
        c.addCircle(layer, point(x, y + radius), point(x - wlow, y + radius), 0)
        c.addCircle(layer, point(xn, yn + radius), point(xn - wlow, yn + radius), 0)
        rval = BoxDraw(c, xgr, ygr, radius, whigh, wlow, layer)

def DrawGuard(c, activeLength, LImp, LCon, LMet, LOxide, LJTE, ngr):
    #
    #  Draw an assembly of guard rings
    # X and Y center of guard ring assembly
    arccgr = activeLength // 2 + 120000
    halfLength = activeLength // 2
    arccx = arccgr
    arccy = arccgr
    # passivation opening parameters
    dlong = 400000
    dshort = 100000
    # extent of implants
    implo = [5000] * ngr
    imphi = [5000] * ngr
    # extent of contacts
    conlo = [3000] * ngr
    conhi = [3000] * ngr
    # JTE
    jtelo = [8000] * ngr
    jtehi = [8000] * ngr
    # metal overhang
    metallo = [11000, 11000, 11000, 11000, 34000]
    metalhi = [15000, 15000, 15000, 18000, 63000]
    # ring radii from xy center
    ringrad = [84000, 115000, 147000, 178000, 244000]
    xguard = [arccx, -arccx, -arccx, arccx]
    yguard = [arccy, arccy, -arccy, -arccy]
    #  layout.drawing.activeLayer=4
    # draw implants
    layer = LImp
    grarray(c, nrings, xguard, yguard, layer, implo, imphi, ringrad)
    # draw jte
    layer = LJTE
    grarray(c, nrings, xguard, yguard, layer, jtelo, jtehi, ringrad)
    # draw contacts
    layer = LCon
    grarray(c, nrings, xguard, yguard, layer, conlo, conhi, ringrad)
    # draw metal
    layer = LMet
    grarray(c, nrings, xguard, yguard, layer, metallo, metalhi, ringrad)
    # Draw inner guard - outer section
    ogrc = 28000
    grdraw(c, xguard, yguard, ogrc, ogrc - 10000, LImp, ogrc)
    grdraw(c, xguard, yguard, 3000, 3000, LCon, ogrc)
    grdraw(c, xguard, yguard, ogrc, ogrc + 10000, LMet, ogrc)
    grdraw(c, xguard, yguard, ogrc + 5000, ogrc - 5000, LJTE, ogrc)
    # inner section
    darc = (arccx - halfLength) // 2
    arccx = activeLength // 2
    arccy = activeLength // 2
    xguard = [arccx, -arccx, -arccx, arccx]
    yguard = [arccy, arccy, -arccy, -arccy]
    grdraw(c, xguard, yguard, darc - 75000, darc, LImp, darc)
    grdraw(c, xguard, yguard, darc - 70000, darc + 5000, LJTE, darc)
    grdraw(c, xguard, yguard, darc - 60000, darc, LMet, darc)
    grdraw(c, xguard, yguard, darc - 85000, darc + 15000, LOxide, darc)
    # fill in arc
    xcorns = [arccgr, arccgr, arccx]
    ycorns = [arccgr, arccy, arccgr]
    drcorner(xcorns, ycorns, 3, LImp)
    drcorner(xcorns, ycorns, 3, LMet)
    drcorner(xcorns, ycorns, 3, LJTE)


#  Edge ring warning - hard wire n type
#  grdraw(c, xguard, yguard, 50000, 200000, Ledg, 800000)
#  grdraw(c, xguard, yguard, 5000, 5000, LCon, 800000)
#  grdraw(c, xguard, yguard, 55000, 200000, LMet, 800000)
#  grdraw(c, xguard, yguard, 100000,100000, LOxide, 920000)

def InnerGuard(c, activeLength, LImp, LCon, LMet, LOxide, Ledg):
    # Draw inner guard - outer section
    arccgr = activeLength // 2 + 120000
    halfLength = activeLength // 2
    arccx = arccgr
    arccy = arccgr
    xguard = [arccx, -arccx, -arccx, arccx]
    yguard = [arccy, arccy, -arccy, -arccy]
    ogrc = 28000
    grdraw(c, xguard, yguard, ogrc, ogrc, LImp, ogrc)
    grdraw(c, xguard, yguard, 3000, 3000, LCon, ogrc)
    grdraw(c, xguard, yguard, ogrc, ogrc + 10000, LMet, ogrc)
    # inner section
    darc = (arccx - halfLength) // 2
    arccx = activeLength // 2
    arccy = activeLength // 2
    xguard = [arccx, -arccx, -arccx, arccx]
    yguard = [arccy, arccy, -arccy, -arccy]
    grdraw(c, xguard, yguard, darc - 70000, darc, LImp, darc)
    grdraw(c, xguard, yguard, darc - 60000, darc, LMet, darc)
    grdraw(c, xguard, yguard, darc - 80000, darc + 20000, LOxide, darc)
    # fill in arc
    xcorns = [arccgr, arccgr, arccx]
    ycorns = [arccgr, arccy, arccgr]
    drcorner(xcorns, ycorns, 3, LImp)
    drcorner(xcorns, ycorns, 3, LMet)


def EdgeRing(c, activeLength, LImp, LCon, LMet, LOxide, Ledg):
    # Draw inner guard - outer section
    arccgr = activeLength // 2 + 120000
    halfLength = activeLength // 2
    arccx = activeLength // 2
    arccy = activeLength // 2
    darc = (arccx - halfLength) // 2
    xguard = [arccx, -arccx, -arccx, arccx]
    yguard = [arccy, arccy, -arccy, -arccy]
    # fill in arc
    xcorns = [arccgr, arccgr, arccx]
    ycorns = [arccgr, arccy, arccgr]
    #  Edge ring warning
    grdraw(c, xguard, yguard, 50000, 200000, Ledg, 800000)
    grdraw(c, xguard, yguard, 5000, 5000, LCon, 800000)
    grdraw(c, xguard, yguard, 60000, 220000, LMet, 800000)
    grdraw(c, xguard, yguard, 45000, 190000, LOxide, 800000)


def TEdgeRing(c, activeLength, LImp, LCon, LMet, LOxide, Ledg):
    # Draw inner guard - outer section
    arccgr = activeLength // 2 + 120000
    halfLength = activeLength // 2
    arccx = activeLength // 2
    arccy = activeLength // 2
    darc = (arccx - halfLength) // 2
    xguard = [arccx, -arccx, -arccx, arccx]
    yguard = [arccy, arccy, -arccy, -arccy]
    # fill in arc
    xcorns = [arccgr, arccgr, arccx]
    ycorns = [arccgr, arccy, arccgr]
    #  Edge ring warning
    grdraw(c, xguard, yguard, 225000, 360000, Ledg, 400000)
    grdraw(c, xguard, yguard, 5000, 5000, LCon, 400000)
    grdraw(c, xguard, yguard, 340000, 420000, LMet, 400000)


#  grdraw(c, xguard, yguard, 45000, 190000, LOxide, 400000)

def makeStrip(cl, activeLength, sgap, siwidth, scwidth, moffs, LP, LM, LC):
    #
    #   make a Strip cell using rounded edges
    #
    slength = activeLength - 2 * sgap - 2 * siwidth + 120000
    mlength = slength - moffs
    smwidth = siwidth + 2000
    if activeLength < 2000000:
        clength = 160000
    else:
        clength = mlength // 2
    # add implant
    cl.addBox(-siwidth // 2, -slength // 2, siwidth, slength, LP)
    cl.addCircle(LP, point(0, -slength // 2), point(-siwidth // 2, -slength // 2), 0)
    cl.addCircle(LP, point(0, slength // 2), point(-siwidth // 2, slength // 2), 0)
    # add metal
    cl.addBox(-smwidth // 2, -mlength // 2, smwidth, mlength, LM)
    cl.addCircle(LM, point(0, -mlength // 2), point(-smwidth // 2, -mlength // 2), 0)
    cl.addCircle(LM, point(0, mlength // 2), point(-smwidth // 2, mlength // 2), 0)
    # add contact
    cl.addBox(-scwidth // 2, -clength // 2, scwidth, clength, LC)
    cl.addCircle(LC, point(0, -clength // 2), point(-scwidth // 2, -clength // 2), 0)
    cl.addCircle(LC, point(0, clength // 2), point(-scwidth // 2, clength // 2), 0)


def makeACStrip(cl, activeLength, sgap, mwidth, LM):
    #
    #	sgap - half gap to active length
    #	mwidth - metal width
    #	LM - metal layer
    #		see about pad locations ...
    #
    slength = activeLength - 2 * sgap + 110000
    mlength = slength
    # add metal
    cl.addBox(-smwidth // 2, -mlength // 2, mwidth, mlength, LM)
    cl.addCircle(LM, point(0, -mlength // 2), point(-mwidth // 2, -mlength // 2), 0)
    cl.addCircle(LM, point(0, mlength // 2), point(-smwidth // 2, mlength // 2), 0)


