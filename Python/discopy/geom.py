import numpy as np


def getXYZ(x1, x2, x3, opts, pars):

    x1, x2, x3 = np.broadcast_arrays(x1, x2, x3)

    if opts['GEOMETRY'] == "cylindrical":
        return getXYZcyl(x1, x2, x3, pars)
    elif opts['GEOMETRY'] == "spherical":
        return getXYZsph(x1, x2, x3, pars)
    else:
        return getXYZcart(x1, x2, x3, pars)


def getCellCentroids(dat, opts, pars):
    x1jph = dat[0]
    x2iph = dat[3]
    x3kph = dat[1]
    nx2 = dat[7]
    index = dat[8]
    x2max = pars['Phi_Max']
    N = x2iph.shape[0]

    nx1 = nx2.shape[1]
    nx3 = nx2.shape[0]

    x11d = getCentroid(x1jph[:-1], x1jph[1:], 1, opts)
    x31d = getCentroid(x3kph[:-1], x3kph[1:], 3, opts)

    x1 = np.empty(N)
    x2 = np.empty(N)
    x3 = np.empty(N)

    for k in range(nx3):
        for j in range(nx1):
            a = index[k, j]
            b = index[k, j] + nx2[k, j]

            x1[a:b] = x11d[j]
            x3[a:b] = x31d[k]

            x2p = x2iph[a:b].copy()
            x2p[x2p < 0.0] += x2max
            x2p[x2p >= x2max] -= x2max
            x2m = np.roll(x2p, 1)
            x2m[x2m > x2p] -= x2max

            x2[a:b] = getCentroid(x2m, x2p, 2, opts)

    return x1, x2, x3


def getCentroid(xm, xp, dim, opts):
    if opts['GEOMETRY'] == "cylindrical":
        return getCentroidcyl(xm, xp, dim)
    elif opts['GEOMETRY'] == "spherical":
        return getCentroidsph(xm, xp, dim)
    else:
        return getCentroidcart(xm, xp, dim)


def getDV(dat, opts, pars):

    x1jph = dat[0]
    x2iph = dat[3]
    x3kph = dat[1]
    nx2 = dat[7]
    index = dat[8]

    if opts['GEOMETRY'] == "cylindrical":
        return getDVcyl(x1jph, x2iph, x3kph, nx2, index, pars)
    elif opts['GEOMETRY'] == "spherical":
        return getDVsph(x1jph, x2iph, x3kph, nx2, index, pars)
    else:
        return getDVcart(x1jph, x2iph, x3kph, nx2, index, pars)


def getVXYZ(x1, x2, x3, v1, v2, v3, opts):

    h1, h2, h3 = getScaleFactors(x1, x2, x3, opts)
    return getVecXYZ(x1, x2, x3, h1*v1, h2*v2, h3*v3, opts)


def getBXYZ(x1, x2, x3, B1, B2, B3, opts):

    return getVecXYZ(x1, x2, x3, B1, B2, B3, opts)


def getVecXYZ(x1, x2, x3, v1, v2, v3, opts):

    if opts['GEOMETRY'] == "cylindrical":
        return getVecXYZcyl(x1, x2, x3, v1, v2, v3)
    elif opts['GEOMETRY'] == "spherical":
        return getVecXYZsph(x1, x2, x3, v1, v2, v3)
    else:
        return getVecXYZcart(x1, x2, x3, v1, v2, v3)


def getScaleFactors(x1, x2, x3, opts):
    if opts['GEOMETRY'] == "cylindrical":
        return getScaleFactorscyl(np.atleast_1d(x1), np.atleast_1d(x2),
                                  np.atleast_1d(x3))
    elif opts['GEOMETRY'] == "spherical":
        return getScaleFactorssph(np.atleast_1d(x1), np.atleast_1d(x2),
                                  np.atleast_1d(x3))
    else:
        return getScaleFactorscart(np.atleast_1d(x1), np.atleast_1d(x2),
                                  np.atleast_1d(x3))


def getXYZcart(x, y, z, pars):
    return x, y, z


def getXYZcyl(r, phi, z, pars):
    x = r*np.cos(phi)
    y = r*np.sin(phi)
    return x, y, z


def getXYZsph(r, phi, th, pars):
    x = r*np.cos(phi)*np.sin(th)
    y = r*np.sin(phi)*np.sin(th)
    z = r*np.cos(th)
    return x, y, z


def getCentroidcart(xm, xp, dim):
    return 0.5*(xm+xp)


def getCentroidcyl(xm, xp, dim):
    if dim == 1:
        xc = 2.0*(xm*xm + xm*xp + xp*xp) / (3.0*(xm+xp))
    else:
        xc = 0.5*(xm+xp)
    return xc


def getCentroidsph(xm, xp, dim):
    if dim == 1:
        xc = 3.0*(xm*xm*xm + xm*xm*xp + xm*xp*xp + xp*xp*xp
                  ) / (4.0*(xm*xm+xm*xp+xp*xp))
    elif dim == 3:
        xc = np.empty(xm.shape)
        xc[xm == xp] = xm[xm == xp]
        dif = xm != xp
        cp = np.cos(xp)
        cm = np.cos(xm)
        sp = np.sin(xp)
        sm = np.sin(xm)
        xc[dif] = ((xm*cm-xp*cp+sp-sm) / (cm-cp))[dif]
    else:
        xc = 0.5*(xm+xp)
    return xc


def getDVcart(xjph, yiph, zkph, ny, index, pars):

    nz = ny.shape[0]
    nx = ny.shape[1]

    dx1d = xjph[1:] - xjph[:-1]
    dz1d = zkph[1:] - zkph[:-1]

    ymax = pars['Phi_Max']

    dx = np.empty(yiph.shape)
    dy = np.empty(yiph.shape)
    dz = np.empty(yiph.shape)

    for k in range(nz):
        for j in range(nx):
            a = index[k, j]
            b = index[k, j]+ny[k, j]

            dx[a:b] = dx1d[j]
            dz[a:b] = dz1d[k]

            yp = yiph[a:b]
            ym = np.roll(yp, 1)
            ym[ym > yp] -= ymax
            dy[a:b] = yp-ym

    if nz > 1:
        dV = dx*dy*dz
    else:
        dV = dx*dy

    return dV


def getDVcyl(rjph, piph, zkph, nphi, index, pars):

    nz = nphi.shape[0]
    nr = nphi.shape[1]

    r1d = 0.5*(rjph[1:] + rjph[:-1])
    dr1d = rjph[1:] - rjph[:-1]
    dz1d = zkph[1:] - zkph[:-1]

    pmax = pars['Phi_Max']

    r = np.empty(piph.shape)
    dr = np.empty(piph.shape)
    dp = np.empty(piph.shape)
    dz = np.empty(piph.shape)

    for k in range(nz):
        for j in range(nr):
            a = index[k, j]
            b = index[k, j]+nphi[k, j]

            r[a:b] = r1d[j]
            dr[a:b] = dr1d[j]
            dz[a:b] = dz1d[k]

            pp = piph[a:b]
            pm = np.roll(pp, 1)
            pm[pm > pp] -= pmax
            dp[a:b] = pp-pm

    if nz > 1:
        dV = r*dr*dp*dz
    else:
        dV = r*dr*dp

    return dV


def getDVsph(rjph, piph, tkph, nphi, index, pars):

    nt = nphi.shape[0]
    nr = nphi.shape[1]

    r21d = (rjph[1:]*rjph[1:] + rjph[1:]*rjph[:-1] + rjph[:-1]*rjph[:-1])/3.0
    sinth1d = np.sin(0.5*(tkph[:-1]+tkph[1:]))
    dr1d = rjph[1:] - rjph[:-1]
    sindth1d = 2.0*np.sin(0.5*(tkph[1:] - tkph[:-1]))

    pmax = pars['Phi_Max']

    r2 = np.empty(piph.shape)
    sinth = np.empty(piph.shape)
    dr = np.empty(piph.shape)
    dp = np.empty(piph.shape)
    sindth = np.empty(piph.shape)

    for k in range(nt):
        for j in range(nr):
            a = index[k, j]
            b = index[k, j]+nphi[k, j]

            r2[a:b] = r21d[j]
            dr[a:b] = dr1d[j]
            sinth[a:b] = sinth1d[k]
            sindth[a:b] = sindth1d[k]

            pp = piph[a:b]
            pm = np.roll(pp, 1)
            pm[pm > pp] -= pmax
            dp[a:b] = pp-pm

    dV = r2*sinth*dr*dp*sindth

    return dV


def integrate(f, dat, opts, pars, dV=None):
    if dV is None:
        dV = getDV(dat, opts, pars)

    nphi = dat[7]
    index = dat[8]

    mask = np.empty(f.shape, dtype=bool)
    mask[:] = True

    nz = nphi.shape[0]
    nr = nphi.shape[1]

    for k in range(nz):
        for j in range(nr):
            a = index[k, j]
            b = index[k, j] + nphi[k, j]
            if (j < 2 and pars["NoBC_Rmin"] == 0)\
                    or (j >= nr-2 and pars["NoBC_Rmax"] == 0)\
                    or (k < 2 and nz > 1 and pars["NoBC_Zmin"] == 0)\
                    or (k >= nz-2 and nz > 1 and pars["NoBC_Zmax"] == 0):
                mask[a:b] = False

    return (f[mask]*dV[mask]).sum()


def getScaleFactorscart(x, y, z):

    hx = np.ones(y.shape)
    hy = np.ones(y.shape)
    hz = np.ones(y.shape)

    return hx, hy, hz


def getScaleFactorscyl(r, phi, z):

    hr = np.ones(phi.shape)
    hp = np.empty(phi.shape)
    hz = np.ones(phi.shape)

    hp[:] = r

    return hr, hp, hz


def getScaleFactorssph(r, phi, th):

    hr = np.ones(phi.shape)
    hp = np.empty(phi.shape)
    ht = np.empty(phi.shape)

    hp[:] = r * np.sin(th)
    ht[:] = r

    return hr, hp, ht


def getVecXYZcart(x, y, z, vx, vy, vz):
    return vx, vy, vz


def getVecXYZcyl(r, phi, z, vr, vp, vz):

    cp = np.cos(phi)
    sp = np.sin(phi)

    vx = cp*vr - sp*vp
    vy = sp*vr + cp*vp

    return vx, vy, vz


def getVecXYZsph(r, phi, th, vr, vp, vt):

    cp = np.cos(phi)
    sp = np.sin(phi)
    ct = np.cos(th)
    st = np.sin(th)

    vx = cp*(st*vr + ct*vt) - sp*vp
    vy = sp*(st*vr + ct*vt) + cp*vp
    vz = ct*vr - st*vt

    return vx, vy, vz


def getSignedDphi(phip, phim, pars):

    pmax = pars['Phi_Max']

    dphi = phip-phim
    while dphi > 0.5*pmax:
        dphi -= pmax
    while dphi < -0.5*pmax:
        dphi += pmax

    return dphi


def calcDivB(dat, opts, pars, dV=None):
    if dV is None:
        dV = getDV(dat, opts, pars)

    Phi = dat[5]
    nphi = dat[7]
    index = dat[8]
    piph = dat[3]
    idPhi0 = dat[6]

    nz = nphi.shape[0]
    nr = nphi.shape[1]

    divB = np.zeros(piph.shape)

    PhiP = Phi[:, 0]
    PhiRm = Phi[:, 1]
    PhiRp = Phi[:, 2]
    PhiZm = Phi[:, 3]
    PhiZp = Phi[:, 4]

    # phi direction
    for k in range(nz):
        for j in range(nr):
            a = index[k, j]
            b = index[k, j] + nphi[k, j]

            divB[a:b] += PhiP[a:b] - np.roll(PhiP[a:b], 1)

    # r direction
    for k in range(nz):
        for jL in range(nr-1):
            jR = jL+1
            aL = index[k, jL]
            bL = index[k, jL] + nphi[k, jL]
            aR = index[k, jR]
            bR = index[k, jR] + nphi[k, jR]
            piphL = piph[aL:bL]
            piphR = piph[aR:bR]
            PhiRL = PhiRp[aL:bL]
            PhiRR = PhiRm[aR:bR]

            iL0 = idPhi0[k, jL] - aL
            iR0 = idPhi0[k, jR] - aR
            iL = iL0
            iR = iR0
            done = False

            while not done:
                phiL = piphL[iL]
                phiR = piphR[iR]

                dp = getSignedDphi(phiL, phiR, pars)

                if dp < 0.0:
                    phi_back = np.roll(piphR, 1)
                    dp1 = getSignedDphi(phiL, phi_back[iR], pars)
                    if dp1 <= 0.0:
                        print("Uh Oh")
                if dp > 0.0:
                    phi_back = np.roll(piphL, 1)
                    dp1 = getSignedDphi(phiR, phi_back[iL], pars)
                    if dp1 <= 0.0:
                        print("Oh Uh")

                if dp == 0.0:
                    divB[aL+iL] += PhiRL[iL]
                    divB[aR+iR] -= PhiRL[iL]
                    iL += 1
                    iR += 1
                elif dp < 0.0:
                    divB[aL+iL] += PhiRL[iL]
                    divB[aR+iR] -= PhiRL[iL]
                    iL += 1
                else:
                    divB[aL+iL] += PhiRR[iR]
                    divB[aR+iR] -= PhiRR[iR]
                    iR += 1

                if iL == nphi[k, jL]:
                    iL = 0
                if iR == nphi[k, jR]:
                    iR = 0

                if iL == iL0 and iR == iR0:
                    done = True
    # z direction

    for kL in range(nz-1):
        for j in range(nr):
            kR = kL+1
            aL = index[kL, j]
            bL = index[kL, j] + nphi[kL, j]
            aR = index[kR, j]
            bR = index[kR, j] + nphi[kR, j]
            piphL = piph[aL:bL]
            piphR = piph[aR:bR]
            PhiZL = PhiZp[aL:bL]
            PhiZR = PhiZm[aR:bR]

            pimhL = np.roll(piphL, 1)
            pimhR = np.roll(piphR, 1)

            iL0 = idPhi0[kL, j] - aL
            iR0 = idPhi0[kR, j] - aR
            iL = iL0
            iR = iR0
            done = False

            while not done:
                phiL = piphL[iL]
                phiR = piphR[iR]

                dp = getSignedDphi(phiL, phiR, pars)

                if dp < 0.0:
                    dp1 = getSignedDphi(phiL, pimhR[iR], pars)
                    if dp1 <= 0.0:
                        print("Uh Oh")
                if dp > 0.0:
                    dp1 = getSignedDphi(phiR, pimhL[iL], pars)
                    if dp1 <= 0.0:
                        print("Oh Uh")

                if dp == 0.0:
                    print("quad?")
                    divB[aL+iL] += PhiZL[iL]
                    divB[aR+iR] -= PhiZL[iL]
                    iL += 1
                    iR += 1
                elif dp < 0.0:
                    divB[aL+iL] += PhiZL[iL]
                    divB[aR+iR] -= PhiZL[iL]
                    iL += 1
                else:
                    divB[aL+iL] += PhiZR[iR]
                    divB[aR+iR] -= PhiZR[iR]
                    iR += 1

                if iL == nphi[kL, j]:
                    iL = 0
                if iR == nphi[kR, j]:
                    iR = 0

                if iL == iL0 and iR == iR0:
                    done = True

    divB /= dV

    return divB


def getDA1(x1, x2ph, x3f, dat, opts, pars):

    x2max = pars['Phi_Max']
    index = dat[8]
    nx2 = dat[7]
    x2mh = np.empty(x2ph.shape)
    x3mh = np.empty(x2ph.shape)
    x3ph = np.empty(x2ph.shape)
    for k in range(index.shape[0]):
        for j in range(index.shape[1]):
            a = index[k, j]
            b = index[k, j] + nx2[k, j]
            x2mh[a:b] = np.roll(x2ph[a:b], 1)
            x3mh[a:b] = x3f[k]
            x3ph[a:b] = x3f[k+1]
    x2mh[x2mh > x2ph] -= x2max

    return getDA1grid(x1, x2ph, x2mh, x3ph, x3mh, opts)


def getDA2(x1f, x2, x3f, dat, opts, pars):

    index = dat[8]
    nx2 = dat[7]
    x3ph = np.empty(x2.shape)
    x3mh = np.empty(x2.shape)
    x1mh = np.empty(x2.shape)
    x1ph = np.empty(x2.shape)
    for k in range(index.shape[0]):
        for j in range(index.shape[1]):
            a = index[k, j]
            b = index[k, j] + nx2[k, j]
            x1mh[a:b] = x1f[j]
            x1ph[a:b] = x1f[j+1]
            x3mh[a:b] = x3f[k]
            x3ph[a:b] = x3f[k+1]

    return getDA2grid(x1ph, x1mh, x2, x3ph, x3mh, opts)


def getDA3(x1f, x2ph, x3, dat, opts, pars):

    x2max = pars['Phi_Max']
    index = dat[8]
    nx2 = dat[7]
    x2mh = np.empty(x2ph.shape)
    x1mh = np.empty(x2ph.shape)
    x1ph = np.empty(x2ph.shape)
    for k in range(index.shape[0]):
        for j in range(index.shape[1]):
            a = index[k, j]
            b = index[k, j] + nx2[k, j]
            x2mh[a:b] = np.roll(x2ph[a:b], 1)
            x1mh[a:b] = x1f[j]
            x1ph[a:b] = x1f[j+1]
    x2mh[x2mh > x2ph] -= x2max

    return DA3grid(x1ph, x1mh, x2ph, x2mh, x3, opts)


def getDA1grid(x1, x2ph, x2mh, x3ph, x3mh, opts):
    if opts['GEOMETRY'] == "cylindrical":
        return getDA1cyl(x1, x2ph, x2mh, x3ph, x3mh)
    elif opts['GEOMETRY'] == "spherical":
        return getDA1sph(x1, x2ph, x2mh, x3ph, x3mh)
    else:
        return getDA1cart(x1, x2ph, x2mh, x3ph, x3mh)


def getDA2grid(x1ph, x1mh, x2, x3ph, x3mh, opts):
    if opts['GEOMETRY'] == "cylindrical":
        return getDA2cyl(x1ph, x1mh, x2, x3ph, x3mh)
    elif opts['GEOMETRY'] == "spherical":
        return getDA2sph(x1ph, x1mh, x2, x3ph, x3mh)
    else:
        return getDA2cart(x1ph, x1mh, x2, x3ph, x3mh)


def getDA3grid(x1ph, x1mh, x2ph, x2mh, x3, opts):
    if opts['GEOMETRY'] == "cylindrical":
        return getDA3cyl(x1ph, x1mh, x2ph, x2mh, x3)
    elif opts['GEOMETRY'] == "spherical":
        return getDA3sph(x1ph, x1mh, x2ph, x2mh, x3)
    else:
        return getDA3cart(x1ph, x1mh, x2ph, x2mh, x3)


def getDA1cart(x1, x2ph, x2mh, x3ph, x3mh):
    return (x2ph-x2mh)*(x3ph-x3mh)


def getDA1cyl(x1, x2ph, x2mh, x3ph, x3mh):
    return x1*(x2ph-x2mh)*(x3ph-x3mh)


def getDA1sph(x1, x2ph, x2mh, x3ph, x3mh):
    sinth = np.sin(0.5*(x3ph+x3mh))
    sindth = 2*np.sin(0.5*(x3ph-x3mh))
    return x1*x1*sinth*sindth*(x2ph-x2mh)


def getDA2cart(xph, xmh, y, zph, zmh):
    return (xph-xmh)*(zph-zmh)

def getDA2cyl(rph, rmh, phi, zph, zmh):
    return (rph-rmh)*(zph-zmh)

def getDA2sph(rph, rmh, phi, tph, tmh):
    return 0.5*(rph+rmh)*(rph-rmh)*(tph-tmh)


def getDA3cart(xph, xmh, yph, ymh, z):
    return (xph-xmh)*(yph-ymh)


def getDA3cyl(rph, rmh, pph, pmh, z):
    return 0.5*(rph+rmh)*(rph-rmh)*(pph-pmh)


def getDA3sph(rph, rmh, pph, pmh, th):
    sinth = np.sin(th)
    return 0.5*(rph+rmh)*(rph-rmh)*(pph-pmh)


def integrate2(f, dat, opts, pars):

    index = dat[8]
    n1 = index.shape[1]
    n2 = dat[7]
    n3 = index.shape[0]
    x2ph = dat[3]
    x2max = pars['Phi_Max']

    integral = np.empty((n3, n1), dtype=f.dtype)

    for k in range(n3):
        for j in range(n1):
            a = index[k, j]
            b = index[k, j] + n2[k, j]
            x2p = x2ph[a:b]
            x2m = np.roll(x2p, 1)
            x2m[x2m > x2p] -= x2max
            dx2 = x2p-x2m
            integral[k, j] = (f[a:b]*dx2).sum()

    if n3 == 1:
        integral = integral[0, :]

    return integral


def integrateTrans1(x1, f, dat, opts, pars, dA=None):

    if dA is None:
        dA = getDA1(x1, dat[3], dat[1], dat, opts, pars)

    mask = np.empty(f.shape, dtype=bool)
    mask[:] = True
    x1ph = dat[0]
    index = dat[8]
    n1 = index.shape[1]
    n2 = dat[7]
    n3 = index.shape[0]

    for k in range(n3):
        for j in range(n1):
            a = index[k, j]
            b = index[k, j] + n2[k, j]
            if (k < 2 and n3 > 1 and pars["NoBC_Zmin"] == 0)\
                    or (k >= n3-2 and n3 > 1 and pars["NoBC_Zmax"] == 0):
                mask[a:b] = False

    integral = np.empty(n1, dtype=f.dtype)
    for j in range(n1):
        ind = (x1 < x1ph[j+1]) & (x1 > x1ph[j])
        integral[j] = (f[mask*ind]*dA[mask*ind]).sum()

    return integral


def calculateGrad(x1, x2, x3, f, dat, opts, pars):

    dfdx2 = calculateGradX2(x2, f, dat, pars)
    dfdx1 = calculateGradX1(x1, x2, f, dat, opts, pars, dfdx2)
    dfdx3 = calculateGradX3(x2, x3, f, dat, opts, pars, dfdx2)

    return dfdx1, dfdx2, dfdx3


def calculateGradX2(x2, f, dat, pars):

    index = dat[8]
    n1 = index.shape[1]
    n2 = dat[7]
    n3 = index.shape[0]
    x2_max = pars['Phi_Max']
    
    dfdx2 = np.empty(f.shape)

    multiD = True if len(f.shape) > 1 else False
    extraDims = (len(f.shape)-1) * (1,)

    for k in range(n3):
        for j in range(n1):
            a = index[k, j]
            b = index[k, j] + n2[k, j]

            fjk = f[a:b]
            x2jk = x2[a:b]
            if multiD:
                x2jk = x2jk.reshape(x2jk.shape+extraDims)

            x2jkR = np.roll(x2jk, -1)
            x2jkL = np.roll(x2jk, 1)

            dx2 = np.roll(x2jk, -1) - np.roll(x2jk, 1)
            dx2[dx2 < 0] += x2_max
            dx2[dx2 > x2_max] -= x2_max

            dfdx2[a:b] = (np.roll(fjk, -1, axis=0)
                          - np.roll(fjk, 1, axis=0)) / dx2

    return dfdx2


def calculateGradX1(x1, x2, f, dat, opts, pars, dfdx2=None):

    index = dat[8]
    n1 = index.shape[1]
    n2 = dat[7]
    n3 = index.shape[0]
    x2_max = pars['Phi_Max']
    x1f = dat[0]
    x3f = dat[1]
    x2ph = dat[3]
    
    dfdx1 = np.zeros(f.shape)

    if n1 <= 1:
        return dfdx1

    dAtot = np.zeros(f.shape)
    if dfdx2 is None:
        dfdx2 = calculateGradX2(x2, f, dat, pars)
    
    multiD = True if len(f.shape) > 1 else False
    extraDims = (len(f.shape)-1) * (1,)

    for k in range(n3):
        for jf in range(1, n1):
            jL = jf-1
            jR = jf
            aL = index[k, jL]
            bL = index[k, jL] + n2[k, jL]
            aR = index[k, jR]
            bR = index[k, jR] + n2[k, jR]

            x2phL = unwrap(x2ph[aL:bL], x2_max)
            x2phR = unwrap(x2ph[aR:bR], x2_max)
            fL = f[aL:bL]
            fR = f[aR:bR]
            dfdx2L = dfdx2[aL:bL]
            dfdx2R = dfdx2[aR:bR]
            x1L = x1[aL:bL]
            x1R = x1[aR:bR]

            x20 = x2phL[0]
            while x20 > x2phR[-1]:
                x2phR += x2_max
            while x20 <= x2phR[-1]-x2_max:
                x2phR -= x2_max
            x2mhL = np.roll(x2phL, 1)
            x2mhL[0] -= x2_max
            x2mhR = np.roll(x2phR, 1)
            x2mhR[0] -= x2_max

            iL0 = 0
            iR0 = np.searchsorted(x2phR, x20)
            if iR0 >= n2[k, jR]:
                iR0 = 0
            # iL = iL0
            # iR = iR0

            n2f = n2[k, jL] + n2[k, jR]

            iL = np.empty(n2f, dtype=int)
            iR = np.empty(n2f, dtype=int)

            iLi = iL0
            iRi = iR0

            adjLi = 0.0
            adjRi = 0.0

            for i in range(n2f):
                iL[i] = iLi
                iR[i] = iRi
                if x2phL[iLi]+adjLi < x2phR[iRi]+adjRi:
                    iLi += 1
                    if iLi >= n2[k, jL]:
                        iLi = 0
                        adjLi += x2_max
                else:
                    iRi += 1
                    if iRi >= n2[k, jR]:
                        iRi = 0
                        adjRi += x2_max

            x2p = np.minimum(unwrap(x2phL[iL], x2_max),
                             unwrap(x2phR[iR], x2_max))
            x2m = np.maximum(unwrap(x2mhL[iL], x2_max),
                             unwrap(x2mhR[iR], x2_max))
            

            x2i = 0.5*(x2p + x2m)
            x2L = unwrap(x2[aL:bL][iL], x2_max)
            x2R = unwrap(x2[aR:bR][iR], x2_max)

            dx2L = x2i-x2L
            dx2L[dx2L > 0.5*x2_max] -= x2_max
            dx2L[dx2L < -0.5*x2_max] += x2_max
            dx2R = x2i-x2R
            dx2R[dx2R > 0.5*x2_max] -= x2_max
            dx2R[dx2R < -0.5*x2_max] += x2_max

            dAi = getDA1grid(x1f[jf], x2p, x2m, x3f[k+1], x3f[k], opts)
            dx1i = x1R[iR] - x1L[iL]
            if multiD:
                dx2L = dx2L.reshape(dx2L.shape + extraDims)
                dx2R = dx2R.reshape(dx2R.shape + extraDims)
                dx1i = dx1i.reshape(dx1i.shape + extraDims)
                dAi = dAi.reshape(dAi.shape + extraDims)

            fLi = fL[iL] + dfdx2L[iL] * dx2L
            fRi = fR[iR] + dfdx2R[iR] * dx2R

            dfdxi = (fRi - fLi) / dx1i

            np.add.at(dfdx1, aL+iL, dfdxi * dAi)
            np.add.at(dAtot, aL+iL, dAi)
            np.add.at(dfdx1, aR+iR, dfdxi * dAi)
            np.add.at(dAtot, aR+iR, dAi)

            if (x2p <= x2m).any():
                print("PROBLEM", k, jf)
            
    dfdx1 /= dAtot

    return dfdx1


def calculateGradX3(x2, x3, f, dat, opts, pars, dfdx2=None):

    index = dat[8]
    n1 = index.shape[1]
    n2 = dat[7]
    n3 = index.shape[0]
    x2_max = pars['Phi_Max']
    x1f = dat[0]
    x3f = dat[1]
    x2ph = dat[3]
    
    dfdx3 = np.zeros(f.shape)

    if n3 <= 1:
        return dfdx3

    dAtot = np.zeros(f.shape)
    if dfdx2 is None:
        dfdx2 = calculateGradX2(x2, f, dat, pars)
    
    multiD = True if len(f.shape) > 1 else False
    extraDims = (len(f.shape)-1) * (1,)

    for kf in range(1, n3):
        for j in range(n1):
            kL = kf-1
            kR = kf
            aL = index[kL, j]
            bL = index[kL, j] + n2[kL, j]
            aR = index[kR, j]
            bR = index[kR, j] + n2[kR, j]

            x2phL = unwrap(x2ph[aL:bL], x2_max)
            x2phR = unwrap(x2ph[aR:bR], x2_max)
            fL = f[aL:bL]
            fR = f[aR:bR]
            dfdx2L = dfdx2[aL:bL]
            dfdx2R = dfdx2[aR:bR]
            x3L = x3[aL:bL]
            x3R = x3[aR:bR]

            x20 = x2phL[0]
            while x20 > x2phR[-1]:
                x2phR += x2_max
            while x20 <= x2phR[-1]-x2_max:
                x2phR -= x2_max
            x2mhL = np.roll(x2phL, 1)
            x2mhL[0] -= x2_max
            x2mhR = np.roll(x2phR, 1)
            x2mhR[0] -= x2_max

            iL0 = 0
            iR0 = np.searchsorted(x2phR, x20)
            if iR0 >= n2[k, jR]:
                iR0 = 0
            # iL = iL0
            # iR = iR0

            n2f = n2[kL, j] + n2[kR, j]

            iL = np.empty(n2f, dtype=int)
            iR = np.empty(n2f, dtype=int)

            iLi = iL0
            iRi = iR0

            adjLi = 0.0
            adjRi = 0.0

            for i in range(n2f):
                iL[i] = iLi
                iR[i] = iRi
                if x2phL[iLi]+adjLi < x2phR[iRi]+adjRi:
                    iLi += 1
                    if iLi >= n2[k, jL]:
                        iLi = 0
                        adjLi += x2_max
                else:
                    iRi += 1
                    if iRi >= n2[k, jR]:
                        iRi = 0
                        adjRi += x2_max

            x2p = np.minimum(unwrap(x2phL[iL], x2_max),
                             unwrap(x2phR[iR], x2_max))
            x2m = np.maximum(unwrap(x2mhL[iL], x2_max),
                             unwrap(x2mhR[iR], x2_max))
            

            x2i = 0.5*(x2p + x2m)
            x2L = unwrap(x2[aL:bL][iL], x2_max)
            x2R = unwrap(x2[aR:bR][iR], x2_max)

            dx2L = x2i-x2L
            dx2L[dx2L > 0.5*x2_max] -= x2_max
            dx2L[dx2L < -0.5*x2_max] += x2_max
            dx2R = x2i-x2R
            dx2R[dx2R > 0.5*x2_max] -= x2_max
            dx2R[dx2R < -0.5*x2_max] += x2_max

            dAi = getDA3grid(x1f[j+1], x1f[j], x2p, x2m, x3f[kf], opts)
            dx3i = x3R[iR] - x3L[iL]
            if multiD:
                dx2L = dx2L.reshape(dx2L.shape + extraDims)
                dx2R = dx2R.reshape(dx2R.shape + extraDims)
                dx3i = dx3i.reshape(dx3i.shape + extraDims)
                dAi = dAi.reshape(dAi.shape + extraDims)

            fLi = fL[iL] + dfdx2L[iL] * dx2L
            fRi = fR[iR] + dfdx2R[iR] * dx2R

            dfdxi = (fRi - fLi) / dx3i

            np.add.at(dfdx3, aL+iL, dfdxi * dAi)
            np.add.at(dAtot, aL+iL, dAi)
            np.add.at(dfdx3, aR+iR, dfdxi * dAi)
            np.add.at(dAtot, aR+iR, dAi)

            if (x2p <= x2m).any():
                print("PROBLEM", k, jf)
    
    dfdx3 /= dAtot

    return dfdx3


def unwrap(x, maxval):
    # 'unwraps' arr x with periodicity maxval
    # result is increasing arr

    diff = x[1:] - x[:-1]

    moddiff = np.mod(diff, maxval)
    corr = np.cumsum(moddiff-diff)

    x2 = np.empty(x.shape)
    x2[0] = x[0]
    x2[1:] = x[1:] + corr

    return x2


def calculateDivV(x1, x2, x3, v1, v2, v3, dat, opts, pars, dV=None):

    index = dat[8]
    n1 = index.shape[1]
    n2 = dat[7]
    n3 = index.shape[0]
    x2_max = pars['Phi_Max']
    x1f = dat[0]
    x2ph = dat[3]
    x3f = dat[1]

    if dV is None:
        dV = getDV(dat, opts, pars)
    
    dv1dx2 = calculateGradX2(x2, v1, dat, pars)
    dv2dx2 = calculateGradX2(x2, v2, dat, pars)
    dv3dx2 = calculateGradX2(x2, v3, dat, pars)

    divV = np.zeros(v1.shape)

    # x2 part
   
    print("x2 div")

    for k in range(n3):
        for j in range(n1):
            a = index[k, j]
            b = index[k, j] + n2[k, j]

            dAph = getDA2grid(x1f[j+1], x1f[j], x2ph[a:b], x3f[k+1], x3f[k],
                              opts)

            v = v2[a:b] 
            vR = np.roll(v, -1)
            vph = 0.5*(v+vR)

            _, hph, _ = getScaleFactors(x1[a:b], x2ph[a:b], x3[a:b], opts)

            fph = vph*hph*dAph
            fmh = np.roll(fph, 1)

            divV[a:b] += fph-fmh


    # x1 part

    print("x1 div")

    for k in range(n3):
        for jf in range(1, n1):
            jL = jf-1
            jR = jf
            aL = index[k, jL]
            bL = index[k, jL] + n2[k, jL]
            aR = index[k, jR]
            bR = index[k, jR] + n2[k, jR]

            x2phL = unwrap(x2ph[aL:bL], x2_max)
            x2phR = unwrap(x2ph[aR:bR], x2_max)
            v1L = v1[aL:bL]
            v1R = v1[aR:bR]
            dv1dx2L = dv1dx2[aL:bL]
            dv1dx2R = dv1dx2[aR:bR]

            x20 = x2phL[0]
            while x20 > x2phR[-1]:
                x2phR += x2_max
            while x20 <= x2phR[-1]-x2_max:
                x2phR -= x2_max
            x2mhL = np.roll(x2phL, 1)
            x2mhL[0] -= x2_max
            x2mhR = np.roll(x2phR, 1)
            x2mhR[0] -= x2_max

            iL0 = 0
            iR0 = np.searchsorted(x2phR, x20)
            if iR0 >= n2[k, jR]:
                iR0 = 0
            iL = iL0
            iR = iR0

            for i in range(n2[k, jL] + n2[k, jR]):
                x2p = min(x2phL[iL], x2phR[iR])
                x2m = max(x2mhL[iL], x2mhR[iR])
                x2i = 0.5*(x2p + x2m)
                x2L = 0.5*(x2phL[iL] + x2mhL[iL])
                x2R = 0.5*(x2phR[iR] + x2mhR[iR])
                vL = v1L[iL] + dv1dx2L[iL]*(x2i-x2L)
                vR = v1R[iR] + dv1dx2R[iR]*(x2i-x2R)

                dA = getDA1grid(x1f[jf], x2p, x2m, x3f[k+1], x3f[k], opts)
                h, _, _ = getScaleFactors(x1f[jf], x2i,
                                          0.5*(x3[aL+iL]+x3[aR+iR]), opts)
                divV[aL+iL] += 0.5*(vL+vR)*h*dA
                divV[aR+iR] -= 0.5*(vL+vR)*h*dA

                if x2phL[iL] < x2phR[iR]:
                    iL += 1
                    if iL >= n2[k, jL]:
                        iL = 0
                        x2phL += x2_max
                        x2mhL += x2_max
                else:
                    iR += 1
                    if iR >= n2[k, jR]:
                        iR = 0
                        x2phR += x2_max
                        x2mhR += x2_max
            if iL != iL0 or iR != iR0:
                print("PROBLEM", aL, bL-aL, iL0, iL,"|", aR, bR-aR, iR0, iR)

    # x3 part
    print("x3 div")

    for k in range(1, n3):
        for jf in range(n1):
            kL = kf-1
            kR = kf
            aL = index[kL, j]
            bL = index[kL, j] + n2[kL, j]
            aR = index[kR, j]
            bR = index[kR, j] + n2[kR, j]

            x2phL = unwrap(x2ph[aL:bL], x2_max)
            x2phR = unwrap(x2ph[aR:bR], x2_max)
            v3L = v3[aL:bL]
            v3R = v3[aR:bR]
            dv3dx2L = dv3dx2[aL:bL]
            dv3dx2R = dv3dx2[aR:bR]

            x20 = x2phL[0]
            while x20 > x2phR[-1]:
                x2phR += x2_max
            while x20 <= x2phR[0]:
                x2phR -= x2_max
            x2mhL = np.roll(x2phL, 1)
            x2mhL[0] -= x2_max
            x2mhR = np.roll(x2phR, 1)
            x2mhR[0] -= x2_max

            iL = 0
            iR = np.searchsorted(x2phR, x20)

            for i in range(n2[kL, j] + n2[kR, j]):
                x2p = min(x2phL[iL], x2phR[iR])
                x2m = max(x2mhL[iL], x2mhR[iR])
                x2i = 0.5*(x2p + x2m)
                x2L = 0.5*(x2phL[iL] + x2mhL[iL])
                x2R = 0.5*(x2phR[iR] + x2mhR[iR])
                vL = v3L[iL] + dv3dx2L[iL]*(x2i-x2L)
                vR = v3R[iR] + dv3dx2R[iR]*(x2i-x2R)

                dA = getDA3grid(x1f[j+1], x1f[j], x2p, x2m, x3f[kf], opts)
                _, _, h = getScaleFactors(0.5*(x1[aL+iL]+x1[aR+iR]),
                                          x2i, x3f[kf], opts)
                
                divV[aL+iL] += 0.5*(vL+vR)*h*dA
                divV[aR+iR] -= 0.5*(vL+vR)*h*dA

                if x2phL[iL] < x2phR[iR]:
                    iL += 1
                    if iL >= n2[k, jL]:
                        iL = 0
                else:
                    iR += 1
                    if iR >= n2[k, jR]:
                        iR = 0

    divV /= dV

    return divV

def calculateCurlV(x1, x2, x3, v1, v2, v3, dat, opts, pars, dV=None):

    index = dat[8]
    n1 = index.shape[1]
    n2 = dat[7]
    n3 = index.shape[0]
    x2_max = pars['Phi_Max']
    x1f = dat[0]
    x2ph = dat[3]
    x3f = dat[1]

    if dV is None:
        dV = getDV(dat, opts, pars)
    dv1dx2 = calculateGradX2(x2, v1, dat, pars)
    dv2dx2 = calculateGradX2(x2, v2, dat, pars)
    dv3dx2 = calculateGradX2(x2, v3, dat, pars)

    vortZ = np.zeros(v1.shape)

    # x2 part
    print("x2 div")

    for k in range(n3):
        for j in range(n1):
            a = index[k, j]
            b = index[k, j] + n2[k, j]

            dAph = getDA2grid(x1f[j+1], x1f[j], x2ph[a:b], x3f[k+1], x3f[k],
                              opts)

            #v = v2[a:b]
            v = v1[a:b]
            vR = np.roll(v, -1)
            vph = 0.5*(v+vR)

            #_, hph, _ = getScaleFactors(x1[a:b], x2ph[a:b], x3[a:b], opts)
            hph, _, _ = getScaleFactors(x1[a:b], x2ph[a:b], x3[a:b], opts)

            fph = vph*hph*dAph
            fmh = np.roll(fph, 1)

            #divV[a:b] += fph-fmh
            vortZ[a:b] -= fph-fmh


    # x1 part

    print("x1 div")

    for k in range(n3):
        #for jf in range(1, 3):
        for jf in range(1, n1):
            #print(jf)
            jL = jf-1
            jR = jf
            aL = index[k, jL]
            bL = index[k, jL] + n2[k, jL]
            aR = index[k, jR]
            bR = index[k, jR] + n2[k, jR]

            x2phL = unwrap(x2ph[aL:bL], x2_max)
            x2phR = unwrap(x2ph[aR:bR], x2_max)
            #v1L = v1[aL:bL]
            #v1R = v1[aR:bR]
            v2L = v2[aL:bL]
            v2R = v2[aR:bR]
            #dv1dx2L = dv1dx2[aL:bL]
            #dv1dx2R = dv1dx2[aR:bR]
            dv2dx2L = dv2dx2[aL:bL]
            dv2dx2R = dv2dx2[aR:bR]

            x20 = x2phL[0]
            while x20 > x2phR[-1]:
                x2phR += x2_max
            while x20 <= x2phR[-1] - x2_max:
                x2phR -= x2_max
            x2mhL = np.roll(x2phL, 1)
            x2mhL[0] -= x2_max
            x2mhR = np.roll(x2phR, 1)
            x2mhR[0] -= x2_max

            iL0 = 0
            iR0 = np.searchsorted(x2phR, x20)

            #start changing things here

            maxL = n2[k, jL]
            maxR = n2[k, jR]

            #not vectorized, but apparently not a bottleneck
            iL = [iL0]
            iR = [iR0]
            diffL = np.zeros(maxL+maxR)
            diffR = np.zeros(maxL+maxR)
            for i in range( maxL + maxR-1):
              if x2phL[iL[i]]+diffL[i] < x2phR[iR[i]]+diffR[i]:
                iR.append(iR[-1])
                iL.append((iL[-1]+1)%maxL)
              else:
                iR.append((iR[-1]+1)%maxR)
                iL.append(iL[-1])
              if x2phL[iL[-2]]>x2phL[iL[-1]]: diffL[i+1:]+=x2_max
              if x2phR[iR[-2]]>x2phR[iR[-1]]: diffR[i+1:]+=x2_max

            x2p = np.minimum(x2phL[iL]+diffL, x2phR[iR]+diffR)
            x2m = np.maximum(x2mhL[iL]+diffL, x2mhR[iR]+diffR)

            x2i = 0.5*(x2p + x2m)

            x2L = 0.5*(x2phL[iL] + x2mhL[iL])+diffL
            x2R = 0.5*(x2phR[iR] + x2mhR[iR])+diffR

            vL = v2L[iL] + dv2dx2L[iL]*(x2i-x2L)
            vR = v2R[iR] + dv2dx2R[iR]*(x2i-x2R)

            dA = getDA1grid(x1f[jf], x2p, x2m, x3f[k+1], x3f[k], opts)
            _, h, _ = getScaleFactors(x1f[jf], x2i, 0.0, opts)

            np.add.at(vortZ, aL+iL, 0.5*(vL+vR)*h*dA)
            np.subtract.at(vortZ, aR+iR, 0.5*(vL+vR)*h*dA)


            '''
            #old version. VERY slow, especially for large N
            iL = iL0
            iR = iR0
            for i in range(n2[k, jL] + n2[k, jR]):
                x2p = min(x2phL[iL], x2phR[iR])
                x2m = max(x2mhL[iL], x2mhR[iR])
                x2i = 0.5*(x2p + x2m)
                x2L = 0.5*(x2phL[iL] + x2mhL[iL])
                x2R = 0.5*(x2phR[iR] + x2mhR[iR])
                vL = v2L[iL] + dv2dx2L[iL]*(x2i-x2L)
                vR = v2R[iR] + dv2dx2R[iR]*(x2i-x2R)

                dA = getDA1grid(x1f[jf], x2p, x2m, x3f[k+1], x3f[k], opts)
                _, h, _ = getScaleFactors(x1f[jf], x2i, 0.5*(x3[aL+iL]+x3[aR+iR]), opts)
                vortZ[aL+iL] += 0.5*(vL+vR)*h*dA
                vortZ[aR+iR] -= 0.5*(vL+vR)*h*dA

                if x2phL[iL] < x2phR[iR]:
                    iL += 1
                    if iL >= n2[k, jL]:
                        iL = 0
                        x2phL += x2_max
                        x2mhL += x2_max
                else:
                    iR += 1
                    if iR >= n2[k, jR]:
                        iR = 0
                        x2phR += x2_max
                        x2mhR += x2_max
            if iL != iL0 or iR != iR0:
                print("PROBLEM", aL, bL-aL, iL0, iL,"|", aR, bR-aR, iR0, iR)
            '''
    vortZ /= dV

    return vortZ
