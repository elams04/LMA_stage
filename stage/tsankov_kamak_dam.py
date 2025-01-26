import math
import os
from math import tan, atan
from typing import List

import numpy as np
from scipy.interpolate import CloughTocher2DInterpolator as CT

from Meshes.msh_transformer.instance import MSHInstance
from Meshes.msh_transformer.transformer import Transformer, BoundingBox
import Meshes.gmesher as m

dirname = os.path.dirname(__file__)
foundation = MSHInstance(os.path.join(dirname, "tsankov_kamak_dam", "Foundation.msh"))

foundation.apply_transform(
    [Transformer.material_filter(2)],
)

foundation.save("tsankov_kamak_dam/Foundation_transformed.msh")
m.gmsh.open("tsankov_kamak_dam/Foundation_transformed.msh")
m.gmsh.fltk.run()
m.gmsh.finalize()

outer = foundation.get_outer_points()

xp_faces = []
xm_faces = []

for p in outer:
    dirs = p.get_faces_directions()

    # if the point is on the x+ and x- side at the same time
    # it is not an outer point this test helps to get rid of
    # outer points that lay between the two faces
    if "x+" in dirs and "x-" in dirs:
        continue

    if "x+" in dirs:
        xp_faces.append(p)

    elif "x-" in dirs:
        xm_faces.append(p)

xp_zy, xp_x = [], []
xm_zy, xm_x = [], []

for p in xp_faces:
    xp_zy.append((p.z, p.y))
    xp_x.append(p.x)

for p in xm_faces:
    xm_zy.append((p.z, p.y))
    xm_x.append(p.x)


def interpolate(xy, z):
    """CT interpolator + nearest-neighbor extrapolation.

    Parameters
    ----------
    xy : ndarray, shape (npoints, ndim)
        Coordinates of data points
    z : ndarray, shape (npoints)
        Values at data points

    Returns
    -------
    func : callable
        A callable object which mirrors the CT behavior,
        with an additional neareast-neighbor extrapolation
        outside of the data range.
    """
    x = xy[:, 0]
    y = xy[:, 1]
    f = CT(xy, z)

    # this inner function will be returned to a user
    def new_f(xx, yy):
        # evaluate the CT interpolator. Out-of-bounds values are nan.
        zz = f(xx, yy)
        nans = np.isnan(zz)

        if nans.any():
            # for each nan point, find its nearest neighbor
            inds = np.argmin(
                (x[:, None] - xx[nans]) ** 2 + (y[:, None] - yy[nans]) ** 2, axis=0
            )
            # ... and use its value
            zz[nans] = z[inds]
        return zz

    return new_f


xp_zy = np.array(xp_zy)
xp_zy[:, 0] -= np.average(xp_zy[:, 0])
xp_zy[:, 1] -= np.average(xp_zy[:, 1])

xp_x = np.array(xp_x)

offset_x = np.average(xp_x)

xp_x -= offset_x

xm_zy = np.array(xm_zy)
xm_zy[:, 0] -= np.average(xm_zy[:, 0])
xm_zy[:, 1] -= np.average(xm_zy[:, 1])

xm_x = np.array(xm_x)
xm_x -= offset_x

forward_deformation = interpolate(xp_zy, xp_x)
backward_deformation = interpolate(xm_zy, xm_x)


if __name__ == "__main__":
    Z = np.linspace(-75, 75, 30)
    Y = np.linspace(-150, 150, 30)

    Z, Y = np.meshgrid(Z, Y)

    Xp = forward_deformation(Z, Y)
    Xm = backward_deformation(Z, Y)

    for i in range(30):
        for j in range(30):
            print(Y[i, j], Z[i, j], "|", Xp[i, j], Xm[i, j])
            m.Point(Y[i, j], Z[i, j], Xp[i, j])
            m.Point(Y[i, j], Z[i, j], (Xp[i, j] + Xm[i, j]) / 2)
            m.Point(Y[i, j], Z[i, j], Xm[i, j])

    m.geo.synchronize()

    m.gmsh.fltk.run()

    m.gmsh.finalize()
