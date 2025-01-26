import Meshes.main_model.gmesher as m
import sys

# ===========
#  Constants
# ===========
scale = 10  # 10m
elements = 2  # per unit of length

# ===========
# Parameters
# ===========
dam_height = 10
wall_height = 11
dam_size = 20
water_length = 7
water_depth = 9
dam_thickness = 1
wall_thickness = 3
ground_height = 1
air_length = 10
# ===========

model_length = water_length + dam_thickness + air_length

assert (wall_height >= dam_height)
assert (wall_height >= water_depth)

geo = m.geo
model = m.model

# bottom center
bc_g_points = m.get_rectangle_from_size(dam_size * scale, ground_height * scale, wall_thickness * scale)
bg_g_curve = m.get_curve_from_points(bc_g_points)

# transfinite curves
for i in m.get_lines_from_points(bc_g_points):
    geo.mesh.set_transfinite_curve(i, 8)

bg_g_plane = geo.add_plane_surface([bg_g_curve])
# transfinite surfaces
geo.mesh.set_transfinite_surface(bg_g_plane)
geo.mesh.set_recombine(2, 1)

bc_ground_volume = geo.extrude([(2, bg_g_plane)],
                               water_length * scale, 0, 0,
                               [water_length * elements],
                               recombine=True)
br_ground_volume = geo.extrude([bc_ground_volume[3]],
                               0, wall_thickness * scale, 0,
                               [wall_thickness * elements],
                               recombine=True)
bl_ground_volume = geo.extrude([bc_ground_volume[5]],
                               0, -wall_thickness * scale, 0,
                               [wall_thickness * elements],
                               recombine=True)

r_bw_ground_volume = geo.extrude([br_ground_volume[3]],
                                 0, 0, water_depth * scale,
                                 [water_depth * elements],
                                 recombine=True)
l_bw_ground_volume = geo.extrude([bl_ground_volume[5]],
                                 0, 0, water_depth * scale,
                                 [water_depth * elements],
                                 recombine=True)

# print('bc_ground_volume', m.get_volumes(bc_ground_volume, 1))
# print('br_ground_volume', m.get_volumes(br_ground_volume, 1))
# print('bl_ground_volume', m.get_volumes(bl_ground_volume, 1))
# print('r_bw_ground_volume', m.get_volumes(r_bw_ground_volume, 1))
# print('l_bw_ground_volume', m.get_volumes(l_bw_ground_volume, 1))

around_dam_extrusion = [r_bw_ground_volume[3],
                        l_bw_ground_volume[5],
                        bc_ground_volume[0],
                        br_ground_volume[4],
                        bl_ground_volume[4]]

ground_around_dam_volume = geo.extrude(
    around_dam_extrusion,
    dam_thickness * scale, 0, 0,
    [dam_thickness * elements],
    recombine=True)

# print('ground_around_dam_volume', m.get_volumes(ground_around_dam_volume, len(around_dam_extrusion)))

ground_after_dam_volume = geo.extrude(
    m.get_faces_in_extrusion_direction(ground_around_dam_volume, len(around_dam_extrusion)),
    air_length * scale, 0, 0,
    [air_length * elements],
    recombine=True)

# print('ground_after_dam_volume', m.get_volumes(ground_after_dam_volume, len(around_dam_extrusion)))

above_water_below_dam_extrusion = [
    ground_after_dam_volume[10],
    ground_after_dam_volume[4],
    ground_around_dam_volume[4],
    ground_around_dam_volume[10],
    r_bw_ground_volume[0],
    l_bw_ground_volume[0],
]

ground_above_water_below_dam_volume = geo.extrude(
    above_water_below_dam_extrusion,
    0, 0, (dam_height - water_depth) * scale,
    [(dam_height - water_depth) * elements],
    recombine=True)

# print('ground_above_water_below_dam_volume', m.get_volumes(ground_above_water_below_dam_volume,
#                                                            len(above_water_below_dam_extrusion)))

ground_above_water_above_dam_volume = geo.extrude(
    m.get_faces_in_extrusion_direction(ground_above_water_below_dam_volume, len(above_water_below_dam_extrusion) - 1),
    0, 0, (wall_height - dam_height) * scale,
    [(wall_height - dam_height) * elements],
    recombine=True)

# print('ground_above_water_above_dam_volume', m.get_volumes(ground_above_water_above_dam_volume,
#                                                            len(above_water_below_dam_extrusion)))

dam_bw_extrusion = [ground_around_dam_volume[16]]

dam_volume_bw = geo.extrude(
    dam_bw_extrusion,
    0, 0, water_depth * scale,
    [water_depth * elements],
    recombine=True)

print('dam_volume_bw', m.get_volumes(dam_volume_bw, len(dam_bw_extrusion)))

dam_volume_aw = geo.extrude(
    m.get_faces_in_extrusion_direction(dam_volume_bw, len(dam_bw_extrusion)),
    0, 0, (dam_height - water_depth) * scale,
    [(dam_height - water_depth) * elements],
    recombine=True)

print('dam_volume_aw', m.get_volumes(dam_volume_aw, len(dam_bw_extrusion)))

water_volume = geo.extrude(
    [bc_ground_volume[4]],
    0, 0, water_depth * scale,
    [water_depth * elements],
    recombine=True)

print('water_volume', water_volume[1])

# physical groups
ground_volumes = [bc_ground_volume[1],
                  br_ground_volume[1],
                  bl_ground_volume[1],
                  r_bw_ground_volume[1],
                  l_bw_ground_volume[1],
                  *m.get_volumes(ground_around_dam_volume, len(around_dam_extrusion) + 1),
                  *m.get_volumes(ground_after_dam_volume, len(around_dam_extrusion) + 1),
                  *m.get_volumes(ground_above_water_below_dam_volume, len(above_water_below_dam_extrusion)),
                  # *get_volumes(ground_above_water_above_dam_volume, len(above_water_below_dam_extrusion)),
                  ]

m.add_volumes_to_phy_group(ground_volumes, 'Ground')
m.add_volumes_to_phy_group([dam_volume_bw[1], dam_volume_aw[1]], 'Dam')
m.add_volumes_to_phy_group([water_volume[1]], 'Water')

geo.synchronize()

m.gmsh.model.mesh.generate(3)

# gmsh.write("barrage_simple3.msh")

if '-nopopup' not in sys.argv:
    m.gmsh.fltk.run()

m.gmsh.finalize()
