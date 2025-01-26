from typing import List, Tuple, Dict

import gmsh
import sys
from math import cos, sin

physical_groups = {}

options = {"transfinite": 3}


class Material:
    def __init__(
        self,
        _type: str,
        vel_p: float,
        vel_s: float,
        dens: float,
        q_p: float = 0,
        q_s: float = 0,
    ):
        self.type = _type
        self.vel_p = vel_p
        self.vel_s = vel_s
        self.dens = dens
        self.q_p = q_p
        self.q_s = q_s

    def copy_as_pml(self):
        if self.type == "S":
            new_type = "P"
        elif self.type == "F":
            new_type = "L"
        else:
            new_type = self.type

        return Material(new_type, self.vel_p, self.vel_s, self.dens, self.q_p, self.q_s)


class PhysicalGroup:
    def __init__(
        self,
        volumes: List["Volume"],
        name: str,
        material: Material,
        is_pml: bool = False,
    ):
        self.volumes = volumes
        self.name = name
        self.material = material.copy_as_pml() if is_pml else material

        group = len(physical_groups) + 1
        physical_groups[name] = group

        volume_ids = []
        for volume in volumes:
            volume_ids += volume.ids

        model.addPhysicalGroup(3, volume_ids, group, name)

    def format(self):
        m = self.material
        return (
            f"{m.type} {m.vel_p:.8f} {m.vel_s:.8f} {m.dens:.8f} {m.q_p:.8f} {m.q_s:.8f}"
        )


class PML:
    def __init__(
        self,
        x: Tuple[float, float],
        y: Tuple[float, float],
        z: Tuple[float, float],
        neighbor_material_offset: int,
    ):
        self.pos_x = x[0]
        self.w_x = x[1]
        self.pos_y = y[0]
        self.w_y = y[1]
        self.pos_z = z[0]
        self.w_z = z[1]
        self.neighbor_material_offset = neighbor_material_offset

    def format(self, n_pow: int, a_pow: float):
        return (
            " ".join(
                [
                    f"{e:.8f}"
                    for e in [
                        n_pow,
                        a_pow,
                        self.pos_x,
                        self.w_x,
                        self.pos_y,
                        self.w_y,
                        self.pos_z,
                        self.w_z,
                    ]
                ]
            )
            + f" {self.neighbor_material_offset}"
        )

    def __repr__(self):
        return f"PML({self.pos_x}, {self.w_x}, {self.pos_y}, {self.w_y}, {self.pos_z}, {self.w_z}, {self.neighbor_material_offset})"


class PhysicalGroupCollection:
    groups: Dict[str, Tuple[int, PhysicalGroup]]
    pml: Dict[str, PML]

    def __init__(self):
        self.groups = {}
        self.pml = {}

    def add(self, group: PhysicalGroup):
        if group.name not in self.groups:
            self.groups[group.name] = (len(self.groups), group)
            return group
        else:
            assert False, "Physical Group with name {} already exists".format(
                group.name
            )

    def add_pml(self, group: PhysicalGroup, pml: PML):
        g = self.add(group)
        self.pml[group.name] = pml

        return g, pml

    def write_material_input(self, path: str, n_pow: int = 2, a_pow: float = 10.0):
        with open(path, "w+") as f:
            f.write(f"{len(self.groups)}\n")

            for m in self.groups.values():
                f.write(f"{m[1].format()}\n")

            f.write("#\n#\n")

            for m in self.pml.values():
                f.write(f"{m.format(n_pow, a_pow)}\n")


class Point:
    def __init__(self, x: float, y: float, z: float = 0, _id: int = -1):
        self.x = x
        self.y = y
        self.z = z
        if _id != -1:
            self.id = _id
        else:
            self.id = self.__init_point()

    def get_id(self):
        return self.id

    @staticmethod
    def from_id(_id: int):
        geo.synchronize()
        [x, y, z] = model.getValue(0, _id, [])
        return Point(x, y, z, _id)

    def __str__(self):
        return f"({self.x}, {self.y})"

    def __repr__(self):
        return str(self)

    def offset(self, ox: float, oy: float, oz: float = 0):
        return Point(self.x + ox, self.y + oy, self.z + oz)

    def offset_x(self, offset: float):
        return Point(self.x + offset, self.y, self.z)

    def offset_y(self, offset: float):
        return Point(self.x, self.y + offset, self.z)

    def offset_z(self, offset: float):
        return Point(self.x, self.y, self.z + offset)

    def offset_angle_in_xy(self, angle: float, distance: float):
        return Point(
            self.x + distance * cos(angle), self.y + distance * sin(angle), self.z
        )

    def __init_point(self):
        return geo.add_point(self.x, self.y, self.z)


class Quad:
    def __init__(self, p1: Point, p2: Point, p3: Point, p4: Point, **kwargs):
        if "transfinite" in kwargs:
            transfinite = kwargs["transfinite"]
        else:
            transfinite = options["transfinite"]

        if "init_lines" in kwargs:
            init_lines = kwargs["init_lines"]
        else:
            init_lines = True

        self.p1 = p1
        self.p2 = p2
        self.p3 = p3
        self.p4 = p4

        if not init_lines:
            assert False, "not implemented"
        else:
            self.id, self.lines = self.__init_shape(transfinite)

    def __str__(self):
        return f"Quad(\n\t{self.p1},\n\t{self.p2},\n\t{self.p3},\n\t{self.p4}\n)"

    def __repr__(self):
        return str(self)

    def get_points(self):
        return [self.p1, self.p2, self.p3, self.p4]

    def __init_shape(self, transfinite: int) -> (int, [int, int, int, int]):
        points = [self.p1.id, self.p2.id, self.p3.id, self.p4.id]

        lines = [
            geo.add_line(points[0], points[1]),
            geo.add_line(points[1], points[2]),
            geo.add_line(points[2], points[3]),
            geo.add_line(points[3], points[0]),
        ]

        for i in lines:
            geo.mesh.set_transfinite_curve(i, transfinite)

        curve_loop_id = geo.add_curve_loop(lines)

        return curve_loop_id, lines

    @staticmethod
    def from_ids(points: [int, int, int, int], **kwargs):
        if "init_lines" in kwargs:
            init_lines = kwargs["init_lines"]
        else:
            init_lines = True

        return Quad(
            Point.from_id(points[0]),
            Point.from_id(points[1]),
            Point.from_id(points[2]),
            Point.from_id(points[3]),
            init_lines=init_lines,
        )


class QuadGroup:
    def __init__(self, quads: List[Quad], id: int = None):
        self.quads = quads
        if id is not None:
            self.id = id
            self.surfaces = []
        else:
            self.id, self.surfaces = self.__init_group()

    @staticmethod
    def from_gmsh(data: Tuple[int, int]):
        return QuadGroup([], data[1])

    def gmsh_surfaces(self):
        return [(2, i) for i in self.surfaces]

    def __init_group(self):
        surfaces = []
        for quad in self.quads:
            surface = geo.add_plane_surface([quad.id])
            surfaces.append(surface)
            geo.mesh.set_transfinite_surface(surface)
            geo.mesh.set_recombine(2, surface)

        return geo.add_surface_loop(surfaces), surfaces

    def __str__(self):
        return "QuadGroup({})".format(self.id)

    def __repr__(self):
        return str(self)

    def extrude_z(self, dz: float, elements=None):
        if elements is None:
            elements = [1]
        to_extrude = [(2, s) for s in self.surfaces]
        return Volume(
            geo.extrude(to_extrude, 0, 0, dz, elements, recombine=True),
            self,
            "z-" if dz < 0 else "z+",
        )


class Volume:
    def __init__(
        self,
        extrusion: List[Tuple[int, int]],
        extruded_quad_group: QuadGroup = None,
        direction: str = None,
    ):
        self.extrusion = extrusion
        self.extruded_quad_group = extruded_quad_group
        self.direction = direction
        self.faces = self.__get_faces_from_extrusion()
        self.ids = [extrusion[i][1] for i in range(1, len(extrusion), 5 + 1)]
        self.extruded_surfaces_count = self.__get_extruded_surfaces_count()

    def __get_faces_from_extrusion(self):
        if self.extruded_quad_group is None:
            return []

        faces = []
        for quad in self.extruded_quad_group.quads:
            if self.direction[0] == "z":
                extrusion_faces = {
                    "x+": QuadGroup.from_gmsh(self.extrusion[3]),
                    "x-": QuadGroup.from_gmsh(self.extrusion[5]),
                    "y+": QuadGroup.from_gmsh(self.extrusion[4]),
                    "y-": QuadGroup.from_gmsh(self.extrusion[2]),
                    "z+": (
                        QuadGroup.from_gmsh(self.extrusion[0])
                        if self.direction[1] == "+"
                        else QuadGroup.from_gmsh((2, quad.id))
                    ),
                    "z-": (
                        QuadGroup.from_gmsh((2, quad.id))
                        if self.direction[1] == "+"
                        else QuadGroup.from_gmsh(self.extrusion[0])
                    ),
                }

    def __str__(self):
        return f"Volume({self.ids})"

    def __get_extruded_surfaces_count(self):
        s = 0
        for i in range(1, len(self.extrusion)):
            if self.extrusion[i][0] == 3:
                s += 1

        return s

    def get_surfaces_in_extrusion_direction(self):
        return [self.extrusion[x] for x in range(0, len(self.extrusion), 5 + 1)]


class VolumeGroup:
    group: PhysicalGroup
    pml: PML | None

    def __init__(self, volumes: list[Volume]):
        self.volumes = volumes

    def set_phy_group(
        self, collection: PhysicalGroupCollection, name: str, material: Material
    ) -> PhysicalGroup:
        self.group = collection.add(
            PhysicalGroup(
                self.volumes,
                name,
                material,
            )
        )
        return self.group

    def set_pml(
        self,
        collection: PhysicalGroupCollection,
        name: str,
        material: Material,
        pml: PML,
    ):
        self.group, self.pml = collection.add_pml(
            PhysicalGroup(self.volumes, name, material, True), pml
        )


gmsh.initialize(sys.argv)

gmsh.model.add("BarrageSimple")

model = gmsh.model
geo = gmsh.model.geo


# def get_rectangle_from_size(dy: float, dz: float, oy: float = 0.0, oz: float = 0.0) -> [int, int, int, int]:
#     return [
#         geo.add_point(0, oy, oz),
#         geo.add_point(0, oy + dy, oz),
#         geo.add_point(0, oy + dy, oz + dz),
#         geo.add_point(0, oy, oz + dz)
#     ]


# def get_lines_from_points(points: [int, int, int, int]) -> [int, int, int, int]:
#     return [
#         geo.add_line(points[0], points[1]),
#         geo.add_line(points[1], points[2]),
#         geo.add_line(points[2], points[3]),
#         geo.add_line(points[3], points[0])
#     ]


# def get_curve_from_points(points: [int, int, int, int]) -> int:
#     return geo.add_curve_loop(get_lines_from_points(points))


# def add_volumes_to_phy_group(volumes, name: str) -> int:
#     if name not in physical_groups:
#         group = len(physical_groups) + 1
#         physical_groups[name] = group
#         model.addPhysicalGroup(3, [volume[1] for volume in volumes], group, name)
#         return group
#     else:
#         assert False, 'Cannot create a Physical Group with name {}: Group exists.'.format(name)


# def get_faces_in_extrusion_direction(created_volume, previous_extrusion_size):
#     if previous_extrusion_size == 1:
#         return [created_volume[0]]
#     # + 1 to account for the created volume reference
#     extrusions = [created_volume[x] for x in range(0, len(created_volume), previous_extrusion_size + 1)]
#     return extrusions


# def get_volumes(volumes, extrusions_count):
#     if extrusions_count == 1:
#         return [volumes[1]]
#
#     v = [volumes[x] for x in range(1, len(volumes), extrusions_count)]
#     return v
