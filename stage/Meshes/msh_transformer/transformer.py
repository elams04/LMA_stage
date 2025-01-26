import math
from typing import Protocol, List, Dict


class Node:
    materials: List[int] | None

    def __init__(self, _id: int, x: float, y: float, z: float, removed: bool = False):
        self.id = _id
        self.x = x
        self.y = y
        self.z = z
        self.is_removed = removed
        self.materials = None
        self.elements = []

    def add_element(self, element: "Element") -> None:
        self.elements.append(element)

    def get_faces_directions(self) -> List[str]:
        faces = []
        for elem in self.elements:
            for key, value in elem.surfaces.items():
                if self in value:
                    faces.append(key)

        return faces

    def transform(self, fn: "TransformerFunction", mat: List[int]):
        if self.is_removed:
            return
        self.x, self.y, self.z, self.is_removed = fn(
            self.x, self.y, self.z, mat, self.get_faces_directions()
        )

    def format(self):
        return f"{self.id} {self.x} {self.y} {self.z}"

    def apply_transform(self, tr: "Transformer") -> "Node":
        if tr.bounding_box.contains(self):
            self.x, self.y, self.z, self.is_removed = tr.fn(
                self.x, self.y, self.z, self.materials, self.get_faces_directions()
            )
        return self

    def is_outer(self):
        if self.materials is None:
            return False
        else:
            return len(self.materials) < 5

    def __repr__(self):
        if self.is_removed:
            return f"Point(..., {self.x}, {self.y}, {self.z})"

        return f"Point({self.id}, {self.x}, {self.y}, {self.z})"


class Materials:
    points: Dict[int, List[int]]

    def __init__(self):
        self.points = {}
        self.materials = []

    def add_element(self, element: "Element"):
        if element.material not in self.materials:
            self.materials.append(element.material)

        for point in element.points:
            if point.id in self.points:
                self.points[point.id].append(element.material)
            else:
                self.points[point.id] = [element.material]

            point.materials = self.points[point.id]


class Element:
    def __init__(
        self,
        _id: int,
        _type: int,
        prop1: int,
        material: int,
        elementary: int,
        points: List[Node],
    ):
        self.id = _id
        self.type = _type
        self.prop1 = prop1
        self.elementary = elementary
        self.material = material
        self.points = points
        self.is_removed = False
        self.surfaces = Element._init_surfaces(points)

    @staticmethod
    def _init_surfaces(points: List[Node]):
        surface_points = {
            "x+": [1, 5, 6, 2],
            "x-": [3, 7, 8, 4],
            "y+": [5, 6, 7, 8],
            "y-": [1, 2, 3, 4],
            "z+": [2, 6, 7, 3],
            "z-": [1, 5, 8, 4],
        }

        surfaces = {}

        for key, value in surface_points.items():
            surfaces[key] = [points[i - 1] for i in value]

        return surfaces

    def apply_transform(self, tr: "Transformer") -> "Element":
        if tr.element_filter_fn:
            if tr.element_filter_fn(self.material):
                self.is_removed = True
                for p in self.points:
                    if all(mat == self.material for mat in p.materials):
                        p.is_removed = True

        return self

    def format(self):
        return " ".join(
            [
                str(self.id),
                str(self.type),
                str(self.prop1),
                str(self.material),
                str(self.elementary),
            ]
            + [str(point.id) for point in self.points]
        )


class BoundingBox:
    def __init__(
        self,
        min_x: float,
        max_x: float,
        min_y: float,
        max_y: float,
        min_z: float,
        max_z: float,
        threshold=0.0,
    ):
        self.min_x = min_x
        self.max_x = max_x
        self.min_y = min_y
        self.max_y = max_y
        self.min_z = min_z
        self.max_z = max_z
        self.threshold = threshold
        self.is_infinite = math.isinf(min_x + min_y + min_z + max_x + max_y + max_z)

    @staticmethod
    def infinite():
        return BoundingBox(
            -math.inf,
            math.inf,
            -math.inf,
            math.inf,
            -math.inf,
            math.inf,
        )

    def render(self, m, offset: float = 0.0):
        if self.infinite:
            return

        b_xm_ym_zm = m.Point(
            self.min_x - offset, self.min_y - offset, self.min_z - offset
        )
        b_xm_ym_zp = m.Point(
            self.min_x - offset, self.min_y - offset, self.max_z + offset
        )
        b_xm_yp_zm = m.Point(
            self.min_x - offset, self.max_y + offset, self.min_z - offset
        )
        b_xm_yp_zp = m.Point(
            self.min_x - offset, self.max_y + offset, self.max_z + offset
        )
        b_xp_ym_zm = m.Point(
            self.max_x + offset, self.min_y - offset, self.min_z - offset
        )
        b_xp_ym_zp = m.Point(
            self.max_x + offset, self.min_y - offset, self.max_z + offset
        )
        b_xp_yp_zm = m.Point(
            self.max_x + offset, self.max_y + offset, self.min_z - offset
        )
        b_xp_yp_zp = m.Point(
            self.max_x + offset, self.max_y + offset, self.max_z + offset
        )

        m.Quad(b_xm_ym_zm, b_xm_ym_zp, b_xm_yp_zp, b_xm_yp_zm)
        m.Quad(b_xm_ym_zm, b_xm_ym_zp, b_xp_ym_zp, b_xp_ym_zm)
        m.Quad(b_xm_ym_zm, b_xp_ym_zm, b_xp_yp_zm, b_xm_yp_zm)
        m.Quad(b_xp_yp_zp, b_xp_yp_zm, b_xp_ym_zm, b_xp_ym_zp)
        m.Quad(b_xp_yp_zp, b_xp_yp_zm, b_xm_yp_zm, b_xm_yp_zp)
        m.Quad(b_xp_yp_zp, b_xm_yp_zp, b_xm_ym_zp, b_xp_ym_zp)

    def contains(self, point: Node) -> bool:
        return (
            self.min_x - self.threshold <= point.x <= self.max_x + self.threshold
            and self.min_y - self.threshold <= point.y <= self.max_y + self.threshold
            and self.min_z - self.threshold <= point.z <= self.max_z + self.threshold
        )

    @staticmethod
    def from_nodes(nodes: List[Node]):
        min_x, min_y, min_z = math.inf, math.inf, math.inf
        max_x, max_y, max_z = -math.inf, -math.inf, -math.inf

        for node in nodes:
            min_x = min(min_x, node.x)
            max_x = max(max_x, node.x)
            min_y = min(min_y, node.y)
            max_y = max(max_y, node.y)
            min_z = min(min_z, node.z)
            max_z = max(max_z, node.z)

        return BoundingBox(min_x, max_x, min_y, max_y, min_z, max_z)

    @staticmethod
    def get_dampening_coefficient(x: float, x_min: float, x_max: float, d=10.0):
        return max(0.0, 1 - (2 * (1 / (x_max - x_min)) * (x - x_min) - 1) ** d)

    def __repr__(self):
        return f"BoundingBox(\n\t{self.min_x}-{self.max_x}\n\t{self.min_y}-{self.max_y}\n\t{self.min_z}-{self.max_z}\n)"


class TransformerFunction(Protocol):
    def __call__(
        self, x: float, y: float, z: float, mat: List[int], on_faces: List[str]
    ) -> (float, float, float, bool): ...

    @staticmethod
    def identity(
        x: float, y: float, z: float, mat: List[int], on_faces: List[str]
    ) -> (float, float, float, bool):
        return x, y, z, False


class ElementFilterFunction(Protocol):
    def __call__(self, mat: int) -> bool: ...


class Transformer:
    def __init__(
        self,
        bounding_box: BoundingBox,
        fn: TransformerFunction,
        element_filter_fn: ElementFilterFunction | None = None,
    ):
        self.bounding_box = bounding_box
        self.fn = fn
        self.element_filter_fn = element_filter_fn

    @staticmethod
    def material_filter(mat: int) -> "Transformer":
        def filter_material(material: int):
            if mat < 0:
                return abs(mat) == material

            return mat != material

        return Transformer(
            BoundingBox.infinite(),
            TransformerFunction.identity,
            filter_material,
        )
