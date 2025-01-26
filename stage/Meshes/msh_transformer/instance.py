from typing import List

from Meshes.msh_transformer.transformer import (
    Materials,
    Node,
    Element,
    Transformer,
    BoundingBox,
)


class MSHInstance:
    _bounding_box: BoundingBox | None

    def __init__(self, input_path: str):
        self.input_path = input_path
        self._bounding_box = None
        (
            self.file_content,
            self.output_template,
            self.nodes,
            self.elements,
            self.materials,
        ) = MSHInstance._read_mesh_file(input_path)

    @staticmethod
    def _normalize(content: str) -> str:
        return content.replace("\r\n", "\n").replace("\r", "\n").replace("  ", " ")

    @staticmethod
    def _read_mesh_file(file_path: str):
        with open(file_path, "r") as f:
            content = f.read()

            content = MSHInstance._normalize(content)

            if not content.startswith("$MeshFormat\n2.2"):
                raise ValueError("Invalid mesh file format, should be .msh 2.2")

            output_template = ""
            materials = Materials()

            def parse_node(line: str) -> Node:
                parts = line.split(" ")
                return Node(
                    int(parts[0]), float(parts[1]), float(parts[2]), float(parts[3])
                )

            nodes_dict = {}
            elements = []
            is_parsing_nodes = False
            is_parsing_elements = False

            lines = content.split("\n")

            def parse_element(line: str) -> Element:
                try:
                    parts = line.split(" ")
                    element = Element(
                        int(parts[0]),
                        int(parts[1]),
                        int(parts[2]),
                        int(parts[3]),
                        int(parts[4]),
                        [nodes_dict[int(part)] for part in parts[5:]],
                    )
                    for part in parts[5:]:
                        nodes_dict[int(part)].add_element(element)
                    materials.add_element(element)
                    return element
                except ValueError:
                    raise ValueError(f'Failed to parse element: "{line}"')

            i = 0
            while i < len(lines):
                if not is_parsing_nodes and not is_parsing_elements:
                    output_template += lines[i] + "\n"

                if lines[i].startswith("$Nodes"):
                    is_parsing_nodes = True
                    i += 1  # skip the node count
                elif lines[i].startswith("$EndNodes"):
                    is_parsing_nodes = False
                    output_template += "{nodes}\n"
                    output_template += lines[i] + "\n"
                elif is_parsing_nodes:
                    node = parse_node(lines[i])
                    nodes_dict[node.id] = node

                if lines[i].startswith("$Elements"):
                    is_parsing_elements = True
                    i += 1  # skip the elements count
                elif lines[i].startswith("$EndElements"):
                    is_parsing_elements = False
                    output_template += "{elements}\n"
                    output_template += lines[i] + "\n"
                elif is_parsing_elements:
                    elements.append(parse_element(lines[i]))

                i += 1

            return (
                content,
                output_template,
                list(nodes_dict.values()),
                elements,
                materials,
            )

    def get_bounding_box(self):
        return BoundingBox.from_nodes(self.nodes)

    def save(
        self, output_path: str, save_nodes: bool = True, save_elements: bool = True
    ):
        if save_nodes or save_elements:
            remaining_nodes = [node.format() for node in self.nodes]
        else:
            remaining_nodes = []

        output_nodes = f"{len(remaining_nodes)}\n" + "\n".join(remaining_nodes)

        if save_elements:
            remaining_elements = [
                element.format() for element in self.elements if not element.is_removed
            ]
        else:
            remaining_elements = []

        output_elements = f"{len(remaining_elements)}\n" + "\n".join(remaining_elements)

        output_content = self.output_template.format(
            nodes=output_nodes, elements=output_elements
        )

        with open(output_path, "w") as f:
            f.write(output_content)

        return output_path

    def apply_transform(self, transformers: List[Transformer | None], output_path: str = None, progress_steps=10):
        for i in range(len(transformers)):
            if transformers[i] is None:
                continue

            print(f'\n# applying transformer {i + 1}')
            print('# transforming nodes')
            new_nodes = []
            N = len(self.nodes)
            print(f'{N} nodes')
            steps = N // progress_steps
            for j in range(N):
                if j % steps == 0:
                    print(f'{j / N * 100:.0f}%')
                new_nodes.append(self.nodes[j].apply_transform(transformers[i]))

            self.nodes = new_nodes

            print('# transforming elements')
            new_elements = []
            E = len(self.elements)
            print(f'{E} elements')
            steps = E // progress_steps
            for j in range(E):
                if j % steps == 0:
                    print(f'{j / E * 100:.0f}%')
                new_elements.append(self.elements[j].apply_transform(transformers[i]))

        if output_path:
            self.save(output_path, True, True)

        return self

    def get_outer_points(self):

        return [node for node in self.nodes if not node.is_removed and node.is_outer()]
