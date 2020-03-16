import os

from felab.mesh import Mesh
from felab.error import UserInputError


def VTU2Genesis(nodtab=None, eletab=None, filename=None):
    if filename is None:
        assert nodtab is not None and eletab is not None
        outfile = "mesh.g"
    elif not os.path.isfile(filename):
        assert nodtab is not None and eletab is not None
        assert filename.endswith(".g")
        outfile = filename
        filename = None
    else:
        assert nodtab is None and eletab is None
        outfile = os.path.splitext(filename)[0] + ".g"
    try:
        mesh = Mesh(nodtab=nodtab, eletab=eletab, filename=filename)
    except KeyError:
        return
    mesh.to_genesis(outfile)


def INP2Genesis(filename):
    lines = open(filename).readlines()
    kw, name = None, None
    nodtab = []
    eletab = []
    nodesets = {}
    elemsets = {}
    element_blocks = {}
    for line in lines:
        line = ",".join([x.strip() for x in line.split(",")])
        if line.startswith("**"):
            continue
        if not line.split():
            continue
        if line.startswith("*"):
            name = None
            line = line.split(",")
            kw = line[0][1:].lower()
            opts = {}
            for opt in line[1:]:
                k, v = opt.split("=")
                opts[k.strip().lower()] = v.strip()
            if kw != "element":
                name = None
            if kw == "element":
                name = opts.get("elset")
                if name is None:
                    raise UserInputError("requires elements be put in elset")
                element_blocks[name.upper()] = []
            elif kw == "nset":
                name = opts["nset"]
                nodesets[name.upper()] = []
            elif kw == "elset":
                elemsets[name.upper()] = []
            continue
        if kw is None:
            continue
        if kw == "node":
            line = line.split(",")
            nodtab.append([int(line[0])] + [float(x) for x in line[1:]])
            continue
        elif kw == "element":
            eledef = [int(n) for n in line.split(",")]
            eletab.append(eledef)
            element_blocks[name].append(eledef[0])
            continue
        elif kw == "nset":
            nodesets[name.upper()].extend(
                [int(n) for n in line.split(",") if n.split()]
            )
            continue
        elif kw == "elset":
            elemsets[name.upper()].extend(
                [int(n) for n in line.split(",") if n.split()]
            )
            continue
    mesh = Mesh(nodtab=nodtab, eletab=eletab)
    for (name, elems) in element_blocks.items():
        mesh.element_block(name, elems)
    for (name, nodes) in nodesets.items():
        mesh.node_set(name, nodes)
    for (name, elems) in elemsets.items():
        mesh.element_set(name, elems)
    outfile = os.path.splitext(filename)[0] + ".g"
    mesh.to_genesis(outfile)
    return
