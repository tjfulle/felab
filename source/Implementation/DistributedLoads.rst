
Source terms and Distributed Loads
==================================

Source terms
------------

Overview
~~~~~~~~

Source terms:

-  are forces per unit length, area, or volume acting on the system in
   solid mechanics problems, or

-  heat sources in heat transfer problems.

Defining source terms
~~~~~~~~~~~~~~~~~~~~~

Source terms are defined in a two-dimensional list where is a sublist
defining the source term type, the region on which the source is active,
and the arguments required by the source term type. Source term types
can be one of

+-------------------+---------------------+-----------------------+------------------------------------------------------+
| Source type       | Symbolic constant   | Valid element types   | Arguments                                            |
+===================+=====================+=======================+======================================================+
| Heat generation   |                     |                       | Magnitude of the heat generation :math:`s`           |
+-------------------+---------------------+-----------------------+------------------------------------------------------+
| Gravity           |                     | All solid elements    | Material density, components of the gravity vector   |
+-------------------+---------------------+-----------------------+------------------------------------------------------+

The specifier is valid only on rectilinear problem domains and is one of
the following symbolic constants:

-  , , , , , or . ,, correspond to the :math:`x`, :math:`y`, and
   :math:`z` coordinate directions and the identifiers and to the
   corresponding low and high boundaries.

Examples
^^^^^^^^

Heat source of magnitude :math:`10` to all nodes in the domain

.. code:: python

    [S, ALL, 10]

Heat source given by to nodes on the bottom edge of a rectangular
domain:

.. code:: python

    def fun(x):
        return 2. * (1. + x[:,1]) / ((3. + x[:,0])**2 + (1 + x[:,1])**2)
    [S, JLO, fun]

Gravity load in the :math:`[1,0,0]` direction for a material having
density :math:`\rho=0`.

.. code:: python

    [GRAV, ALL, 10, 1, 0, 0]

Internal representation of source terms
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Internally, source terms for every degree of freedom in the model are
stored in a single array . is generated from the user defined
containersby the function , shown in Listing [lst:src.src.1]. is defined
in and is invoked as

.. code:: python

    f = FormatSourceTermsUDOF(coord, elemap, eletyp, elecon, src)

The arguments to have been described in previous sections. The output
from is:

+---------+-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| ``f``   | Nodal source contributions as an array of length nd\*n where nd is the number of degrees of freedom per node and n the total number of nodes. The interpretation of the source term magnitude is dependent on the element type.   |
+---------+-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+

.. code:: python

    def FormatSourceTermsUDOF(coord, elemap, eletyp, elecon, src):
        numele = elecon.shape[0]
        numnod = coord.shape[0]
        ndof = NumberOfDOFPerNode(eletyp[0])
        f = zeros(numnod*ndof)
        for fi in src:
            ftype = fi[0]
            if ftype == SOURCE:
                # Heat generation
                region, magnitude = fi[1:3]
                try:
                    magnitude(array([[1,1,1]]))
                    fun = magnitude
                except TypeError:
                    fun = lambda x: magnitude
                if region == ALL:
                    els = range(numele)
                else:
                    els = [elemap[el] for el in region]
                nodes = unique(elecon[els])
                f[nodes] += fun(coord[nodes])

            elif ftype == GRAV:
                region, density = fi[1:3]
                components = fi[3:]
                if len(components) != ndof:
                    raise ValueError('Incorrect number of gravity load components')
                if region == ALL:
                    els = range(numele)
                else:
                    els = [elemap[el] for el in region]
                nodes = unique(elecon[els])
                g = array(components) * density
                for node in nodes:
                    dofs = [ndof*node + i for i in range(ndof)]
                    f[dofs] += g

        return f

Distributed loads
-----------------

Overview
~~~~~~~~

Distributed loads:

-  apply forces per unit length/area/volume on element edges, faces, or
   volumes,

-  require that an appropriate distributed load type be specified.

Defining distributed loads
~~~~~~~~~~~~~~~~~~~~~~~~~~

Distributed loads are defined in a two-dimensional list where is a
sublist defining the distributed load type, the region on which the
distributed load is prescribed, and the arguments required by the
distributed load type. Distributed load types can be one of

+----------------------+---------------------+-----------------------+---------------------------------------------------------------------------------------------------------------+
| Load type            | Symbolic constant   | Valid element types   | Arguments                                                                                                     |
+======================+=====================+=======================+===============================================================================================================+
| Heat conduction      |                     |                       | Magnitude of the heat flux normal :math:`q_n`                                                                 |
+----------------------+---------------------+-----------------------+---------------------------------------------------------------------------------------------------------------+
| Surface convection   |                     |                       | Convection heat transfer constant :math:`h` and the temperature of the surrounding fluid :math:`T_{\infty}`   |
+----------------------+---------------------+-----------------------+---------------------------------------------------------------------------------------------------------------+
| Surface tractions    |                     | ,                     | magnitude of the surface traction.                                                                            |
+----------------------+---------------------+-----------------------+---------------------------------------------------------------------------------------------------------------+

The specifier defines the surfaces on which the distributed load is
applied. Define a surface by an element number and an edge identifier
label. The edge identifier is one of , , :math:`\ldots` , where is the
j\ :math:`^{\rm th}` element edge. Element edge numbers are defined in
Appendix [app:elem.lib]. For example, to define a surface composed of
edge 1 of elements 1, 2, and 3:

Alternatively, the regions

-  to represent all nodes on the boundary

-  , , , , , or . ,, correspond to the :math:`x`, :math:`y`, and
   :math:`z` coordinate directions and the identifiers and to the
   corresponding low and high boundaries.

can be specified to automatically generate the surface.

Examples
^^^^^^^^

Heat conduction with :math:`q_n=2000` on the top boundary of a
rectangular domain:

.. code:: python

    [QCOND, JHI, 2000]

Heat convection with :math:`h=250` and :math:`T_{\infty}=25` on the
bottom of a rectangular domain:

.. code:: python

    [QCONV, JLO, 250, 25]

Surface traction on edge 1 of elements 1, 2, and 5 with magnitude
:math:`qn=2000`.

.. code:: python

    [TRAC, [(1, S1), (2, S1), (5, S1)], 2000]

Internal representation of distributed loads
--------------------------------------------

Internally, distributed load data is stored in two arrays:

-  specifies surfaces on which distributed loads are applied.

-  specifies the magnitudes of distributed loads.

and are generated from the user defined container by the function ,
shown in Listing [lst:src.fmtdl]. is defined in and is invoked as

.. code:: python

    dltags, dlvals = FormatDistributedLoadsUDOF(coord, elemap, eletyp,
                                                elecon, dload)

The arguments to have described in previous sections.

The outputs from are the containers and . and are best explained by
example. Consider the triangular element with distributed load in Figure
[fig:src.dload]. The element connectivity is

.. code:: python

    elecon = [[10, 1, 2, 5, 6],
              [20, 7, 5, 2],
              [30, 7, 4, 5],
              [40, 7, 3, 4],
              [50, 7, 2, 3]]

The list is:

.. code:: python

    # Distributed load
    dload = [[T, [(40, S2)], qn]]

.. figure:: DLoadEx.png
   :align: center

   Distributed load example.

The corresponding and containers are:

.. code:: python

    dltags = [[0, 0, 0],
              [0, 0, 0],
              [0, 0, 0],
              [0, T, 0],
              [0, 0, 0]]
    dlvals = [[[0., 0., 0.], [0., 0., 0.], [0., 0., 0.]],
              [[0., 0., 0.], [0., 0., 0.], [0., 0., 0.]],
              [[0., 0., 0.], [0., 0., 0.], [0., 0., 0.]],
              [[0., 0., 0.], [qn, 0., 0.], [0., 0., 0.]],
              [[0., 0., 0.], [0., 0., 0.], [0., 0., 0.]]]

is the j\ :math:`^{\rm th}` edge of the i\ :math:`^{\rm th}` element.
indicates that no distributed load is applied over that edge and a
distributed is applied, otherwise.

is the k\ :math:`^{\rm th}` distributed load specifier for edge j of
element i.

.. code:: python

    def FormatDistributedLoadsUDOF(coord, elemap, eletyp, elecon, dload):
        from mesh import FindSurface
        numele, numnod = elecon.shape[0], coord.shape[0]
        # Surface load tags and vals
        numedge = max([EdgesPerElement(et) for et in eletyp])
        dltags = zeros((numele, numedge), dtype=int)
        dlvals = zeros((numele, numedge, 3))
        errors = 0
        for load in dload:
            loadtyp, region = load[:2]
            args = load[2:]
            nargs = len(args)
            if loadtyp in (QCOND, QCONV):
                # Surface heat flux
                if loadtyp == QCOND:
                    if nargs != 1:
                        logging.error('Expected 1 QCOND dload component')
                        errors += 1
                        continue
                if loadtyp == QCONV:
                    if nargs != 2:
                        logging.error('Expected 2 QCONV dload components')
                        errors += 1
                        continue
                surf = FindSurface(coord, elemap, eletyp, elecon, region)
            elif loadtyp == TRACTION:
                if nargs != 1:
                    logging.error('Expected 1 TRACTION dload component')
                    errors += 1
                    continue
                surf = FindSurface(coord, elemap, eletyp, elecon, region)
            else:
                raise ValueError('Unrecognized dload type')
            for (el, edge) in surf:
                dltags[el, edge] = loadtyp
                dlvals[el, edge][:nargs] = args
        if errors:
            raise Exception('Stopping due to previous errors')
        return dltags, dlvals
