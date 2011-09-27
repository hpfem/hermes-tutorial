Essential and Natural Boundary Conditions
-----------------------------------------

Hermes distinguishes between **essential** and **natural** boundary conditions. 
The former type eliminates degrees of freedom from the domain's boundary
(the solution is prescribed) while the latter does not. 
Examples of essential boundary conditions are Dirichlet conditions for 
H1-problems and perfect conductor conditions in the space H(curl).
Examples of natural boundary conditions are Neumann or Newton (Robin) 
conditions for H1-problems and impedance conditions in the space 
H(curl). Only essential conditions are treated explicitly in Hermes, 
while the natural ones are defined using surface integrals 
in the weak formulation. Let us start with showing default ways 
to define essential boundary conditions.

Default constant essential boundary conditions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

We already met a default essential boundary condition in the previous example 03-poisson:

.. sourcecode::
    .

    // Initialize essential boundary conditions.
    DefaultEssentialBCConst<double> bc_essential(Hermes::vector<std::string>("Bottom", "Inner", "Outer", "Left"), 
                                                 FIXED_BDY_TEMP);

.. latexcode::
    .

    // Initialize essential boundary conditions.
    DefaultEssentialBCConst<double> bc_essential(Hermes::vector<std::string>("Bottom", "Inner", 
                                                 "Outer", "Left"), FIXED_BDY_TEMP);

This one assigned a constant value FIXED_BDY_TEMP to all boundary edges with the markers 
"Bottom", "Inner", "Outer" or "Left". 

After creating one or more essential boundary conditions, they are passed into the container 
class EssentialBCs::

    EssentialBCs<double> bcs(&bc_essential);

The purpose of this container is to collect all essential boundary conditions to be passed into a Space, 
as we shall see shortly. The constructor of the EssentialBCs class can accept a Hermes::vector of
essential boundary conditions. 

Default nonconstant essential boundary conditions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Another default essential boundary condition reads boundary values from a given 
function. This is useful, for example, in benchmarks with known exact solution. 
Typical usage is as follows::

    // Set exact solution.
    CustomExactSolution exact(&mesh, EXACT_SOL_P);

    // Initialize boundary conditions.
    DefaultEssentialBCNonConst<double> bc_essential("Bdy", &exact);
    EssentialBCs<double> bcs(&bc_essential);

This technique is used in the tutorial example P01/07-general and in several 
benchmarks that are part of the repository hermes-examples.

Custom essential boundary conditions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Custom essential conditions can be created by subclassing the abstract class
EssentialBoundaryCondition.

This class can represent arbitrary essential boundary conditions that depend 
on space and time. Every descendant of this class must redefine the purely 
virtual functions get_value_type() and value(). This will be explained in
more detail in the following example 04-bc-dirichlet.

