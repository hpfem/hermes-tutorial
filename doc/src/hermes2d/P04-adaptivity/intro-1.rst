Adaptive low-order FEM and hp-FEM
---------------------------------

Our team has a passion for adaptive higher-order finite element methods
because they are amazingly efficient compared to conventional low-order 
FEM. These methods have been an integral part of Hermes since its origin. 
Let us begin this tutorial chapter by explaining what gains 
the user can expect when using adaptive *hp*-FEM instead of conventional
adaptive low-order FEM.

Low-order FEM is inefficient
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Many practitioners are skeptical about using adaptive FEM because adaptivity can 
slow down and after some time virtually freeze a computation. The reason is that 
in conventional, algebraically-convergent FEM such as linear or quadratic, one needs to add more 
and more degrees of freedom to decrease the error. In fact the number of DOF needs 
to grow **exponentially** if the error is to decrease linearly. This phenomenon is 
problem-independent and it is due to poor approximation properties of low-order elements.
In contrast to this, adaptive *hp*-FEM achieves linear decrease of error by adding degrees of freedom 
**linearly**. The different convergence rates of adaptive FEM with linear and quadratic 
elements, and adaptive *hp*-FEM is shown below. Note that 
on the log-log scale, the slope of the convergence curve tells the convergence rate 
of the method. If the slope is constant then the method converges algebraically.


.. figure:: intro/conv_dof.png
   :align: center
   :scale: 60% 
   :figclass: align-center
   :alt: Typical convergence curves for adaptive linear FEM, quadratic FEM, and hp-FEM.

The reader can see that the 
linear FEM would need in the order of 1,000,000,000,000,000,000 degrees of freedom 
(DOF) to reach a level of accuracy where the hp-FEM is with less than 10,000 DOF. 
A similar effect can be observed in the CPU-time convergence graph:

.. figure:: intro/conv_cpu.png
   :align: center
   :scale: 60% 
   :figclass: align-center
   :alt: Typical convergence curves for adaptive linear FEM, quadratic FEM, and hp-FEM.

These convergence curves are typical representative examples, confirmed with
many numerical experiments of independent researchers, and supported with
theory. The low-order FEM is doomed by the underlying math -- its poor convergence cannot 
be fixed by designing smarter adaptivity algorithms.

What makes adaptive hp-FEM so fast
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In order to obtain fast, usable adaptivity (the red curve), one
has to resort to adaptive hp-FEM. The hp-FEM takes advantage of 
the following facts:

* Large high-order elements approximate smooth parts of the solution much more efficiently 
  than small linear ones. 
  We have couple of benchmarks such as `smooth-iso <http://hpfem.org/hermes/doc/src/hermes2d/benchmarks-general/smooth-iso.html>`_ 
  that illustrate this in an impressive way - please spend a few minutes to check them out. If this is your
  first interaction with adaptive hp-FEM, don't skip this and run the benchmark on your computer, since the 
  results are worth it. 
* This holds the other way round for singularities, steep gradients, oscillations and other "bad" spots: 
  These are approximated best using small low-order elements.
* In order to capture efficiently anisotropic solution behavior, one needs adaptivity algorithms 
  that can refine meshes anisotropically both in $h$ and $p$. This is illustrated 
  in several benchmarks such as 
  `smooth-aniso-x <http://hpfem.org/hermes/doc/src/hermes2d/benchmarks-general/smooth-aniso-x.html>`_  
  or `nist-07 (line singularity) <http://hpfem.org/hermes/doc/src/hermes2d/benchmarks-nist/nist-07.html>`_.

What does it take to do adaptive hp-FEM
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Automatic adaptivity in the hp-FEM is substantially different from adaptivity
in low-order FEM, since every element can be refined in many different ways.
The following figure shows several illustrative refinement candidates for 
a fourth-order element.

.. figure:: intro/refinements.png
   :align: center
   :scale: 50% 
   :figclass: align-center
   :alt: Examples of hp-refinements.

The number of allowed element refinements is implementation-dependent, but in general
it is very low in $h$ or $p$ adaptivity, much higher in hp-adaptivity, 
and it rises even more when quadrilateral elements are used and their anisotropic 
refinements are enabled. 

Eight hp-adaptivity modes in Hermes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Hermes has eight different adaptivity modes P_ISO, P_ANISO, H_ISO, H_ANISO,
HP_ISO, HP_ANISO_P, HP_ANISO_H, HP_ANISO. In this order, usually P_ISO yields the 
worst results and HP_ANISO the best. However, even P_ISO can be very efficient 
with a good a-priori locally refined starting mesh. 

The most general mode HP_ANISO considers around 100 refinement candidates 
for each element. The difference between the next best mode HP_ANISO_H
and HP_ANISO is only significant for problems that exhibit strong 
anisotropic behavior. The selection of the hp-refinement mode is 
where the user can use his a-priori knowledge of the problem to make 
the computation faster. 

Why do we need more than standard error estimates
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Due to the large number of refinement options in each element, classical error estimators that
provide just one number per element are not enough. To guide hp-adaptivity, one really needs 
to know the **shape** of the approximation error, not only its magnitude.

In analogy to the most successful adaptive ODE solvers,
Hermes uses a pair of approximations with different orders of accuracy 
to obtain this information. *fine mesh solution* and *orthogonal projection on 
a coarse submesh*. The initial coarse submesh mesh is read from the mesh 
file, and the initial fine mesh is created through its global refinement 
both in $h$ and $p$. The fine mesh solution is the approximation of interest 
both during the adaptive process and at the end of computation. 

Robustness of the reference solution approach
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Note that the reference solution approach is PDE independent, which is truly great 
for multiphysics coupled problems. Hermes does not use a single analytical error 
estimate or any other technique that would narrow down its applicability to selected 
equations or low-order FEM. 

Room for improvement
~~~~~~~~~~~~~~~~~~~~

An obvious disadvantage of the reference solution approach to automatic adaptivity is its higher 
computational cost, especially in 3D. We are aware of this fact and would not mind 
at all replacing the current paradigm with some cheaper technique -- as long as it is 
PDE-independent, works for elements of high orders, and handles anisotropy in both 
'h' and 'p'. Seemingly, however, no such alternatives exist. If you have any ideas, let 
us know.
