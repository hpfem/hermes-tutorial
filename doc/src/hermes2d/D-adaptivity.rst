Tutorial Part D (Automatic Adaptivity)
=======================================

So far we have not paid any attention to the accuracy of results. In general, 
computations on a fixed mesh are not very accurate and they should not be trusted. 
Controlled accuracy can be achieved through *adaptive mesh refinement (AMR)*. 
Advanced adaptive higher-order FEM (*hp*-FEM) is one of the main 
strengths of Hermes. With eight modes of automatic *hp*-adaptivity, Hermes 
excels at delivering highly accurate results with much lower numbers of degrees 
of freedom than conventional FEM codes.

.. toctree::
    :maxdepth: 2

    D-adaptivity/intro-1   
    D-adaptivity/intro-2
    D-adaptivity/01-intro
    D-adaptivity/01-intro-matrix-free
    D-adaptivity/02-kelly
    D-adaptivity/intro-3
    D-adaptivity/03-system
    D-adaptivity/04-complex
    D-adaptivity/05-hcurl
    D-adaptivity/06-exact
    D-adaptivity/07-nonlinear
    D-adaptivity/08-transient-space-only
    D-adaptivity/09-transient-time-only
    D-adaptivity/10-transient-space-and-time








