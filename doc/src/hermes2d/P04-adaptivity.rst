Tutorial Part IV (Automatic Adaptivity)
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

    P04-adaptivity/intro-1   
    P04-adaptivity/intro-2
    P04-adaptivity/01-intro
    P04-adaptivity/intro-3
    P04-adaptivity/02-system
    P04-adaptivity/03-complex
    P04-adaptivity/04-hcurl
    P04-adaptivity/05-exact
    P04-adaptivity/06-nonlinear
    P04-adaptivity/07-transient-space-only
    P04-adaptivity/08-transient-time-only
    P04-adaptivity/09-transient-space-and-time

 








