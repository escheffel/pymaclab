.. index:: cloud; sphinx theme, sphinx theme; cloud

========================
PyMacLab Tutorial Series
========================

The Python DSGE instance
========================

*Introduction*

  As already stated in the introduction of the introductory basic tutorial, PyMacLab's strength or orginal design goal has been that of providing
  users with a rich and flexible DSGE data structure (also called `Class` in object-oriented programming speak) which allows them to do lots of
  interesting things with DSGE models and to treat them as if they were some kind of primitive `data type` in their own right.
  While the previous tutorial described some basics as well as the all-important DSGE model file structure and syntax conventions,
  in this section I am going to stress some of the object-oriented programming features of PyMacLab, in particular the
  structure of a PyMacLab DSGE model `instance` or data structure.

  Readers with a background in modern programming languages supporting the object-oriented programming (OOP) paradigm will easily relate to
  concepts in this sections, while for others it may appear more cryptic at first sight. But to convey these concepts to researchers is
  important, as it stresses many particular advantages of PyMacLab over other programs, and in particular its `flexibility`, `transparency`,
  `consistency`, `persistence` and enormous scope for `extensibility`. All example code fragments provided here assume that you are replicating
  them within an IPyton interactive session, but they could also be called from a Python program "batch" file.

