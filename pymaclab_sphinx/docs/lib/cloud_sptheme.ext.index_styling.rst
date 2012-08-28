==========================================================================
:mod:`cloud_sptheme.ext.index_styling` - improves css styling for genindex
==========================================================================

.. module:: cloud_sptheme.ext.index_styling
    :synopsis: adds additional css styling to general index

This Sphinx extension intercepts & modifies the general index data
before it is rendered to html, adding some additional css classes
to help Sphinx themes (e.g. :doc:`/cloud_theme`)
provide addtional per-type styling for index entries.

This extension adds the following css classes to ``genindex``:

* For all entries referencing an ``attribute``, ``method``, ``class``,
  ``function``, or ``module``:

  The type of the entry (e.g. ``attribute``) is wrapped in a
  :samp:`<span class="category {type}">` element.
  If the element has a location (e.g. ``myclass in module myapp``),
  the ``myapp`` portion is wrapped in a ``<span class="location">`` element.

* Entries which don't fit into one of the above categories are not modified.
