# The [[scuff-em]] documentation system

The documentation for [[scuff-em]] is maintained in the form
of plain-text (markdown) files stored in the ``doc`` subdirectory
of the [[scuff-em]] repository. These files are processed by
the wonderful open-source [MkDocs](http://www.mkdocs.org) system
to build the web-based documentation hierarchy.

A major reason for this choice of documentation system is to
make it easy for [[scuff-em]] users to edit the existing 
documentation and to contribute new documentation. If you 
discover incorrect or incomplete portions of the documentation,
or if you would like to add new documentation, please consider
contributing to [[scuff-em]] by making the relevant changes
in the ``doc`` subdirectory of your [[scuff-em]] repository,
then [submitting a pull request on GitHub.][PullRequest]

## Building and serving your own local copy of the documentation

It's easy to build your own local copy of the entire [[scuff-em]]
documentation hierarchy, which you can then view offline.
This allows you to access the documentation without Internet 
access, and also to preview any changes you might make to
the documentation before submitting them in the form of a pull
request.

### One-time only setup operations

To build the [[scuff-em]] documentation, you will need a [[python]]
installation on your system, and you will need the ``mkdocs`` and 
``python-markdown-math`` packages. On my system I was able to 
install these using the following commands:

````bash
% sudo pip install mkdocs
% sudo pip install --upgrade pyinotify
% git clone https://github.com/mitya57/python-markdown-math.git
% cd python-markdown-math 
% python setup.py build
% sudo python setup.py install
````

### Building the documentation

To build the documentation, starting from the top-level
directory of your [[scuff-em]] repository, simply say

````bash
% cd doc
% mkdocs build
````

This will create the HTML hierarchy in the subdirectory
``doc/site.``

### Serving the documentation

A wonderfully convenient feature provided by MkDocs is the
ability to serve your local version of the documentation 
locally to a web browser running on your machine, without
having to mess around with configuring apache or any other
webserver software. To do this, starting from the ``doc`` 
subdirectory of the [[scuff-em]] repository you simply go
like this:

````bash
% cd doc
% mkdocs serve
````

Then direct your favorite web browser to the site
``127.0.0.1:8000``, i.e.

````bash
% google-chrome 127.0.0.1:8000
````

This should pull up the top-level page of the [[scuff-em]]
documentation tree, with internal links pointing to 
your local copies of the various pages.

The great thing about this is that, whenever you save changes 
to a file in the ``doc`` subdirectory of your [[scuff-em]] 
repository, the documentation is automatically rebuilt and
the webserver automatically refreshed. This allows you to
edit the [[scuff-em]] documentation in WYSIWIG fashion, 
simply by working on an ``.md`` file in a text editor in one window,
while having a web browser open to the corresponding subpage
of ``127.0.0.1:8000`` in another window. Then, every time
you save changes to the text file, the web page is automatically
updated! I love this brilliant system.

[PullRequest]: https://help.github.com/articles/creating-a-pull-request/
