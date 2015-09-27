# Notes on SCUFF-EM documentation system

## Homer Reid 20150603

The SCUFF-EM documentation consists of a bunch of
individual files written in plain-text format,
which are eventually interpreted as 'markdown'
language to produce HTML syntax for the 
web-based documentation.

However, to contribute to the documentation you
do not need to know anything about markdown or 
HTML or anything! In general you can just edit
the plain-text files, and your changes will
be reflected on the web-based documentation
after the documentation build process has been
run.

On the other hand, if you wish to build your own
local version of the documentation and view it 
on your own machine, this is easy to do, by 
following these steps:

**(a)** * One-time only setup operations: *

````bash
% sudo pip install mkdocs
% sudo pip install --upgrade pyinotify
% git clone https://github.com/mitya57/python-markdown-math.git
% cd python-markdown-math 
% python setup.py build
% sudo python setup.py install
````

**(b)** *To launch the local documentation server:*

(from the top-level SCUFF-EM repository)

````bash
% cd doc
% mkdocs serve
````bash

**(c)** Then launch your favorite web browser and
navigate to the following URL:
 http://127.0.0.1:8000

## Other notes

* I couldn't figure out how to get small-caps text
in markdown, so I developed a kluge, described in 
`scuff-em-mkdocs-theme/css/theme_extra.css.`

* To get inline math working, I needed to install

  1. Install python-markdown-math:

````bash
 % git clone https://github.com/mitya57/python-markdown-math.git
 % cd python-markdown-math 
 % python setup.py build
 % sudo python setup.py install
````

 % git clone https://github.com/mitya57/python-markdown-math.git

  2. Add this line to doc/mkdocs.yml

extra_javascript: ['https://cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-AMS_HTML','mathjaxhelper.js']

  3. Copy the following files from the python-markdown-math package to `doc/docs:`
    -- `setup.py`
    -- `mdx_math.py`

  4. Create the file `doc/docs/mathjaxhelper.js` with the following content:

````
MathJax.Hub.Config({
  config: ["MMLorHTML.js"],
  jax: ["input/TeX", "output/HTML-CSS", "output/NativeMML"],
  extensions: ["MathMenu.js", "MathZoom.js"]
});
````

* The mkdocs system puts the file hierarchy it builds in 
the subdirectory doc/site off the top-level scuff-em 
folder. The files in this directory are considered the
results of a build operation and thus are not included
among the version-controlled files in the scuff-em git
repository.

However, there is a *separate* git repository, rooted in
doc/site, in which the built files *are* version controlled.
This repository corresponds to the GitHub repository
homerreid/scuff-em-documentation, which is a *separate
repository* from the main scuff-em repository 
(homerreid/scuff-em). 

The reason for this is that the "GitHub Pages" service
wants a git repository for which (a) the root directory 
is the root directory of the documentation tree, and 
(b) the repository is committed to the special branch
"gh-pages." I don't want items (a) or (b) to be true 
for the main scuff-em repository, so it's easier to
have a separate repository for documentation alone.
