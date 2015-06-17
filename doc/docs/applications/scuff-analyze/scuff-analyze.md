# scuff-analyze

The documentation for [[scuff-analyze]] has not yet been
ported from its earlier version. For the time being, please
[access the earlier version of the documentation.][EarlierVersion]

<a name="RenameMesh"></a>
### `bash` script to rename a mesh file according to its number of edges

I call this script `RenameMesh.` Usage example:

````bash
 % RenameMesh Paraboloid.msh
 Moved Paraboloid.msh to Paraboloid_750.msh.
````

````bash
#!/bin/bash

if [ $# -eq  0 ]
then
  echo "usage: $0 Object.msh"
fi

NUMEDGES=`scuff-analyze --mesh $1  | grep 'interior edges' | head -1 | cut -f2 -d' '`
if [ "x${NUMEDGES}" == "x" ]
then 
  echo "Something's gone wrong."
  exit
fi
  
/bin/mv $1 ${1%%.msh}_${NUMEDGES}.msh

echo "Moved $1 to ${1%%.msh}_${NUMEDGES}.msh."
````

[EarlierVersion]: http://homerreid.com/scuff-em/scuff-analyze
