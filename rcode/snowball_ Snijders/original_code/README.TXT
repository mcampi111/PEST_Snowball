SNOWBALL is a program to compute estimates according to the methods
in

O. Frank & T.A.B. Snijders,
"Estimating the size of hidden populations using snowball sampling",
 Journal of Official Statistics, 10 (1994), 53 - 67.

An ASCII input file is required that gives each sampled vertex `i'
 and the vertices to which this vertex `i' has an arc.
Vertices are coded by numbers in the range 1 ... vmax
 (vmax depending on implementation; e.g., 1,000).
Each line must have the vertex code and then, separated by blanks,
 the codes of the other vertices to which an arrow from this vertex points.
This means that, for vertex `i', the number of codes on the line
 for this vertex is the outdegree of `i' plus 1.
Blanks after the last code on a line lead to input errors, so take care
 that there are no trailing blanks.
Empty lines also lead to input errors, so take care that there are
 no empty lines in the data file.
Small example data files are SN1.DAT and SN2.DAT.

One of the estimates (the Horvitz-Thompson estimate) makes use
 of the indegrees. Therefore it can be computed only if the
 indegrees are known (which is rarely the case).
IF indegrees are also available, these can be given in a separate input file
 (format: each line has first the vertex code and then the
 corresponding indegree),
 and then the Horvitz-Thompson estimators also will be computed.
An example data files for indegrees to be used with SN2.DAT is SN2.DIN.
