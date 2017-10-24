# Discrete regression curves on rotation group SO(n)

Code accompanying the paper:

Interpolation and Regression of Rotation Matrices,

Nicolas Boumal, 2013

In: Nielsen F., Barbaresco F. (eds) Geometric Science of Information.
Lecture Notes in Computer Science, vol 8085. Springer, Berlin, Heidelberg

https://link.springer.com/chapter/10.1007/978-3-642-40020-9_37

PDF included on this GIT repository.

Code was adapted from the original in October 2017.

Requires Manopt, freely available at http://www.manopt.org


# About the parameters

Large w's (weights on individual control points) favor interpolation over smoothness (though you would need to let w go to infinity to get actual interpolation with the current code; there are better ways to interpolate.)

Large lambda (velocity regularization) favors shorter curves. If mu (acceleration regularization) is zero, then the result will be piecewise geodesic, but not exactly passing through the control points (it will take shortcuts.)

Large mu (acceleration regularization) will favor an almost geodesic regression curve. For mu > 0 but not necessarily large, "corners" will be penalized strongly, so that the curve will feel overall smoother (not piecewise geodesic).

Something to keep in mind is: the larger mu, the slower the solver is. This is because a large mu tends to translate into poor conditioning, unfortunately. Designing a good preconditioner for this is an open question.
