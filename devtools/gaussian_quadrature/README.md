These files were downloaded from https://pomax.github.io/bezierinfo/legendre-gauss.html on 9-3-2025. They are the abscissae and weights for doing Gaussian Quadrature integration on the integral from -1 to 1. 

lgvalues-abscissa.php

lgvalues-weights.php

To convert to interval 0 to 1:

Abscissa01_0 = 0.5 * (Abscissa11 + 1)

Abscissa01_1 = 1 - Abscissa01_0

Weight01 = 0.5 * Weight11
