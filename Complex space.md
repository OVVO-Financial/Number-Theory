# Factors in the complex space

Below is the complex space associated with `n=309`.  Each triangle represents a complex number whose corresponding real numbers lie within the interval `[3,sqrt(n)]`.

Blue triangles are complex numbers whose real product is less than `n`.

Red triangles are complex numbers whose real product is greater than `n`.

The green triangle represents the complex number associated with the factors of `n=309`, `53+50i` which has corresponding real numbers `[3,103]`.

![Complex Space](https://github.com/OVVO-Financial/Number-Theory/blob/master/Images/Complex%20plane.png)

## Possible complex factors
All possible factors reside along the intersection of red & blue triangles.  We can highlight this strip in green.  The question becomes, **_how do we best navigate this strip?**_

![Factor Strip](https://github.com/OVVO-Financial/Number-Theory/blob/master/Images/Factor%20Strip%20in%20Green.jpeg)

### Trial Division
Trial division diagonally checks all complex points as illuastrated below, where all of the highlighted complex numbers have a lower real number = 3.

![Trial Div](https://github.com/OVVO-Financial/Number-Theory/blob/master/Images/Trial%20Division%20by%203.jpeg)

### Fermat
Fermat's difference of squares `N = a^2 - b^2` vertically checks all complex points illustrated below, where all of the highlighted complex numbers are equal to `a` in `b^2 = a^2 - N`.

![Fermat](https://github.com/OVVO-Financial/Number-Theory/blob/master/Images/Vertical%20Fermat.jpeg)
